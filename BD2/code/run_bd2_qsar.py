import os,csv,json,warnings,zipfile, math
import numpy as np,pandas as pd
from rdkit import Chem,DataStructs,RDLogger
from rdkit.Chem import AllChem,rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.model_selection import GroupShuffleSplit,GroupKFold,GridSearchCV
from sklearn.svm import SVR
from sklearn.metrics import r2_score,mean_squared_error,mean_absolute_error
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.base import clone
from sklearn.dummy import DummyRegressor
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
RDLogger.DisableLog('rdApp.*');warnings.filterwarnings('ignore')
F='/mnt/data/BRD4_BD2_STRICT_MODELING_775.csv'; out='/mnt/data/BRD4_BD2_QSAR_RESULTS';os.makedirs(out,exist_ok=True)
df=pd.read_csv(F).copy(); mols=[Chem.MolFromSmiles(s) for s in df.canonical_isomeric_smiles];assert all(m is not None for m in mols)
fpgen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
fps=[fpgen.GetFingerprint(m) for m in mols]
X=np.zeros((len(mols),2048),dtype=np.float32)
for i,fp in enumerate(fps):DataStructs.ConvertToNumpyArray(fp,X[i])
y=df.pIC50.to_numpy(dtype=float)
scaff=[MurckoScaffold.MurckoScaffoldSmiles(mol=m,includeChirality=False) or 'NO_SCAFFOLD_'+df.canonical_isomeric_smiles.iloc[i] for i,m in enumerate(mols)]
df['murcko_scaffold']=scaff
holdout=GroupShuffleSplit(n_splits=1,test_size=.20,random_state=2024)
train_idx,test_idx=next(holdout.split(X,y,scaff))
assert not (set(np.array(scaff)[train_idx]) & set(np.array(scaff)[test_idx]))
pipe=Pipeline([('svr',SVR(kernel='rbf',epsilon=.1))]); params={'svr__C':[1,10,100], 'svr__gamma':['scale',0.01,0.1], 'svr__epsilon':[0.1,0.2]}
gcv=GroupKFold(n_splits=5)
grid=GridSearchCV(pipe,params,cv=gcv.split(X[train_idx],y[train_idx],np.array(scaff)[train_idx]),scoring='neg_root_mean_squared_error',n_jobs=-1,refit=True)
grid.fit(X[train_idx],y[train_idx]); model=grid.best_estimator_
def metric(yt,yp):return {'n':len(yt),'r2':float(r2_score(yt,yp)),'rmse':float(np.sqrt(mean_squared_error(yt,yp))),'mae':float(mean_absolute_error(yt,yp)),'bias':float(np.mean(yp-yt))}
pred_train=model.predict(X[train_idx]); pred_test=model.predict(X[test_idx]);train_metrics=metric(y[train_idx],pred_train); test_metrics=metric(y[test_idx],pred_test)
# repeated scaffold splits exclusively within development subset; fixed external holdout not reused
rows=[];dev_scaff=np.array(scaff)[train_idx]
for seed in range(20):
 a,b=next(GroupShuffleSplit(n_splits=1,test_size=.20,random_state=1000+seed).split(X[train_idx],y[train_idx],dev_scaff))
 m=clone(model).fit(X[train_idx][a],y[train_idx][a]); pred=m.predict(X[train_idx][b]);d=metric(y[train_idx][b],pred)
 rows.append({'repeat':seed+1,'train_n':len(a),'test_n':len(b),'train_scaffolds':len(set(dev_scaff[a])),'test_scaffolds':len(set(dev_scaff[b])),**d})
rep=pd.DataFrame(rows);rep.to_csv(out+'/repeated_scaffold_20x_metrics.csv',index=False)
# y randomization across development; evaluate on same heldout for negative control, no tuning
randoms=[]
for seed in range(30):
 shuffled=np.random.default_rng(seed+991).permutation(y[train_idx]);m=clone(model).fit(X[train_idx],shuffled);p=m.predict(X[test_idx]);randoms.append({'seed':seed+991,**metric(y[test_idx],p)})
pd.DataFrame(randoms).to_csv(out+'/y_randomization_30x.csv',index=False)
# nearest-neighbor Tanimoto to dev train for all external holdout compounds
sim=[]
for ix in test_idx:
 similarities=DataStructs.BulkTanimotoSimilarity(fps[ix],[fps[j] for j in train_idx])
 sim.append(max(similarities))
predictions=df.iloc[test_idx].copy();predictions['observed_pIC50']=y[test_idx];predictions['predicted_pIC50']=pred_test
predictions['predicted_IC50_nM']=10**(9-pred_test);predictions['residual_pIC50']=pred_test-y[test_idx]
predictions['nearest_training_tanimoto']=sim;predictions['AD_ge_0_4']=np.array(sim)>=.4;predictions['AD_ge_0_5']=np.array(sim)>=.5
predictions.to_csv(out+'/heldout_predictions_with_AD.csv',index=False)
summary={'input_rows':int(len(df)),'unique_scaffolds':int(len(set(scaff))),'n_train':int(len(train_idx)),'n_holdout':int(len(test_idx)),'train_scaffolds':len(set(np.array(scaff)[train_idx])),'holdout_scaffolds':len(set(np.array(scaff)[test_idx])),'best_params':grid.best_params_,'cross_validation_best_rmse':-float(grid.best_score_),'train':train_metrics,'fixed_scaffold_holdout':test_metrics,'repeated_scaffold_20x':{'r2_mean':float(rep.r2.mean()),'r2_sd':float(rep.r2.std()),'rmse_mean':float(rep.rmse.mean()),'rmse_sd':float(rep.rmse.std()),'mae_mean':float(rep.mae.mean()),'r2_min':float(rep.r2.min()),'r2_max':float(rep.r2.max())},'AD_holdout':{'n_ge_0_4':int((np.array(sim)>=.4).sum()),'n_ge_0_5':int((np.array(sim)>=.5).sum()),'r2_ge_0_4':metric(y[test_idx][np.array(sim)>=.4],pred_test[np.array(sim)>=.4]) if (np.array(sim)>=.4).sum()>1 else None,'r2_ge_0_5':metric(y[test_idx][np.array(sim)>=.5],pred_test[np.array(sim)>=.5]) if (np.array(sim)>=.5).sum()>1 else None},'y_randomization':{'holdout_r2_mean':float(np.mean([d['r2'] for d in randoms])),'holdout_r2_max':float(np.max([d['r2'] for d in randoms]))},'limitations':['Preliminary structure-level computational audit; not all source publications manually confirmed','Median IC50 aggregation across potentially different assays','Scaffold split and tuning on development only; external holdout evaluated once','AD similarity is an operational threshold, not calibrated confidence']}
with open(out+'/summary.json','w') as f:json.dump(summary,f,indent=2)
with open(out+'/README.txt','w') as f:f.write('BRD4-BD2 QSAR: ECFP4 2048-bit + RBF-SVR; fixed scaffold holdout (20%, seed 2024); 5-fold group CV grid search on train; 20 repeated scaffold splits strictly in development set; 30 shuffled-target negative controls; maximum Tanimoto similarity AD at 0.4/0.5. The input remains provisional and is not independently experimental-source verified.\n')
fig,ax=plt.subplots(figsize=(6,5));ax.scatter(y[test_idx],pred_test,s=22,alpha=.7);mmin=min(y[test_idx].min(),pred_test.min());mmax=max(y[test_idx].max(),pred_test.max());ax.plot([mmin,mmax],[mmin,mmax],linestyle='--',color='black');ax.set(xlabel='Observed pIC50',ylabel='Predicted pIC50',title='BRD4-BD2 fixed scaffold holdout');ax.text(.04,.96,f"R² = {test_metrics['r2']:.3f}\nRMSE = {test_metrics['rmse']:.3f}\nn = {len(test_idx)}",transform=ax.transAxes,va='top');fig.tight_layout();fig.savefig(out+'/holdout_observed_vs_predicted.png',dpi=220);plt.close(fig)
fig,ax=plt.subplots(figsize=(6,4));ax.hist(rep.r2,bins=10);ax.set(xlabel='R²',ylabel='Number of scaffold splits',title='20 repeated scaffold validations');fig.tight_layout();fig.savefig(out+'/repeated_scaffold_r2_distribution.png',dpi=220);plt.close(fig)
import joblib
joblib.dump({'pipeline':model,'fingerprint':{'type':'Morgan ECFP4','radius':2,'bits':2048},'training_smiles':df.iloc[train_idx].canonical_isomeric_smiles.tolist()},out+'/svr_model.joblib')
with zipfile.ZipFile('/mnt/data/BRD4_BD2_QSAR_RESULTS.zip','w',zipfile.ZIP_DEFLATED) as z:
 for name in os.listdir(out):z.write(out+'/'+name,arcname=name)
print(json.dumps(summary,indent=2))