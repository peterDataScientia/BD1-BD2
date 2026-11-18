import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from pathlib import Path
from zipfile import ZipFile
src=Path('/mnt/data/BRD4_BD2_AD_CALIBRATION/development_repeated_scaffold_oof_predictions.csv')
out=Path('/mnt/data/BRD4_BD2_VALIDATION_ERRORS');out.mkdir(exist_ok=True)
d=pd.read_csv(src)
s='nearest_fold_training_tanimoto';e='absolute_error_pIC50'
d['signed_error_pIC50']=d.predicted_pIC50-d.observed_pIC50
d['sq_error']=d.signed_error_pIC50**2
bins=np.array([0,.2,.3,.4,.5,.6,.7,.8,.9,1.000001]);labels=['0.00–0.20','0.20–0.30','0.30–0.40','0.40–0.50','0.50–0.60','0.60–0.70','0.70–0.80','0.80–0.90','0.90–1.00']
d['similarity_bin']=pd.cut(d[s],bins=bins,labels=labels,right=False,include_lowest=True)
g=d.groupby('similarity_bin',observed=True)
stats=g.agg(n_predictions=(e,'size'),n_unique_compounds=('molecule_chembl_ids','nunique'),mean_similarity=(s,'mean'),MAE=(e,'mean'),mean_signed_error=('signed_error_pIC50','mean'),median_abs_error=(e,'median'),p90_abs_error=(e,lambda x:np.percentile(x,90)),RMSE=('sq_error',lambda x:np.sqrt(x.mean()))).reset_index()
# repeat-block bootstrap (accounts for dependencies within each scaffold split, but not shared compounds across repeats)
rng=np.random.default_rng(2024);reps=np.sort(d['repeat'].unique());rb={r:d[d['repeat']==r] for r in reps};boots=[]
for j in range(1200):
 sampled=rng.choice(reps,len(reps),replace=True)
 dat=pd.concat([rb[r] for r in sampled],ignore_index=True)
 boots.append(dat.groupby('similarity_bin',observed=True)[e].mean())
b=pd.DataFrame(boots)
for _,r in stats.iterrows():
 col=r.similarity_bin;v=b[col].dropna().to_numpy(); stats.loc[stats.similarity_bin==col,'MAE_CI_low']=np.percentile(v,2.5) if len(v) else np.nan;stats.loc[stats.similarity_bin==col,'MAE_CI_high']=np.percentile(v,97.5) if len(v) else np.nan
stats.to_csv(out/'similarity_bin_validation_errors.csv',index=False)
# threshold comparisons and uncertainties by bootstrapping repeats
thrs=np.round(np.arange(.2,.81,.05),2);rows=[]
for t in thrs:
 z=d[d[s]>=t];vals=[]
 for j in range(500):
  rr=rng.choice(reps,len(reps),replace=True);bb=pd.concat([rb[v] for v in rr],ignore_index=True);zz=bb[bb[s]>=t]
  vals.append((len(zz)/len(bb), zz[e].mean()))
 v=np.array(vals)
 rows.append(dict(threshold=t,n_predictions=len(z),coverage=len(z)/len(d),MAE=z[e].mean(),RMSE=np.sqrt(z.sq_error.mean()),MAE_CI_low=np.percentile(v[:,1],2.5),MAE_CI_high=np.percentile(v[:,1],97.5),coverage_CI_low=np.percentile(v[:,0],2.5),coverage_CI_high=np.percentile(v[:,0],97.5)))
tab=pd.DataFrame(rows);tab.to_csv(out/'threshold_reliability_coverage.csv',index=False)
rho,pval=spearmanr(d[s],d[e]);perrepeat=[]
for rr,gg in d.groupby('repeat'):
 rrho,rp=spearmanr(gg[s],gg[e]);perrepeat.append((rr, len(gg),rrho,rp))
pr=pd.DataFrame(perrepeat,columns=['repeat','n_predictions','spearman_rho','p_value']);pr.to_csv(out/'per_repeat_similarity_error_spearman.csv',index=False)
# colors bright and matched approved chart theme
plt.rcParams.update({'font.family':'DejaVu Sans','axes.labelsize':15,'axes.labelweight':'bold','xtick.labelsize':11.5,'ytick.labelsize':11.5,'axes.linewidth':2,'pdf.fonttype':42,'savefig.dpi':600})
def clean(ax):
 ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False)
 ax.spines['bottom'].set_linewidth(2);ax.spines['left'].set_linewidth(2)
 ax.tick_params(axis='both',width=2,length=5)
 for lbl in list(ax.get_xticklabels())+list(ax.get_yticklabels()):lbl.set_weight('bold')
 ax.grid(False)
fig,axes=plt.subplots(1,2,figsize=(14.5,5.9),constrained_layout=True)
a,c=axes
color='#1686CE'; orange='#EE7F20';green='#159C89'
xx=np.arange(len(stats)); a.bar(xx,stats.MAE,color=[orange if str(t)=='0.30–0.40' else ('#57B3E1' if i%2==0 else color) for i,t in enumerate(stats.similarity_bin)],width=.76,edgecolor='#17293B',linewidth=.6)
a.errorbar(xx,stats.MAE,yerr=[stats.MAE-stats.MAE_CI_low,stats.MAE_CI_high-stats.MAE],fmt='none',capsize=3,color='#17293B',lw=1.4)
a.set_xticks(xx,stats.similarity_bin.astype(str),rotation=44,ha='right');a.set_ylabel('Validation-derived MAE (pIC$_{50}$)');a.set_xlabel('Nearest-training ECFP4 Tanimoto bin');a.set_ylim(0,max(stats.MAE_CI_high)*1.18)
for j,r in stats.iterrows():a.text(j, r.MAE_CI_high+.025,f'n={int(r.n_predictions)}',ha='center',va='bottom',fontsize=9,fontweight='bold')
clean(a)
c.plot(tab.threshold,tab.MAE,color=orange,lw=3,marker='o',markersize=6,label='MAE')
c.fill_between(tab.threshold,tab.MAE_CI_low,tab.MAE_CI_high,color=orange,alpha=.12)
c.set_xlabel('Minimum ECFP4 Tanimoto threshold');c.set_ylabel('Validation-derived MAE (pIC$_{50}$)');c.axvline(.3,color='#4D55B5',lw=2,ls='--',label='Cutoff = 0.30');c.set_xlim(.19,.81);c.set_ylim(.30,.65)
c2=c.twinx();c2.plot(tab.threshold,tab.coverage*100,color=green,lw=3,marker='s',markersize=5,label='Coverage');c2.set_ylabel('Validation coverage (%)',size=15,weight='bold');c2.set_ylim(0,110);c2.spines['top'].set_visible(False);c2.spines['right'].set_linewidth(2);c2.tick_params(width=2,length=5)
for label in c2.get_yticklabels():label.set_weight('bold')
clean(c)
c.legend([c.lines[0],c2.lines[0],c.lines[1]],['MAE','Coverage','AD = 0.30'],frameon=False,fontsize=10,loc='lower left')
for ext in ['png','pdf','svg']:fig.savefig(out/f'BRD4_BD2_Validation_Derived_Errors.{ext}',bbox_inches='tight',facecolor='white',dpi=600)
plt.close(fig)
# plot scatter with binned points in same style
fig,ax=plt.subplots(figsize=(7.7,5.5),layout='constrained')
ax.scatter(d[s],d[e],s=16,alpha=.18,color=color,edgecolors='none',rasterized=True)
ax.errorbar(stats.mean_similarity,stats.MAE,yerr=[stats.MAE-stats.MAE_CI_low,stats.MAE_CI_high-stats.MAE],color=orange,marker='o',lw=2.7,capsize=3,label='Binned MAE ± bootstrap interval')
ax.axvline(.3,color='#4D55B5',ls='--',lw=2,label='AD threshold 0.30')
ax.set_xlabel('Nearest-training ECFP4 Tanimoto similarity');ax.set_ylabel('Absolute validation error (pIC$_{50}$)');ax.set_xlim(0,1.01);ax.set_ylim(0,min(4,d[e].max()+.2));clean(ax);ax.legend(frameon=False,fontsize=10)
for ext in ['png','pdf']:fig.savefig(out/f'BRD4_BD2_Validation_Error_Scatter.{ext}',bbox_inches='tight',facecolor='white',dpi=600)
plt.close(fig)
readme=f'''BRD4-BD2 validation-derived error audit\nSource: {src.name}\n20 repeated scaffold validation partitions, 2,520 validation predictions from 587 unique compounds.\nAbsolute error = abs(predicted_pIC50 - observed_pIC50). Nearest-training ECFP4 Tanimoto similarity computed against the fold-specific training set.\nOverall Spearman correlation (exploratory; repeated predictions are not independent): {rho:.4f}, nominal p={pval:.4g}.\nRepeat-wise Spearman median {pr.spearman_rho.median():.4f}; range {pr.spearman_rho.min():.4f} to {pr.spearman_rho.max():.4f}.\nInterval bars are percentile intervals from resampling 20 repeated splits. Repeats share molecules/scaffolds: these intervals do not fully account for dependence and are descriptive, not independent-sample uncertainty estimates.\nModel parameters/fingerprint were not retuned here; this script audits existing out-of-fold predictions. No FDA drugs are used to determine performance.\nNo claim that threshold 0.30 is a statistically optimized reliability boundary.\n'''
(out/'README.txt').write_text(readme)
with ZipFile('/mnt/data/BRD4_BD2_VALIDATION_ERRORS.zip','w') as z:
 for f in out.iterdir():z.write(f,arcname=f.name)
print(readme)
print(stats.to_string(index=False,float_format=lambda x:f'{x:.3f}'))
print('thresholds:\n',tab[tab.threshold.isin([.2,.3,.4,.5,.6,.7])].to_string(index=False,float_format=lambda x:f'{x:.3f}'))