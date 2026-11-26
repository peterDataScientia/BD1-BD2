from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, MultipleLocator
import csv

OUT = Path(__file__).resolve().parent
plt.rcParams.update({
    'font.family':'DejaVu Sans', 'font.size':12, 'font.weight':'bold',
    'axes.labelsize':15, 'axes.labelweight':'bold','axes.linewidth':2,
    'xtick.labelsize':12,'ytick.labelsize':12,'xtick.major.size':6,'ytick.major.size':6,
    'xtick.major.width':2,'ytick.major.width':2,'xtick.direction':'out','ytick.direction':'out',
    'legend.fontsize':11,'legend.frameon':False,'lines.linewidth':2.2,
    'figure.facecolor':'white','axes.facecolor':'white','savefig.facecolor':'white',
    'axes.grid':False,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none',
    'mathtext.default':'regular',
})
SIZE=(130/25.4,100/25.4)

def finish(fig, ax, stem):
    ax.minorticks_off()
    ax.tick_params(axis='both', pad=5)
    for t in ax.get_xticklabels()+ax.get_yticklabels(): t.set_fontweight('bold')
    for ext in ['png','pdf','svg']:
        fig.savefig(OUT/f'{stem}.{ext}', dpi=600)
    plt.close(fig)

pH=np.array([7.0,7.4,7.5])
residues=[('BD1 His77',5.77,'#0072E5','o'),('BD2 His369',6.45,'#F07C00','s'),('BD2 His388',6.40,'#009E73','^'),('BD2 His396',6.24,'#B43DBA','D'),('BD2 His437',6.27,'#C63B30','v')]
fig,ax=plt.subplots(figsize=SIZE)
fig.subplots_adjust(left=.20,right=.965,bottom=.20,top=.745)
for name,pka,color,marker in residues:
    y=100/(1+10**(pH-pka))
    ax.plot(pH,y,label=name,color=color,marker=marker,markersize=5.2,markeredgewidth=.8)
ax.set(xlabel='pH',ylabel='Protonated fraction [%]',xlim=(6.96,7.54),ylim=(0,25))
ax.set_xticks([7.0,7.1,7.2,7.3,7.4,7.5]); ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.set_yticks([0,5,10,15,20,25])
ax.legend(ncol=2,loc='lower center',bbox_to_anchor=(.5,1.055),columnspacing=1.1,handlelength=1.7,handletextpad=.45,labelspacing=.4,borderaxespad=0,prop={'size':11,'weight':'bold'})
finish(fig,ax,'Fig1_Histidine_protonation')

pH_full=np.arange(15)
bd1=[7.47,7.44,7.18,5.93,3.83,3.43,3.31,2.90,2.82,3.61,7.09,14.81,23.33,28.95,32.03]
bd2=[.57,.54,.25,-1.56,-5.71,-8.46,-9.25,-9.79,-9.93,-8.71,-4.23,3.37,9.60,13.74,16.69]
fig,ax=plt.subplots(figsize=SIZE)
fig.subplots_adjust(left=.22,right=.965,bottom=.20,top=.92)
ax.plot(pH_full,bd1,marker='o',markersize=4.7,color='#0072E5',label='BD1')
ax.plot(pH_full,bd2,marker='s',markersize=4.7,color='#F07C00',label='BD2')
ax.set(xlabel='pH',ylabel='Relative folding free energy\n[kcal/mol]',xlim=(-.35,14.35),ylim=(-12,35))
ax.set_xticks(np.arange(0,15,2)); ax.set_yticks([-10,0,10,20,30])
ax.legend(loc='upper left',prop={'size':11,'weight':'bold'},handlelength=1.8)
finish(fig,ax,'Fig2_Folding_free_energy')

with open(OUT/'histidine_protonation.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['Residue','Predicted_pKa','pH','Estimated_protonated_fraction_percent'])
    for name,pka,_,_ in residues:
        for ph in pH: w.writerow([name,pka,ph,100/(1+10**(ph-pka))])
with open(OUT/'folding_free_energy.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['pH','BD1_kcal_mol','BD2_kcal_mol']); w.writerows(zip(pH_full,bd1,bd2))
