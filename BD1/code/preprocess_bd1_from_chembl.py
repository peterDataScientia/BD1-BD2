#!/usr/bin/env python3
"""Conservative BRD4-BD1 biochemical IC50 curation from ChEMBL exact-IC50 export.
Source: BRD4_IC50_curated_preliminary.csv, 10,919 prefiltered ChEMBL target records.
Output: auditable assay decisions, eligibility records, standardized compound aggregates.
Important: automated screening; construct/endpoint assignment requires manual assay review.
"""
import argparse,json,re
from pathlib import Path
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

BD1 = re.compile(r'(?<!\w)BD[\s\-_]?1(?!\w)|\bBDI\b|\bfirst bromodomain\b|\bbromodomain\s*(?:domain\s*)?1\b|\bBRD4\s*\(\s*1\s*\)',re.I)
BD2 = re.compile(r'(?<!\w)BD[\s\-_]?2(?!\w)|\bBDII\b|\bsecond bromodomain\b|\bbromodomain\s*(?:domain\s*)?2\b|\bBRD4\s*\(\s*2\s*\)',re.I)
DSF=re.compile(r'\bDSF\b|thermal shift|melting temperature|protein stability shift|differential scanning fluorimetry',re.I)
CELL=re.compile(r'\bcell(?:ular)? (?:viability|proliferation|growth|death|inhibition|IC50|assay|line)s?\b|cytotoxicity|anti[- ]?proliferat|tumou?r cell|\bMTT\b|\bCCK[- ]?8\b',re.I)
MUT=re.compile(r'\bmutant\b|\bmutation\b|\bY390A\b|\bY97A\b|\bvariant\b',re.I)

def decision(row):
    s=str(row['Assay Description'] or '')
    if not BD1.search(s):return 'review_BD1_not_explicit'
    if row['Assay Type']!='B':return 'review_non_biochemical'
    if DSF.search(s):return 'exclude_DSF_thermal_shift'
    if CELL.search(s) or str(row['BAO Label']).lower() in ('cell-based format','organism-based format','tissue-based format'):return 'review_cellular_or_nonbiochemical_format'
    if MUT.search(s):return 'review_mutant_or_multi_construct'
    if BD2.search(s):return 'review_BD1_BD2_mixed_description'
    return 'accept_BD1_biochemical'

def standardize(s):
    try:
        m=Chem.MolFromSmiles(str(s));
        if m is None:return None
        m=rdMolStandardize.Cleanup(m)
        m=rdMolStandardize.FragmentParent(m)
        return Chem.MolToSmiles(m,isomericSmiles=True)
    except Exception:return None

def join_unique(x):return ';'.join(sorted(set(str(v) for v in x if pd.notna(v))))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',required=True);ap.add_argument('--out',required=True);args=ap.parse_args()
    out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    raw=pd.read_csv(args.input,low_memory=False)
    assays=raw.sort_values('Assay ChEMBL ID').drop_duplicates('Assay ChEMBL ID').copy()
    assays['decision']=assays.apply(decision,axis=1)
    assays_audit=assays[['Assay ChEMBL ID','Assay Description','Assay Type','BAO Label','Document ChEMBL ID','decision']].copy()
    assays_audit['n_records_in_input']=assays_audit['Assay ChEMBL ID'].map(raw['Assay ChEMBL ID'].value_counts())
    assays_audit.to_csv(out/'BD1_assay_audit.csv',index=False)
    raw['assay_decision']=raw['Assay ChEMBL ID'].map(dict(zip(assays['Assay ChEMBL ID'],assays['decision'])))
    raw['standard_value_nM']=pd.to_numeric(raw['Standard Value'],errors='coerce')
    quality=(raw['Standard Type'].eq('IC50') & raw['Standard Units'].eq('nM') & raw['Standard Relation'].astype(str).str.strip().isin(['=','\'=\'']) & raw['standard_value_nM'].gt(0) & raw['Data Validity Comment'].isna() & raw['Potential Duplicate'].isin([0,False,np.nan]))
    raw['eligibility']=np.where(raw['assay_decision']=='accept_BD1_biochemical',np.where(quality,'eligible','exclude_measurement_quality'),raw['assay_decision'])
    raw[['Molecule ChEMBL ID','Smiles','Assay ChEMBL ID','Document ChEMBL ID','Assay Description','Assay Type','BAO Label','Standard Value','Standard Units','Standard Relation','Data Validity Comment','Potential Duplicate','eligibility']].to_csv(out/'BD1_record_decisions_all.csv',index=False)
    eligible=raw.loc[raw['eligibility']=='eligible'].copy()
    eligible['standardized_smiles']=eligible['Smiles'].map(standardize)
    eligible['pIC50']=9-np.log10(eligible['standard_value_nM'])
    eligible[['Molecule ChEMBL ID','Smiles','standardized_smiles','Assay ChEMBL ID','Document ChEMBL ID','standard_value_nM','pIC50','Assay Description']].to_csv(out/'BD1_eligible_activity_records.csv',index=False)
    good=eligible.loc[eligible['standardized_smiles'].notna()].copy()
    per_assay=good.groupby(['standardized_smiles','Assay ChEMBL ID'],as_index=False).agg(assay_median_pIC50=('pIC50','median'),n_measurements=('pIC50','size'))
    agg=per_assay.groupby('standardized_smiles',as_index=False).agg(pIC50=('assay_median_pIC50','median'),n_assays=('Assay ChEMBL ID','nunique'),min_assay_pIC50=('assay_median_pIC50','min'),max_assay_pIC50=('assay_median_pIC50','max'))
    agg['interassay_pIC50_range']=agg['max_assay_pIC50']-agg['min_assay_pIC50']
    agg['IC50_median_nM']=10**(9-agg['pIC50'])
    metadata=good.groupby('standardized_smiles',as_index=False).agg(molecule_chembl_ids=('Molecule ChEMBL ID',join_unique),assay_chembl_ids=('Assay ChEMBL ID',join_unique),document_chembl_ids=('Document ChEMBL ID',join_unique),n_raw_records=('pIC50','size'))
    agg=agg.merge(metadata,on='standardized_smiles',how='left')
    agg['decision']=np.where(agg['interassay_pIC50_range']>1,'review_assay_disagreement_gt_1log','model_ready_automated')
    agg.to_csv(out/'BD1_compound_aggregation_audit.csv',index=False)
    final=agg[agg.decision=='model_ready_automated'].copy()
    final.to_csv(out/'BRD4_BD1_QSAR_model_ready.csv',index=False)
    summary={'source_file':Path(args.input).name,'source_note':'Preliminary curated ChEMBL export, not original unfiltered 58,225 activity ZIP','input_exact_IC50_records':len(raw),'input_assays':len(assays),'input_unique_compound_ids':raw['Molecule ChEMBL ID'].nunique(),'assay_decisions':assays_audit['decision'].value_counts().to_dict(),'record_decisions':raw['eligibility'].value_counts().to_dict(),'eligible_records':len(eligible),'eligible_valid_structures':len(good),'unique_standardized_compounds':len(agg),'compounds_disagreement_review':int((agg['decision']!='model_ready_automated').sum()),'model_ready_compounds':len(final),'domain_review':'Automated conservative filtering; assay descriptions should be manually verified for actual BD1 construct prior to publication.','aggregation':'Within-assay median pIC50, then median of assay medians; >1 pIC50 unit across assays held for review'}
    (out/'summary.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))
if __name__=='__main__':main()