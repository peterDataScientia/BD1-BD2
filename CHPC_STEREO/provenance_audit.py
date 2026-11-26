"""Trace consensus candidates to COCONUT source metadata; never infer confirmed origin."""
import argparse, hashlib, json, os, re, zipfile
from pathlib import Path
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import inchi

BASE=Path("/mnt/lustre/users/apeter1/BRD4_NP_SCREENING/STEREO_SECONDARY/060c6dff792a34ad")
ZIP=Path("/mnt/lustre/users/apeter1/CDK_NP_SCREENING/coconut_csv-11-2026.zip")
CATEGORIES={"BD1_selective":("consensus_BD1_selective","BD1_conservative_pIC50","BD1_conservative_gap"),
            "BD2_selective":("consensus_BD2_selective","BD2_conservative_pIC50","BD2_conservative_gap"),
            "dual":("consensus_dual","dual_conservative_pIC50",None)}
def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(4*1024*1024),b""): h.update(b)
    return h.hexdigest()
def canonical(s):
    m=Chem.MolFromSmiles(s)
    if m is None: return "", ""
    return Chem.MolToSmiles(m,isomericSmiles=True), inchi.MolToInchiKey(m)
def meaningful(s):
    return str(s).strip().lower() not in ("","[]","{}","nan","none","null","n/a","na","unknown","not available","-")
def valid_doi(s):
    return bool(re.search(r"10\.\d{4,9}/\S+",str(s),re.I))
def write(p,d):
    tmp=p.with_name(p.name+".tmp")
    d.to_csv(tmp,index=False)
    os.replace(tmp,p)
def main():
    a=argparse.ArgumentParser()
    a.add_argument("--base",type=Path,default=BASE)
    a.add_argument("--coconut",type=Path,default=ZIP)
    args=a.parse_args()
    root=args.base
    u=pd.read_csv(root/"unique_scoring_structures.csv.gz",dtype={"identifier":str},keep_default_na=False)
    aliases=pd.read_csv(root/"source_aliases.csv",dtype={"identifier":str},keep_default_na=False)
    assert u.chemical_key.is_unique
    clean=u[(u.joint_min_AD>=.3)&~u.known_connectivity_either&u.manual_exclusion_reason.eq("")].copy()
    for name,(flag,_,__) in CATEGORIES.items():
        old=pd.read_csv(root/(name+"_AD_0.3.csv"),dtype={"identifier":str})
        assert set(old.chemical_key)==set(clean.loc[clean[flag],"chemical_key"]),name
    d=clean.loc[clean[[v[0] for v in CATEGORIES.values()]].any(axis=1)].copy()
    assert len(d)>0
    al=aliases[aliases.chemical_key.isin(d.chemical_key)].copy()
    assert set(al.chemical_key)==set(d.chemical_key)
    print("Unique consensus candidates:",len(d),"source aliases:",len(al),flush=True)
    al[["raw_canonical_SMILES","raw_InChIKey"]]=pd.DataFrame([canonical(s) for s in al.smiles],index=al.index)
    assert al.raw_canonical_SMILES.ne("").all()
    assert al.raw_InChIKey.ne("").any(),"RDKit InChI support required"
    ids=set(al.loc[al.source.eq("COCONUT"),"identifier"])
    keys=set(al.raw_InChIKey)-{""}
    matches=[]
    with zipfile.ZipFile(args.coconut) as z:
        names=[x for x in z.namelist() if x.lower().endswith(".csv")]
        assert len(names)==1,names
        with z.open(names[0]) as f:
            for chunk in pd.read_csv(f,chunksize=20000,dtype=str,keep_default_na=False):
                mask=chunk.identifier.isin(ids)|chunk.standard_inchi_key.isin(keys)
                if mask.any(): matches.append(chunk[mask].copy())
    if not matches: raise RuntimeError("No candidate metadata matches")
    m=pd.concat(matches,ignore_index=True)
    m.insert(0,"metadata_row",np.arange(len(m)))
    m[["metadata_canonical_SMILES","recomputed_InChIKey"]]=pd.DataFrame([canonical(s) for s in m.canonical_smiles],index=m.index)
    # Join via source ID and independent exact full InChIKey. Confirm exact
    # canonical raw structures; ID or InChIKey alone never qualifies.
    left=al[al.source.eq("COCONUT")].merge(m,left_on="identifier",right_on="identifier",how="inner",suffixes=("_alias","_meta"))
    left["match_method"]="COCONUT_ID"
    left=left.rename(columns={"identifier":"alias_identifier"})
    left["coconut_identifier"]=left.alias_identifier
    right=al[al.raw_InChIKey.ne("")].merge(m,left_on="raw_InChIKey",right_on="standard_inchi_key",how="inner",suffixes=("_alias","_meta"))
    right=right.rename(columns={"identifier_alias":"alias_identifier","identifier_meta":"coconut_identifier"})
    right["match_method"]="FULL_INCHIKEY"
    joined=pd.concat([left,right],ignore_index=True,sort=False)
    assert len(joined)>0
    joined["exact_structure_match"]=joined.raw_canonical_SMILES.eq(joined.metadata_canonical_SMILES)&joined.metadata_canonical_SMILES.ne("")
    joined["inchi_key_consistent"]=joined.raw_InChIKey.eq(joined.recomputed_InChIKey)&joined.raw_InChIKey.eq(joined.standard_inchi_key)&joined.raw_InChIKey.ne("")
    joined["identity_verified_in_snapshot"]=joined.exact_structure_match&joined.inchi_key_consistent
    joined["organism_present"]=joined.organisms.map(meaningful)
    joined["doi_present"]=joined.dois.map(valid_doi)
    joined["organism_and_doi_same_record"]=joined.identity_verified_in_snapshot&joined.organism_present&joined.doi_present
    # These are metadata flags, not verified isolation evidence.
    def join_values(s):
        return " || ".join(sorted(set(str(x) for x in s if meaningful(x))))
    verified=joined[joined.identity_verified_in_snapshot].copy()
    evidence=verified[verified.organism_and_doi_same_record].copy()
    records=[]
    for key in d.chemical_key:
        j=joined[joined.chemical_key.eq(key)]
        v=verified[verified.chemical_key.eq(key)]
        e=evidence[evidence.chemical_key.eq(key)]
        if len(e): tier="ORGANISM_AND_DOI"
        elif len(v) and v.organism_present.any(): tier="ORGANISM_ONLY"
        elif len(v) and v.doi_present.any(): tier="DOI_ONLY"
        elif len(v): tier="NO_ORIGIN_METADATA"
        elif len(j): tier="STRUCTURE_MISMATCH"
        else: tier="NO_COCONUT_MATCH"
        records.append(dict(chemical_key=key,metadata_tier=tier,
             origin_status="UNVERIFIED",
             coconut_ids=join_values(v.coconut_identifier),
             names=join_values(v.name),organisms=join_values(v.organisms),
             dois=join_values(v.dois),collections=join_values(v.collections),
             evidence_coconut_ids=join_values(e.coconut_identifier),
             evidence_organisms=join_values(e.organisms),evidence_dois=join_values(e.dois),
             matched_metadata_rows=int(v.metadata_row.nunique()),
             unmatched_identity_rows=int(j.loc[~j.identity_verified_in_snapshot,"metadata_row"].nunique())))
    d=d.merge(pd.DataFrame(records),on="chemical_key",validate="one_to_one")
    tier_order={"ORGANISM_AND_DOI":0,"ORGANISM_ONLY":1,"DOI_ONLY":2,"NO_ORIGIN_METADATA":3,"STRUCTURE_MISMATCH":4,"NO_COCONUT_MATCH":5}
    d["metadata_rank"]=d.metadata_tier.map(tier_order)
    # Store complete, nonaggregated evidence to retain organism-reference associations.
    out=root/"PROVENANCE"
    out.mkdir(exist_ok=True)
    write(out/"all_consensus_candidates.csv",d.sort_values(["metadata_rank","joint_min_AD","max_representation_disagreement","chemical_key"],ascending=[True,False,True,True]))
    write(out/"candidate_aliases.csv",al)
    write(out/"coconut_metadata_matches.csv",joined)
    write(out/"isolation_review_queue.csv",evidence.drop_duplicates(["chemical_key","metadata_row"]))
    write(out/"identity_review.csv",joined.loc[~joined.identity_verified_in_snapshot])
    # Empty fields for human verification, explicitly not automatically approved.
    manual=d.loc[d.metadata_tier.eq("ORGANISM_AND_DOI"),["chemical_key","identifier","evidence_coconut_ids","evidence_organisms","evidence_dois","BD1_SMILES","BD2_SMILES"]].copy()
    for col in ("verified_name","verified_organism","isolation_reference","structure_verification","origin_decision","evidence_location","known_BRD4_activity","review_notes"):
        manual[col]=""
    write(out/"manual_review_template.csv",manual)
    report={"unique_candidates":len(d),"source_aliases":len(al),"metadata_tiers":{k:int(v) for k,v in d.metadata_tier.value_counts().items()},
            "confirmed_natural_products":None,"note":"Organism and DOI are metadata leads, not confirmation. Exact raw structure plus full InChIKey required. LOTUS-only structures without COCONUT match require independent LOTUS provenance. No AD or potency thresholds changed.","categories":{}}
    for name,(flag,pot,gap) in CATEGORIES.items():
        x=d[d[flag]].copy()
        # Prioritize provenance evidence, then stereo completeness, joint AD, conservative potency and selectivity.
        sort=["metadata_rank","unspecified_stereo_either","joint_min_AD",pot]+([gap] if gap else [])+["max_representation_disagreement","chemical_key"]
        asc=[True,True,False,False]+([False] if gap else [])+[True,True]
        x=x.sort_values(sort,ascending=asc)
        x.insert(0,"provenance_rank",np.arange(1,len(x)+1))
        write(out/(name+"_all_ranked.csv"),x)
        y=x[x.metadata_tier.eq("ORGANISM_AND_DOI")]
        write(out/(name+"_organism_and_doi.csv"),y)
        report["categories"][name]={"total":len(x),"organism_and_doi":len(y),
           "organism_and_doi_AD_0.4":int((y.joint_min_AD>=.4).sum()),
           "organism_and_doi_AD_0.5":int((y.joint_min_AD>=.5).sum())}
    report["input_sha256"]={"coconut_zip":sha(args.coconut),"unique_scoring_structures":sha(root/"unique_scoring_structures.csv.gz"),"source_aliases":sha(root/"source_aliases.csv")}
    report["script_sha256"]=sha(Path(__file__))
    temp=out/"summary.tmp"
    temp.write_text(json.dumps(report,indent=2))
    os.replace(temp,out/"summary.json")
    print(json.dumps(report,indent=2),flush=True)
    print("Output:",out,flush=True)
if __name__=="__main__": main()
