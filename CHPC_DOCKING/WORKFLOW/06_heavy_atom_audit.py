#!/usr/bin/env python3
"""Consolidated version of the coordinate and displacement audits."""
import math
from pathlib import Path
import os
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--stage", choices=["pdb", "pdbqt"], default="pdb")
args = parser.parse_args()
os.chdir(Path.home() / "BRD4_NP_DOCKING/08_PROTONATION")

def atoms(path):
    data = {}
    for l in open(path):
        if not l.startswith("ATOM  "):
            continue
        name = l[12:16].strip()
        if name.lstrip("0123456789").startswith("H"):
            continue
        key = (l[21], l[22:27].strip(), name)
        if key in data:
            raise ValueError("Duplicate atom: " + str(key))
        data[key] = (l[17:20].strip(), tuple(float(l[i:i+8]) for i in (30,38,46)))
    return data

for d in ("BD1","BD2"):
    if args.stage == "pdb":
        a, b = atoms(f"{d}_apo.pdb"), atoms(f"PH75/{d}_pH75.pdb")
    else:
        a, b = atoms(f"PH75/{d}_pH75.pdb"), atoms(f"PH75/PDBQT/{d}_pH75_dry.pdbqt")
    ka, kb = set(a), set(b)
    common = sorted(ka & kb)
    assert common, "No matched atoms"
    ds = [math.dist(a[k][1],b[k][1]) for k in common]
    print(f"\n===== {d} =====")
    print("Before heavy atoms:",len(a))
    print("After heavy atoms:",len(b))
    print("Matched:",len(ds))
    print("Removed:",sorted(ka-kb))
    print("Added:",sorted(kb-ka))
    print("RMSD: %.6f A" % math.sqrt(sum(x*x for x in ds)/len(ds)))
    print("Maximum displacement: %.6f A" % max(ds))
    print("Atoms displaced > 0.1 A:",sum(x>0.1 for x in ds))
    for k, dist in zip(common,ds):
        if dist > 0.1:
            print(a[k][0],k,round(dist,6),"A")
            print("Before:",a[k][1])
            print("After:",b[k][1])
    if args.stage == "pdbqt":
        assert ka == kb and max(ds) <= 0.001, "Conversion geometry changed"
    else:
        expected = {("B","459","OXT")} if d == "BD2" else set()
        assert not (ka-kb) and (kb-ka)==expected, "Unexpected heavy atom changes"
