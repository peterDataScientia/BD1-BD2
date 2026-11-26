#!/usr/bin/env python3
import math
from pathlib import Path
import os
os.chdir(Path.home() / "BRD4_NP_DOCKING/08_PROTONATION")
def atoms(path):
    return [(l[12:16].strip(), l[17:20].strip(), l[22:27].strip(),
             tuple(float(l[i:i+8]) for i in (30,38,46)))
            for l in open(path)
            if l.startswith(("ATOM  ","HETATM"))
            and not l[12:16].strip().lstrip("0123456789").startswith("H")]
for d,res in [("BD1","116"),("BD2","353")]:
    lig = atoms(f"../04_REDOCKING/{d}_WSH_crystal.pdb")
    waters = [a for a in atoms(f"../02_RECEPTORS/{d}_wet_matched_receptor.pdbqt")
              if a[1] in ("HOH","WAT")]
    assert len(waters) == 5
    print(f"\n===== {d} =====")
    for state,path in [("Before",f"{d}_apo.pdb"),("After",f"PH75/{d}_pH75.pdb")]:
        residue = [a for a in atoms(path) if a[2] == res]
        print(state)
        print("Minimum ligand distance:",round(min(math.dist(a[3],b[3]) for a in residue for b in lig),3),"A")
        print("Minimum retained water distance:",round(min(math.dist(a[3],b[3]) for a in residue for b in waters),3),"A")
