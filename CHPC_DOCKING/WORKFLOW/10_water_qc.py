#!/usr/bin/env python3
from pathlib import Path
import os
os.chdir(Path.home() / "BRD4_NP_DOCKING/08_PROTONATION")
def records(path):
    return [l for l in Path(path).read_text().splitlines() if l.startswith(("ATOM  ","HETATM"))]
for d in ("BD1","BD2"):
    dry = records(f"PH75/PDBQT/{d}_pH75_dry.pdbqt")
    wet = records(f"PH75/PDBQT/{d}_pH75_wet.pdbqt")
    old = records(f"../02_RECEPTORS/{d}_wet_matched_receptor.pdbqt")
    waters = [l for l in wet if l[17:20].strip() in ("HOH","WAT")]
    original = [l for l in old if l[17:20].strip() in ("HOH","WAT")]
    protein = [l for l in wet if l[17:20].strip() not in ("HOH","WAT")]
    assert dry == protein
    assert len(waters)==len(original)==5
    assert sorted(l[11:] for l in waters)==sorted(l[11:] for l in original)
    assert len({int(l[6:11]) for l in wet})==len(wet), "Duplicate atom serials"
    print(f"\n===== {d} FINAL QC =====")
    print("Protein records unchanged: PASS")
    print("Water coordinates unchanged: PASS")
    print("Water charges and atom types unchanged: PASS")
    print("Water count:",len(waters))
    print("FINAL WATER RESTORATION: PASS")
