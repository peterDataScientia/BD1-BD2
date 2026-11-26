#!/usr/bin/env python3
from pathlib import Path
import os
os.chdir(Path.home() / "BRD4_NP_DOCKING/08_PROTONATION")
for d in ("BD1","BD2"):
    dry = Path(f"PH75/PDBQT/{d}_pH75_dry.pdbqt")
    old = Path(f"../02_RECEPTORS/{d}_wet_matched_receptor.pdbqt")
    out = Path(f"PH75/PDBQT/{d}_pH75_wet.pdbqt")
    protein = [l for l in dry.read_text().splitlines() if l.startswith("ATOM  ")]
    waters = [l for l in old.read_text().splitlines()
              if l.startswith(("ATOM  ","HETATM")) and l[17:20].strip() in ("HOH","WAT")]
    assert len(waters)==5, (d,len(waters))
    serial = max(int(l[6:11]) for l in protein)
    restored = [l[:6]+f"{serial+i:5d}"+l[11:] for i,l in enumerate(waters,1)]
    out.write_text("\n".join(protein+restored)+"\n")
    print(f"\n===== {d} =====")
    print("Protein atoms:",len(protein))
    print("Restored waters:",len(restored))
    print("Total atoms:",len(protein)+len(restored))
    print("Output:",out)
