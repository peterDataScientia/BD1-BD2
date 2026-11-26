"""Prepare the 15 manifest structures with RDKit ETKDGv3 and MMFF94.
This step preserves the supplied molecular identity; it is not a pH
microstate prediction or final PDBQT preparation.
Run from the directory containing the manifest.
"""
import csv
from pathlib import Path
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

rows = list(csv.DictReader(open("BRD4_NP15_ligand_manifest.csv", encoding="utf-8")))
assert len(rows) == 15
assert len({r["inchikey"] for r in rows}) == 15
for folder in ("SDF", "MOL2", "PDBQT"):
    Path(folder).mkdir(exist_ok=True)

results = []
for r in rows:
    lid, name = r["ligand_id"], r["name"]
    mol = Chem.MolFromSmiles(r["smiles"])
    assert mol is not None, lid
    assert Chem.MolToInchiKey(mol) == r["inchikey"], lid
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 2024
    params.numThreads = 1
    params.enforceChirality = True
    ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=10, params=params))
    assert ids, lid + ": embedding failed"
    assert AllChem.MMFFHasAllMoleculeParams(mol), lid
    opt = AllChem.MMFFOptimizeMoleculeConfs(
        mol, numThreads=1, maxIters=2000, mmffVariant="MMFF94"
    )
    valid = [(energy, cid) for cid, (status, energy) in zip(ids, opt)
             if status == 0 and np.isfinite(energy)]
    assert valid, lid + ": No converged conformer"
    energy, best = min(valid)
    final = Chem.Mol(mol)
    final.RemoveAllConformers()
    final.AddConformer(mol.GetConformer(best), assignId=True)
    Chem.AssignStereochemistryFrom3D(final)
    key = Chem.MolToInchiKey(Chem.RemoveHs(final))
    assert key == r["inchikey"], lid + ": Identity mismatch"
    final.SetProp("_Name", lid + "_" + name)
    final.SetProp("InChIKey", key)
    final.SetProp("MMFF94_energy", str(energy))
    writer = Chem.SDWriter("SDF/" + lid + ".sdf")
    writer.write(final)
    writer.close()
    # Verify the written intermediate file, not just the in-memory molecule.
    check = Chem.SDMolSupplier("SDF/" + lid + ".sdf", removeHs=True)[0]
    assert check is not None and Chem.MolToInchiKey(check) == key, lid
    results.append([lid, name, len(ids), len(valid), energy,
                    Chem.GetFormalCharge(final), key, "PASS"])
    print("%s | %s | %d/%d converged | %.3f | PASS" %
          (lid, name, len(valid), len(ids), energy), flush=True)

with open("ligand_preparation_audit.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["Ligand_ID","Name","Embedded_conformers","Converged_conformers",
                "MMFF94_energy_kcal_mol","Formal_charge","InChIKey","Status"])
    w.writerows(results)
print("\nPrepared:", len(results), "/ 15")
