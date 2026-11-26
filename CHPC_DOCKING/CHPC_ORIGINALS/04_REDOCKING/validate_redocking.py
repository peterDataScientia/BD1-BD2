import numpy as np
import csv
from rdkit import Chem
from scipy.optimize import linear_sum_assignment

def read_atoms(lines):
    atoms = []
    for l in lines:
        if l.startswith(("ATOM  ", "HETATM")):
            t = l.split()[-1]
            if t in ("H", "HD", "HS"):
                continue

            atoms.append({
                "key": (l[6:11], l[12:16], t),
                "element": {
                    "A":"C", "NA":"N", "OA":"O",
                    "SA":"S"
                }.get(t,t),
                "xyz": [float(l[i:i+8])
                        for i in (30,38,46)]
            })
    return atoms

results = []

for d in ["BD1", "BD2"]:

    ref = Chem.SDMolSupplier(
        d+"_WSH_prepared.sdf",
        removeHs=True
    )[0]

    assert ref is not None
    assert ref.GetNumAtoms() == 32

    ref_xyz = ref.GetConformer().GetPositions()
    ref_elements = [a.GetSymbol() for a in ref.GetAtoms()]

    inp = read_atoms(open(d+"_WSH.pdbqt"))
    inp_xyz = np.array([a["xyz"] for a in inp])

    assert len(inp) == 32

    D = np.linalg.norm(
        ref_xyz[:,None,:]-inp_xyz[None,:,:],
        axis=2
    )

    for i,e in enumerate(ref_elements):
        for j,a in enumerate(inp):
            if e != a["element"]:
                D[i,j] = 1e6

    a,b = linear_sum_assignment(D)

    assert D[a,b].max() < 0.02

    automorphisms = ref.GetSubstructMatches(
        ref,
        uniquify=False,
        useChirality=True,
        maxMatches=100000
    )

    assert 0 < len(automorphisms) < 100000

    for state in ["dry","wet"]:

        print("\n===== %s_%s =====" % (d,state))
        print("Mode | Vina score | Crystal RMSD | Status")

        path = d+"_"+state+"_WSH_redocked.pdbqt"

        lines = open(path).readlines()

        models = []
        current = None

        for l in lines:
            if l.startswith("MODEL"):
                current = []
            elif l.startswith("ENDMDL"):
                assert current is not None
                models.append(current)
                current = None
            elif current is not None:
                current.append(l)

        assert models

        for mode,model in enumerate(models,1):

            pose = read_atoms(model)

            assert len(pose) == 32
            assert [x["key"] for x in pose] == \
                   [x["key"] for x in inp]

            pose_xyz = np.array([x["xyz"] for x in pose])

            mapped = np.empty_like(ref_xyz)
            mapped[a] = pose_xyz[b]

            rmsd = min(
                np.sqrt(np.mean(np.sum(
                    (ref_xyz-mapped[list(p)])**2,
                    axis=1
                )))
                for p in automorphisms
            )

            score = next(
                float(l.split()[3])
                for l in model
                if l.startswith("REMARK VINA RESULT:")
            )

            status = "PASS" if rmsd < 2.0 else "FAIL"

            print("%4d | %10.3f | %12.3f | %s" %
                  (mode,score,rmsd,status))

            results.append([
                d,state,mode,score,rmsd,status
            ])

with open("redocking_validation.csv","w") as f:
    w = csv.writer(f)
    w.writerow([
        "Domain","Receptor","Mode",
        "Vina_score_kcal_mol",
        "Crystal_RMSD_A","Status"
    ])
    w.writerows(results)

print("\n===== TOP-RANKED POSE VALIDATION =====")

for d,s,m,score,rmsd,status in results:
    if m == 1:
        print("%s_%s | %.3f | %.3f | %s" %
              (d,s,score,rmsd,status))

print("\nSaved: redocking_validation.csv")
