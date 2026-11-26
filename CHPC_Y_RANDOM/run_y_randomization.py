"""Nested scaffold CV permutation test. Run from any directory."""
import argparse
import hashlib
import json
import os
import time
from pathlib import Path

import numpy as np
import pandas as pd
import rdkit
import sklearn
from joblib import Parallel, delayed
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.model_selection import GroupKFold
from sklearn.svm import SVR

FILES = {
    "BD1": ("BD1/data/BRD4_BD1_QSAR_model_ready_2816.csv",
            "BD1/results/heldout_predictions.csv", "standardized_smiles", 2816, 541),
    "BD2": ("BD2/data/BRD4_BD2_STRICT_MODELING_775.csv",
            "BD2/validation/heldout_predictions_with_AD.csv",
            "canonical_isomeric_smiles", 775, 187),
}
GAMMAS = (0.005, 0.01)
CS = (1, 10, 100)
EPSILON = 0.1
SEED = 20261126
N_PERM = 200

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def atomic_json(path, obj):
    tmp = path.with_name(path.name + "." + str(os.getpid()) + ".tmp")
    tmp.write_text(json.dumps(obj, indent=2, allow_nan=False))
    os.replace(tmp, path)

def atomic_csv(path, df):
    tmp = path.with_name(path.name + "." + str(os.getpid()) + ".tmp")
    df.to_csv(tmp, index=False)
    os.replace(tmp, path)

def scaffold(m, s):
    return MurckoScaffold.MurckoScaffoldSmiles(mol=m, includeChirality=False) or "ACYCLIC_" + s

def prepare(root, domain):
    data, hold, col, n, h = FILES[domain]
    d = pd.read_csv(root / data)
    test = pd.read_csv(root / hold, usecols=[col])
    assert len(d) == n and len(test) == h
    assert d[col].is_unique and test[col].is_unique
    assert set(test[col]) <= set(d[col])
    d = d.loc[~d[col].isin(test[col])].reset_index(drop=True)
    assert len(d) == n-h
    y = d.pIC50.to_numpy(float)
    assert np.isfinite(y).all() and np.var(y) > 0
    mols = [Chem.MolFromSmiles(s) for s in d[col]]
    hold_mols = [Chem.MolFromSmiles(s) for s in test[col]]
    assert all(m is not None for m in mols + hold_mols)
    groups = np.array([scaffold(m, s) for m, s in zip(mols, d[col])])
    assert set(groups).isdisjoint({scaffold(m, s) for m, s in zip(hold_mols, test[col])})
    idx = np.arange(len(y))
    outer = []
    fold_ids = np.full(len(y), -1, dtype=int)
    for f, (tr, va) in enumerate(GroupKFold(5).split(idx, groups=groups)):
        assert set(groups[tr]).isdisjoint(groups[va])
        inner = [(tr[t], tr[v]) for t, v in GroupKFold(4).split(tr, groups=groups[tr])]
        for t, v in inner:
            assert set(groups[t]).isdisjoint(groups[v])
        outer.append((tr, va, inner))
        assert (fold_ids[va] == -1).all()
        fold_ids[va] = f
    assert (fold_ids >= 0).all()
    # Fold identities are fixed using structures only and reused for all permutations.
    meta = {
        "domain": domain, "n_development": len(d), "n_holdout": h,
        "n_permutations": N_PERM, "seed": SEED,
        "fingerprint": "Morgan radius 2, 2048 binary bits",
        "outer": "GroupKFold(5)", "inner": "GroupKFold(4)",
        "scaffold_rule": "Murcko no chirality; acyclic unique by SMILES (legacy LOOCV rule)",
        "grid": {"C": list(CS), "gamma": list(GAMMAS), "epsilon": EPSILON},
        "selection": "minimum mean inner-fold RMSE; ties gamma then C",
        "score": "pooled out-of-fold Q2 = 1-SSE/SST of evaluated labels",
        "inputs": {data: sha(root / data), hold: sha(root / hold)},
        "fold_sha256": hashlib.sha256(fold_ids.astype("<i8").tobytes()).hexdigest(),
        "script_sha256": sha(Path(__file__)),
        "versions": {"numpy": np.__version__, "pandas": pd.__version__,
                     "rdkit": rdkit.__version__, "sklearn": sklearn.__version__},
    }
    key = hashlib.sha256(json.dumps(meta, sort_keys=True).encode()).hexdigest()[:16]
    return d, col, mols, y, groups, fold_ids, outer, meta, key

def fit(K, C, tr, va, y):
    return SVR(kernel="precomputed", C=C, epsilon=EPSILON, cache_size=100).fit(
        K[np.ix_(tr, tr)], y[tr]).predict(K[np.ix_(va, tr)])

def evaluate(y, kernels, outer):
    pred = np.full(len(y), np.nan)
    choices = []
    for f, (tr, va, inner) in enumerate(outer):
        scores = []
        for g in GAMMAS:
            K = kernels[g]
            for C in CS:
                score = np.mean([np.sqrt(np.mean((fit(K, C, t, v, y)-y[v])**2))
                                 for t, v in inner])
                scores.append((float(score), g, C))
        score, g, C = min(scores)
        pred[va] = fit(kernels[g], C, tr, va, y)
        choices.append({"fold": f, "C": C, "gamma": g, "inner_RMSE": score})
    assert np.isfinite(pred).all()
    err = pred-y
    return {"Q2": float(1-np.sum(err**2)/np.sum((y-y.mean())**2)),
            "RMSE": float(np.sqrt(np.mean(err**2))),
            "MAE": float(np.mean(np.abs(err))), "bias": float(np.mean(err)),
            "fold_choices": choices}, pred

def summarize(out, meta):
    ids = range(N_PERM+1)
    missing = [i for i in ids if not (out / ("perm_%03d.json" % i)).is_file()]
    print("Completed %d/%d; missing: %s" % (N_PERM+1-len(missing), N_PERM+1, missing), flush=True)
    if missing:
        raise SystemExit("Incomplete test; no final p-value calculated.")
    rows = []
    for i in ids:
        r = json.loads((out / ("perm_%03d.json" % i)).read_text())
        assert r["permutation"] == i and r["n"] == meta["n_development"]
        assert np.isfinite([r[k] for k in ("Q2","RMSE","MAE","bias")]).all()
        rows.append({k: r[k] for k in ("permutation","n","Q2","RMSE","MAE","bias","seconds")})
    frame = pd.DataFrame(rows)
    original = rows[0]
    null = frame.iloc[1:]
    b = int((null.Q2 >= original["Q2"]).sum())
    summary = {
        "domain": meta["domain"], "n_development": meta["n_development"],
        "n_permutations": N_PERM, "original_Q2": original["Q2"],
        "original_RMSE": original["RMSE"], "original_MAE": original["MAE"],
        "random_Q2_mean": float(null.Q2.mean()), "random_Q2_sd": float(null.Q2.std(ddof=1)),
        "random_Q2_min": float(null.Q2.min()), "random_Q2_max": float(null.Q2.max()),
        "random_Q2_95th_percentile": float(null.Q2.quantile(.95)),
        "random_RMSE_mean": float(null.RMSE.mean()), "b_ge_original": b,
        "p_value_one_sided": (b+1)/(N_PERM+1),
    }
    atomic_csv(out / "Y_randomization_results.csv", frame)
    atomic_json(out / "summary.json", summary)
    print(json.dumps(summary, indent=2), flush=True)

def main():
    p = argparse.ArgumentParser()
    p.add_argument("domain", choices=FILES)
    p.add_argument("--start", type=int, default=0)
    p.add_argument("--end", type=int, default=N_PERM)
    p.add_argument("--workers", type=int, default=1)
    p.add_argument("--summarize", action="store_true")
    a = p.parse_args()
    if not (0 <= a.start <= a.end <= N_PERM and a.workers >= 1):
        p.error("Require 0 <= start <= end <= 200 and workers >= 1.")
    root = Path(__file__).resolve().parents[1]
    d, col, mols, y, groups, folds, outer, meta, key = prepare(root, a.domain)
    out = root / "CHPC_Y_RANDOM" / "outputs" / a.domain / key
    out.mkdir(parents=True, exist_ok=True)
    atomic_json(out / "protocol.json", meta)
    print("Output:", out, flush=True)
    if a.summarize:
        summarize(out, meta)
        return
    pending = [i for i in range(a.start, a.end+1) if not (out / ("perm_%03d.json" % i)).exists()]
    if pending:
        gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
        X = np.asarray([gen.GetFingerprintAsNumPy(m) for m in mols], dtype=np.float64)
        dist = euclidean_distances(X, X, squared=True)
        kernels = {g: np.exp(-g*dist) for g in GAMMAS}
        del X, dist

        def one(i):
            t0 = time.monotonic()
            # Permutation zero is the unshuffled baseline. Target-specific deterministic streams.
            rng = np.random.default_rng(np.random.SeedSequence([SEED, 1 if a.domain == "BD1" else 2, i]))
            yp = y if i == 0 else rng.permutation(y)
            row, pred = evaluate(yp, kernels, outer)
            row.update({"permutation": i, "n": len(y), "seconds": time.monotonic()-t0})
            if i == 0:
                atomic_csv(out / "original_OOF_predictions.csv", pd.DataFrame({
                    "smiles": d[col], "scaffold": groups, "outer_fold": folds,
                    "observed_pIC50": y, "predicted_pIC50": pred}))
            atomic_json(out / ("perm_%03d.json" % i), row)
            print("%s perm %03d Q2=%+.5f RMSE=%.5f seconds=%.1f" %
                  (a.domain, i, row["Q2"], row["RMSE"], row["seconds"]), flush=True)

        # Threads share the read-only kernels instead of replicating them in processes.
        Parallel(n_jobs=min(a.workers, len(pending)), prefer="threads", batch_size=1)(
            delayed(one)(i) for i in pending)
    count = len(list(out.glob("perm_*.json")))
    print("Completed %d/%d" % (count, N_PERM+1), flush=True)
    if count == N_PERM+1:
        summarize(out, meta)

if __name__ == "__main__":
    main()
