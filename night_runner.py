#!/usr/bin/env python3
import os, sys, json, math, subprocess
from pathlib import Path
import numpy as np
import pandas as pd

def sh(cmd):
    print("RUN:", " ".join(cmd))
    subprocess.run(cmd, check=True)

def pick_best(csv_path):
    df = pd.read_csv(csv_path)
    # drop rows with NaNs in metrics we need
    df = df.dropna(subset=["recon_err", "H_nonzero_frac", "stability_mean"])
    if df.empty: raise RuntimeError("No valid rows in grid_summary.csv")

    # normalize metrics (lower is better for recon_err & H_nonzero_frac; higher is better for stability)
    def z(x): return (x - x.mean()) / (x.std() + 1e-12)
    z_re = z(df["recon_err"])
    z_hn = z(df["H_nonzero_frac"])
    stab = df["stability_mean"]
    z_stab = (stab - stab.mean()) / (stab.std() + 1e-12)

    # composite score (lower is better): recon + sparsity - stability
    score = z_re + 0.5*z_hn - 1.5*z_stab
    df["score"] = score

    df = df.sort_values("score")
    best = df.iloc[0].to_dict()
    return best, df

def load_duration(times_path):
    t = np.load(times_path)
    if t.size == 0: raise RuntimeError("Empty times file")
    return float(np.max(t))

def main():
    # ---- user paths (edit here) ----
    PY  = "/opt/anaconda3/bin/python"
    SCRIPT = "/Volumes/Extreme SSD/Neuropixels/Python/seqnmf_like.py"
    CLUS = "/Volumes/Extreme SSD/Neuropixels/9153_01302025_tagging_g0/9153_01302025_tagging_g0_imec0/kilosort4/spike_clusters.npy"
    TIMES= "/Volumes/Extreme SSD/Neuropixels/9153_01302025_tagging_g0/9153_01302025_tagging_g0_imec0/kilosort4/spike_seconds_adj.npy"
    LAB  = "/Volumes/Extreme SSD/Neuropixels/9153_01302025_tagging_g0/9153_01302025_tagging_g0_imec0/kilosort4qMetrics/templates._bc_unit_labels.tsv"
    CLASSCSV = "/Volumes/Extreme SSD/Neuropixels/9153_01302025_tagging_g0/9153_01302025_tagging_g0_imec0/kilosort4/unit_classification_rulebased.csv"
    FR_MIN = 0.2
    FR_MAX = 12

    GRID = "/Volumes/Extreme SSD/Neuropixels/seqnmf_grid"
    FINAL= "/Volumes/Extreme SSD/Neuropixels/seqnmf_final"

    # ---- grid hyperparams (short slice) ----
    tmin = 0.0
    tmax_grid = 600.0           # 10 min grid slice
    binsize = "0.005,0.01,0.02"
    topn    = "200,300"
    L       = "30,60"
    hop     = "5"
    k       = "4,8,16"
    alphaH  = "0.2,0.4,0.8"
    seeds   = "0,1"

    # 1) run grid
    Path(GRID).mkdir(parents=True, exist_ok=True)
    sh([PY, "/Volumes/Extreme SSD/Neuropixels/Python/run_seqnmf_grid.py",
        "--python", PY, "--script", SCRIPT,
        "--clusters", CLUS, "--times", TIMES, "--labels", LAB,
        "--classcsv", CLASSCSV,   # <-- add this
        "--tmin", str(tmin), "--tmax", str(tmax_grid),
        "--binsize", binsize, "--topn", topn, "--L", L, "--hop", hop,
        "--k", k, "--alphaH", alphaH, "--seeds", seeds,
        "--fr_min", str(FR_MIN), "--fr_max", str(FR_MAX),
        "--outbase", GRID, "--resume"])

    # 2) pick best
    csv_path = str(Path(GRID) / "grid_summary.csv")
    best, df = pick_best(csv_path)
    print("\nBEST COMBO (grid slice):")
    print(json.dumps(best, indent=2))

    # save ranking
    df.to_csv(str(Path(GRID) / "grid_summary_ranked.csv"), index=False)

    # 3) run full duration with best hyperparams (both seeds)
    T_full = load_duration(TIMES)
    print(f"\nFull duration detected: {T_full:.3f} s")

    Path(FINAL).mkdir(parents=True, exist_ok=True)
    params = dict(
        binsize=str(best["binsize"]),
        topn=str(int(best["topn"])),
        L=str(int(best["L"])),
        hop=str(int(best["hop"])),
        k=str(int(best["k"])),
        alphaH=str(best["alphaH"]),
    )

    for seed in ["0","1"]:
        outdir = str(Path(FINAL) / f"FULL_bs{params['binsize']}_top{params['topn']}_L{params['L']}_hop{params['hop']}_k{params['k']}_a{params['alphaH']}_seed{seed}")
        sh([PY, SCRIPT, "--clusters", CLUS, "--times", TIMES, "--labels", LAB,
            "--classcsv", CLASSCSV,   # <-- add this
            "--tmin", "0.0", "--tmax", f"{T_full}",
            "--binsize", params["binsize"], "--topn", params["topn"],
            "--L", params["L"], "--hop", params["hop"],
            "--k", params["k"], "--alphaH", params["alphaH"],
            "--fr_min", str(FR_MIN), "--fr_max", str(FR_MAX),
            "--seed", seed, "--outdir", outdir])

    # 4) write final choice
    with open(str(Path(FINAL) / "FINAL_CHOICE.txt"), "w") as f:
        f.write(json.dumps(best, indent=2))
        f.write("\nFull duration seconds: %.3f\n" % T_full)
        f.write("Two seeds were fit and saved under: %s\n" % FINAL)
    print("\nWrote FINAL_CHOICE.txt")
    print("Done.")
if __name__ == "__main__":
    main()