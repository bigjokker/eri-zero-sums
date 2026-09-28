"""Compute 2 Re ERi(ρ log x) for zeros with γ ≤ T_max, save a table.

Uses Odlyzko zeros1 (first 100k, 3e-9) and Gram's series (Grobner (1.1)).
"""
from __future__ import annotations

import csv
import sys
import time
from pathlib import Path

from mpmath import mp, mpc, mpf, log, zeta

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from compute.eri import eri_gram  # noqa: E402

ZEROS = ROOT / "data" / "zeros1.txt"
OUT = ROOT / "data" / "eri_list.csv"
X = 1.1
N_ZEROS = 10142  # numerics notes: 10142 zeros to T = 9998.85


def load_gammas(n: int = N_ZEROS) -> list[float]:
    gs = []
    with ZEROS.open() as f:
        for line in f:
            gs.append(float(line))
            if len(gs) >= n:
                break
    return gs


def precompute_zeta(max_k: int, dps: int) -> None:
    prev = mp.dps
    mp.dps = dps
    # Warm the cache inside eri_gram via direct zeta calls it will reuse
    # by filling eri._zeta_cache.
    from compute import eri as eri_mod

    for k in range(2, max_k + 2):
        eri_mod._zeta_cache[(dps, k)] = zeta(k)
    mp.dps = prev
    print(f"precomputed ζ(2..{max_k+1}) at dps={dps}", flush=True)


def main() -> None:
    gammas = load_gammas()
    print(f"{len(gammas)} zeros, γ1={gammas[0]}, γN={gammas[-1]}", flush=True)
    logx = log(mpf(str(X)))
    t0 = time.time()
    rows = []
    # Highest |z| ≈ 953; dps = |z|/ln10+35 ≈ 448. Fix 460 so ζ(k) is cached once.
    DPS = 460
    mp.dps = DPS
    precompute_zeta(max_k=4200, dps=DPS)

    for i, g in enumerate(gammas, 1):
        z = mpc(mpf("0.5"), mpf(str(g))) * logx
        val = eri_gram(z, dps=DPS)
        v = float(2 * val.real)  # ρ and conjugate
        rows.append(
            {
                "n": i,
                "gamma": g,
                "eri_re": float(val.real),
                "eri_im": float(val.imag),
                "pair": v,
            }
        )
        if i == 1 or i % 200 == 0 or i == len(gammas):
            elapsed = time.time() - t0
            print(
                f"n={i:5d}/{len(gammas)}  γ={g:10.4f}  pair={v:12.6f}  "
                f"{elapsed:.1f}s",
                flush=True,
            )

    with OUT.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["n", "gamma", "eri_re", "eri_im", "pair"])
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {OUT}  ({len(rows)} rows, {time.time()-t0:.1f}s)", flush=True)


if __name__ == "__main__":
    main()
