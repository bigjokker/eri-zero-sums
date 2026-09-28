"""Fejér means, residue-model check, and Bessel-product fit from eri_list.csv."""
from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import j0

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
FIGS = ROOT / "data" / "figures"
FIGS.mkdir(parents=True, exist_ok=True)

X = 1.1
LOGX = math.log(X)
SQRT_LOGX = math.sqrt(LOGX)
PI = math.pi
L = -8.941836473520125  # R(x) - (0.4), π0=0; from compute.eri.L_regularized


def load_eri():
    gammas, pairs = [], []
    with (DATA / "eri_list.csv").open() as f:
        for row in csv.DictReader(f):
            gammas.append(float(row["gamma"]))
            pairs.append(float(row["pair"]))
    return np.asarray(gammas), np.asarray(pairs)


def load_residues():
    rows = []
    with (DATA / "residues_60.csv").open() as f:
        for row in csv.DictReader(f):
            rows.append(row)
    return rows


def prefix_sums(pairs):
    return np.cumsum(pairs)


def F_of_T(S, T):
    """(ΣR_T - L) / (√T log(T/2π))."""
    return (S - L) / (np.sqrt(T) * np.log(T / (2 * PI)))


def fejer(gammas, pairs, T, k=1):
    """Σ (1 - γ/T)_+^k pair."""
    w = np.clip(1.0 - gammas / T, 0.0, None) ** k
    return float(np.dot(w, pairs))


def residue_model(gammas_r, r_abs, r_re, r_im, T):
    """L + (1/(π log x)) · 2 Re Σ r_ρ (Y^ρ/ρ) (log(Y/C) - 1/ρ), Y=T log x, C=2π log x."""
    Y = T * LOGX
    C = 2 * PI * LOGX
    logYC = math.log(Y / C)
    acc = 0j
    for g, ra, rr, ri in zip(gammas_r, r_abs, r_re, r_im):
        rho = 0.5 + 1j * g
        r = rr + 1j * ri
        term = r * (Y ** rho) / rho * (logYC - 1 / rho)
        acc += term
    return L + (1 / (PI * LOGX)) * 2 * acc.real


def bessel_density(amps, x_grid, xi_max=80.0, n_xi=4000):
    """p(x) = (1/2π) ∫ e^{-i ξ x} ∏ J_0(2 a_n ξ) dξ, amps = |a_n|."""
    xi = np.linspace(0.0, xi_max, n_xi)
    dxi = xi[1] - xi[0]
    log_phi = np.zeros_like(xi)
    for a in amps:
        # J0(0)=1; skip xi=0 in log
        v = j0(2 * a * xi)
        # clip to avoid log(0)
        v = np.clip(np.abs(v), 1e-16, None) * np.sign(np.where(v == 0, 1, v))
        # product of possibly negative J0: track sign separately
    phi = np.ones_like(xi)
    for a in amps:
        phi *= j0(2 * a * xi)
    phi[0] = 1.0
    # even function: p(x) = (1/π) ∫_0^∞ cos(ξ x) φ(ξ) dξ
    p = np.zeros_like(x_grid)
    for i, x in enumerate(x_grid):
        p[i] = (1 / PI) * np.trapezoid(phi * np.cos(xi * x), dx=dxi)
    return p


def main():
    g, v = load_eri()
    S = prefix_sums(v)
    T = g  # evaluate just after including that zero
    F = F_of_T(S, T)

    # --- Fejér ---
    marks = np.array([64.0, 128, 256, 512, 1024, 2048, 4096, 8192, g[-1]])
    T_grid = np.unique(np.concatenate([marks, np.geomspace(64, g[-1], 80)]))
    fej1 = np.array([fejer(g, v, t, 1) for t in T_grid])
    fej2 = np.array([fejer(g, v, t, 2) for t in T_grid])
    fej3 = np.array([fejer(g, v, t, 3) for t in T_grid])
    raw = np.interp(T_grid, T, S)

    print("L =", L)
    print("raw ΣR at Tmax =", S[-1], "  vs L  delta=", S[-1] - L)
    print("Riesz at Tmax: k=1", fej1[-1], "Δ", fej1[-1] - L)
    print("               k=2", fej2[-1], "Δ", fej2[-1] - L)
    print("               k=3", fej3[-1], "Δ", fej3[-1] - L)
    print(f"{'T':>10} {'ΣR':>10} {'k=1':>10} {'k=2':>10} {'k=3':>10}  |k1-L| |k2-L| |k3-L|")
    for t in marks:
        r = float(np.interp(t, T, S))
        f1, f2, f3 = fejer(g, v, t, 1), fejer(g, v, t, 2), fejer(g, v, t, 3)
        print(
            f"{t:10.2f} {r:10.4f} {f1:10.4f} {f2:10.4f} {f3:10.4f}  "
            f"{abs(f1-L):6.3f} {abs(f2-L):6.3f} {abs(f3-L):6.3f}"
        )

    # --- F(u) stats ---
    mask = T >= 64
    Fm = F[mask]
    print()
    print("F(u) for T≥64:")
    print(f"  n={mask.sum()}  min={Fm.min():.4f}  max={Fm.max():.4f}  "
          f"rms={np.sqrt(np.mean(Fm**2)):.4f}  mean={Fm.mean():.4f}")

    res = load_residues()
    gammas_r = np.array([float(r["gamma"]) for r in res])
    r_abs = np.array([float(r["r_gamma_abs"]) for r in res])
    r_re = np.array([float(r["r_gamma_re"]) for r in res])
    r_im = np.array([float(r["r_gamma_im"]) for r in res])
    rho_abs = np.sqrt(0.25 + gammas_r ** 2)
    amps = (1 / (PI * SQRT_LOGX)) * r_abs / rho_abs  # |a_ρ|
    Renv = float(np.sum(2 * amps))
    print(f"  first-60 amplitude sum = {Renv:.6f} (partial envelope)")
    print(f"  observed max|F| / first-60 sum = {max(abs(Fm.min()), Fm.max()) / Renv:.3f}")
    imin, imax = int(np.argmin(Fm)), int(np.argmax(Fm))
    Tm = T[mask]
    print(f"  min F at T={Tm[imin]:.1f}, max F at T={Tm[imax]:.1f}")

    print()
    print("per-octave max |F| (compare numerics 0.014–0.026):")
    print(f"{'T range':>14} {'max|ΣR-L|':>12} {'max|F|':>10} {'max/√T':>10}")
    edges = [64, 128, 256, 512, 1024, 2048, 4096, 8192, 10000]
    for a, b in zip(edges, edges[1:]):
        m = (T >= a) & (T < b)
        if not np.any(m):
            continue
        d = np.max(np.abs(S[m] - L))
        fmax = np.max(np.abs(F[m]))
        mx = np.max(np.abs(S[m] - L) / np.sqrt(T[m]))
        suffix = ' (incomplete octave)' if T[-1] < b else ''
        print(f"  {a:5d}–{b:<5d} {d:12.3f} {fmax:10.4f} {mx:10.4f}{suffix}")

    # model vs data at the numerics checkpoints
    checkpoints = [1238, 2188, 3068, 4722, 5514, 7051, 7802, 8543, 9274, 9998]
    print()
    print("residue model (60 ρ) vs full sum:")
    print(f"{'T':>8} {'ΣR-L':>10} {'model':>10} {'ratio':>8}")
    for t0 in checkpoints:
        idx = int(np.searchsorted(T, t0, side="right") - 1)
        if idx < 0:
            continue
        t = T[idx]
        data = S[idx] - L
        model = residue_model(gammas_r, r_abs, r_re, r_im, t) - L
        ratio = data / model if model != 0 else float("nan")
        print(f"{t:8.1f} {data:10.3f} {model:10.3f} {ratio:8.3f}")

    # --- plots ---
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(T, S, lw=0.7, color="0.6", label=r"$\Sigma R_T$")
    ax.plot(T_grid, fej1, lw=1.3, label=r"Riesz $k=1$")
    ax.plot(T_grid, fej2, lw=1.3, label=r"Riesz $k=2$")
    ax.plot(T_grid, fej3, lw=1.1, ls="--", label=r"Riesz $k=3$")
    ax.axhline(L, color="k", ls="--", lw=0.8, label=f"numerical candidate L = {L:.4f}")
    ax.set_xscale("log")
    ax.set_xlabel("T")
    ax.set_ylabel("sum")
    ax.legend()
    ax.set_title(r"Raw $\Sigma R_T$ vs Riesz means at $x=1.1$")
    fig.tight_layout()
    fig.savefig(FIGS / "fejer.png", dpi=140)
    plt.close()

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(T[mask], Fm, lw=0.7)
    ax.axhline(Renv, color="C1", ls="--", lw=0.8, label=f"first-60 amplitude {Renv:.4f}")
    ax.axhline(-Renv, color="C1", ls="--", lw=0.8)
    ax.set_xscale("log")
    ax.set_xlabel("T")
    ax.set_ylabel("F(u)")
    ax.set_title(r"$(\Sigma R_T-L)/(\sqrt{T}\log(T/2\pi))$")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGS / "F_of_T.png", dpi=140)
    plt.close()

    # histogram vs Bessel density
    x_grid = np.linspace(-Renv * 1.05, Renv * 1.05, 400)
    pdf = bessel_density(amps, x_grid)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.hist(Fm, bins=40, density=True, alpha=0.45, label="empirical F (T≥64)")
    ax.plot(x_grid, np.clip(pdf, 0, None), lw=1.6, label=r"Bessel product density (60 zeros)")
    ax.axvline(Renv, color="C1", ls="--", lw=0.8)
    ax.axvline(-Renv, color="C1", ls="--", lw=0.8)
    ax.set_xlabel("F")
    ax.set_ylabel("density")
    ax.legend()
    ax.set_title("Empirical F vs Bessel-product law (short window)")
    fig.tight_layout()
    fig.savefig(FIGS / "bessel_fit.png", dpi=140)
    plt.close()

    print()
    print("wrote", FIGS / "fejer.png")
    print("wrote", FIGS / "F_of_T.png")
    print("wrote", FIGS / "bessel_fit.png")


if __name__ == "__main__":
    main()
