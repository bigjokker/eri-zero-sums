"""Workstream A: r_ρ from the closed form, not from a numerical Mellin pole.

    r_ρ = Γ(1-ρ) sin(πρ/2) / ((ρ-1) ζ'(ρ))
        = π (2π)^{-ρ} χ(ρ) / ((ρ-1) ζ'(ρ))

On Re ρ = 1/2, |Γ(1-ρ) sin(πρ/2)| = √(π/2), so
    |r_ρ| = √(π/2) / (|ρ-1| |ζ'(ρ)|).

Target (numerics notes, ρ1..ρ8):
    0.1117  0.05243  0.03652  0.03159
    0.02753 0.01722  0.02055  0.01578
"""
from __future__ import annotations

import csv
import argparse
from pathlib import Path

from mpmath import mp, gamma, sin, pi, sqrt, zeta, zetazero, power, conj, exp, log

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data"
OUT.mkdir(exist_ok=True)

TARGET = [
    0.1117,
    0.05243,
    0.03652,
    0.03159,
    0.02753,
    0.01722,
    0.02055,
    0.01578,
]


def chi(s):
    """χ(s) = 2^s π^{s-1} sin(πs/2) Γ(1-s)."""
    return power(2, s) * power(pi, s - 1) * sin(pi * s / 2) * gamma(1 - s)


def residue_pair(rho, zpp):
    """Return (r_from_gamma, r_from_chi, abs_on_line, |ζ'(ρ)|)."""
    r_gamma = gamma(1 - rho) * sin(pi * rho / 2) / ((rho - 1) * zpp)
    r_chi = pi * power(2 * pi, -rho) * chi(rho) / ((rho - 1) * zpp)
    abs_line = sqrt(pi / 2) / (abs(rho - 1) * abs(zpp))
    return r_gamma, r_chi, abs_line, abs(zpp)


def compute(n_zeros: int = 60, dps: int = 40) -> list[dict]:
    mp.dps = dps
    rows = []
    for n in range(1, n_zeros + 1):
        rho = zetazero(n)
        # ζ' at a simple zero: mpmath's zeta(s, derivative=1)
        zpp = zeta(rho, derivative=1)
        r_g, r_c, abs_line, abs_zpp = residue_pair(rho, zpp)
        gamma_n = float(rho.imag)
        row = {
            "n": n,
            "gamma": gamma_n,
            "re_rho": float(rho.real),
            "zeta_prime_abs": float(abs_zpp),
            "zeta_prime_re": float(zpp.real),
            "zeta_prime_im": float(zpp.imag),
            "r_gamma_re": float(r_g.real),
            "r_gamma_im": float(r_g.imag),
            "r_gamma_abs": float(abs(r_g)),
            "r_chi_abs": float(abs(r_c)),
            "abs_on_line": float(abs_line),
            "gamma_chi_rel": float(abs(r_g - r_c) / abs(r_g)),
            "line_vs_gamma_rel": float(abs(abs(r_g) - abs_line) / abs(r_g)),
        }
        if n <= len(TARGET):
            row["target"] = TARGET[n - 1]
            row["abs_minus_target"] = row["r_gamma_abs"] - TARGET[n - 1]
        rows.append(row)
        print(
            f"{n:3d}  γ={gamma_n:12.6f}  |r|={row['r_gamma_abs']:.8f}"
            f"  |r|_line={row['abs_on_line']:.8f}  |ζ'|={row['zeta_prime_abs']:.6f}"
            f"  Δ(γ,χ)={row['gamma_chi_rel']:.1e}  Δ(line)={row['line_vs_gamma_rel']:.1e}",
            flush=True,
        )
    return rows


def write_csv(rows: list[dict], path: Path) -> None:
    fields = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def envelope_constant(rows: list[dict], x: float = 1.1) -> dict:
    """Σ 2|r_ρ|/|ρ| and the F-envelope (1/(π√log x)) Σ 2|r|/|ρ|."""
    from math import log as flog, pi as fpi, sqrt as fsqrt

    s = 0.0
    for row in rows:
        rho_abs = (0.25 + row["gamma"] ** 2) ** 0.5
        s += 2.0 * row["r_gamma_abs"] / rho_abs
    logx = flog(x)
    env = s / (fpi * fsqrt(logx))
    return {
        "n_zeros": len(rows),
        "sum_2_abs_r_over_abs_rho": s,
        "envelope_F": env,
        "logx": logx,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--count', type=int, default=60)
    parser.add_argument('--dps', type=int, default=40)
    parser.add_argument('--output', type=Path, help='CSV or NPY output; NPY columns are ordinate, Re r, Im r')
    args = parser.parse_args()
    if args.count < 1 or args.dps < 25:
        parser.error('require count >= 1 and dps >= 25')
    destination = args.output or OUT / f'residues_{args.count}.csv'
    if destination.suffix.lower() not in ('.csv', '.npy'):
        parser.error('output must end in .csv or .npy')
    rows = compute(args.count, dps=args.dps)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.suffix.lower() == '.npy':
        import numpy as np
        np.save(destination, np.array([[r['gamma'], r['r_gamma_re'], r['r_gamma_im']] for r in rows]))
    else:
        write_csv(rows, destination)
    env = envelope_constant(rows)
    print()
    print("target ρ1..ρ8:", TARGET)
    print("got    ρ1..ρ8:", [round(r["r_gamma_abs"], 5) for r in rows[:8]])
    print(
        f"partial envelope from {env['n_zeros']} zeros: "
        f"Σ 2|r|/|ρ| = {env['sum_2_abs_r_over_abs_rho']:.6f},  "
        f"F partial envelope = {env['envelope_F']:.6f}  (0.0379 is the unnormalised 60-term sum)"
    )
    # sanity: ρ1 arithmetic from the notes
    r1 = rows[0]
    print(
        f"ρ1 check: √(π/2)/(|ρ-1||ζ'|) = {r1['abs_on_line']:.6f}; "
        f"|ζ'|={r1['zeta_prime_abs']:.5f}; |ρ-1|={(0.25+r1['gamma']**2)**0.5:.4f}"
    )
