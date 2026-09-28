"""Gram-series ERi, matching Grobner's entire-function definition (1.1).

    ERi(z) = 1 + Σ_{k≥1} z^k / (k · k! · ζ(k+1))

Avoids Ei branches (Elliott's Ei differs from ExpIntegralEi by πi sgn(Im s)).
Precision and truncation follow the numerics notes:
    dps = |z|/ln(10) + 35,   J = 4|z| + 90.
"""
from __future__ import annotations

import math

from mpmath import mp, mpf, mpc, zeta, log, pi, exp, sin, cos, sqrt

# Cache ζ(k+1) at working dps. Keys are (dps, k).
_zeta_cache: dict[tuple[int, int], object] = {}


def _zeta_int(k: int, dps: int):
    key = (dps, k)
    if key not in _zeta_cache:
        prev = mp.dps
        mp.dps = dps
        _zeta_cache[key] = zeta(k)
        mp.dps = prev
    return _zeta_cache[key]


def eri_gram(z, dps: int | None = None, nterms: int | None = None):
    """ERi(z) via Gram's series at working precision `dps`."""
    z = mpc(z)
    az = abs(z)
    if az == 0:
        return mp.one
    if dps is None:
        dps = int(az / mp.log(10) + 35) + 5
    if nterms is None:
        nterms = int(4 * az + 90) + 5
    dps = max(dps, 25)
    prev = mp.dps
    mp.dps = dps
    try:
        s = mp.one
        # u_k = z^k / k!
        u = mp.one
        for k in range(1, nterms + 1):
            u *= z / k
            term = u / (k * _zeta_int(k + 1, dps))
            s += term
            # Terms peak at k ~ |z| (size ~ e^{|z|}); do not stop near the peak.
            if k > 2 * az + 20 and abs(term) < mpf("10") ** (-(dps - 8)):
                break
        return +s
    finally:
        mp.dps = prev


def R_of_x(x: float, dps: int = 40):
    """R(x) = ERi(log x) for real x > 1."""
    prev = mp.dps
    mp.dps = dps
    try:
        return eri_gram(log(mpf(x)), dps=dps)
    finally:
        mp.dps = prev


def trivial_closed_form(x: float):
    """Grobner (0.4): 1/log x − (1/π) arctan(π/log x)."""
    lx = log(mpf(x))
    return 1 / lx - (1 / pi) * mp.atan(pi / lx)


def averaged_prime_count(x: float) -> float:
    """Count primes through x, assigning half weight when x is prime."""
    x = float(x)
    if not math.isfinite(x) or x <= 1:
        raise ValueError("x must be finite and greater than 1")
    n = math.floor(x)
    primes: list[int] = []
    for candidate in range(2, n + 1):
        if all(candidate % p for p in primes if p * p <= candidate):
            primes.append(candidate)
    return float(len(primes)) - (0.5 if x == n and primes and primes[-1] == n else 0.0)


def L_regularized(x: float = 1.1, dps: int = 40):
    """Numerical candidate L = R(x) − π₀(x) − I(x), not a proved contour centre."""
    pi0 = averaged_prime_count(x)
    prev = mp.dps
    mp.dps = dps
    try:
        Rx = R_of_x(x, dps=dps)
        triv = trivial_closed_form(x)
        L = Rx.real - pi0 - triv
        return {
            "x": x,
            "logx": float(log(mpf(x))),
            "R": float(Rx.real),
            "trivial": float(triv),
            "pi0": pi0,
            "L": float(L),
        }
    finally:
        mp.dps = prev
