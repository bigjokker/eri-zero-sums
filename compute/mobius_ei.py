"""ERi(z) by Mobius-Ei in double precision.

Gram's series needs ~|z|/ln(10) digits (4180 at |z|=9531).  This needs 16.

    ERi(z) = sum_{n<=N} mu(n)/n * Ei(z/n)  +  tail(N),        N = 2^ceil(log2(|z|/TAILARG))

    tail(N) = (gammaE + Log(-z))*A0 - A1 + sum_{k>=1} z^k/(k k!) * Gk
    A0 = sum_{n>N} mu(n)/n      = -m(N)                  [ sum mu(n)/n      =  0 ]
    A1 = sum_{n>N} mu(n)log n/n = -1 - sum_{n<=N} ...    [ sum mu(n)log n/n = -1 ]
    Gk = sum_{n>N} mu(n)/n^{k+1}

Gk MUST be formed in high precision: as 1/zeta(k+1) - partial it is the difference of
two numbers ~1 whose difference is ~N^{-k} (1e-139 at k=60, N=191).  In double that
returns 1e-16 of noise, which the factor z^k/k! ~ 1e96 then amplifies to 1e80.

Branch: Elliott's. For Im w > 0, Ei(w) = -E1(-w) = gammaE + Log(-w)
+ sum_k w^k/(k k!). Thus Ei(w) is the usual principal Ei minus i*pi.
The same branch must be used in the head and tail; their constant shifts cancel.
"""
import math
import numpy as np
from numpy.polynomial.laguerre import laggauss

GAMMA_E = 0.5772156649015328606
TAILARG = 5.0
WSPLIT  = 12.0     # |w| below this -> power series; above -> Gauss-Laguerre
KMAX    = 48      # terms of the tail series; (|z|/N)^k/(k^2 k!) < 1e-23 by k=40
MLAG    = 80

def mobius_sieve(N):
    mu = np.ones(N + 1, dtype=np.int64); primes = []; is_c = np.zeros(N + 1, dtype=bool)
    for i in range(2, N + 1):
        if not is_c[i]:
            primes.append(i); mu[i] = -1
        for p in primes:
            if i * p > N: break
            is_c[i * p] = True
            if i % p == 0: mu[i * p] = 0; break
            mu[i * p] = -mu[i]
    mu[0] = 0
    return mu

_LAGQ = laggauss(MLAG)

def _ei_large(w):
    """Ei on Elliott's branch, Im w > 0, |w| large: Ei(w) = -E1(-w),
    E1(v) = e^{-v} int_0^inf e^{-t}/(t+v) dt by Gauss-Laguerre."""
    t, wq = _LAGQ
    v = -w
    E1 = np.exp(-v) * np.sum(wq[:, None] / (t[:, None] + v[None, :]), axis=0)
    return -E1

def _ei_small(w):
    """Elliott Ei(w) = gammaE + Log(-w) + sum_{k>=1} w^k/(k k!)."""
    out = GAMMA_E + np.log(-w)
    term = np.ones_like(w); acc = np.zeros_like(w)
    for k in range(1, 200):
        term = term * w / k
        acc = acc + term / k
        if k > 8 and np.max(np.abs(term / k)) < 1e-18 * np.max(np.abs(acc)): break
    return out + acc

def ei_elliott(w):
    w = np.atleast_1d(np.asarray(w, dtype=complex))
    if not np.isfinite(w).all() or np.any(w.imag <= 0):
        raise ValueError('fast Elliott Ei requires finite w with Im w > 0')
    out = np.empty_like(w)
    big = np.abs(w) >= WSPLIT
    if big.any():   out[big]  = _ei_large(w[big])
    if (~big).any():out[~big] = _ei_small(w[~big])
    return out

_CACHE = {}
def _tail_data(N):
    if N in _CACHE: return _CACHE[N]
    from mpmath import mp, zeta, mpf
    # G_k = 1/zeta(k+1) - partial is a difference of two numbers ~1 whose value is ~N^-k.
    # Resolving it needs k*log10(N) significant digits; a fixed dps silently returns noise
    # once k*log10(N) exceeds it, and z^k/k! then amplifies that noise without limit.
    dps = int(KMAX * math.log10(N)) + 60
    mu = mobius_sieve(N)
    n = np.arange(1, N + 1); m = mu[1:N + 1].astype(float)
    nz = m != 0; idx = n[nz]; mun = m[nz]
    A0 = -float(np.sum(mun / idx))
    A1 = -1.0 - float(np.sum(mun * np.log(idx) / idx))
    Gk = np.empty(KMAX)
    with mp.workdps(dps):
        for k in range(1, KMAX + 1):
            s = mpf(1) / zeta(k + 1)
            for i, mi in zip(idx, mun):
                s -= mpf(int(mi)) / mpf(int(i)) ** (k + 1)
            Gk[k - 1] = float(s)
            if s != 0 and Gk[k - 1] == 0.0:
                raise RuntimeError('G_%d underflowed double precision at N=%d' % (k, N))
            if Gk[k - 1] != 0.0 and abs(Gk[k - 1]) < 10.0 ** (-(dps - 20)):
                raise RuntimeError('G_%d underflowed working precision at N=%d' % (k, N))
    _CACHE[N] = (idx, mun, A0, A1, Gk)
    return _CACHE[N]

def choose_N(az):
    return int(2 ** max(5, int(np.ceil(np.log2(max(az / TAILARG, 32.0))))))

def eri(z):
    """ERi(z), single complex z with Im z > 0."""
    z = complex(z)
    if not (math.isfinite(z.real) and math.isfinite(z.imag) and z.imag > 0):
        raise ValueError('fast ERi requires finite z with Im z > 0; use eri_gram elsewhere')
    idx, mun, A0, A1, Gk = _tail_data(choose_N(abs(z)))
    head = np.sum(mun / idx * ei_elliott(z / idx))
    acc = 0.0 + 0j; term = 1.0 + 0j
    for k in range(1, KMAX + 1):
        term = term * z / k
        acc += term / k * Gk[k - 1]
    return head + (GAMMA_E + np.log(-z)) * A0 - A1 + acc
