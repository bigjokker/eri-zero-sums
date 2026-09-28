"""Theorem A (Paper A) at q = 4, evaluated on a list of ordinates of zeros of L(s, chi_{-4}).

    python compute/theoremA_q4.py data/chi4_zeros_pari10000_clean.txt [--eri]

Theorem A (manuscripts/real-character-prime-counts.md §2.4).  chi = chi_{-4} is odd, so I(y) = int_y^oo dt/((t^2-1) log t).
    pi*(x;4,1) - pi*(x;4,3) + pi*(sqrt x;4,3) + sum_{|gamma|<=T} m_rho R_N(x^rho)
        = C(x) - log(2)*m(N) + E(x,T), E -> 0; m(N) = sum_{n<=N} mu(n)/n.
    R_N(x^rho) = sum_{n<=N} mu(n)/n * li(x^{rho/n}),   li(y^rho) := Ei(rho log y),   N = floor(log x / log 2)
    C(x)       = sum_{n<=N} mu(n)/n * I(x^{1/n}),       I(y) = sum_{k>=0} E1((2k+1) log y)      (§2.3)
Branch: Elliott's, Ei(w) = -E1(-w) for Im w > 0 (zeta paper App. B); the conjugate zero gives the
conjugate term, so the sum over |gamma|<=T is 2 Re of the sum over 0<gamma<=T.  m_rho = 1 on the list.

Two more quantities on the same zeros:
    B_T        = 2 Re sum_{0<gamma<=T} li(2^rho) - I(2)      the lower-limit bracket of (†); Paper A §2.3 lets
                                                              evaluated limit B_T = -log 2
    Etilde(x,T)= E(x,T) with R = ERi in place of R_N           Paper B: Omega_pm(sqrt T log T)   (--eri, slower)
x is chosen so that no x^{1/n}, n<=N, is a prime power; the nearest distances are printed.
Output: data/theoremA_q4.txt, or data/theoremA_q4_eri.txt with --eri
"""

import sys, os, time, math, argparse
from pathlib import Path
import mpmath as mp

mp.mp.dps = 20
XS = [mp.mpf('30.5'), mp.mpf('100.5'), mp.mpf('1500.5')]
CHECKPOINTS = [50, 100, 200, 300, 500, 700, 1000, 1500, 2000, 3000, 5000, 7000, 10000, 15000, 20000]


def mobius(n):
    if n == 1:
        return 1
    m, r, p = n, 1, 2
    while p * p <= m:
        if m % p == 0:
            m //= p
            if m % p == 0:
                return 0
            r = -r
        p += 1
    if m > 1:
        r = -r
    return r


def primes_upto(n):
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]


def is_prime_power(m):
    if m < 2:
        return False
    for p in primes_upto(int(m ** 0.5) + 1) + [m]:
        if m % p == 0:
            while m % p == 0:
                m //= p
            return m == 1
    return False


def dist_to_prime_power(y):
    y = float(y)
    lo, hi = int(math.floor(y)), int(math.ceil(y))
    while not is_prime_power(lo) and lo > 1:
        lo -= 1
    while not is_prime_power(hi):
        hi += 1
    return min(y - lo, hi - y)


def Ei(w):
    """Elliott's branch: Ei(w) = -E1(-w) for Im w > 0."""
    return -mp.e1(-w)


def I_odd(y):
    L = mp.log(y)
    return mp.nsum(lambda k: mp.e1((2 * k + 1) * L), [0, mp.inf])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('zeros', type=Path)
    parser.add_argument('--eri', action='store_true')
    parser.add_argument('--output', type=Path, help='Explicit output path for staged runs')
    args = parser.parse_args()
    zfile = args.zeros
    use_eri = args.eri
    from validate_zeros import read
    read(zfile)  # Reject invalid inputs before creating any output.
    gammas = [mp.mpf(l.strip()) for l in zfile.read_text(encoding='utf-8-sig').splitlines()
              if l.strip() and not l.strip().startswith('#')]
    Tmax = gammas[-1]
    output = args.output or Path(__file__).resolve().parents[1] / 'data' / ('theoremA_q4' + ('_eri' if use_eri else '') + '.txt')
    out = output.open('w', encoding='utf-8')
    def P(*a):
        s = ' '.join(str(x) for x in a); print(s, flush=True); out.write(s + '\n'); out.flush()
    P('# Theorem A at q=4 on %d zeros of L(s,chi_-4) from %s, gamma <= %s; mpmath dps=%d'
      % (len(gammas), zfile, mp.nstr(Tmax, 8), mp.mp.dps))

    # branch cross-check against the project evaluator, at one w
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    try:
        from mobius_ei import ei_elliott, eri
        w = (mp.mpf('0.5') + 1j * gammas[0]) * mp.log(XS[1])
        d = abs(complex(Ei(w)) - complex(ei_elliott(complex(w))[0]))
        P('# branch check vs compute/mobius_ei.ei_elliott at w=%s: |diff| = %.1e' % (mp.nstr(w, 6), d))
    except Exception as e:
        if use_eri:
            raise RuntimeError('requested --eri evaluator failed') from e
        P('# optional evaluator cross-check failed (%s)' % e)

    primes = primes_upto(int(max(XS)) + 2)
    I2 = I_odd(mp.mpf(2))
    P('# I(2) = %s ; Paper A gives lim B_T = -log 2 = %s (log L(0,chi_-4))' % (mp.nstr(I2, 12), mp.nstr(-mp.log(2), 12)))

    for x in XS:
        N = int(mp.floor(mp.log(x) / mp.log(2)))
        roots = [x ** (mp.mpf(1) / n) for n in range(1, N + 1)]
        dmin = min(dist_to_prime_power(r) for r in roots)
        mu = [mobius(n) for n in range(N + 1)]
        mN = sum(mp.mpf(mu[n]) / n for n in range(1, N + 1))
        p1 = sum(1 for p in primes if p <= x and p % 4 == 1)
        p3 = sum(1 for p in primes if p <= x and p % 4 == 3)
        p3s = sum(1 for p in primes if p <= mp.sqrt(x) and p % 4 == 3)
        C = sum(mp.mpf(mu[n]) / n * I_odd(roots[n - 1]) for n in range(1, N + 1))
        combin = mp.mpf(p1 - p3 + p3s)
        P('')
        P('== x = %s : N = %d, min_n <x^{1/n}> = %.3f, m(N) = sum_{n<=N} mu(n)/n = %s' % (x, N, dmin, mp.nstr(mN, 8)))
        P('   pi*(x;4,1) = %d, pi*(x;4,3) = %d, pi*(sqrt x;4,3) = %d  ->  combinatorial side = %d' % (p1, p3, p3s, int(combin)))
        P('   C(x) = %s ; kappa*m(N) = %s; corrected E(x,T) tends to 0' % (mp.nstr(C, 12), mp.nstr(-mp.log(2) * mN, 8)))
        P('   %8s  %14s  %14s  %24s%s' % ('T', 'E(x,T)', 'B_T', 'E - m(N)(B_T + log2)', '  Etilde(x,T)' if use_eri else ''))
        S = mp.mpf(0); B = mp.mpf(0); Se = mp.mpf(0)
        cps = [c for c in CHECKPOINTS if c <= Tmax] + [Tmax]
        ci = 0
        logx = mp.log(x); log2 = mp.log(2)
        oct_lo, oct_max, oct_min = 100.0, -1e9, 1e9
        t0 = time.time()
        for j, g in enumerate(gammas):
            rho = mp.mpf('0.5') + 1j * g
            RN = sum(mp.mpf(mu[n]) / n * Ei(rho * logx / n) for n in range(1, N + 1))
            S += 2 * RN.real
            B += 2 * Ei(rho * log2).real
            if use_eri:
                Se += 2 * complex(eri(complex(rho * logx))).real
            E = combin + S - C + mp.log(2) * mN
            if g >= oct_lo * 2:
                P('   octave [%g,%g): max E = %+.4f, min E = %+.4f' % (oct_lo, oct_lo * 2, oct_max, oct_min))
                oct_lo *= 2; oct_max, oct_min = -1e9, 1e9
            oct_max = max(oct_max, float(E)); oct_min = min(oct_min, float(E))
            while ci < len(cps) and (g >= cps[ci] or j == len(gammas) - 1):
                BT = B - I2
                line = '   %8s  %+14.6f  %+14.6f  %+24.6f' % (mp.nstr(min(g, cps[ci]), 7), float(E), float(BT), float(E - mN * (BT + mp.log(2))))
                if use_eri:
                    line += '  %+14.4f' % float(combin + Se - C + mp.log(2) * mN)
                P(line)
                ci += 1
                if j == len(gammas) - 1:
                    break
        P('   (%.0f s)' % (time.time() - t0))
    out.close()


if __name__ == '__main__':
    main()
