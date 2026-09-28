import sys, time, math, csv, json
import numpy as np
sys.path.insert(0, 'compute')
from mobius_ei import eri
from eri import R_of_x, averaged_prime_count
from mpmath import mp, mpc, beta as Bt, digamma as psi, re as RE, mpf
mp.dps = 30

g = np.loadtxt('data/zeros1.txt')
res = []
for a, b, c in np.load('data/residues_400.npy'):
    res.append((float(a), complex(float(b), float(c))))

def trivial(lx):                 # closed form for the trivial-zero sum
    return 1.0 / lx - math.atan(math.pi / lx) / math.pi

out = {}
for x in [1.5, 2.0, 3.0]:
    lx = math.log(x)
    Rx = float(R_of_x(x).real)  # Gram series supports real log x
    pi0 = averaged_prime_count(x)
    L = Rx - pi0 - trivial(lx)
    print('x=%s  log x=%.6f  R(x)=%.9f  pi0=%.1f  triv=%.9f  L=%.9f'
          % (x, lx, Rx, pi0, trivial(lx), L), flush=True)

    t0 = time.time(); pair = np.empty(len(g))
    for i, gi in enumerate(g):
        pair[i] = 2 * eri(complex(0.5, gi) * lx).real
        if i % 25000 == 0: print('   %6d  %.0fs' % (i, time.time() - t0), flush=True)
    print('   ERi done %.0fs' % (time.time() - t0), flush=True)
    np.save('data/pair_100k_x%s.npy' % x, pair)

    # --- comparison with the first 60 positive amplitudes ------------
    S60 = sum(2 * abs(r) / abs(complex(0.5, gamma)) for gamma, r in res[:60])
    Rpred = S60 / (math.pi * math.sqrt(lx))  # historical JSON key; finite partial sum
    csum = np.cumsum(pair)
    def SR(T):
        j = np.searchsorted(g, T, 'right')
        return csum[j - 1] if j else 0.0
    rows = []
    lo = 1024.0
    while lo < g[-1]:
        hi = min(2 * lo, g[-1])
        Ts = np.linspace(lo, hi, 400)
        m = max(abs(SR(T) - L) / (math.sqrt(T) * math.log(T / (2 * math.pi))) for T in Ts)
        rows.append((lo, hi, m, m / Rpred)); lo = hi
    out[x] = dict(logx=lx, L=L, Rpred=Rpred,
                  normalization_terms=60, normalization_kind='partial positive amplitude sum',
                  octaves=[[a, b, c, d] for a, b, c, d in rows],
                  octaves_complete=[b == 2 * a for a, b, c, d in rows])
    print('   first-60 partial amplitude R60(x) = %.6f' % Rpred, flush=True)
    for a, b, m, r in rows:
        print('     %7.0f-%7.0f   max|F| = %.5f   ratio = %.3f%s' %
              (a, b, m, r, '' if b == 2 * a else ' (incomplete octave)'), flush=True)

    # --- fitted transient coefficient; not a derived asymptotic ----------
    pre = []
    for gam, r in res:
        rho = mpc(0.5, gam)
        pre.append((rho, (2 / (math.pi * lx)) * mpc(r.real, r.imag) * mp.power(lx, rho) / rho))
    def model(T, k):
        tot = mpf(0)
        for rho, A in pre:
            fac = mp.power(T, rho) * k * Bt(rho + 1, k)
            brk = mp.log(T / (2 * mp.pi)) + psi(rho + 1) - psi(rho + 1 + k) - 1 / rho
            tot += RE(A * fac * brk)
        return float(tot)
    def riesz(T, k):
        j = np.searchsorted(g, T, 'right')
        return float(np.sum((1 - g[:j] / T) ** k * pair[:j]))
    Is = []
    for k in (1, 2, 3):
        for T in (2048., 4096., 8192.):
            Is.append(((riesz(T, k) - L) - model(T, k)) * T / k)
    out[x]['tau_estimates'] = Is
    print('   fitted tau(x) estimates (resid*T/k): ' + ' '.join('%.1f' % v for v in Is), flush=True)
    print(flush=True)

json.dump(out, open('data/x_sweep.json', 'w'), indent=1)
print('written data/x_sweep.json')
