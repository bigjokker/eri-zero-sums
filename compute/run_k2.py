"""Archived Riesz comparison; amplitude guides use 60 terms and are not upper bounds.
Run from the repository root after obtaining the external zeros1 input.
"""
import csv, math, sys, time, numpy as np
sys.path.insert(0, 'compute')
from mobius_ei import eri, choose_N
logx = math.log(1.1); L = -8.94183647352

g = np.loadtxt('data/zeros1.txt')
print('zeros loaded: %d, gamma_max = %.4f' % (len(g), g[-1]), flush=True)
for N in sorted({choose_N(abs(complex(0.5, gi) * logx)) for gi in g[::500]}):
    t0 = time.time(); __import__('mobius_ei')._tail_data(N)
    print('  tail data N=%-5d %.1fs' % (N, time.time() - t0), flush=True)

t0 = time.time(); pair = np.empty(len(g))
for i, gi in enumerate(g):
    pair[i] = 2 * eri(complex(0.5, gi) * logx).real
    if i % 20000 == 0: print('   %6d  %.0fs' % (i, time.time() - t0), flush=True)
print('ERi done in %.0fs' % (time.time() - t0), flush=True)
np.save('data/pair_100k.npy', pair)

def riesz(T, k):
    m = g <= T; return np.sum((1 - g[m] / T) ** k * pair[m])
pref = lambda T: math.sqrt(T) * math.log(T / (2 * math.pi)) / (math.pi * math.sqrt(logx))
# First-60 partial amplitudes; the final octave is incomplete.
S = {0: 0.03793060, 1: 0.00175808, 2: 0.00020064, 3: 0.00003728}

print('\n  T        raw       k=1      k=2      k=3  |  |k2-L|   60-term guide   refuted law')
rows = []
for T in [9999, 16000, 25000, 40000, 60000, 74920]:
    r = [riesz(T, k) for k in range(4)]
    d2 = abs(r[2] - L); mine = pref(T) * S[2]
    old = 0.1525 * math.sqrt(9999 / T) * math.log(T / (2 * math.pi)) / math.log(9999 / (2 * math.pi))
    print('%6d %9.2f %8.3f %8.3f %8.3f | %7.3f  %7.3f     %7.3f' % (T, r[0], r[1], r[2], r[3], d2, mine, old), flush=True)
    rows.append((T, d2, mine, old))

print('\nper-octave max |k=2 - L| (finite amplitude comparison):')
lo = 9999
while lo < 74920:
    hi = min(2 * lo, 74920)
    Ts = np.linspace(lo, hi, 200)
    m = max(abs(riesz(T, 2) - L) for T in Ts)
    print('  %6d-%6d  max %.4f   60-term guide %.4f   refuted %.4f' % (
        lo, hi, m, pref(hi) * S[2],
        0.1525 * math.sqrt(9999 / hi) * math.log(hi / (2 * math.pi)) / math.log(9999 / (2 * math.pi))), flush=True)
    lo = hi
