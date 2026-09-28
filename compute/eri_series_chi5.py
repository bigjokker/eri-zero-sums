"""ERi series over the zeros of L(s, chi) for the order-4 character chi mod 5, both signs of gamma.

    python compute/eri_series_chi5.py CHI_ZEROS CHIBAR_ZEROS [x ...]

Sum_{|gamma|<T} ERi(rho log x): the zeros with gamma > 0 come from the chi list; those with gamma < 0
are the conjugates of the zeros of L(s, chi-bar), so their terms are conj ERi(rho' log x) over the
chi-bar list (ERi(conj z) = conj ERi(z)).  Ordered by |gamma|.  Prints Re and Im at checkpoints and
the per-octave range of each: Paper B extended to complex chi predicts Re = Omega(sqrt T log T) and
Im convergent.  ERi from compute/mobius_ei.eri (Elliott branch, Im z > 0).
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from mobius_ei import eri


def read(f):
    return [float(l) for l in open(f, encoding='utf-8', errors='replace') if l.strip() and l.strip()[0].isdigit()]


def main():
    zc, zb = read(sys.argv[1]), read(sys.argv[2])
    xs = [float(a) for a in sys.argv[3:]] or [1.1, 2.0]
    terms = sorted([(g, +1) for g in zc] + [(g, -1) for g in zb])
    T = min(zc[-1], zb[-1])
    terms = [t for t in terms if t[0] <= T]
    print('# %d zeros of chi (gamma>0) + %d of chi-bar (as gamma<0), |gamma| <= %.3f' % (len(zc), len(zb), T))
    for x in xs:
        lx = np.log(x)
        S = 0j
        chk = {}
        oct_lo, oct = 50.0, []
        rows = []
        for g, sgn in terms:
            v = eri((0.5 + 1j * g) * lx)
            S += v if sgn > 0 else np.conj(v)
            if g >= 2 * oct_lo:
                if oct:
                    rows.append((oct_lo, 2 * oct_lo, min(o.real for o in oct), max(o.real for o in oct), min(o.imag for o in oct), max(o.imag for o in oct)))
                oct_lo, oct = 2 * oct_lo, []
            oct.append(S)
            for c in (100, 300, 1000, 3000):
                if abs(g - c) < 1.0 and c not in chk:
                    chk[c] = S
        if oct:
            rows.append((oct_lo, min(T, 2 * oct_lo), min(o.real for o in oct), max(o.real for o in oct), min(o.imag for o in oct), max(o.imag for o in oct)))
        print('\n== x = %g : checkpoints ' % x + '  '.join('T=%d: %+.4f%+.4fi' % (c, s.real, s.imag) for c, s in sorted(chk.items())))
        print('   octave        Re range                 Im range        (widths)')
        for lo, hi, r0, r1, i0, i1 in rows:
            print('   [%5.0f,%5.0f)  [%+.4f, %+.4f]   [%+.4f, %+.4f]   %.4f  %.4f' % (lo, hi, r0, r1, i0, i1, r1 - r0, i1 - i0))


if __name__ == '__main__':
    main()
