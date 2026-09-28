"""The lower-limit bracket of Paper A's (dagger) for the quadratic character chi_D, on a PARI zero list.

    python compute/bracket_D.py ZEROS_FILE [ZEROS_FILE ...]

Each file is the output of a quadratic-character recipe in compute/gp/: an optional header line
    # D=.. L0=.. L1=.. Lp0=.. h=.. w=.. reg=.. chi2=..
(PARI's values) followed by one ordinate per line.  Computes

    B_T := Pi*(2,chi) + 2 Re sum_{0<gamma<=T} li(2^rho) - I(2,chi)

with Pi*(2,chi) = chi(2)/2 (the star at the jump; 0 when chi(2) = 0), li(y^rho) := Ei(rho log y) in
Elliott's branch, and I(2,chi) the trivial part of Paper A's I(y,chi) at y = 2:
    odd  chi (D<0):  I(2) =              sum_{k>=0} E1((2k+1) log 2)
    even chi (D>0):  I(2) = -log log 2 + sum_{k>=1} E1( 2k    log 2)
Paper A (odd chi): B_T -> log L(0,chi) = log(2h/w).  Even chi: the expected limit is log L'(0,chi) minus Euler's constant; the
residual B_T - log L'(0,chi) is printed as a numerical diagnostic.
Also prints the Riemann-von Mangoldt count and a theta(T)/pi count check with a self-calibrated
constant: a missed pair shows as a persistent step of -2.
"""

import sys, math
import mpmath as mp

mp.mp.dps = 20
Ei = lambda w: -mp.e1(-w)  # Elliott's branch for the positive ordinates used below.
LOG2 = mp.log(2)


def theta(t, q, a):
    t = mp.mpf(t)
    return (t / 2) * mp.log(q / mp.pi) + mp.loggamma((mp.mpf(1) / 2 + a + 1j * t) / 2).imag


def load(f):
    hdr = {}
    g = []
    for l in open(f, encoding='utf-8', errors='replace'):
        l = l.strip()
        if not l:
            continue
        if l.startswith('#'):
            for kv in l[1:].split():
                if '=' in kv:
                    k, v = kv.split('=', 1)
                    hdr[k] = v
        elif l[0].isdigit():
            g.append(mp.mpf(l))
    return hdr, g


def count_check(g, q, a):
    T = float(g[-1])
    rvm = (T / (2 * math.pi)) * math.log(q * T / (2 * math.pi * math.e))
    devs = []
    ts = list(range(500, int(T) + 1, 500)) + [T]
    for t in ts:
        n = sum(1 for x in g if x <= t)
        devs.append((t, n, n - float(theta(t, q, a) / mp.pi)))
    med = sorted(d for _, _, d in devs)[len(devs) // 2]
    c = round(2 * med) / 2
    print('  count: %d zeros to %.3f; RvM (T/2pi)log(qT/2pi e) = %.1f (diff %+.1f); theta/pi + %.1f: '
          % (len(g), T, rvm, len(g) - rvm, c) + ' '.join('%d:%+.2f' % (t, d - c) for t, _, d in devs))
    bad = [t for t, _, d in devs if abs(d - c) > 1.5]
    if bad:
        print('  *** COUNT CHECK FAILS at T =', bad)
    return not bad


def run(f):
    hdr, g = load(f)
    D = int(hdr.get('D', 0)) or int(sys.argv[sys.argv.index('--D') + 1])
    odd = D < 0
    q = abs(D)
    if odd:
        I2 = mp.nsum(lambda k: mp.e1((2 * k + 1) * LOG2), [0, mp.inf])
    else:
        I2 = -mp.log(LOG2) + mp.nsum(lambda k: mp.e1(2 * k * LOG2), [1, mp.inf])
    chi2 = int(hdr['chi2']) if 'chi2' in hdr else None
    Pi2 = mp.mpf(chi2) / 2 if chi2 is not None else mp.mpf(0)
    print('%s: D=%d (%s chi, q=%d)  I(2,chi)=%s  Pi*(2,chi)=%s' % (f, D, 'odd' if odd else 'even', q, mp.nstr(I2, 10), mp.nstr(Pi2, 4)))
    ok = count_check(g, q, 1 if odd else 0)
    target = None
    if hdr:
        L0, L1, Lp0 = (mp.mpf(hdr[k]) for k in ('L0', 'L1', 'Lp0'))
        h, w, reg = int(hdr['h']), int(hdr['w']), mp.mpf(hdr['reg'])
        print('  PARI: L(0)=%s  L(1)=%s  L\'(0)=%s  h=%d  w=%d  reg=%s  chi(2)=%d' % (mp.nstr(L0, 10), mp.nstr(L1, 10), mp.nstr(Lp0, 10), h, w, mp.nstr(reg, 10), chi2))
        if odd:
            target = mp.log(L0)
            print('  class number formula: L(0,chi) - 2h/w = %s;  prediction B_T -> log L(0,chi) = %s' % (mp.nstr(L0 - mp.mpf(2 * h) / w, 3), mp.nstr(target, 8)))
        else:
            target = mp.log(Lp0)
            print('  class number formula: L\'(0,chi) - h log eps = %s;  reference log L\'(0,chi) = %s (no prediction: residual printed)' % (mp.nstr(Lp0 - h * reg, 3), mp.nstr(target, 8)))
    B = mp.mpf(0)
    rows = {}
    for x in g:
        B += 2 * Ei((mp.mpf('0.5') + 1j * x) * LOG2).real
        for T in (500, 1000, 1500, 2000, 2500, 3000):
            if abs(x - T) < 1.0 and T not in rows:
                rows[T] = Pi2 + B - I2
    rows['end'] = Pi2 + B - I2
    line = '  B_T:  ' + '  '.join('T=%s: %+.4f' % (T, float(v)) for T, v in rows.items())
    if target is not None:
        line += '\n  B_T - %s:  ' % ('log L(0)' if odd else 'log L\'(0)') + '  '.join('T=%s: %+.4f' % (T, float(v - target)) for T, v in rows.items())
    print(line)
    if chi2:
        print('  sum-minus-integral alone (no Pi*(2)): %+.4f at the end; Paper A section 2.5 predicts log L(0) - chi(2)/2 = %s for odd chi' % (float(B - I2), mp.nstr(target - Pi2, 6) if (odd and target is not None) else 'n/a'))
    print()
    return ok


if __name__ == '__main__':
    files = [a for a in sys.argv[1:] if not a.startswith('--') and not a.lstrip('-').isdigit()]
    allok = True
    for f in files:
        allok &= run(f)
    sys.exit(0 if allok else 1)
