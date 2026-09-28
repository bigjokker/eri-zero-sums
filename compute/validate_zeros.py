"""Validate a list of ordinates of zeros of L(s, chi_{-4}) on Re s = 1/2.

    python compute/validate_zeros.py LIST [REF ...]

Checks:
  1. Riemann-von Mangoldt count (T/2pi) log(qT/(2 pi e)), q = 4, at the top of the list.
  2. theta(T)/pi count check at every 1000 (and at the top), theta(t) = (t/2) log(4/pi)
     + Im log Gamma((3/2+it)/2), constant calibrated by the 25 LMFDB-verified zeros below 60.
     A missed pair shows as a step of -2; S(T) stays within about +-1.5.
  3. Against each REF list: on the overlap of ranges [max(bottom), min(top)], zeros unmatched in
     either direction (tolerance 1e-7) and the maximum difference over matched pairs.  A REF may be
     a partial list (e.g. a window scan still running); only its covered range is compared.
Exit status 1 if any 1000-window deviation exceeds 1.5 in magnitude or any zero is unmatched.
"""

import sys, math, bisect
import mpmath as mp

mp.mp.dps = 20


def theta(t):
    t = mp.mpf(t)
    return (t / 2) * mp.log(4 / mp.pi) + mp.loggamma((mp.mpf(3) / 2 + 1j * t) / 2).imag


def read(f):
    values = []
    with open(f, encoding='utf-8-sig') as fh:
        for number, line in enumerate(fh, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            try:
                values.append(float(line))
            except ValueError as exc:
                raise ValueError('%s:%d: invalid ordinate %r' % (f, number, line)) from exc
    if not values or any(not math.isfinite(g) or g <= 0 for g in values):
        raise ValueError('%s: expected nonempty, finite, positive ordinates' % f)
    if any(a >= b for a, b in zip(values, values[1:])):
        raise ValueError('%s: ordinates must increase strictly' % f)
    return values


def unmatched(a, b, tol=1e-7):
    out = []
    for g in a:
        i = bisect.bisect_left(b, g - tol)
        if not (i < len(b) and abs(b[i] - g) <= tol):
            out.append(g)
    return out


def main():
    if len(sys.argv) < 2:
        print('usage: python compute/validate_zeros.py LIST [REF ...]', file=sys.stderr)
        sys.exit(2)
    try:
        z = read(sys.argv[1])
    except (OSError, ValueError) as exc:
        print('Input check failed: %s\nINVALID' % exc)
        sys.exit(1)
    T = z[-1]
    ok = True
    print('%s: %d zeros, %.6f .. %.6f' % (sys.argv[1], len(z), z[0], T))
    rvm = (T / (2 * math.pi)) * math.log(4 * T / (2 * math.pi * math.e))
    print('Riemann-von Mangoldt (T/2pi) log(4T/2pi e) = %.2f  (difference %+.2f, O(log T) = %.1f)' % (rvm, len(z) - rvm, math.log(T)))
    c = 25 - theta(60) / mp.pi
    for t in list(range(1000, int(T) + 1, 1000)) + [T]:
        n = sum(1 for g in z if g <= t)
        dev = n - float(theta(t) / mp.pi + c)
        flag = '' if abs(dev) <= 1.5 else '   <-- FAIL'
        if flag:
            ok = False
        print('  theta/pi check T=%9.3f  zeros %6d  dev %+5.2f%s' % (t, n, dev, flag))
    for ref in sys.argv[2:]:
        try:
            r = read(ref)
        except (OSError, ValueError) as exc:
            print('  vs %s: incomplete check (%s)' % (ref, exc))
            ok = False
            continue
        if r[-1] < z[0] or z[-1] < r[0]:
            print('  vs %s: incomplete check (no overlapping ordinates)' % ref)
            ok = False
            continue
        lo = max(z[0], r[0]) - 1e-6
        hi = min(T, r[-1]) + 1e-6
        za = [g for g in z if lo <= g <= hi]
        ra = [g for g in r if lo <= g <= hi]
        um, ur = unmatched(za, ra), unmatched(ra, za)
        md = max(abs(a - b) for a, b in zip(za, ra)) if za and len(za) == len(ra) else float('nan')
        print('  vs %s on [%.3f, %.3f]: %d vs %d zeros; unmatched %d / %d; max|diff| %.2e' % (ref, lo, hi, len(za), len(ra), len(um), len(ur), md))
        if um or ur:
            ok = False
            print('     in LIST not REF:', ['%.6f' % g for g in um][:10])
            print('     in REF not LIST:', ['%.6f' % g for g in ur][:10])
    print('VALID' if ok else 'INVALID')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
