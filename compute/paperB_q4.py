"""Paper B's Theorem at q = 4, numerically: the ERi series over the zeros of L(s, chi_{-4}) at x = 1.1.

    python compute/paperB_q4.py

Paper B says Sigma R_T^chi(x) - ell_chi = (1/(pi log x)) P_A(T log x) + O(log T) with the SAME P_A
as the zeta paper up to C = (2 pi / q) log x in place of 2 pi log x.  So the zeta paper's residue
model (compute/analyze.py, residues_60.csv: the sixty r_rho of the zeta paper's section 2) should
reproduce the chi-sum with log(Y/C) shifted by log q, and the normalised F = (S - ell)/(sqrt T log(T/2pi))
should show the zeta paper's per-octave envelope (0.014-0.026 in its Appendix B).

Outputs data/eri_list_chi4.csv (pair sums per zero) and data/paperB_q4.txt (the tables).
The centre ell_chi is not known in closed form (Paper B: C_{a,x,chi} is an absorbed constant), so it
is fitted once by least squares against the residue model. The same fit on the zeta data is
compared numerically with L = R(x) - pi_0(x) - I(x) = -8.9418; equality of this candidate with
the contour centre is unproved. These fits do not derive either centre.
"""

import csv, math, os, sys, json, hashlib
from pathlib import Path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from mobius_ei import eri
from validate_zeros import read

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, 'data')
X = 1.1
LOGX = math.log(X)
PI = math.pi
L_ZETA = -8.941836473520125          # zeta paper: R(x) - pi_0(x) - I(x) at x = 1.1
EDGES = [64, 128, 256, 512, 1024, 2048, 4096, 8192, 10000]


def zeros_chi4(data=DATA):
    return np.array(read(Path(data) / 'chi4_zeros_pari10000_clean.txt'))


def pair_sums_chi4(data=DATA):
    data = Path(data)
    out = data / 'eri_list_chi4.csv'
    meta = out.with_suffix('.meta.json')
    digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    g = zeros_chi4(data)
    expected = dict(schema=1, x=X, branch='Elliott upper half-plane; pair = 2 Re ERi',
                    source_sha256=digest(data / 'chi4_zeros_pari10000_clean.txt'),
                    evaluator_sha256=digest(Path(__file__).with_name('mobius_ei.py')),
                    generator_sha256=digest(Path(__file__)),
                    rows=len(g), numpy=np.__version__)
    if out.exists() and meta.exists():
        try:
            info = json.loads(meta.read_text(encoding='utf-8'))
            if isinstance(info, dict) and all(info.get(k) == v for k, v in expected.items()) and info.get('csv_sha256') == digest(out):
                with out.open(newline='', encoding='utf-8') as fh:
                    rows = list(csv.DictReader(fh))
                cached_g = np.array([float(r['gamma']) for r in rows])
                pair = np.array([float(r['pair']) for r in rows])
                if len(cached_g) == len(g) and np.isfinite(pair).all() and np.allclose(cached_g, g, atol=2e-12, rtol=0):
                    return g, pair
        except (OSError, ValueError, KeyError, TypeError):
            pass  # Missing, changed or malformed provenance requires recomputation.
    print('Recomputing chi_-4 pair cache: absent or changed provenance.', flush=True)
    vals = np.array([eri((0.5 + 1j * gg) * LOGX) for gg in g])
    pair = 2 * vals.real
    if not np.isfinite(vals).all():
        raise ValueError('nonfinite ERi values; cache not written')
    pending = out.with_suffix('.csv.tmp')
    with pending.open('w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh); w.writerow(['n', 'gamma', 'eri_re', 'eri_im', 'pair'])
        for i, (gg, v, p) in enumerate(zip(g, vals, pair)):
            w.writerow([i + 1, '%.12f' % gg, repr(float(v.real)), repr(float(v.imag)), repr(float(p))])
    pending.replace(out)
    expected['csv_sha256'] = digest(out)
    pending_meta = meta.with_suffix('.json.tmp')
    pending_meta.write_text(json.dumps(expected, indent=2) + '\n', encoding='utf-8')
    pending_meta.replace(meta)
    return g, pair


def pair_sums_zeta():
    rows = list(csv.DictReader(open(os.path.join(DATA, 'eri_list.csv'))))
    return np.array([float(r['gamma']) for r in rows]), np.array([float(r['pair']) for r in rows])


def residues():
    rows = list(csv.DictReader(open(os.path.join(DATA, 'residues_60.csv'))))
    g = np.array([float(r['gamma']) for r in rows])
    r = np.array([float(r['r_gamma_re']) + 1j * float(r['r_gamma_im']) for r in rows])
    return g, r


def model0(T, C, g_r, r_r):
    """(1/(pi log x)) 2 Re sum_rho r_rho (Y^rho/rho)(log(Y/C) - 1/rho), Y = T log x; no centre."""
    Y = T * LOGX
    rho = 0.5 + 1j * g_r
    logYC = np.log(Y / C)
    acc = np.zeros_like(T, dtype=complex)
    for rh, rr in zip(rho, r_r):
        acc += rr * (Y[:, None] ** rh).ravel() / rh * (logYC - 1 / rh)
    return (1 / (PI * LOGX)) * 2 * acc.real


def table(name, T, S, centre, out, q=1):
    """F = (S - c)/(sqrt T log(T/2pi)) as in analyze.py, and F_q with log(qT/2pi) = log(Y/C), the
    logarithm P_A actually carries for the modulus q."""
    F = (S - centre) / (np.sqrt(T) * np.log(T / (2 * PI)))
    Fq = (S - centre) / (np.sqrt(T) * np.log(q * T / (2 * PI)))
    out.append('%s: centre %.4f' % (name, centre))
    out.append('  %12s %12s %10s %12s %10s' % ('T range', 'max|S-c|', 'max|F|', 'max|F_q|', 'max/sqrtT'))
    for a, b in zip(EDGES, EDGES[1:]):
        m = (T >= a) & (T < b)
        if not m.any():
            continue
        out.append('  %5d-%-6d %12.3f %10.4f %12.4f %10.4f' % (a, b, np.max(np.abs(S[m] - centre)), np.max(np.abs(F[m])), np.max(np.abs(Fq[m])), np.max(np.abs(S[m] - centre) / np.sqrt(T[m]))))
    return Fq


def main():
    out = ['Numerical comparison using the first 60 zeta residues; centres fitted unless stated.',
           'Finite model amplitudes are not the full conditional envelope. Final [8192,10000) range is incomplete.']
    g_r, r_r = residues()
    C_ZETA = 2 * PI * LOGX
    C_CHI = (2 * PI / 4) * LOGX
    # zeta: the method check
    Tz, pz = pair_sums_zeta(); Sz = np.cumsum(pz)
    mz = Tz >= 64
    for C, lab in ((C_ZETA, 'C = 2 pi log x (its own)'), (C_CHI, 'C = (pi/2) log x (wrong for zeta)')):
        M = model0(Tz[mz], C, g_r, r_r)
        ell = np.mean(Sz[mz] - M)
        rms = np.sqrt(np.mean((Sz[mz] - M - ell) ** 2))
        out.append('zeta sum, %d zeros to %.1f, %s: fitted centre %.4f (paper: L = %.4f), rms model residual %.3f' % (len(Tz), Tz[-1], lab, ell, L_ZETA, rms))
    # chi_{-4}
    Tc, pc = pair_sums_chi4(); Sc = np.cumsum(pc)
    mc = Tc >= 64
    fits = {}
    for C, lab in ((C_CHI, 'C = (pi/2) log x (Paper B)'), (C_ZETA, 'C = 2 pi log x (zeta value)')):
        M = model0(Tc[mc], C, g_r, r_r)
        ell = np.mean(Sc[mc] - M)
        resid = Sc[mc] - M - ell
        rms = np.sqrt(np.mean(resid ** 2))
        fits[lab] = (ell, rms)
        out.append('chi_-4 sum, %d zeros to %.1f, %s: fitted centre %.4f, rms model residual %.3f, max |resid| %.3f' % (len(Tc), Tc[-1], lab, ell, rms, np.max(np.abs(resid))))
    ell_chi = fits['C = (pi/2) log x (Paper B)'][0]
    out.append('')
    out.append('model vs data at checkpoints (chi_-4, Paper B C), relative error of the oscillation:')
    Mfull = model0(Tc, C_CHI, g_r, r_r)
    for cp in (1238, 2188, 3068, 4722, 5514, 7051, 7802, 8543, 9274, 9998):
        i = np.searchsorted(Tc, cp)
        if i >= len(Tc): continue
        d, m = Sc[i] - ell_chi, Mfull[i]
        out.append('  T=%7.1f  S-ell=%9.3f  model=%9.3f  diff=%7.3f  rel=%6.3f' % (Tc[i], d, m, d - m, abs(d - m) / max(abs(d), 1e-9)))
    out.append('')
    Fz = table('zeta (centre = L of the zeta paper)', Tz, Sz, L_ZETA, out)
    out.append('')
    Fc = table('chi_-4 (centre fitted, Paper B C); F_q uses log(4T/2pi)', Tc, Sc, ell_chi, out, q=4)
    out.append('')
    out.append('the two normalised functions on the common range T in [64, %.0f]:' % min(Tz[-1], Tc[-1]))
    grid = np.geomspace(64, min(Tz[-1], Tc[-1]), 4000)
    fz = np.interp(grid, Tz, Sz) - L_ZETA
    fc = np.interp(grid, Tc, Sc) - ell_chi
    norm = np.sqrt(grid) * np.log(grid / (2 * PI))
    normq = np.sqrt(grid) * np.log(4 * grid / (2 * PI))
    out.append('  rms F_zeta = %.4f, rms F_chi = %.4f (log(T/2pi)) / %.4f (log(4T/2pi)); correlation(F_zeta, F_chi) = %.3f; rms(F_chi[log 4T/2pi] - F_zeta) = %.4f' % (np.sqrt(np.mean((fz / norm) ** 2)), np.sqrt(np.mean((fc / norm) ** 2)), np.sqrt(np.mean((fc / normq) ** 2)), np.corrcoef(fz / norm, fc / normq)[0, 1], np.sqrt(np.mean((fc / normq - fz / norm) ** 2))))
    out.append('  per-octave max |(S_chi - ell) - (S_zeta - L)| / sqrt T   (Paper B: the difference is O(sqrt T), from log q times the integral of f_A):')
    for a, b in zip(EDGES, EDGES[1:]):
        m = (grid >= a) & (grid < b)
        if m.any():
            out.append('    %5d-%-6d %8.4f' % (a, b, np.max(np.abs(fc[m] - fz[m]) / np.sqrt(grid[m]))))
    txt = '\n'.join(out)
    print(txt)
    open(os.path.join(DATA, 'paperB_q4.txt'), 'w').write(txt + '\n')


if __name__ == '__main__':
    main()
