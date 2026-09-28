"""Read-only representative checks for the archived ERi computations."""

from __future__ import annotations

import math
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from mpmath import mp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "compute"))
from eri import L_regularized, averaged_prime_count, eri_gram  # noqa: E402
from mobius_ei import eri, ei_elliott  # noqa: E402


def main() -> None:
    assert [averaged_prime_count(x) for x in (1.1, 2, 2.5, 3, 3.5)] == [0, .5, 1, 1.5, 2]
    assert abs(L_regularized(3)["L"] + 0.0126441817631335) < 1e-12

    old_dps = mp.dps
    for w in (.05 + 1j, .05 + 11.999j, .05 + 12.001j, .05 + 100j):
        assert abs(complex(ei_elliott(w)[0]) + complex(mp.e1(-w))) < 3e-13, w
    for z in (.05 + 1j, .05 + 11.999j, .05 + 12.001j,
              .05 + 100j, .05 + 159.999j, .05 + 160.001j, .05 + 500j):
        assert abs(complex(eri_gram(z)) - eri(z)) < 3e-13, z
    assert mp.dps == old_dps, "fast evaluator changed global mpmath precision"
    for z in (0, .05 - 100j, complex(float('nan'), 1), complex(1, float('inf'))):
        try:
            eri(z)
        except ValueError:
            pass
        else:
            raise AssertionError('unsupported fast ERi input accepted: %s' % z)

    rows = np.load(ROOT / "data" / "residues_400.npy")
    amps = 2 * np.hypot(rows[:, 1], rows[:, 2]) / np.hypot(.5, rows[:, 0])
    amps /= math.pi * math.sqrt(math.log(1.1))
    assert abs(float(amps[:60].sum()) - .039108408609149985) < 1e-12
    assert abs(float(amps.sum()) - .041718458436398695) < 1e-12

    validator = [sys.executable, str(ROOT / "compute" / "validate_zeros.py"),
                 str(ROOT / "data" / "chi4_zeros_pari10000_clean.txt")]
    refs = [str(ROOT / "data" / name) for name in
            ("chi4_zeros_mp_snapshot.txt", "chi4_zeros_mp_6k8k.txt")]
    passed = subprocess.run(validator + refs, cwd=ROOT, capture_output=True, text=True)
    assert passed.returncode == 0 and passed.stdout.rstrip().endswith("VALID"), passed.stdout
    absent = subprocess.run(validator + [str(ROOT / "data" / "no-such-reference.txt")],
                            cwd=ROOT, capture_output=True, text=True)
    assert absent.returncode != 0 and absent.stdout.rstrip().endswith("INVALID")
    with tempfile.TemporaryDirectory() as directory:
        bad = Path(directory) / 'bad.txt'
        for content in ('10\n6\n', '6\n-1\n10\n', '6\nnan\n10\n', '6\nGP error\n10\n', ''):
            bad.write_text(content, encoding='utf-8')
            checked = subprocess.run([*validator[:2], str(bad)], cwd=ROOT, capture_output=True, text=True)
            assert checked.returncode == 1 and checked.stdout.rstrip().endswith('INVALID'), checked.stdout
            assert 'Traceback' not in checked.stderr, checked.stderr
    print("Quick validation passed: centering, evaluator switches, precision, envelopes, zero lists.")


if __name__ == "__main__":
    main()
