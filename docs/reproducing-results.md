# Reproducing the numerical results

Run commands from the repository root. Python 3.11.9 with the pinned packages
was used for release checks. The core checks need only mpmath and NumPy;
`analyze.py` additionally needs SciPy and Matplotlib. PARI/GP 2.17.4 generated
the canonical character-zero lists.

## Fast checks using included data

```shell
python -m pip install -r requirements-core.txt
python compute/quick_validate.py
python compute/check_regeneration.py
python compute/validate_zeros.py data/chi4_zeros_pari10000_clean.txt data/chi4_zeros_mp_snapshot.txt data/chi4_zeros_mp_6k8k.txt data/chi4_zeros_to60.txt
```

The validator's success marker is `VALID`; invalid data or requested missing
references produce `INVALID` and exit status 1. These are numerical checks,
not a rigorous certificate that every zero has been found.

To evaluate the small included zero list through the staged pipeline:

```shell
python compute/regenerate_chi4.py --reuse-zeros data/chi4_zeros_to60.txt --reference data/chi4_zeros_mp_snapshot.txt
```

Each run creates a new directory under `data/regenerated/`, rejects an existing
output directory and writes `COMPLETE.json` only after validation and a new
theorem table succeed. A failed run writes `FAILED.json`. Canonical inputs
are not replaced. These local staging directories are ignored by Git.

## Figures and Dirichlet comparison

```shell
python -m pip install -r requirements.txt
python compute/analyze.py
python compute/paperB_q4.py
```

`analyze.py` regenerates the PNG files under `data/figures/` from the included
10,142-row Gram table and 60-residue CSV. `paperB_q4.py` verifies or rebuilds
the full chi-minus-4 ERi cache and writes `data/paperB_q4.txt`. A cache is
reused only when the zero list, evaluator, generator, NumPy version, row count
and CSV checksum agree with its metadata. A rebuild takes a few seconds on
the reviewed machine; plotting time varies by machine.

The plotted model and Bessel product use finite residues. Finite-window
histograms are illustrative comparisons, not proofs of the limiting distribution.

## Generate new chi-minus-4 zero lists

Install PARI/GP separately, then supply its executable or place `gp` on PATH:

```shell
python compute/regenerate_chi4.py --gp gp --tmax 60 --reference data/chi4_zeros_mp_snapshot.txt
```

Windows example:

```powershell
python compute/regenerate_chi4.py --gp "C:\Program Files (x86)\Pari64-2-17-4\gp.exe" --tmax 60 --reference data/chi4_zeros_mp_snapshot.txt
```

The canonical height-10,000 settings are `--tmax 10000 --precision 38
--subdivision 32`. The generation script requests a 4 GB PARI stack, so it
requires sufficient memory. Full generation is much more expensive than
the small smoke check. Review its manifest and diagnostics before adopting
a new list.

The independent mpmath sign-change finder can be run on a small interval:

```shell
python compute/chi4_zeros.py --tmax 60 --out data/local-checks/chi4-to60.txt --validate
```

It checks the file named by `--out` against the included canonical PARI list.
Sign scans may miss close pairs or even-multiplicity zeros. For general
character lists, the PARI recipes are in `compute/gp/`; they preserve the
parameters of the archived computations. Redirect their output to a new
local file rather than replacing an archived list. The quadratic-character
files include constants in comment headers used by `bracket_D.py`.

```shell
python compute/bracket_D.py data/chi_D-8_zeros_pari3000.txt data/chi_D5_zeros_pari3000.txt
python compute/eri_series_chi5.py data/chi5_c2_zeros_pari3000.txt data/chi5_c3_zeros_pari3000.txt
python compute/trivial_integrals.py
```

The first two commands print diagnostics; the final command also replaces
`data/trivial_integrals.txt`. The complex-character spectra must be kept
separate; their positive zero counts differ by one in the archived run.

## Residues

Recompute 400 residues into a separate file:

```shell
python compute/residues.py --count 400 --dps 40 --output data/local-checks/residues-400.npy
```

The NPY columns are ordinate, real residue and imaginary residue. The default
60-residue calculation writes `data/residues_60.csv`; use `--output` to keep
the archived file. Hundreds of mpmath zeta zeros and derivatives take
substantially longer than the quick checks.

## Larger zeta experiments: external input required

Obtain the [first 100,000 zeta ordinates from Andrew Odlyzko](https://www-users.cse.umn.edu/~odlyzko/zeta_tables/)
and save the linked [text table](https://www-users.cse.umn.edu/~odlyzko/zeta_tables/zeros1)
as `data/zeros1.txt`. The source describes its accuracy as within `3 × 10⁻⁹`.
This third-party input is excluded from the repository and ignored by Git.
Its historical SHA-256 is recorded in the provenance guide.

```shell
python compute/run_k2.py
python compute/x_sweep.py
```

These scripts replace their matching pair arrays and sweep JSON in `data/`.
The historical 100,000-zero fast evaluation at x=1.1 took about one minute;
the three-point sweep takes longer. Its real-axis centre uses the Gram
evaluator and its outputs distinguish finite normalization and incomplete
octaves. Fitted transient coefficients remain numerical estimates.

```shell
python compute/run_eri_list.py
```

This regenerates `data/eri_list.csv` using the first 10,142 external ordinates
and the independent 460-digit Gram series. It is substantially slower than
the fast evaluator and was not rerun as part of packaging this release.
