# Research status

Author: **Open**. Release: **0.1.0**, 27 September 2026.
Zenodo concept DOI: [10.5281/zenodo.23004250](https://doi.org/10.5281/zenodo.23004250).
Version 0.1.0: [10.5281/zenodo.23004251](https://doi.org/10.5281/zenodo.23004251).

This repository is ready to share as a research draft with a numerical companion.
It is not a journal acceptance, independent theorem certification or claim that
priority has been established by an exhaustive literature search. The Zenodo DOI
archives this release. It does not add a peer review.

## What the manuscripts claim

The zeta manuscript develops unconditional two-sided oscillation at scale
`sqrt(T) log(T)` and failure of every fixed nonnegative real Riesz order.
It gives a Bohr expansion under the Weak Mertens Conjecture and an exact
envelope and Bessel-product distribution with an additional linear-independence
hypothesis on the positive zeta ordinates.

The real-character manuscript uses finite Möbius inversion to obtain a formula
for `π*(x,χ) + π_-*(√x)`, with half weights at jumps, an empty-cutoff case and
the parity-dependent constant. It is a formula for that character combination;
an individual residue-class prime count needs additional character formulas.

The Dirichlet manuscript applies the same untwisted ERi function to Dirichlet
zeros. Its leading oscillation is driven by zeta residues. For complex
characters it treats the real and imaginary parts separately. The conditional
envelope hypotheses concern zeta; simplicity of the Dirichlet zeros is not assumed.

## Review and numerical checks

Seven model-assisted audit reports were integrated into the working manuscripts.
Direct source checks covered the Ei convention, averaged endpoints,
parity constants and the local zero-count and height-choice bounds used in the
arguments. The package uses the cleaned manuscripts, with displayed formula
contents retained and GitHub-compatible display delimiters.

The fast Ei evaluator uses Elliott's upper-half-plane value `Ei(w) = -E1(-w)`.
Changing both the Möbius head and tail from the alternative Ei values leaves
the entire ERi sum unchanged. The fast ERi helper is deliberately limited to
finite upper-half-plane inputs; `eri_gram` supplies the independent entire-series
evaluation, including real arguments.

Checks cover endpoint centering, evaluator switches, preservation of mpmath
precision, finite residue amplitudes, malformed zero lists, reference overlap,
staged subprocess failures, rejection of old output directories and cache
invalidation. PARI/GP generation has also been exercised through height 60.
The larger archived searches have not all been rerun after review.

## Remaining qualifications

- The numerical zeta centre `R(x) − π₀(x) − I(x)` has not been identified with
  the absorbed contour constant by a proof.
- The 60- and 400-residue amplitudes are finite partial sums. They do not
  give rigorous upper bounds for the infinite conditional envelope.
- The fitted Dirichlet centre and transient coefficients are observations,
  not derived asymptotics. An incomplete last octave cannot be read as a full
  octave maximum.
- Zero-count, sign-change and reference-overlap checks are numerical diagnostics;
  they do not certify completeness or establish GRH for an infinite spectrum.
- A complete MathSciNet or zbMATH subject search has not been performed.
  The comparison with existing work in the manuscripts does not establish
  novelty solely from an absence of search hits.
- Montgomery–Vaughan Theorem 10.17 was not read in full. The Dirichlet manuscript
  cites the narrower local-count bound recorded in the verified proof of
  Lemma 12.6 and the height choice in Lemma 12.7, rather than importing a broader
  claim from the unread theorem.

Final submission to a journal or arXiv would also require venue-specific
typesetting and publication metadata. Those do not prevent sharing these
Markdown research drafts on GitHub.
