# Notation and hypotheses

This sheet accompanies the three manuscripts. Their numbered definitions and theorem
statements give the full domains and constants.

| Symbol or convention | Meaning |
|---|---|
| `x > 1` | Fixed real argument; `log` is the natural logarithm |
| `ERi(z)` | Entire Gram function `1 + Σ z^k/(k k! ζ(k+1))` |
| `R(x^ρ)` | `ERi(ρ log x)`; the exponent is used before evaluating the function |
| `Ei(w)` | Elliott's branch; `Ei(w)=-E1(-w)` in the upper half-plane |
| `ρ=β+iγ`, `m_ρ` | A non-trivial zero and its multiplicity; the relevant zeta or L-function is specified |
| `Σ R_T` | Zero sum ordered by absolute ordinate; multiplicities are retained |
| `r_ρ` | Residue coefficient from zeta; Paper B's leading oscillation also uses these zeta coefficients |
| `D(a,x)` / `D_{a,x}` | Absorbed contour constant in the zeta transfer |
| `L` in the zeta numerical discussion | Candidate centre `R(x)−π₀(x)−I(x)`; identification with the contour constant is unproved |
| `L(s,χ)` | Dirichlet L-function; distinct from the numerical constant `L` |
| `I(x)` in the zeta manuscript | Regularized trivial-zero value `1/log x−arctan(π/log x)/π` |
| `I(y,χ)` in Paper A | Parity-dependent trivial-zero tail, including `−log log y` for even χ |
| `C(x,χ)` in Paper A | Finite sum `Σ_{n≤N}(μ(n)/n)I(x^{1/n},χ)` |
| `ℓ_χ(x)` | Absorbed centre in Paper B; the numerical q=4 value is fitted |
| `π₀`, `π*`, `Π*` | Average of left/right limits at a jump; an endpoint prime has half weight |
| `χ` in Paper A | Real primitive non-principal character modulo `q>1` |
| `π_-^*(y)` | Averaged count of primes with `χ(p)=-1`; primes dividing q are excluded |
| `N=⌊log x/log 2⌋` | Finite Möbius cutoff; relevant sums are empty when `1<x<2` |
| `R_N(x^ρ)` | `Σ_{n≤N}(μ(n)/n) Ei(ρ log x/n)` |
| `κ(χ)` | `log L(0,χ)` for odd χ; `log L'(0,χ)−γ_E` for even χ |
| `C₀`, `γ_E` | Euler's constant; `γ` without a subscript is a zero ordinate |
| `k≥0` | Fixed real Riesz order; the k=0 multiplier is interpreted separately as 1 |
| `WMC` | Weak Mertens Conjecture for the ordinary Möbius summatory function |
| `LI` | Rational linear independence of the positive zeta ordinates |

The two-sided oscillation and failure of fixed real Riesz orders are unconditional statements
in the audited arguments. The Bohr expansion uses WMC; the Bessel product and exact envelope
also use LI. These are hypotheses on zeta in Paper B, and no simplicity assumption on the
Dirichlet zeros is used.

Numerical sums from 60 or 400 residues are partial amplitudes. They are not the infinite
conditional envelope or rigorous upper bounds for it. The fitted transient and candidate
centres are numerical observations. Zero-count and overlap checks do not establish
completeness of a zero list.
