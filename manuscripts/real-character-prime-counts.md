# A finite Riemann-R formula for real-character prime counts

**Open**

## Abstract

For a real primitive non-principal Dirichlet character, finite Möbius inversion
gives an identity for the averaged character prime count together with a correction at
the square root. Combining this identity with the truncated explicit formula gives a
finite Riemann-R formula, including the parity-dependent constant and a remainder
tending to zero for fixed x. The identity retains half weights at jumps and includes
the empty cutoff case. Replacing the finite Möbius cutoff by the entire ERi function
introduces the divergent zero series studied in the companion paper.

For a real primitive non-principal character, the formula concerns
`π*(x,χ) + π_-^*(√x)`, a combination of residue classes. An individual progression requires
additional character formulas; at modulus 4 the ordinary averaged prime count supplies the
second equation. The analytic companion is [Paper B](dirichlet-zero-sums.md).

---

## 0. Notation

`χ` a **real primitive non-principal** character mod `q` — equivalently the Kronecker symbol
`(D/·)` of a fundamental discriminant `D` with `|D| = q`; `h(D)` is the class number of `ℚ(√D)`,
`w(D)` the number of roots of unity in it, and for `D > 0`, `ε_D > 1` its fundamental unit. At
`s = 0` the class number formula reads `L(0,χ) = 2h(D)/w(D)` for `D < 0` and
`L'(0,χ) = h(D) log ε_D` for `D > 0` (where `w = 2`). `C₀ = 0.5772…` is Euler's constant, in MV's
notation (12.7); `γ` is always the ordinate of a zero, and the `c₀(χ)` of §2.2 is a different
object.

**The star is Riemann's convention**, as in Manzoni: at a jump the value is the mean of the
left- and right-hand limits. Every starred function below is so normalised, and `π_-^*` likewise.
This is not decoration — the explicit formula (†) produces the averaged function, so the left
side of Theorem A must be averaged too. Otherwise the two sides differ by `½χ(p)` at each
prime, where `π*(·,χ)` jumps, and by more at prime powers, where `Π*(·,χ)` does; and
`E_Π → 0` fails on exactly that set. Write

```
S_±      = {a mod q : χ(a) = ±1}
π*(y,χ)  = Σ_{p ≤ y} χ(p),  averaged at jumps
           (= Σ_{a∈S_+}π*(y;q,a) − Σ_{a∈S_-}π*(y;q,a))
π_-^*(y) = #{p ≤ y : χ(p) = −1}, averaged   (= Σ_{a∈S_-} π*(y;q,a))
Π*(x,χ)  = Σ_{k≥1} (1/k) Σ_{p ≤ x^{1/k}} χ(p^k) = Σ_{k≥1} π*(x^{1/k}, χ^k)/k,  averaged
```

Averaging is linear and commutes with the finite sums of §1, so Lemmas 1–2 and Proposition 3
hold verbatim for the starred functions.

`ρ = β+iγ` runs over non-trivial zeros of `L(s,χ)`, and `m_ρ` is the multiplicity of `ρ`. Sums
over `ρ` are written with `m_ρ` explicitly rather than by listing zeros with repetition: **no
simplicity is assumed**, and nothing here needs it — Theorem A is an exact identity in which `m_ρ`
is simply the multiplicity. `R(x^ρ) := ERi(ρ log x)`, and for `M ≥ 1`

```
R_M(x^ρ) := Σ_{n≤M} (μ(n)/n) li(x^{ρ/n})        (so R_M → R as M → ∞)
```

`N := ⌊log x / log 2⌋`, the largest `n` with `x^{1/n} ≥ 2`.

For `1 < x < 2`, `N=0`: `R_0`, `m(0)`, `C(x,χ)` and `E(x,T)` are all empty sums,
hence zero, as are both prime counts in Theorem A.

Throughout, `li(y^ρ) := Ei(ρ log y)`, the exponential integral continued from the negative axis.
The naive principal branch of `li` on `ℂ`, cut along `(−∞,1]`, is **not** what appears here; it is
wrong for complex exponents; the exponent is applied before evaluating Ei.

---

## 1. The obstruction, and why it closes for real `χ`

`Π*` mixes characters: the `k`-th term carries `χ^k`, not `χ`, so Möbius inversion does **not**
give `π*(x,χ) = Σ_n (μ(n)/n) Π*(x^{1/n}, χ)`. That is the reason no `R(x;q,a)` appears in the
monographs.

For **real** `χ` it closes, because `χ^k` takes only two values. For `p ∤ q`,
`χ(p^k) = χ(p)^k = χ(p)` when `k` is odd and `= 1` when `k` is even; and for `p | q` every term
vanishes, so all sums below run over `p ∤ q`.

> **Lemma 1.** For `χ` real primitive non-principal mod `q` and `x > 1`,
> ```
> Π*(x,χ) = Σ_{k≥1} π*(x^{1/k},χ)/k  +  Σ_{j≥1} π_-^*(x^{1/(2j)})/j .
> ```

*Proof.* Split by parity of `k`. Odd `k` contribute `π*(x^{1/k},χ)/k` directly. For even `k`,
`χ^k` is `1` on `p ∤ q`, so `π*(x^{1/k},χ^k) = π^{(q)*}(x^{1/k})`, the count of `p ≤ x^{1/k}` with
`p ∤ q`. Since `π^{(q)*}(y) = π_+^*(y) + π_-^*(y)` and `π*(y,χ) = π_+^*(y) − π_-^*(y)`, we have
`π^{(q)*}(y) − π*(y,χ) = 2π_-^*(y)`. Adding and subtracting `π*(x^{1/k},χ)/k` over even `k` gives
`Σ_{k even} 2π_-^*(x^{1/k})/k = Σ_{j≥1} π_-^*(x^{1/(2j)})/j`. ∎

Write `A(x) := Σ_{k≥1} π*(x^{1/k},χ)/k`, so `Π*(x,χ) = A(x) + Δ(x)` with
`Δ(x) := Σ_{j≥1} π_-^*(x^{1/(2j)})/j`. Both sums are finite: every term with `x^{1/k} < 2` vanishes.

`A` is now a clean `Σ_k (·)/k`, so Möbius inversion applies:

> **Lemma 2.** `π*(x,χ) = Σ_{n≥1} (μ(n)/n) A(x^{1/n})`, a finite sum (terms with `n > N` vanish).

More generally, for any function `f(y)` vanishing on `1 < y < 2`, if
`A_f(x) = Σ_{k≥1} f(x^{1/k})/k`, then
`f(x) = Σ_{n≤N} (μ(n)/n) A_f(x^{1/n})`. Grouping the finite double sum by `kn=m`
gives coefficient `(1/m)Σ_{n|m}μ(n)`, proving this version as well.

The same inversion applied to `Δ` collapses it. Setting `y = √x`, one has
`Δ(x) = Σ_j π_-^*(y^{1/j})/j`, so Lemma 2 for `π_-^*` itself gives
`Σ_n (μ(n)/n) Δ(x^{1/n}) = π_-^*(y)`. The range `n ≤ N` includes every surviving term
(`Δ(x^{1/n})` vanishes for `n > N/2`). Thus:

> **Proposition 3 (the combinatorial identity — exact, finite, elementary).**
> ```
> π*(x,χ) + π_-^*(√x) = Σ_{n≤N} (μ(n)/n) Π*(x^{1/n}, χ) .
> ```

The double sum `Σ_{n,j} (μ(n)/(nj)) π_-^*(x^{1/(2jn)})` is `π_-^*(√x)`: grouping by `m = jn`, the
coefficient of `π_-^*(x^{1/(2m)})` is `(1/m) Σ_{n|m} μ(n)`, which is `1` at `m = 1` and `0`
otherwise. (Checked numerically for `χ₋₄` at `x = 10³, 10⁴, 10⁶, 10⁸, 10¹⁰`: the double sum
equals `π*(√x;4,3)` exactly — 6, 13, 87, 619, 4808.)
Manzoni's correction is this identity, not a truncation of it. For `q = 4`, `S_- = {3}`,
`π_-^*(y) = π*(y;4,3)`, and the left side is `π*(x;4,1) − π*(x;4,3) + π*(√x;4,3)`.

---

## 2. The explicit formula, by parity

`χ` real primitive mod `q` is **odd** or **even** according as `χ(−1) = −1` or `+1`; for
`χ = (D/·)` that is `sign(D)`. Write `𝔞 = 0` (even), `𝔞 = 1` (odd). The completed function is

```
Λ(s,χ) = (q/π)^{(s+𝔞)/2} Γ((s+𝔞)/2) L(s,χ),        Λ(s,χ) = ε(χ) Λ(1−s, χ̄),   |ε(χ)| = 1.
```

The two parities differ in three linked ways, and the split has to be made here rather than
patched later: **`χ₋₄` is odd**, so the odd case is the one to write first; even is the other
family, not the default.

| | odd (`𝔞 = 1`) | even (`𝔞 = 0`) |
|---|---|---|
| Γ-factor | `Γ((s+1)/2)` | `Γ(s/2)` |
| trivial zeros | `s = −1, −3, −5, …` | `s = 0, −2, −4, …` |
| `L(0,χ)` | `= −B_{1,χ} ≠ 0` | `= 0`, a simple zero |
| constant of (†), §2.3 | `log L(0,χ) = log(2h(D)/w(D))` | `log L'(0,χ) − C₀ = log(h(D) log ε_D) − C₀` |
| `s = 0` **in `ψ`** (§2.3: neither survives into `Π*`; the `log L(0,χ)` that does is not a `ψ`-level object) | a **value**: `−(L'/L)(0,χ)` | a **residue**: the trivial zero collides with the `1/s` of Perron |
| example | `χ₋₄`, `χ₋₃`, `χ₋₈` | `χ₅`, `χ₈`, `χ₁₃` |

For even `χ`, `L(0,χ) = 0`, so **`(L'/L)(0,χ)` is not a holomorphic constant and must not be
written as one.** `(−L'/L)(s,χ)·y^s/s` has a *double* pole at `s = 0`, and the residue carries a
`log y`.

### 2.1 Odd `χ` — the `χ₋₄` case

Shifting `(1/2πi)∫(−L'/L)(s,χ)y^s/s\,ds` to the left picks up the non-trivial zeros, each with its
multiplicity `m_ρ`, the simple pole of `1/s` at `s = 0` (where `L(0,χ) ≠ 0`, so `−L'/L` is
regular), and the trivial zeros
`s = −(2m+1)`, each a simple zero of `L` contributing `−y^{s_0}/s_0`:

```
ψ₀(y,χ) = −Σ_ρ m_ρ y^ρ/ρ − (L'/L)(0,χ) + Σ_{m≥0} y^{−(2m+1)}/(2m+1)
        = −Σ_ρ m_ρ y^ρ/ρ − (L'/L)(0,χ) + arctanh(1/y).
```

The left side is the averaged `ψ₀` of §0, as MV state it (§2.5); at a prime power the unaveraged
`ψ` differs by half the jump. No `log y`: `1 − 𝔞 = 0`.

### 2.2 Even `χ`

Now the trivial zero at `s = 0` sits on the Perron pole. With
`L(s,χ) = L'(0,χ)s + ½L''(0,χ)s² + ⋯`, one has `(L'/L)(s,χ) = 1/s + c₀(χ) + O(s)` with
`c₀(χ) = L''(0,χ)/(2L'(0,χ))`, so the residue of `(−L'/L)(s,χ)y^s/s` at `s = 0` is
`−(log y + c₀(χ))`. The remaining trivial zeros are `s = −2m`, `m ≥ 1`:

```
ψ₀(y,χ) = −Σ_ρ m_ρ y^ρ/ρ − log y − c₀(χ) + Σ_{m≥1} y^{−2m}/(2m)
        = −Σ_ρ m_ρ y^ρ/ρ − log y − c₀(χ) − ½log(1 − y^{−2}).
```

Again the formula uses `ψ₀`; the parity term is `−(1−𝔞)log x`.

### 2.3 Passing to `Π*`

`Π*(y,χ) = ∫ dψ₀(t,χ)/log t`, the same partial summation that turns `y^ρ/ρ` into `li(y^ρ)`:
`∫ t^{ρ−1}dt/log t = li(y^ρ)`. Applying it termwise from `2` to `y`, and **evaluating** the
lower-limit bracket this produces (below: `log L(0,χ)` for odd `χ`, `log L'(0,χ) − C₀` for even),
gives (†) with

```
odd  :  I(y,χ) = ∫_y^∞ dt / ((t²−1) log t)                    [ from d/dt arctanh(1/t) = −1/(t²−1) ]
even :  I(y,χ) = −log log y + ∫_y^∞ dt / (t(t²−1) log t)      [ the integral is Riemann's, from ζ ]
```

**The `ψ`-constants do not survive.** `Π* = ∫dψ₀/log t`, and a constant in `ψ₀` has `d(const)=0`;
integrating by parts, `Π*(y) = ψ₀(y)/log y − ψ₀(y₀)/log y₀ + ∫_{y₀}^y ψ₀(t)dt/(t log²t)`, and a
constant `K` contributes `K/log y − K/log y₀ + K[−1/log t]_{y₀}^y = 0` exactly. So
**`−(L'/L)(0,χ)` (odd) and `−c₀(χ)` (even) do not appear in `Π*` at all.** What survives is `I`,
including the `−log log y` that `−log y` becomes. (Checked numerically: a constant contributes
`−3.5·10⁻⁸`, i.e. quadrature error; `−log y ↦ −log log y + log log y₀` to ten digits.)

**The lower-limit constant.** The cancellation shows that no *value of `L'/L`* survives: not `−(L'/L)(0,χ)`, not `−c₀(χ)`, and not MV's `C(χ)` of
(12.7). It does **not** show that (†) has no constant term. Writing `I` as a tail `∫_y^∞` leaves a
lower-limit bracket at `y₀ = 2`, and that bracket is a constant with a value:

```
Π*(2,χ) + Σ_ρ m_ρ li(2^ρ) − I(2,χ)  =  κ(χ)  :=  log L(0,χ)          (odd χ)
                                              =  log L'(0,χ) − C₀     (even χ).
```

For odd `χ` this is the residue of `log L(s,χ)·y^s/s` at `s = 0` in the Mellin form of `Π*`. It is
not new: it is the constant `log L(0)` of Landau's formula (38) — E. Landau, *Nouvelle démonstration
pour la formule de Riemann sur le nombre des nombres premiers inférieurs à une limite donnée et
démonstration d'une formule plus générale pour le cas des nombres premiers d'une progression
arithmétique*, Ann. Sci. École Norm. Sup. (3) **25** (1908) 399–442, Part II §XV, p. 441 — which is
(†) untruncated, for every character of either parity, with `log L(0)` "la valeur obtenue en
prolongeant `log L(s)` suivant l'axe réel en évitant les points singuliers par de petits détours
dans le demi-plan supérieur", and with the trivial zeros as the integral
`∫_x^∞ dy/(y^{1−a} log y (y²−1))`, `a = 𝔞` [La, pp. 441–442].

For even `χ` there is no residue to take: `L(0,χ) = 0`, and `log L(s,χ) = log L'(0,χ) + log s + O(s)`
is a logarithmic branch point of `log L` sitting on Perron's pole. Landau's convention for that
case, p. 442 — *"Si `L(0) = 0`, le résultat (38) subsiste encore si l'on fait la convention
d'entendre par `−λ Li(1) + log L(0)` (cas du zéro `𝔯 = 0` de multiplicité `λ ≥ 1`) la limite
`lim_{h→+0} {−λ Li(x^h) + log L(−h)}`"* — evaluates in three lines, with `λ = 1` since the zero of
a primitive even `χ` at `0` is simple. `Li(x^h) = Ei(h log x) = C₀ + log h + log log x + O(h)`.
`L(−h,χ) = −h L'(0,χ) + O(h²)`, so on Landau's path — the small detour above `s = 0` —
`log L(−h,χ) = log L'(0,χ) + log h + iπ + O(h)`, with `L'(0,χ) = h(D) log ε_D > 0`. The `log h`
cancel:

```
lim_{h→0+} {−Li(x^h) + log L(−h,χ)}  =  log L'(0,χ) − C₀ − log log x + iπ .
```

The `−log log x` is the term `I(y,χ)` already carries — §2.2's `−log y`, integrated; the `iπ` is
absorbed by Landau's `−(α − c)πi`, the zero at `0` being counted among his real zeros `𝔯`; and what
remains is the constant `log L'(0,χ) − C₀`. In Landau's form the trivial zero at `s = 0` does
produce a `li` term, `−Li(x⁰)`, regularised jointly with `log L(0)`; in the `ψ₀`-route above the
same quantity arrives through the double pole. The two bookkeepings agree, as the numbers show.

Thus (†) carries the constant `κ(χ)`, while constants arising from `L'/L`
are removed by partial summation. The class-number expressions for `κ` are stated
in Theorem A and compared numerically in §4.

`I(y,χ) → 0` as `y → ∞` in the odd case. In the even case it does **not**: `−log log y → −∞`.
It is bounded on `2 ≤ y ≤ x` with `x` fixed, which is all Theorem A needs, but it is not `O(1)`
as a statement about `I`.

**Lower limit: `y₀ = 2`, fixed.** The termwise transform cannot be run from `t = 1` — the odd
trivial-zero series `arctanh(1/t)` diverges as `t → 1⁺`, and `∫dt/((t²−1)log t)` blows up there
like `(t−1)^{−2}`. It does not need to be. `Π*(·,χ)` vanishes on `[1,2)`, and on `[2,∞)` both
`arctanh(1/t)` and the even series are regular; Theorem A calls `I(x^{1/n},χ)` only for `n ≤ N`,
so only at `y = x^{1/n} ≥ 2`. **The singularity at `t = 1` lives in the `n > N` region that §3
already refuses to enter** — and Riemann's grouping exists to continue `li(y^ρ)` through `y = 1`,
which is that same region. So integrate `∫dψ₀/log t` from `2` to `y`, take `I` as the tail
`∫_y^∞` plus the even `−log log y`, and **evaluate** the lower-limit bracket that this produces —
it is `κ(χ)`, as displayed above — rather than letting `E_Π` absorb it. `E_Π` is then
`∫_2^y dR/log t` minus the tail `Σ_{|γ|>T} m_ρ li(2^ρ)` of the lower-limit bracket (§2.5),
which does tend to `0`. No free `y₀`, and no unevaluated one.

**Both integrals are exponential-integral series.** On `t > 1` expand the trivial-zero factor
geometrically — `1/(t²−1) = Σ_{k≥0} t^{−2k−2}` and `1/(t(t²−1)) = Σ_{k≥0} t^{−2k−3}` — and
substitute `u = log t`, so that `∫_{log y}^∞ e^{−mu} du/u = E₁(m log y)`. Then

```
odd  :  ∫_y^∞ dt/((t²−1) log t)   =  Σ_{k≥0} E₁((2k+1) log y)
even :  ∫_y^∞ dt/(t(t²−1) log t)  =  Σ_{k≥0} E₁((2k+2) log y)
```

Here `E₁(z) = ∫_z^∞ e^{−u}du/u` is the exponential integral, written throughout with a single
argument. (MV's `E₁(x,χ)` of §2.5 is an unrelated object and always carries two.)

The two families of exponents are the two families of trivial zeros — odd `χ` has them at the
negative odd integers, even `χ` at the negative even integers — and this is not an algebraic
accident. In the branch convention of §0, `E₁(m log y) = −li(y^{−m})`, so each series is

```
Σ_{k≥0} E₁(m_k log y)  =  −Σ_{ρ trivial} li(y^ρ),
```

which is (†)'s nontrivial sum in the same form and with the same sign. The trivial zeros are not a
correction bolted onto (†); they are the rest of it. For odd `χ` that is the whole of `I(y,χ)`. For
even `χ` it is not:

```
I(y,χ)  =  −log log y  −  Σ_{ρ = −2, −4, …} li(y^ρ)                 (even χ)
```

and the sum runs over `−2, −4, …` **only**. By §2's table the trivial zeros of an even `χ` are
`0, −2, −4, …`, so `s = 0` is one of them — but it yields no `li` term. It descends instead through
the **double pole of `(−L'/L)(s,χ)·y^s/s` at `s = 0`**, whose residue at the `ψ` level is `−log y`
(§2.2), and that is the `−log log y`. The zero of `L(s,χ)` there is simple; the pole is double
because that zero collides with Perron's `1/s`. Odd `χ` has `L(0,χ) ≠ 0` and no trivial zero at the
origin, so no collision and no such term — which is why the two cases differ by more than a shift
of exponents.

Since `E₁(m log y) ≍ y^{−m}/(m log y)`, both series converge geometrically for `y ≥ 2` with ratio
`y^{−2}`, so `C(x,χ)` is evaluable to any precision from the definition. Checked against numerical
quadrature at `y = 2, 3, 10, 100` to 30 digits, worst discrepancy `4.93·10⁻³²`
(`compute/trivial_integrals.py`; output `data/trivial_integrals.txt`).

Two things this does **not** give. It is not a named entry in any table we have found — Landau
writes it as the integral, `∫_x^∞ dy/(y^{1−a} log y (y²−1))`, and the series is its closed form.
And there is no odd analogue of the folklore `ζ` formula
`Σ_{k≥1} R(y^{−2k}) = 1/log y − (1/π)arctan(π/log y)` — which is not an identity: the series
diverges ([G] Prop. 4.1) and the right side is its regularised value. The formula is
Möbius-smoothed, i.e. it belongs to the `n`-sum run to infinity, which is precisely what the
paragraph below refuses.

**Do not Möbius-invert `I` into a closed package.** In Theorem A, `C(x,χ)` is the *finite* sum
`Σ_{n≤N}(μ(n)/n) I(x^{1/n},χ)` and stays that way. Running the `n`-sum to infinity and
collapsing it — the `ζ` analogue is `Σ_k R(x^{−2k})`, whose folklore value is
`1/log x − (1/π)arctan(π/log x)` — is the
`n > N` interchange of §3 in another costume, and [G] Prop. 4.1 shows the resulting series
diverges for `ζ`. The arctan package is not available here either.

The parity terms agree with [MV, (12.6)]. Its constant `C(χ)` is removed by
partial summation, while the evaluated lower-limit bracket `κ(χ)` remains.
All real zeros are retained in the indexing `|γ| ≤ T`; see §2.5.

### 2.4 Truncation

For non-principal `χ`, the chain `ψ₀(x,χ) → Π*(x,χ)` — MV Theorem 12.10, then partial summation
(§2.5) — gives, for `T ≥ 2` not an ordinate of a zero of `L(s,χ)`,

```
Π*(y,χ) = −Σ_{|γ|≤T} m_ρ li(y^ρ) + I(y,χ) + κ(χ) + E_Π(y,T)          (†)
```

for `χ` of either parity, with `I(y,χ)` as in §2.3 — the trivial-zero tail, plus `−log log y` when
`χ` is even — and `E_Π(y,T)` the truncation remainder of §2.5 (the integral of `dR/log t` together
with the tail of the lower-limit bracket), `E_Π(y,T) → 0` as `T → ∞` for each fixed `y > 1`. The
constant is the evaluated lower-limit bracket of §2.3,
`κ(χ) = log L(0,χ)` for odd `χ` and `log L'(0,χ) − C₀` for even; there is no `c(χ)` of `L'/L` type,
since by §2.3 the `ψ`-constants die under `∫dψ₀/log t`. Untruncated, (†) is Landau's (38).

Substituting (†) into Proposition 3 and interchanging the **finite** `n`-sum with the
**finite** `γ`-sum:

> **Theorem A.** Let `χ = (D/·)` be real primitive non-principal mod `q = |D| > 1`, let `x > 1` be fixed, put
> `N = ⌊log x/log 2⌋` and `m(N) = Σ_{n≤N} μ(n)/n`, and let `T ≥ 2` not be an ordinate of a zero of
> `L(s,χ)`. Then
> ```
> π*(x,χ) + π_-^*(√x) + Σ_{|γ|≤T} m_ρ R_N(x^ρ) = C(x,χ) + κ(χ)·m(N) + E(x,T),
> ```
> where `C(x,χ) = Σ_{n≤N}(μ(n)/n) I(x^{1/n},χ)`, `E(x,T) = Σ_{n≤N}(μ(n)/n) E_Π(x^{1/n},T) → 0` as
> `T → ∞`, and
> ```
> κ(χ) = log L(0,χ)       = log(2h(D)/w(D))              (D < 0, χ odd),
> κ(χ) = log L'(0,χ) − C₀ = log(h(D) log ε_D) − C₀        (D > 0, χ even),
> ```
> `C₀` Euler's constant.

Everything here converges. Both interchanges are between finite sums. The hypothesis `T ≥ 2` is
MV 12.10's, and it is not vacuous: for large `q` the lowest zero of `L(s,χ)` sits below height `2`.

The accepted answer to MathOverflow question 386213 (2021; [source discussion](https://mathoverflow.net/questions/386213))
is the `ζ`-twin of this statement: the Möbius sum cut at `m ≤ X` with `X > log₂ x`, the constant
`−log 2·Σ_{m≤X} μ(m)/m`, and no remainder, the zero-sum being left untruncated. At `X = N` the
constant is `log|ζ(0)|·m(N)`; for larger `X` the extra terms are `(μ(m)/m)Π*(x^{1/m}) = 0` each, their
share of the constant absorbed by their zero-sums. There it is followed by the passage to the infinite
form that §3 refuses. For `L(s,χ)` the untruncated identity with its constant is Landau's (38) of
1908 (§2.3); Theorem A is its truncation in `T` and in `n`, with the remainder made explicit.

The zero-sum is over `|γ| ≤ T`, matching MV (12.6). Real zeros are therefore **inside** it, and no
hypothesis is imposed. This is not a convenience: if `χ` is real and `β ∈ (0,1)` is a zero of
`L(s,χ)`, so is `1−β`, and both have `γ = 0`. Writing `Σ_{|γ|≤T}` carries the pair without naming
either.

The alternative packagings are both worse here. Granville's (11.8.4) pulls the pair out
explicitly as `−x^β/β − (x^{1−β} − 1)/(1−β)`; that is honest but names a pair the identity has no
need to name. And it must be the *pair*: at the `ψ` level the constant `1/(1−β)` can be absorbed
into `c_χ`, but §2.3 kills constants in `ψ`, so what reaches `Π*` is `R_N(x^β)` and `R_N(x^{1−β})` with
nothing cancelling either. Dropping `(x^{1−β} − 1)/(1−β)` into an error term is Granville's
(11.8.5), a coarsening calibrated to `x → ∞` at PNT scale. The objection is not that the discarded
term is large — it is not: with `β` the zero near `1` it is `∼ log x` for a Siegel zero and
`x^{1−β}`-sized otherwise, against the `x^β/β` that (11.8.5) keeps. The objection is that it does
not vanish. Its constant part dies in §2.3, and what would then sit inside `E(x,T)` is
`R_N(x^{1−β})`, which is independent of `T`; `E(x,T)` is required to tend to `0` as `T → ∞`, and
with that term present it does not.

---

### 2.5 The remainder, and where `q` enters

**The input.** Montgomery–Vaughan, *Multiplicative Number Theory I*, Theorem 12.10, book pp.
403–405. Let `c > 1`, `x ≥ c`, `T ≥ 2`, and let `χ` be primitive mod `q`, `q > 1`.
Then

```
ψ₀(x,χ) = − Σ_{ρ, |γ|≤T} x^ρ/ρ − ½log(x−1) − (χ(−1)/2) log(x+1) + C(χ) + R(x,T;χ)     (12.6)

C(χ)     = (L'/L)(1,χ̄) + log(q/2π) − C₀                                               (12.7)

R(x,T;χ) ≪ (log x)·min(1, x/(T⟨x⟩)) + (x/T)(log qxT)²                                 (12.8)
```

where `⟨x⟩` is the distance from `x` to the nearest prime power **other than `x` itself**, `C₀` is
Euler's constant, and the sum over `ρ` lists zeros with repetition — MV's convention, our `Σ m_ρ`.
Four things follow.

**It is stated on `ψ₀`.** MV's left side is the averaged function, `ψ₀(x,χ) = ½(ψ(x+,χ) + ψ(x−,χ))`.
That is §0's convention, so the passage to the averaged `Π*` needs no separate argument.

**The parity split is (12.6).** Setting `χ(−1) = ±1` in the two log terms:

| | (12.6) | §2.1 / §2.2 |
|---|---|---|
| even, `χ(−1) = +1` | `−½log(x−1) − ½log(x+1) = −½log(x²−1)` | `−log x − ½log(1 − x^{−2})` |
| odd, `χ(−1) = −1` | `−½log(x−1) + ½log(x+1) = ½log((x+1)/(x−1))` | `arctanh(1/x)` |

Both rows are equalities. Granville's notes on the same formula give the mechanism: for
`χ(−1) = −1` the trivial zeros sit at the negative odd integers, and for `χ(−1) = +1` there is
— in his phrasing — *"a double zero at `s = 0`, leaving a residue of `log x`"*. The residue is
right; the object is not. `L(s,χ)` has a **simple** zero at `s = 0` for even `χ`, and what is
double is the pole of `(−L'/L)(s,χ)·y^s/s` there, the zero having met Perron's `1/s`. §2.2 is
canonical in this file.

**The `T`-independent `log x` is MV's own coarsening, and MV prints both forms.** Theorem 12.12
gives `E₂(x,T,χ) ≪ log x + (x/T)(log xTq)²`, equation (12.16), and its proof reads *"It is also
clear that (12.8) gives (12.16)"* — that is, (12.16) is (12.8) with `min(·) ≤ 1`. So the lumped
shape is a deliberate weakening of the sharp one, and it is the sharp one Paper A needs: at fixed
`y` the quantity `⟨y⟩` is positive — the distance is to the nearest prime power *other than `y`
itself*, so it does not vanish even when `y` is a prime power — hence

```
min(1, y/(T⟨y⟩)) → 0     and     R(y,T;χ) = O_y(1/T)     as T → ∞.
```

MV make the same observation after Theorem 12.5 in order to pass to the untruncated formula.

**The constant.** MV's `C(χ)` is explicit and does not depend on parity; §2.1 and §2.2 carry the
derived constants `−(L'/L)(0,χ)` and `−c₀(χ)`, which must agree with it through the functional
equation. That check has not been done and is not needed: what §2.3 removes under `∫dψ₀/log t` is
every constant *in `ψ`* — values of `L'/L`, including this one — whatever its value. The constant
that survives into `Π*`, `log L(0,χ)`, is not of that kind; it is the lower-limit bracket, and
`C(χ)` has nothing to do with it. If (12.6) is preferred, `C(χ)` may simply be substituted.

**Naming.** MV write `E₁(x,χ)` in Theorem 12.12 for the imprimitivity correction, a different
object carrying a total-variation bound `∫_c^x |dE₁(u,χ)| ≪ (log xq)²`, and it is unrelated to the
exponential integral `E₁` of §2.3. This paper's transfer remainder is written `E_Π(y,T)` to keep
all three apart; MV's is written only with its two arguments.

**Transfer to `Π*`.** Integrate (12.6) against `dt/log t` from `2` to `y`. With `ψ₀` on the left
and the two log terms of (12.6) on the right — `arctanh(1/t)` for odd `χ`, `−log t − ½log(1 − t⁻²)`
for even, the `−log t` integrating to the `−log log y` of `I` —

the identity with the *truncated* lower-limit bracket
`Π*(2,χ) + Σ_{|γ|≤T} m_ρ li(2^ρ) − I(2,χ)`, plus `∫_2^y dR/log t`. `C(χ)` has died as in §2.3.
The truncated bracket converges as `T → ∞` (its summand is `2^ρ/(ρ log 2) + O(|ρ|^{−2})`) to
`κ(χ)`, and that limit is what (†) carries. The difference is `−Σ_{|γ|>T} m_ρ li(2^ρ)`,
independent of `y`, and is put in `E_Π`. When `χ(2) = 0`, as at `q = 4`, `Π*(2,χ) = 0`; when
`χ(2) ≠ 0` the star gives `Π*(2,χ) = χ(2)/2` and the sum-minus-integral part of the bracket
tends to `κ(χ) − χ(2)/2` instead — the same constant, differently split (tested in both signs,
§4). The remainder in (†) is
```
E_Π(y,T) := ∫_2^y dR(t,T;χ)/log t  −  Σ_{|γ|>T} m_ρ li(2^ρ).
```
The tail is `O_y(log²(qT)/T)` by Theorem 12.10 at `y = 2` (the jump there is `O(1/T)`, since
`⟨2⟩ = 1`). `R` is a remainder, not a differentiable
function, so integrate the first piece by parts rather than differentiating it:

```
∫_2^y dR/log t  =  R(y,T;χ)/log y − R(2,T;χ)/log 2 + ∫_2^y R(t,T;χ) dt/(t log²t).
```

For `y ≥ 2` we have `log y ≥ log 2`, and `log²(qtT) ≤ log²(qyT)` on `[2,y]` with
`∫_2^y dt/log²t ≪ y/log²y`, so the second term of (12.8) contributes `≪ y log²(qyT)/(T log²y)` to
the integral — smaller than its boundary term by a factor `log y`. The first term of (12.8), the
jump-correction `(log t)·min(1, t/(T⟨t⟩))`, needs more care, because `⟨t⟩ → 0` as `t` approaches a
prime power, so it cannot be bounded by a minimum over `[2,y]`. At the two boundary points it is
`O(1/(⟨y⟩T))` and `O(1/T)`, `⟨y⟩` being positive at a point. Inside the integral the `min` saturates
at `1` only in a window of width `≍ p^k/T` about each of the finitely many prime powers `p^k ≤ y`,
each window contributing `O(1/T)` against `dt/(t log t)`; outside the windows it is `t/(T⟨t⟩)`, and
`∫ dt/⟨t⟩` over `⟨t⟩ > t/T` is `O_y(log T)`. The jump-correction therefore contributes
`O_y((log T)/T)` in all. Hence

```
E_Π(y,T)  ≪  y log²(qyT)/(T log y)  +  O_y(log²(qT)/T),
```

and for fixed `y`, `E_Π(y,T) → 0` as `T → ∞`, which is all that (†) claims of its remainder. The
constant `κ(χ)` in (†) is the bracket above, not this bound.

**The `n`-sum.** `E(x,T) = Σ_{n≤N}(μ(n)/n) E_Π(x^{1/n},T)`. The factor
`1/log(x^{1/n}) = n/log x` grows with `n` and looks dangerous. Two things stop it:

1. `n ≤ N = ⌊log x/log 2⌋` forces `x^{1/n} ≥ 2`, so `log(x^{1/n}) ≥ log 2`. The denominator
   **bottoms out at a constant**; it does not tend to `0`.
2. The `n` cancels against the `μ(n)/n`:
   ```
   |(μ(n)/n) E_Π(x^{1/n},T)|  ≪  x^{1/n} log²(qTx^{1/n}) / (T log x)  +  O_x((log T)/(nT)).
   ```

The `n = 1` term dominates, so

```
E(x,T)  =  O_x( log²(qxT)/T )  →  0    (T → ∞, x fixed), uniformly in q.
```

The `O_x` absorbs `1/⟨x^{1/n}⟩` from the boundary terms — finite for each `x`, since `⟨·⟩ > 0` at
every point — and the prime-power windows below each `x^{1/n}` from the integrals. Neither is
uniform in `x`; the first degrades as any `x^{1/n}` approaches a prime power. Only finitely many `n`
occur, so this is a statement about the implied constant, not about the limit.

**Dependence on the modulus.** In (12.8) the modulus appears **only** inside `log(qxT)`, with no factor `q`,
`q^{1/2}` or `q^ε` outside, and MV's implied constant is absolute under the stated hypotheses
`x ≥ c`, `T ≥ 2`, `χ` primitive with `q > 1`. Therefore

```
E(x,T)  ≪_x  log²(qxT)/T      uniformly in q; the constant depends on x.
```

The one place `q` does escape a logarithm in this circle of results is `q^{1/2}(log q)²`, arising
from `(L'/L)(1,χ̄) ≪ 1/(1−β₁)` when `χ` has a Landau–Siegel zero. By (12.7) it sits in the constant
`C(χ)`, a constant in `ψ`, and §2.3 removes those. It never reaches `E(x,T)`.

**Real zeros.** MV sum over `|γ| ≤ T`, which **includes** any real zero of `L(s,χ)`; none is
written separately in (12.6). Theorem A adopts the same indexing, so nothing is added and nothing
is dropped. For real `χ` a real zero comes in a pair `β, 1−β`, and after §2.3 both survive into
`Π*` as `R_N(x^β)` and `R_N(x^{1−β})` — the constant `1/(1−β)` that lets MV and Granville tidy the
pair at the `ψ` level is exactly what `∫dψ₀/log t` annihilates. Granville's notes give the
alternative packaging that names the pair, and a further coarsening that discards `x^{1−β}`; the
latter is calibrated to `x → ∞` and is not available at fixed `x` with `T → ∞`. For `q = 4` the
question is empty for an elementary reason; see §4.

---

## 3. Replacing the finite cutoff by the entire function

The form one wants to write — and the form Manzoni and Hutama do write — carries Riemann's `R`
itself:

```
π*(x,χ) + π_-^*(√x) + Σ_{|γ|≤T} m_ρ R(x^ρ)  =  C(x,χ) + Ẽ(x,T).             (3.1)
```

Passing from `R_N` to `R = ERi` means adjoining the terms `n > N`. Those are legitimate *before*
truncation in `T`: for `n > N` we have `x^{1/n} < 2`, so `Π*(x^{1/n},χ) = 0` identically, and one
is adding zero. **After** truncation they are not zero. What (†) says at `y = x^{1/n} < 2` is

```
0 = −Σ_{|γ|≤T} m_ρ li(x^{ρ/n}) + I(x^{1/n},χ) + κ(χ) + E_Π(x^{1/n},T),
```

so `E_Π(x^{1/n},T)` is exactly the amount by which the truncated zero-sum fails to reproduce the
whole of `Π*` at that scale, and the tail

```
Ẽ(x,T) − E(x,T) = Σ_{n>N} (μ(n)/n) [ I(x^{1/n},χ) + E_Π(x^{1/n},T) ]
```

is an infinite sum of such terms — the constant has already cancelled, since
`κ(χ)·[m(N) + Σ_{n>N} μ(n)/n] = κ(χ)·Σ_{n≥1} μ(n)/n = 0`. That is why (3.1), the form
with `R` and the `n`-sum run to infinity, shows no constant while Theorem A does: the two are
reconciled by `Σ μ(n)/n = 0`, and `m(N) → 0` as `N → ∞`. The `I`-tail is the `t → 1` singularity of
§2.3 in the `n > N` region (`I(x^{1/n},χ) ∼ n/(2 log x)`); the `E_Π`-tail's behaviour as `T → ∞`
**is** the convergence question.

> **This is why (3.1) must not be stated with `Ẽ(x,T) = o_x(1)`, even along good ordinates.**
> `Ẽ(x,T) → 0` is *equivalent* to convergence of `Σ_{|γ|≤T} m_ρ R(x^ρ)`. For `ζ` that is exactly
> what Grobner disproves, and for real primitive `χ` it is what Paper B disproves: the `L`-zero-sum
> is `Ω_±(√T log T)`, unconditionally, so `Ẽ(x,T)` does not tend to `0` and (3.1) with `=` is
> **false**. Good ordinates do not rescue it: they sit in every `[T,T+1]`, the increment across a
> unit interval is `O(log T)` at most, against an `Ω(√T log T)` oscillation, so unboundedness
> along all `T` persists along good `T`.

Hence Theorem A is stated with `R_N`, which is what the argument actually delivers, and (3.1) is
recorded as the form whose remainder is not small. Before Paper B the question was open, and not
unnoticed. Manzoni's own comment under his answer says the passage from his (9') to (11') "should
be (more than!) double checked", citing Ingham's warning that Riemann's route to `π` does not
generalise as the route to `ψ` does (MSE 269997, answer 282848, comment of 2013-01-21;
[source discussion](https://math.stackexchange.com/questions/269997)). For `ζ` the convergence of `Σ_ρ R(x^ρ)` was
asked outright on MathOverflow in 2021 (question 386213, [source discussion](https://mathoverflow.net/questions/386213)),
and the accepted answer passes to the infinite form without proof — see the `ζ` paper's §1. For the `β`-zeros, Paper B proves divergence by
running Grobner's contour with `−L'/L` in place of `−ζ'/ζ`: the summand `ERi` is untwisted, so
`f_A`, `H_A` and `P_A` are the `ζ` paper's own with only the constant `C` changing; the horizontals
are MV Lemma 12.7 in place of Lemma 12.2; and the `Ω` comes from `ζ`'s first zero through `P_A`.
Theorems 1–4 of the `ζ` paper are not citations for this — the contour has to be run — but once
it is, the divergence is unconditional. See [Paper B](dirichlet-zero-sums.md).

**Manzoni's identity is stated as exact.** The correction term is the right one (`π_-^*(√x)`,
not a truncation of a longer `Corr`). What he writes is (3.1) with `R_N` replaced by `R`,
`T = ∞`, and `=` in place of a remainder — the same claim the `ζ` paper exists to correct,
transplanted to `q = 4` and never examined. By Paper B it is false.

---

## 4. `q = 4`, and the constant at other `D`

To recover the individual progressions, combine Theorem A's difference with
`π*(x;4,1) + π*(x;4,3) = π*(x) − π₂*(x)`, where `π₂*(x)` is the averaged contribution
of the prime 2: zero below 2, one half at 2, and one above 2. Adding or subtracting the
two equations and dividing by 2 gives the two progression counts. The half weight at
`x=2` is essential.

`D = −4`, `S_+ = {1}`, `S_- = {3}`, `π*(x,χ₋₄) = π*(x;4,1) − π*(x;4,3)`, `π_-^*(y) = π*(y;4,3)`,
zeros those of `L(s,χ₋₄) = β(s)`. Theorem A reads

```
π*(x;4,1) − π*(x;4,3) + π*(√x;4,3) + Σ_{|γ|≤T} m_ρ R_N(x^ρ) = C − (log 2)·m(N) + E(x,T),
```

the constant being `log L(0,χ₋₄) = log ½`: `h(−4) = 1`, `w(−4) = 4`.

The sum runs over `|γ| ≤ T` as everywhere else, but for this modulus it has no real zeros to
carry: `β(s) = Σ_{k≥0}(−1)^k(2k+1)^{−s}` is alternating with terms strictly decreasing to `0` for
`s > 0`, so `β(s) > 1 − 3^{−s} > 0` on `(0,1)`. No zero-free region is invoked and no hypothesis is
made — for `q = 4`, `Σ_{|γ|≤T}` and `Σ_{0<|γ|≤T}` are the same sum. For a general real primitive
`χ` the two differ by the Landau–Siegel pair, which the indexing above carries; see §2.5. (In this
section `β(s)` is the Dirichlet beta function; elsewhere `β` is the real part of a zero. The letter
is standard for both, so it is written here only with its argument.)

**Numerically.** The zeros of `L(s,χ₋₄)`: PARI 2.17.4, `lfunzeros(lfuncreate(-4), 10000, 32)` at
`realprecision 38` — 12349 ordinates to `T = 9999.88`, Riemann–von Mangoldt predicting `12348.3`
(`data/chi4_zeros_pari10000_clean.txt`). Validated by the count `θ(T)/π`,
`θ(t) = (t/2)log(4/π) + Im log Γ((3/2+it)/2)`, calibrated on LMFDB's 25 ordinates below 60, which
stays within `±1` at every thousand (a missed pair would show as `−2`), and zero for zero to
`5.5·10⁻¹²` against an independent mpmath sign-change finder on the completed function over
`[0, 3000]` and `[6000, 8000]` — 5803 of the 12349, including the four closest pairs in that range,
with gaps down to `0.025` (`compute/chi4_zeros.py`, `compute/validate_zeros.py`). Theorem A is
evaluated by `compute/theoremA_q4.py` in Elliott's branch, checked against the `ζ` paper's evaluator
to `1.8·10⁻¹⁵`; results `data/theoremA_q4_pari10000_clean.txt`.

At `x = 30.5`, `N = 4`, `m(N) = 1/6`, the combinatorial side is `4 − 5 + 1 = 0`. The left side
minus `C(x,χ)` reads `−0.1128, −0.1111, −0.1131, −0.1143` at `T = 3000, 5000, 7000, 10000`,
converging to `−(log 2)/6 = −0.1155` and not to `0` — which is how the constant was found. The
bracket `Σ_{|γ|≤T} li(2^ρ) − I(2)` reads `−0.6955, −0.6926, −0.6926, −0.6934` against
`−log 2 = −0.6931`. The archived diagnostic, the left side minus `C` minus the **finite**
bracket `B_T` times `m(N)`, reads `+0.0031, +0.0043, +0.0024, +0.0013`. At `x = 100.5`
(`m(N) = 2/15`, combinatorial side `11 − 13 + 2 = 0`) that diagnostic is
`+0.016, −0.012, +0.006, −0.005`; at `x = 1500.5` (`m(N) = 19/210`, `116 − 122 + 6 = 0`),
which needs `T ≫ x`, it is `+0.072, +0.002, −0.065, +0.050`. Theorem A's remainder subtracts
the **constant** `κ=−log 2`, so `E = diagnostic + m(N)(B_T−κ)`. Current staged outputs label
these separately. The left side minus `C` — `E + κ(χ)·m(N)`, with
`κ(χ)·m(N) = −(log 2)·19/210 = −0.063` — has octave band
`[−0.13, +0.05]` on `[3200, 6400)` about that `−0.063`.

The class-number identification was then tested where it predicts different numbers, on PARI lists
to `T = 3000` (`lfunzeros`, subdivision 16; PARI/GP recipes in `compute/gp/`, `compute/bracket_D.py`, output
`data/bracket_D_pari3000.txt`) for `D = −8` (3461 zeros) and `D = −24` (3985 zeros), both with
`χ(2) = 0` so that the lower limit `2` is clean exactly as at `q = 4`: the bracket reads `−0.002` at
`D = −8` (`h = 1`, `w = 2`, `L(0,χ) = 1`) against `0`, and `+0.691` at `D = −24` (`h = 2`,
`L(0,χ) = 2`) against `+0.693`, both within `±0.005` of their targets at every `T` from `1000` up.
It is `log L(0,χ)`, not `log Λ(0,χ)`, which would be `log √8 = 1.04` at `D = −8`. (A first
`T = 1000` list for `D = −8`, 977 zeros, read `+0.009`; it was missing the pair `948.031, 948.079`,
which the `θ(T)/π` step check caught when the list was redone.)

With `χ(2) ≠ 0` the same lists for `D = −3, −7, −11, −15` (`χ(2) = −1, +1, −1, +1`;
`2h/w = 1/3, 1, 1, 2`) give brackets `−1.101, −0.003, −0.002, +0.690` against
`log L(0,χ) = −1.099, 0, 0, +0.693`, and the sum-minus-integral parts alone read
`−0.601, −0.503, +0.498, +0.190` against §2.5's `κ(χ) − χ(2)/2 = −0.599, −0.5, +0.5, +0.193`: the
split `Π*(2,χ) = χ(2)/2` is what moves, in both signs, and the constant does not.

**Even `χ`.** The same bracket, with `I(2,χ) = −log log 2 + Σ_{k≥1} E₁(2k log 2) = 0.50652` and
`Π*(2,χ) = χ(2)/2`, on PARI lists to `T = 3000` for ten even discriminants
(`data/chi_D{D}_zeros_pari3000.txt`, count-checked as above; PARI's `L'(0,χ)` agrees with
`h(D) log ε_D` to all printed digits):

| `D` | `h` | `log ε_D` | `log L'(0,χ)` | `χ(2)` | bracket at `T = 3000` | bracket `− log L'(0,χ)` |
|---|---|---|---|---|---|---|
| 5 | 1 | 0.4812 | −0.7314 | −1 | −1.3113 | −0.5798 |
| 8 | 1 | 0.8814 | −0.1263 | 0 | −0.7061 | −0.5798 |
| 12 | 1 | 1.3170 | +0.2753 | 0 | −0.3046 | −0.5799 |
| 13 | 1 | 1.1948 | +0.1779 | −1 | −0.4020 | −0.5799 |
| 17 | 1 | 2.0947 | +0.7394 | +1 | +0.1589 | −0.5805 |
| 21 | 1 | 1.5668 | +0.4490 | −1 | −0.1308 | −0.5798 |
| 24 | 1 | 2.2924 | +0.8296 | 0 | +0.2499 | −0.5797 |
| 28 | 1 | 2.7687 | +1.0184 | 0 | +0.4383 | −0.5801 |
| 40 | 2 | 1.8184 | +1.2911 | 0 | +0.7113 | −0.5798 |
| 65 | 2 | 2.7765 | +1.7143 | +1 | +1.1336 | −0.5807 |

`log L'(0,χ)` moves by `2.4` across the table and the residual does not move: it reads
`−0.582, −0.575, −0.572, −0.574, −0.578, −0.580` at `T = 500, 1000, …, 3000` — the same six
numbers to two or three decimals at every `D`, the tail of the bracket sum being nearly
`q`-independent — oscillating about `−C₀ = −0.5772`, with mean `−0.576` over `T ≥ 1000`. The odd
runs show the same tail: at `T = 3000` all six sit `0.002`–`0.003` below their targets, and
subtracting that common offset from the even residual leaves `−0.5772`. That is
`κ(χ) = log L'(0,χ) − C₀`, as Landau's convention says (§2.3); the number was measured before the
page was read. `h = 2` at `D = 40, 65` exercises the class number, and `χ(2) = −1, 0, +1` all
occur, so the split `Π*(2,χ) = χ(2)/2` of §2.5 is exercised in both signs — the residual does not
see it, as it must not.

On the same 12349 zeros, replacing `R_N` by `R` — the substitution §3 refuses — gives at `x = 30.5`
`−0.80, −0.43, −1.96, −2.44` at the same four `T`, against `E` fixed at `−0.11`: the growth doubles
as `√T log T` doubles between `T = 3000` and `10000`, at the scale Paper B proves
(`data/theoremA_q4_eri_pari10000_clean.txt`).

Manzoni's display is this with `R_N` replaced by `R`, `T = ∞`, `=` in place of a remainder, and no
constant — which for him is not an omission but an artefact of the infinite form: with the `n`-sum
run to infinity the constant multiplies `Σ_{n≥1} μ(n)/n = 0` (§3). It is false by Paper B.

---

## 5. Prime ideals

Counting prime **ideals** by norm removes the obstruction of §1 outright: the Möbius inversion of
`J_K` runs over ideal norms, so no character powers arise and no `Corr` is needed. That is the
cleaner inversion, and it is the right route if the object is `π_K(x)`; but its combinatorial
content is then the ordinary Möbius inversion of `J_K`, i.e. a remark rather than a theorem. It is
worth stating in that spirit, together with the observation that Gaussian primes are the `q = 4`
case of Theorem A combined with the `ζ`-identity, via `ζ_{Q(i)} = ζ·L(s,χ₋₄)` — which is Hutama's
`π^G_0` formula, and inherits from `ζ` the same failure of (3.1). The analytic side of `Q(i)` —
`Γ_C(s) = (2π)^{−s}Γ(s)`, Mellin partner `e^{−t}`, reflection `sin(πs)` — belongs to Paper B.

## References

- [G] H. Grobner, *On divergence related to Riemann–von Mangoldt's explicit formula of the
  prime-counting function*, arXiv:2609.02713v1 (2026), https://arxiv.org/abs/2609.02713.
- [MV] H. L. Montgomery and R. C. Vaughan, *Multiplicative Number Theory I: Classical Theory*,
  Cambridge Studies in Advanced Mathematics 97, Cambridge University Press (2007),
  Theorem 12.10, pp. 403–405; ISBN 978-0-521-84903-6.
- [La] E. Landau, *Nouvelle démonstration pour la formule de Riemann sur le nombre des nombres
  premiers inférieurs à une limite donnée et démonstration d'une formule plus générale pour le
  cas des nombres premiers d'une progression arithmétique*, Ann. Sci. École Norm. Sup. (3)
  25 (1908), 399–442, Part II §XV, (38), pp. 441–442; DOI 10.24033/asens.595.
- A. Granville, *The prime number theorem for arithmetic progressions*, course notes,
  Chapter 11, especially (11.8.4)–(11.8.5); archived extract in the research workspace.
- H. Riesel and G. Göhl, *Some calculations related to Riemann's prime number formula*,
  Math. Comp. 24 (1970), 969–983.
- D. J. Hutama, *Implementation of Riemann's Explicit Formula for Rational and Gaussian
  Primes in Sage*, manuscript dated 14 August 2017; archived in the research workspace.
- J. Elliott, *Analytic Number Theory and Algebraic Asymptotic Analysis*, World Scientific
  (2025); arXiv:2407.17820, §4.4.
- [MO] *Is π(x)=R(x)−Σρ R(xρ) correct at all?*, MathOverflow question 386213 (2021),
  https://mathoverflow.net/questions/386213/.
- G. Manzoni, answer and comments to Math StackExchange question 269997,
  https://math.stackexchange.com/a/282848.
- Companions by Open: [the zeta series](zeta-zero-sums.md) and [the Dirichlet zero series](dirichlet-zero-sums.md).
