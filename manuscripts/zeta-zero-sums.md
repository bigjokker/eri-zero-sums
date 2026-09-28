# The sum $\sum_\rho \mathrm{ERi}(\rho\log x)$: a Bohr series, its envelope, and the failure of Riesz summation

**Open**

## Abstract

The partial sums of the Riemann–von Mangoldt ERi series over non-trivial zeta
zeros oscillate unconditionally in both directions at scale $\sqrt T\log T$, with an
explicit lower-bound constant obtained from the first zero. Every fixed nonnegative real
Riesz order retains that oscillation. Under the Weak Mertens Conjecture, the normalized
sums have an absolutely convergent Bohr expansion in $\log T$. Linear independence of
the positive zeta ordinates identifies its Bessel-product distribution and exact envelope.
We also treat the trivial-zero series and distinguish finite numerical amplitude sums
from the infinite conditional envelope.

## 1. Introduction

Riemann's function
$$
R(x)=\sum_{n\ge1}\frac{\mu(n)}{n}\,\mathrm{li}\bigl(x^{1/n}\bigr)
$$
is usually extended to complex arguments through the entire function
$$
\mathrm{ERi}(z):=\sum_{n\ge1}\frac{\mu(n)}{n}\,\mathrm{Ei}\Bigl(\frac zn\Bigr),
\qquad \mathrm{ERi}(0):=1,
\qquad R(x^\rho)=\mathrm{ERi}(\rho\log x),
$$
equivalently by Gram's series $\mathrm{ERi}(z)=1+\sum_{k\ge1}z^k/(k\,k!\,\zeta(k+1))$. The
identity
$$
\pi_0(x)=R(x)-\sum_\rho R(x^\rho)
$$
appears in Zagier's *The first 50 million prime numbers*, in Elliott's monograph, and in the
textbook literature, the sum being taken over all zeros of $\zeta$ with the non-trivial ones
ordered by $|\Im\rho|$. Grobner [G] has shown that the sum over the non-trivial zeros does not
converge: writing
$$
\Sigma R_T(x):=\sum_{\substack{\zeta(\rho)=0\\ 0<|\Im\rho|\le T}} m_\rho\,\mathrm{ERi}(\rho\log x),
\qquad
\Theta:=\sup\{\Re\rho:\zeta(\rho)=0,\ 0<\Re\rho<1\},
$$
he proves that for every fixed $x>1$ and every $\theta<\Theta$ the partial sums $\Sigma R_T(x)$
are not $O(T^\theta)$, so that $\limsup_T|\Sigma R_T(x)|=\infty$.

The identity's standing before [G] deserves a sentence. At the level of $\Pi(x)$, with
$\mathrm{li}(x^\rho)$ in place of $R(x^\rho)$, it is Riemann's formula as proved by von Mangoldt
and, by a new route, by Landau [La] — his (25)–(26), p. 425: the zero-sum in conjugate pairs, in
order of increasing ordinate, the limit along a sequence of heights and then along all; its
$R$-form asks for an interchange of that sum with the Möbius sum over $n$, which no source
proves. On MathOverflow the convergence was asked outright in 2021 [MO]; the accepted answer
writes the identity with the Möbius sum cut at $m\le X$, $X>\log_2x$ — where it is exact, with the
constant $-\log2\sum_{m\le X}\mu(m)/m$ — and then passes to the infinite form without proof; an
answer of 2024 reports that the sum over the non-trivial zeros diverges numerically; Grobner's
own answer, posted 2026-09-05, points to [G]. The $\arctan$ form of the identity that stood on
Wikipedia until July 2021, discussed in the same thread, is Riesel–Göhl's truncation read as an
identity (§1.1). Landau's formula is an identity because it sums $\mathrm{li}(x^\rho)$ and keeps
the trivial zeros as an integral; the $R$-form is not, because it Möbius-inverts first and
truncates in $T$ afterwards. The two are not rearrangements of one series: they differ by the tail
$\sum_{n>N}(\mu(n)/n)[\cdots]$ displayed in Paper A §3, and it is that tail, not the ordering of the
zeros, whose limit fails to exist.

This paper asks what the sum *is*. Unconditionally, it diverges in both directions at the
$\sqrt T\log T$ scale with an explicit lower-bound constant, and no fixed Riesz order repairs the
divergence. Under WMC, after dividing by $\sqrt T\log(T/2\pi)$, it equals a Bohr almost periodic
function of $\log T$ with an absolutely convergent Fourier expansion, plus a term tending to zero;
the limiting distribution then has compact support. Under WMC and LI, its characteristic function
is a Bessel product and its exact envelope is the infinite sum in Theorem 2.

### 1.1 Results

Throughout, $x>1$ is fixed, $\rho=\beta+i\gamma$ runs over non-trivial zeros, and
$$
L:=R(x)-\pi_0(x)-I(x),\qquad I(x):=\frac1{\log x}-\frac1\pi\arctan\frac{\pi}{\log x},
$$
is a numerical candidate for the constant about which $\Sigma R_T$ oscillates. Two cautions about
$L$. First, $I(x)$ is the
value the folklore assigns to the trivial-zero sum $\sum_{k\ge1}R(x^{-2k})$ (Riesel–Göhl; Elliott
(5.1.8)); that series in fact **diverges** ([G] Prop. 4.1, and §5.6 at the resolution of Theorem
1), and $I(x)$ is the value at $s=0$ of the continuation of its Dirichlet series $D_x(s)$ (§5.6) —
an analytic home, not a sum. Second, the contour of §3.1 produces a constant $D_{a,x}$, and its
identification with $L$ is numerical (§3.1, Appendix B), not derived. Neither caution touches
Theorems 1–4: each is invariant under replacing $L$ by any other constant, a constant being
$o(\sqrt T\log T)$. WMC denotes the Weak Mertens
conjecture $\int_1^X(M(u)/u)^2\,du=O(\log X)$, which implies RH and the simplicity of the
critical-line zeros (the unconditional integral $1/\zeta(s)=s\int_1^\infty M(u)u^{-s-1}\,du$
for $\sigma>1$ extends holomorphically to $\sigma>\tfrac12$ by Cauchy–Schwarz and WMC, giving RH;
Titchmarsh Theorem 14.29(A) gives simplicity). The RH-conditional display (14.26.3) is not used
to deduce RH. LI denotes linear independence of the positive
ordinates over $\mathbb{Q}$.

| § | statement | hypotheses |
|---|---|---|
| 2 | $r_\rho=\pi(2\pi)^{-\rho}\chi(\rho)/((\rho-1)\zeta'(\rho))$; on the critical line $\lvert r_\rho\rvert=\sqrt{\pi/2}\,/(\lvert\rho-1\rvert\lvert\zeta'(\rho)\rvert)$; with $N$ zeros the error is $\Omega(\sqrt{N\log N})$ | **none** |
| 3 | $\Sigma R_T(x)=\Omega_\pm(\sqrt T\log T)$, constant $0.008144$ at $x=1.1$ | **none** |
| 5 | $\Sigma R_T^{(k)}(x)=\Omega_\pm(\sqrt T\log T)$ for every fixed $k\ge0$: no Riesz mean sums the series | **none** |
| 4 | $F=\phi+o(1)$ with $\phi$ Bohr u.a.p.; $\hat\mu=\prod_{\gamma>0}J_0(2\lvert a_\rho\rvert\xi)$; $\limsup\lvert F\rvert=\sum_{\gamma>0}2\lvert a_\rho\rvert$ ($0.039108$ from $60$ zeros) | WMC, LI |
| 5 | exact $\limsup$ of the weighted sums, with $W_k(\gamma)$; first-sixty-term ratios $21.6$, $189.0$, $1017.5$ for $k=1,2,3$ | WMC, LI |

Against [G]'s $\ne O(T^\theta)$ for $\theta<\Theta$: the $\varepsilon$ is gone, both signs are
present, the constant is a number, and the exponent no longer waits on whether $\Theta$ is
attained — which is open, and would be needed even for $\Omega_\pm(x^\Theta)$ (Elliott
Prop. 9.1.3(2)).

Section 3 runs on [G]'s auxiliary functions $H_A$, $P_A$ and his rectangle $\Delta_T(a,c)$; that
contour is used, not reproved, and his Lemma 2.5 supplies its good ordinates. Section 4 uses the
same $H_A$ and $P_A$ but a different rectangle — the Mellin shift for $P_A$, written out in §4.4 —
whose height is chosen by Titchmarsh's Theorem 14.16, not by Lemma 2.5; §4.5 keeps the two
notions of good ordinate apart.

### 1.2 One extra factor of $1/\rho$

The organising fact is a comparison with $\psi(x)-x$. There the explicit formula has amplitudes
$1/\rho$; here, by §2, the residue carries an extra $1/\rho$ and the amplitudes are
$$
a_\rho\;\asymp\;\frac1{\rho^2\zeta'(\rho)} .
$$
That single factor does three separate things, and each is a section of this paper.

1. **$\sum_\rho|a_\rho|<\infty$** (§4.1, by Cauchy–Schwarz against $\sum1/|\rho|^2$ under WMC).
   So the frequency series is Bohr *uniformly* almost periodic, not merely $B^2$ as
   $\psi(x)-x$ and $M(x)$ are; $F$ equals that series plus $o(1)$, the limiting law has compact
   support, and the $\limsup$ is an equality rather than an $\Omega$-bound.
2. **$\int|\Phi(\sigma_0+it)|\,dt<\infty$** (§4.3): the factor $1/s$ in $\Phi=G/s$ buys a second
   power of $t$, and without it the contour shift on which §4 rests does not close.
3. **A pointwise truncation tail** (§4.1): the same Cauchy–Schwarz tail as (1), now read as a
   pointwise bound rather than as absolute convergence. Akbary–Ng–Shahabi's condition (1.7) is a
   mean-square condition because for $M(x)$ the tail is large on a thin set; ours is never large,
   so (1.7) is free.

The same factor is what makes the Riesz weight useless in §5: it is already responsible for the
convergence that a summability method is supposed to supply, and the divergence that remains is
of a kind no kernel of fixed width can touch.

### 1.3 Conventions

We keep [G]'s notation. With $B(y)=-\sum_{n\le y}\mu(n)/n$, $Q(u)=u^{-1}B(u^{-1})\in L^1(0,1)$,
and $q_A(u)=e^{-Au}Q(u)$ for a fixed $A>0$:

| object | definition | source |
|---|---|---|
| $\mathrm{ERi}(z)$ | $=-\int_0^1 Q(u)e^{zu}\,du$ | [G] Prop. 1.4 |
| $M_A(z)$ | $=\int_0^1 q_A(u)u^z\,du=\sum_{j\ge0}\dfrac{(-1)^{j+1}A^j}{j!\,(z+j)\zeta(z+j+1)}$ | [G] Lem. 1.5, (1.6) |
| $f_A(Y)$ | $=\Re\,\mathrm{ERi}(-A+iY)=-\int_0^1q_A(u)\cos(Yu)\,du$; real and bounded | [G] §1.3 |
| $H_A(z)$ | $=\displaystyle\int_1^\infty f_A(Y)\,Y^{-z}\,dY$ | [G] §1.3 |
| | $=-\Gamma(1-z)\sin\dfrac{\pi z}2\,M_A(z-1)+\displaystyle\sum_{k\ge0}\frac{(-1)^kM_A(2k)}{(2k)!\,(2k+1-z)}$ | [G] Prop. 1.7, (1.8) |
| $P_A(Y)$ | $=\displaystyle\int_1^Y\log(y/C)\,f_A(y)\,dy$ | [G] §1.4 |
| $C$ | $=2\pi\log x$, forced by [G] Prop. 2.15 — not a free parameter | §2 |

**A clash of conventions, stated in the open.** $H_A$ is defined with $Y^{-z}\,dY$; the Mellin
transforms of §3 and §4 are taken with $Y^{-s-1}\,dY$, which is the convention in which Landau's
theorem and Lemma A are stated. A pole of $H_A$ at $z=\rho$ therefore sits at $s=\rho-1$ in the
second convention. Applying Lemma A to the *bounded* function $\log(\cdot/C)f_A$ gives
$\sigma_0=\beta-1=-\tfrac12$ and the conclusion $\Omega_\pm(Y^{-1/2}\log Y)$ — decaying, true,
and useless. The growing object is the primitive $P_A$, whose Mellin transform in the
$Y^{-s-1}$ convention has a *double* pole at $s=\rho$. This is set up once, in §3.2, and used
again in §4.

Finally, $\Omega_\pm(g)$ means $\limsup f/|g|>k$ and $\liminf f/|g|<-k$ for some $k>0$; and near
a pole of order $r$ at $s_0$ we write $F(s)=\sum_{m=1}^r c_m(s-s_0)^{-m}+\phi(s)$, so that $c_r$
is the leading Laurent coefficient.

---

## 2. The residue is the functional equation

[G]'s meromorphic continuation (Prop. 1.7, eq. (1.8)) reads
$$
H_A(z)
  =-\Gamma(1-z)\sin\Bigl(\frac{\pi z}{2}\Bigr)\,M_A(z-1)
   +\sum_{k\ge0}\frac{(-1)^k M_A(2k)}{(2k)!\,(2k+1-z)} .
$$
The second sum is holomorphic on a neighbourhood of every non-trivial zero ([G], p. 8). The
$j\ge1$ terms of $M_A(z-1)$ are holomorphic at a non-trivial zero $\rho$ as well
($\Re(\rho+j)>1$ for $j\ge1$, unconditionally). The only singularity is the $j=0$ term of $M_A(z-1)$,
$$
M_A(z-1)=-\frac1{(z-1)\zeta(z)}+\text{holomorphic at }\rho .
$$
Thus $A$ drops out of the residue: at a **simple** zero,
$$
\operatorname{Res}_{z=\rho}H_A=\frac{\Gamma(1-\rho)\sin(\pi\rho/2)}{(\rho-1)\zeta'(\rho)} .
$$
Write $\chi(s)=2^s\pi^{s-1}\sin(\pi s/2)\Gamma(1-s)$, so that $\zeta(s)=\chi(s)\zeta(1-s)$
(Titchmarsh (2.1.9)). Then
$$
\pi\,(2\pi)^{-\rho}\chi(\rho)=\Gamma(1-\rho)\sin\Bigl(\frac{\pi\rho}{2}\Bigr),
\qquad\text{because}\qquad
\pi\cdot(2\pi)^{-\rho}\cdot 2^\rho\pi^{\rho-1}=\pi\cdot\pi^{-\rho}\cdot\pi^{\rho-1}=1 .
$$
Hence, at a simple zero,
$$
\tag{$\ast$}
r_\rho:=\operatorname{Res}_{z=\rho}H_A=\frac{\pi\,(2\pi)^{-\rho}\chi(\rho)}{(\rho-1)\zeta'(\rho)} .
$$
This is an identity, not an asymptotic.

On the critical line $|\chi(\tfrac12+it)|=1$ (Titchmarsh (2.1.9); Montgomery–Vaughan §10.1) and
$|(2\pi)^{-\rho}|=(2\pi)^{-1/2}$, so
$$
\bigl|\Gamma(1-\rho)\sin(\pi\rho/2)\bigr|=\frac\pi{\sqrt{2\pi}}=\sqrt{\frac\pi2},
\qquad
|r_\rho|=\frac{\sqrt{\pi/2}}{|\rho-1|\,|\zeta'(\rho)|} .
$$
The cancellation one observes numerically to 25 digits —
$|\Gamma(\tfrac12-it)|^2=\pi/\cosh(\pi t)$ against
$|\sin(\pi/4+i\pi t/2)|^2=\cosh(\pi t)/2$ — is this identity and nothing more. No $t$-dependence
survives except through the denominator, and no $A$.

Off the line, Cao–Tanigawa–Zhai (arXiv:2003.11349, (1.3)) give, uniformly for $\sigma$ in a fixed
strip,
$$
\chi(\sigma+it)=\Bigl(\frac{|t|}{2\pi}\Bigr)^{1/2-\sigma-it}e^{i(t\pm\pi/4)}\bigl(1+O(|t|^{-1})\bigr),
\qquad
|\chi(\sigma+it)|\asymp\Bigl(\frac{|t|}{2\pi}\Bigr)^{1/2-\sigma} .
$$
The residue at a zero with $\Re\rho=\beta$ picks up a factor $|t|^{1/2-\beta}$ — polynomial, not
exponential. This is the one thing that must not be carried over from [G] §4, where the residues
of $D_x(s)$ carry a bare $\Gamma(1-\rho)$ and $|\Gamma(1-\rho_1)|=5.71\cdot10^{-10}$; the
divergence there is exponentially suppressed, and here it is of order one.

The same $2\pi$ is why $C=2\pi\log x$ is forced: on $\Re s=-a$ the functional equation gives
$-\zeta'/\zeta(-a+it)=\log(t/2\pi)+\cdots$ ([G] (2.14)), and the substitution $y=t\log x$ in
Prop. 2.15 turns $\log(t/2\pi)$ into $\log(y/C)$ with exactly this $C$. The natural variable is
$T/2\pi$, the Riemann–von Mangoldt counting variable, and $\log(Y/C)=\log(T/2\pi)$ carries no
leftover $x$. In terms of the number $N$ of zeros used, $T\sim2\pi N/\log N$, so
$\sqrt T\log T\asymp\sqrt{N\log N}$: **using $N$ zeros in the $\mathrm{ERi}$-form incurs an error
of at least that order** (Theorem 1), and of exactly that order under WMC and LI (Theorem 2).

**Multiple zeros.** Platt, Math. Comp. 86 (2017), Theorem 5.1: every zero with
$|\Im\rho|\le3.0610046\cdot10^{10}$ is simple and has $\Re\rho=\tfrac12$. Every computation in
this paper sits six orders of magnitude inside that certificate, so $m_\rho=1$ throughout is a
theorem. In the general statement, a zero of order $m$ makes $1/\zeta$ a pole of order $m$;
replace $1/\zeta'(\rho)$ by the leading principal part $m!/\zeta^{(m)}(\rho)$. The logarithmic
derivative $-\zeta'/\zeta$ always has a simple pole of residue $-m_\rho$, which is why [G]'s
contour (Prop. 2.1) still produces $m_\rho\,\mathrm{ERi}(\rho\log x)$ and not a derivative of
$\mathrm{ERi}$.

**Check.** The closed form $(\ast)$ reproduces
$$
|r_{\rho_1}|,\ldots,|r_{\rho_8}|=0.11172,\;0.05243,\;0.03652,\;0.03159,\;0.02753,\;0.01722,\;
0.02055,\;0.01578
$$
to all printed digits; the $\Gamma$-form and the $\chi$-form agree to $10^{-40}$ (Appendix B).

$(\ast)$ is also the growth of $H_A$ *off* the zeros, and in that role it is what makes the
contour shift of §4 possible: the object to bound on a vertical line is not a bare
$\Gamma(1-z)\sin(\pi z/2)$ but $\chi(z)/((z-1)\zeta(z))$, of known modulus on every line. See
§4.2.

---

## 3. From $P_A$ to $\Sigma R_T$: an unconditional $\Omega_\pm$

### 3.1 What [G] supplies

Fix $x>1$. For $T>0$ not a zero-ordinate, [G]'s rectangle $\Delta_T(a,c)$ (left side
$\Re s=-a\in(-1,0)$, right side $\Re s=c\in(1,2]$) and the residue theorem (Prop. 2.1) give
$$
\mathrm{ERi}(\log x)-\Sigma R_T(x)=I_c(T)-I_{-a}(T)+H_{a,c}(T).
$$
Three evaluations, as $T\to\infty$ along **good ordinates** (Lem. 2.5: every large interval
$[T,T+1]$ contains a $T'$ off the zeros on which the horizontal integrals are $o(1)$):

- **Right vertical** (Prop. 2.6): $I_c(T)\to\pi_0(x)$.
- **Left vertical** (Prop. 2.15), with the specific choice $A=a\log x$ and $C=2\pi\log x$:
  $$
  I_{-a}(T)=\frac1{\pi\log x}\,P_A\bigl(T\log x\bigr)+C_{a,x}+o(1).
  \]
  Every term furnished by the functional equation of $\zeta$ on $\Re s=-a$, except the
  $\log(t/2\pi)$ summand, contributes a constant plus $o(1)$. The logarithmic summand *is* $P_A$,
  after the change of variables $y=t\log x$.
- **Horizontals:** $o(1)$ along the same good ordinates.

Collecting constants into $D_{a,x}$, [G]'s §3 yields, along good $T'\in[T,T+1]$,
$$
\tag{$\dagger$}
\Sigma R_{T'}(x)=\frac1{\pi\log x}\,P_A(T'\log x)+D_{a,x}+o(1),
$$
the $o(1)$ uniform for $T'$ in the window. Explicitly $D_{a,x}=R(x)-\pi_0(x)+C_{a,x}$ with
$C_{a,x}$ the constant of Prop. 2.15; the left side does not depend on $a$, so the $a$-dependence
of $C_{a,x}$ cancels against that of $P_A$ through $A=a\log x$. Numerically $D_{a,x}$ appears to
equal $L$ as defined in §1.1, i.e. $C_{a,x}=-I(x)$; this is not derived. The normalised theorems
below hold with $D_{a,x}$ or any other fixed centre. We write $L$ in their normalisations. Since
$P_A'(Y)=\log(Y/C)f_A(Y)=O(\log Y)$ ($f_A$ bounded, [G] Prop. 1.4), a shift $T'=T+O(1)$ changes
$P_A(T\log x)$ by $O(\log T)$ and does not affect any $\Omega(\sqrt T\log T)$ statement.

### 3.2 The Mellin transform of $P_A$

Let
$$
F(s)=\int_1^\infty P_A(Y)\,Y^{-s-1}\,dY,
\qquad
G(z)=\int_1^\infty\log(Y/C)\,f_A(Y)\,Y^{-z}\,dY=-H_A'(z)-(\log C)H_A(z),
$$
both initially in $\Re>1$. [G]'s integration by parts (proof of Prop. 1.13) gives $G(s)=s\,F(s)$
there, hence as meromorphic functions:
$$
\tag{3.1}
\Phi(s):=F(s)=\frac{G(s)}{s}.
$$
This is the only setup of $\Phi=G/s$ in the paper; §4 refers back to (3.1).

By $(\ast)$, $H_A$ has a simple pole of residue $r_\rho$ at each simple zero, so
$$
-H_A'(z)=\frac{r_\rho}{(z-\rho)^2}+O(1)
$$
and $G$ has a **double** pole there with leading coefficient $r_\rho$ — the $(\log C)H_A$ term is
only simple and cannot cancel it. Hence $F$ has a double pole at $\rho\ne0$ with
$$
c_2:=c_{-2}(F)=\frac{r_\rho}{\rho}.
$$
(As warned in §1.3, this is the point at which one must work with $P_A$ and not with
$\log(\cdot/C)f_A$.)

### 3.3 Lemma A

> **Lemma A.** Let $f\colon[1,\infty)\to\mathbb{R}$ be locally bounded and locally integrable, and
> suppose $F(s)=\int_1^\infty f(x)x^{-s-1}\,dx$ converges in some half-plane
> and continues meromorphically. Let $s_0=\sigma_0+it_0$ with $t_0\ne0$ be a pole of $F$ of order
> $r\ge1$, and assume $F$ has no real pole with real part $\ge\sigma_0$ (otherwise subtract the
> residues of $F(s)X^s$ at those poles as a main term $\mathrm{MT}$ and replace $f$ by
> $f-\mathrm{MT}$). Let $c_r$ be the leading Laurent coefficient of $F$ at $s_0$. Then
> $$
> \limsup_{x\to\infty}\frac{f(x)}{x^{\sigma_0}(\log x)^{r-1}}\ge\frac{|c_r|}{(r-1)!},
> \qquad
> \liminf_{x\to\infty}\frac{f(x)}{x^{\sigma_0}(\log x)^{r-1}}\le-\frac{|c_r|}{(r-1)!} .
> \]
> In particular $f(x)=\Omega_\pm\bigl(x^{\sigma_0}(\log x)^{r-1}\bigr)$.

*Mellin identity used.* For $\Re s>\sigma_0$,
$\int_1^\infty x^{\sigma_0}(\log x)^{r-1}x^{-s-1}\,dx=(r-1)!\,(s-\sigma_0)^{-r}$.

*Proof.* Let $k>0$ and
$$
g_k(x)=f(x)+k\,x^{\sigma_0}(\log x)^{r-1},
\qquad
G_k(s)=\int_1^\infty g_k(x)x^{-s-1}\,dx=F(s)+\frac{k\,(r-1)!}{(s-\sigma_0)^r}
$$
(the subscripted $G_k$ is an auxiliary transform internal to this proof, not the $G$ of (3.1)).
Suppose $g_k(x)>0$ for all large $x$. Landau's theorem for Dirichlet integrals (Ingham,
*The Distribution of Prime Numbers*, Theorem H, pp. 88–89) puts the abscissa of convergence
$\sigma_c$ of $G_k$ on the real axis: if $c(x)$ is real and of constant sign for all sufficiently
large $x$, then the real point of the line of convergence of the integral is a singularity of
the function it represents. (No global boundedness is required. Lowry-Duda's Theorem 3 adds
that extra hypothesis; the comparison integral defining $K$ below uses only local
boundedness of $P_A$ on $[1,Y]$.)
The uncancelled pole of $F$ at $s_0$ forces $\sigma_c\ge\sigma_0$; $F$ has
no real pole with real part $\ge\sigma_0$ and the added term is singular only at $\sigma_0$, so
$G_k$ has no real singularity with real part $>\sigma_0$, hence $\sigma_c\le\sigma_0$. Thus
$\sigma_c=\sigma_0$, and $G_k(\sigma)$ is given by the integral for $\sigma>\sigma_0$. (If $F$ had
a pole with real part $>\sigma_0$, this already contradicts, before any Laurent comparison.) For
$\sigma>\sigma_0$ one has $|G_k(\sigma+it_0)|\le K+G_k(\sigma)$ with
$K=2\int_1^Y|g_k|x^{-\sigma_0-1}\,dx<\infty$ independent of $\sigma$. Multiply by
$(\sigma-\sigma_0)^r$ and let $\sigma\downarrow\sigma_0$: the left side tends to $|c_r|$, the
right side to $k\,(r-1)!$. Thus $k\ge|c_r|/(r-1)!$, and choosing $0<k<|c_r|/(r-1)!$ is a
contradiction. Hence $g_k$ is not eventually positive, i.e.
$f(x)<-k\,x^{\sigma_0}(\log x)^{r-1}$ infinitely often, which is the $\Omega_-$ statement; the
$\Omega_+$ argument is the same with $-f$. $\square$

**Remarks.** (i) The method is Lowry-Duda, arXiv:1910.09969v3, proof of Theorem 1, after deleting
the one sum-specific line: his discrete data $(a(n),\lambda_n,v)$ serve only to produce
$F=D_\lambda V$ in the half-plane of absolute convergence, and we start from $F$. His Theorem 1
is stated for the weighted sum $\sum a(n)v(\lambda_n/x)$ and does not literally cover a Dirichlet
integral of a continuous function. (ii) His printed $G_k$ carries $k/(r-1)!$ rather than
$k\,(r-1)!$, which is inconsistent with his own Mellin identity, and a minus in $g_k$ where the
$\Omega_-$ argument needs a plus; we take the plus and the identity. For a double pole
$(r-1)!=1$, so the factorial discrepancy is invisible; the sign is not. (iii) This is a bound
from **one** non-real pole. A conjugate pair can be twice as large, and Pintz and Révész
extract $\pi/2$ rather than $2$ by a kernel argument; that sharpening does not transfer here,
for reasons given in Appendix A.

### 3.4 Application at $\rho_1$

Take $\rho=\rho_1=\tfrac12+i\gamma_1$, $\gamma_1=14.1347251417346\ldots$, simple and on the
critical line unconditionally (Platt, Theorem 5.1). Then $r=2$, $\sigma_0=\tfrac12$,
$t_0=\gamma_1\ne0$ and $(r-1)!=1$. There is no real pole of $F$ with real part $\ge\tfrac12$.
For $\Re s>1$, the integral defining $H_A$ is holomorphic. At $s=1$, the zero of $1/\zeta$
cancels the $j=0$ denominator in $M_A(s-1)$, while the remaining apparent poles cancel in
Grobner's integral kernel. On $(\tfrac12,1)$, $\zeta(s)$ has no real zero, and the $j\ge1$
terms have $\Re(s+j)>1$. Division by $s$ adds no pole on this interval. Thus $\mathrm{MT}=0$.
If $\zeta$ had a zero with $\beta>\tfrac12$,
the pole of $F$ there would force $\sigma_c>\tfrac12$ with no real singularity in that
half-plane, so Landau would already contradict $g_k>0$ and the Laurent comparison would not be
reached; Theorem 1 is unconditional either way. Lemma A gives
$$
\limsup_{Y\to\infty}\frac{P_A(Y)}{Y^{1/2}\log Y}\ \ge\ \Bigl|\frac{r_{\rho_1}}{\rho_1}\Bigr|,
\qquad
\liminf_{Y\to\infty}\frac{P_A(Y)}{Y^{1/2}\log Y}\ \le\ -\Bigl|\frac{r_{\rho_1}}{\rho_1}\Bigr| .
$$
Numerically (40 digits, Appendix B),
$$
|r_{\rho_1}|=0.11172232897,\qquad |\rho_1|=14.1435658457,\qquad
\Bigl|\frac{r_{\rho_1}}{\rho_1}\Bigr|=0.00789916278 .
$$

### 3.5 Transfer

Set $Y=T\log x$ in $(\dagger)$. Then $Y^{1/2}=\sqrt T\sqrt{\log x}$ and $\log Y\sim\log T$, so
$$
\limsup_{T\to\infty}\frac{\bigl|\Sigma R_T(x)-L\bigr|}{\sqrt T\,\log T}
  \ \ge\ \frac1{\pi\sqrt{\log x}}\Bigl|\frac{r_{\rho_1}}{\rho_1}\Bigr| ,
$$
which at $x=1.1$ is $0.00814444597$. Both signs survive, hence

> **Theorem 1.** For every fixed $x>1$,
> $\Sigma R_T(x)=\Omega_\pm\bigl(\sqrt T\,\log T\bigr)$, unconditionally, with constant
> $\bigl|r_{\rho_1}/\rho_1\bigr|/(\pi\sqrt{\log x})$; at $x=1.1$ this is $0.008144$.

This is [G]'s Main Theorem with the $\varepsilon$ removed, signs added, and an explicit constant,
at the cost of naming a zero rather than $\Theta$. Naming $\Theta$ would require attainment,
which is open off RH (Elliott Prop. 9.1.3(2)); the argument above needs only that $\rho_1$ is a
zero, which is a computation.

The constant is smaller than the observed peaks $0.014$–$0.026$ of $F$, as it must be: those see
many zeros. It is also smaller than the first sixty terms of the conditional envelope, $0.039108$
at $x=1.1$. The full envelope is the unevaluated positive series of Theorem 2 under WMC and LI.
Two further rows of that ladder, and why the natural sharpening
does not reach them, are in Appendix A.

---

## 4. The Bohr series: the truncation error is pointwise

Akbary–Ng–Shahabi [ANS] (1.6)–(1.7) ask for
$$
\phi(y)=c+\Re\sum_{\lambda_n\le X}r_ne^{i\lambda_ny}+E(y,X),
\qquad
\frac1Y\int_{y_0}^Y\bigl|E(y,e^Y)\bigr|^2dy\to0 .
$$
Ng needs a mean-square tail because $\sum1/|\rho\zeta'(\rho)|$ diverges. Ours converges
absolutely, which is the same fact that puts the frequency series in the Bohr class rather than
$B^2$, and it collapses (1.7) once a pointwise expansion
$F(u)=\sum_\rho a_\rho e^{i\gamma_\rho u}+o(1)$ is in hand. That expansion is a contour shift
for the inverse Mellin of $P_A$, written out in §4.4. Throughout this section, WMC is in force
(hence RH, and the critical-line zeros simple: Titchmarsh §14.29 and Thm 14.29(A)), $x>1$ is fixed, and $\delta=1/4$,
$\sigma_1=\tfrac98$ are fixed once and for all. No limit $\delta\to0$ is taken.

### 4.1 Why (1.7) is free

Under WMC, $\sum_\rho1/|\rho\zeta'(\rho)|^2<\infty$ (Titchmarsh Theorem 14.29(B), eq. (14.29.4)). With
$|a_\rho|\asymp1/(|\rho|^2|\zeta'(\rho)|)$ and $\sum1/|\rho|^2<\infty$ — under RH this is
$2+\gamma-\log4\pi=0.0461914\ldots$ — Cauchy–Schwarz on the tail gives
$$
\tag{4.1}
\sum_{\gamma>X}|a_\rho|
  \ \ll_x\ \Bigl(\sum_{\gamma>X}\frac1{|\rho|^2}\Bigr)^{1/2}
          \Bigl(\sum_{\rho}\frac1{|\rho\zeta'(\rho)|^2}\Bigr)^{1/2}
  \ \ll\ \sqrt{\frac{\log X}{X}}\qquad(X\to\infty),
$$
the second factor bounded by the full convergent sum, the square of the first $\ll(\log X)/X$ by
Riemann–von Mangoldt. Once
$F(u)=\sum_{\text{all }\rho}a_\rho e^{i\gamma u}+\varepsilon(u)$ with $\varepsilon(u)\to0$, the
truncation error in (1.6) satisfies $|E(u,X)|\le\sum_{\gamma>X}|a_\rho|+|\varepsilon(u)|$. For
the diagonal $X=e^Y$ the first summand is $o(1)$ **uniformly in $u$**. The second tends to $0$
as $u\to\infty$, hence also in Cesàro $L^2$ mean along $0\le u\le Y$. This is (1.7). The
coupling that makes (1.7) awkward for $M(x)$ never bites: Ng's tail is large on a thin set and
must be averaged, ours is never large.

The content of (1.6) is therefore the expansion with $\varepsilon(u)\to0$, which is the contour
shift. Absolute convergence is used here, and again in §4.4 to drop residues above any
$X\to\infty$; it is not a substitute for the shift.

### 4.2 $H_A$ on a vertical line

In a region free of the poles $z=\rho$ and of $\{1,3,5,\ldots\}$, use (1.8) and, by §2,
$\Gamma(1-z)\sin(\pi z/2)=\pi(2\pi)^{-z}\chi(z)$. Split
$$
M_A(z-1)=-\frac1{(z-1)\zeta(z)}+\sum_{j\ge1}\frac{(-1)^{j+1}A^j}{j!\,(z-1+j)\zeta(z+j)} .
$$

*Leading term:* $\pi(2\pi)^{-z}\chi(z)/((z-1)\zeta(z))$.

*The $j\ge1$ terms.* For $j\ge1$ and $\sigma=\tfrac12-\delta=\tfrac14$, one has
$\Re(z+j)\ge\tfrac54>1$, where $|\zeta(z+j)|\ge1/\zeta(\tfrac54)$ and $\sum A^j/j!\le e^A$. Each
denominator $z-1+j$ is $\asymp|t|$ or larger, so this piece is
$O\bigl(e^A|\chi(z)|/|t|\bigr)$. Here $A=a\log x$ is **fixed**, so $e^A$ is a constant.

*The sum over $k$.* [G] p. 8: it converges absolutely and locally uniformly on
$\mathbb{C}\setminus\{2m+1\}$ and is holomorphic throughout the critical strip. On a vertical
line, $|2k+1-z|\ge c(1+|t|)$ for $k=0$ and $\ge k$ for large $k$, while
$|M_A(2k)|\le\|q_A\|_1$; so the sum is $O(1/|t|)$. Differentiated, $O(1/t^2)$.

*Odd integers.* The apparent poles of (1.8) at $z=1,3,5,\ldots$ cancel ([G] pp. 7–8: the
kernel $\int_1^\infty Y^{-z}\cos(uY)\,dY$ continues to an entire function of $z$). Thus $H_A$ is
holomorphic at $s=1$. The line $\Re z=\tfrac14$ meets none of them, nor does a horizontal at
height $U\to\infty$.

*The derivative.* The three pieces must be differentiated and bounded absolutely: they may
cancel in $H_A$, so no bound relative to $|H_A|$ follows. Differentiating the leading term costs
$\log|t|$ from the phases of $\chi$ and $(2\pi)^{-z}$, plus $\zeta'/\zeta$. Titchmarsh
Theorem 14.2 gives, under RH,
$$
\log\zeta(\sigma+it)=O\bigl((\log t)^{2-2\sigma+\varepsilon}\bigr)
$$
uniformly for $\tfrac12<\sigma_0\le\sigma\le1$. The formula cannot be used for $\sigma>1$:
the exponent $2-2\sigma+\varepsilon$ is then negative, so it would force $\log\zeta=o(1)$,
which is false. For $\sigma>1$ the Euler product gives $|\zeta(\sigma+it)|\ge\zeta(2\sigma)/\zeta(\sigma)>0$
outright, hence $1/\zeta\ll_\sigma 1$ with no $t$ at all (already used on the initial line in
§4.3). In particular $1/\zeta(\sigma+it)\ll t^\varepsilon$
for $\tfrac12+\delta\le\sigma\le1$ by (14.2.6), and for $\sigma>1$ by the Euler product. Cauchy's formula in the disk of radius $\delta/2$
about $\tfrac12+\delta+it$ — empty of zeros under RH — then yields
$\zeta'/\zeta\ll_\delta(\log t)^{1-\delta+\varepsilon}\ll t^\varepsilon$ (the disk extends left
to $\sigma=\tfrac12+\delta/2$, where the exponent in Theorem 14.2 is $1-\delta+\varepsilon$).
On $\sigma\le\tfrac12-\delta$ the functional equation is $1/\zeta(s)=1/(\chi(s)\zeta(1-s))$,
so $1/\zeta$ *gains* a factor $|\chi|^{-1}$ on reflecting to $\Re(1-s)\ge\tfrac12+\delta$; the
same Cauchy bound on $\zeta'/\zeta$ comes back unchanged. Termwise differentiation of the
absolutely convergent $j\ge1$ series gives
$O_A(|\chi(z)|\log|t|/|t|)$; the differentiated $k$-series is $O(1/t^2)$.
Thus on the initial line $\sigma=9/8$, $|H_A'|\ll t^{-13/8}\log t$ and
$|H_A|\ll t^{-1}$; on the left line $\sigma=1/4$, $|H_A'|\ll t^{-3/4}\log t$ and
$|H_A|\ll t^{-3/4}$. These are bounds on each summand followed by the triangle inequality,
not bounds relative to $H_A$. Consequently $G=-H_A'-(\log C)H_A$ is $O(t^{-1})$ on the
initial line and $O(t^{-3/4}\log t)$ on the left line, which is what the contour estimates use.
(Theorem 14.2 is Littlewood's three-circles bound, stated for
$\tfrac12<\sigma_0\le\sigma\le1$; the weaker $t^\varepsilon$ is all the vertical line uses.)

Collecting, on any vertical line in a compact subinterval of $(0,1)$ at positive distance from
the zeros,
$$
H_A(\sigma+it)=\pi(2\pi)^{-\sigma-it}\frac{\chi(\sigma+it)}{(\sigma-1+it)\zeta(\sigma+it)}
  +O_A\Bigl(\frac{|\chi(\sigma+it)|+1}{|t|}\Bigr),
$$
and $G$ obeys the absolute line bounds just established; $\Phi=G/s$ of (3.1) costs one more
$|t|$.

### 4.3 The table, and the second power

Under RH, $1/\zeta(\sigma+it)\ll t^\varepsilon$ for $\sigma\ge\tfrac12+\delta$, by Theorem 14.2
as just cited; the implied constant depends on $\delta$ and $\varepsilon$ and is otherwise
absolute. For $\sigma<\tfrac12$ the identity is $1/\zeta(s)=1/(\chi(s)\zeta(1-s))$ with
$\Re(1-s)>\tfrac12$ — not $\chi(s)/\zeta(1-s)$, which would *cost* $|\chi|$ rather than gain it —
so $|1/\zeta(\sigma+it)|\ll|t|^{\sigma-1/2+\varepsilon}$. The leading term of $H_A$ is
$\pi(2\pi)^{-z}/((z-1)\zeta(1-z))$, hence $\ll t^{-1+\varepsilon}$ on the left. The remainder
$O_A((|\chi|+1)/|t|)$ has two pieces: $|\chi|/t$ from the $j\ge1$ terms, and a floor of order
$1/t$ from the sum over $k$. One integration by parts on the definition of $H_A$ produces the
boundary term $f_A(1)/(1-z)$, and $f_A(1)=\Re\mathrm{ERi}(-A+i)\ne0$, so the floor is sharp:
$|H_A(\sigma+it)|\asymp t^{-1}$ for $\sigma>\tfrac12$. (At $x=1.1$, $a=\tfrac12$,
$|f_A(1)|=0.780\ldots$ and $|\Phi(2+it)|\,t^2\to|\log C\cdot f_A(1)|=0.400\ldots$.) On the left,
$|\chi|/t\asymp t^{-1/2-\sigma}$ sits above that floor. Hence:

| line | $1/\zeta$ | $\lvert H_A(\sigma+it)\rvert$ | $\lvert\Phi(\sigma+it)\rvert$ |
|---|---|---|---|
| $\sigma>\tfrac12$ | $\ll t^\varepsilon$ (RH) | $\ll t^{-1}$ | $\ll t^{-2}$ |
| $\sigma<\tfrac12$ | $\ll t^{\sigma-1/2+\varepsilon}$ | $\ll t^{-1/2-\sigma+\varepsilon}$ | $\ll t^{-3/2-\sigma+\varepsilon}$ |

(The $j\ge1$ remainder, size $|\chi|/t$, dominates the true leading term on the left: on
$\sigma=\tfrac12-\delta$ the leading term is $t^{-1+\varepsilon}$ and the $j\ge1$ piece is
$t^{-1+\delta}$.) At $\sigma_0=\tfrac12-\delta=\tfrac14$,
$$
|\Phi(\sigma_0+it)|\ \ll\ t^{-3/2-\sigma_0+\varepsilon}=t^{-7/4+\varepsilon}.
$$
Take $\varepsilon=\tfrac18$. Then $\int_1^\infty|\Phi(\sigma_0+it)|\,dt<\infty$, the integral a
constant depending on $\delta$, $\varepsilon$, $A$, $x$ only. **The factor $1/s$ in (3.1) is
what buys the second power**; without it the exponent would be $-1+\delta+\varepsilon>-1$ and
the integral would diverge. It is the same $1/\rho$ that makes $\sum|a_\rho|<\infty$ in (4.1).

On the initial line $\sigma_1=\tfrac98>1$, no RH is used: $1/\zeta$ is bounded and
$|\Phi(\tfrac98+it)|\ll t^{-2}$, so $\Phi\in L^1$ there. The line $\sigma=2$ is unavailable:
the same floor would make the tails $Y^2/U$ and the right wing $Y^2U^{-2}$, neither
$o(Y^{1/2})$ at $U=Y^{2/3}$. Both pieces are $o(Y^{1/2})$ for $\sigma_1<7/6$.

### 4.4 The shift

The trivial bound $P_A(Y)\ll Y\log Y$ ($f_A$ bounded) gives absolute convergence of
$\Phi(s)=\int_1^\infty P_A(Y)\,Y^{-s-1}\,dY$ only for $\Re s>1$. Identification with $G/s$ is
made there ([G], proof of Prop. 1.13). Since $P_A$ is $C^1$ and $\Phi\in L^1$ on $\Re s=\tfrac98$,
Mellin inversion holds:
$$
\tag{4.2}
P_A(Y)=\frac1{2\pi i}\int_{\tfrac98-i\infty}^{\tfrac98+i\infty}\Phi(s)\,Y^s\,ds .
$$
(The abscissa of $\Phi$ is not known *a priori* to be $\tfrac12$; that would be circular.)

Truncate at height $U$, shift the rectangle with corners $\tfrac98\pm iU$ and $\sigma_0\pm iU$, and
take $U$ as in Titchmarsh Theorem 14.16: every interval $[V,V+1]$ of large $V$ contains a $U$
with
$$
\bigl|\zeta(\sigma+iU)\bigr|>\exp\Bigl(-\kappa\frac{\log U}{\log\log U}\Bigr)
\qquad\bigl(\tfrac12\le\sigma\le2\bigr),
$$
for an absolute constant $\kappa$ — Titchmarsh writes $A$, which we rename here, since $A=a\log x$
is already the smoothing parameter of $f_A$, $M_A$, $H_A$, $P_A$ throughout this paper —
hence $1/\zeta(\sigma+iU)\ll U^\varepsilon$ uniformly on that segment, and
$\zeta(\sigma+iU)\ne0$, so the horizontal does not hit a pole. Bounding $\Phi=G/s$ on the same
segment also uses $\zeta'/\zeta$, because $G=-H_A'-(\log C)H_A$. Set
$\eta=c/\log\log U$ with fixed $c>0$ small enough for disks of radius $\eta$ centred in
$|\sigma-\tfrac12|\le\eta$ to remain inside the band of Titchmarsh (14.14.2).
The upper bound in Theorems 14.14(A) and (B), the functional equation on the left, and
Cauchy's formula on those disks give $\zeta'(s)\ll U^{o(1)}$ in this band. At the selected
height, Theorem 14.16 gives $1/\zeta(s)\ll U^{o(1)}$ on the right half of the horizontal;
reflection gives the corresponding left-half bound. Thus $\zeta'/\zeta\ll U^{o(1)}$
in the band. Away from the critical line, use the fixed-distance bounds of §4.2;
on the intervening short strip, (14.14.4) and Cauchy's formula give the same bound.
The upper bounds are applied on neighbourhoods before dividing by the lower bound on
the chosen horizontal; a lower bound on that segment alone does not justify Cauchy estimates.
Choose $V=Y^{2/3}$, so
$U=Y^{2/3}+O(1)$. ([G] Lemma 2.5 still supplies the good ordinates for $(\dagger)$ in §3; it is
a different integral.) Cauchy's theorem gives
$$
P_A(Y)
  =\sum_{|\gamma|<U}\operatorname{Res}_{s=\rho}\bigl(\Phi(s)Y^s\bigr)
   +V(Y)+H_\pm(Y)+T_2(Y),
$$
the four terms being residues, the left vertical, the two horizontals, and the tails of (4.2)
above height $U$.

*Poles inside the rectangle.* Under RH the only poles of $H_A$ with $\Re s>0$ are the
non-trivial zeros. $H_A$ is holomorphic at $s=1$ (§4.2). $\Phi=G/s$ has a possible simple pole
at $s=0$, which lies to the left of $\sigma_0=\tfrac14$ and is not picked up. At a simple zero,
$G(s)=r_\rho(s-\rho)^{-2}+b_\rho(s-\rho)^{-1}+O(1)$, so
$$
\Phi(s)=\frac{r_\rho}{\rho}\,(s-\rho)^{-2}+\Bigl(\frac{b_\rho}{\rho}-\frac{r_\rho}{\rho^2}\Bigr)(s-\rho)^{-1}+O(1).
$$
Thus $c_2=r_\rho/\rho$, and
$$
\operatorname{Res}_{s=\rho}\bigl(\Phi(s)Y^s\bigr)=Y^\rho\bigl(c_2\log Y+c_1\bigr),
\qquad
c_1=\frac{b_\rho}{\rho}-\frac{r_\rho}{\rho^2}.
$$
Here $b_\rho=-(\log C)\,r_\rho$ exactly — the simple part of $G=-H_A'-(\log C)H_A$ at $\rho$ comes
only from the second term — so $c_1=-r_\rho\bigl(\log C/\rho+1/\rho^2\bigr)$ and
$\sum|c_1|\le(1+|\log C|)\sum|r_\rho/\rho|<\infty$ under WMC by (4.1). The $c_1$ terms are then
$O(Y^{1/2})$ and die after dividing by $Y^{1/2}\log Y$.

*Tails $T_2$.* $|\Phi(\tfrac98+it)|\ll t^{-2}$, so
$$
|T_2(Y)|\ll Y^{9/8}\int_U^\infty t^{-2}\,dt\ll Y^{9/8}U^{-1}.
$$
With $U=Y^{2/3}$ this is $Y^{9/8-2/3}=Y^{11/24}=o(Y^{1/2})$.

*Left vertical $V(Y)$.* The integrand is $\Phi(\sigma_0+it)Y^{\sigma_0+it}$. By the table and
the choice $\varepsilon=\tfrac18$,
$$
|V(Y)|\ll Y^{1/4}\int_{-\infty}^\infty\bigl(1+|t|\bigr)^{-7/4+1/8}\,dt\ll Y^{1/4}.
$$
After dividing by $Y^{1/2}\log Y$ this is $O(Y^{-1/4}/\log Y)\to0$, for every $Y$, with no
restriction to special ordinates. The implied constant depends on $\delta=1/4$ and is otherwise
absolute: $\delta$ is fixed.

*Dropped residues $\gamma>U$.* Each is $\ll Y^{1/2}\log Y\cdot|a_\rho|$, so the sum over
$\gamma>U=Y^{2/3}$ is $\ll Y^{1/2}\log Y\cdot\sum_{\gamma>Y^{2/3}}|a_\rho|
\ll Y^{1/6}(\log Y)^{3/2}=o(Y^{1/2})$ as $Y\to\infty$, uniformly, by (4.1).

*Horizontals.* The contribution of one side is
$\ll\int_{\sigma_0}^{9/8}|\Phi(\sigma+iU)|\,Y^\sigma\,d\sigma$. Do not obtain $1/\zeta\ll t^\varepsilon$
on the crossing by integrating $\zeta'/\zeta$ from the right-hand edge: the unconditional
$\zeta'/\zeta\ll\log^2t$ gives $\exp(C\log^2t)$, which is not $\ll t^\varepsilon$ (at $t=10^6$ it
is $10^{82.9}$ against $t^{0.1}=10^{0.6}$), and even an RH bound $\zeta'/\zeta\ll t^\varepsilon$
across a segment of fixed length produces a positive power of $t$. Theorem 14.16 supplies the
bound at our chosen $U$, uniformly for $\tfrac12\le\sigma\le2$. Do not treat
$Y^\sigma|\Phi(\sigma+iU)|$ as dominated by its left endpoint either: $Y^\sigma$ grows with
$\sigma$, so if $U=o(\sqrt Y)$ the maximum on the left half sits at $\sigma=\tfrac12$, not at
$\sigma_0$. Split at $\eta=c/\log\log U$ as above, retaining fixed constants in the
Titchmarsh band estimates.

(i) *Right wing,* $\sigma\in[\tfrac12+\eta,\tfrac98]$. On any compact subinterval
$\tfrac12+\delta_0\le\sigma\le1$ with $\delta_0>0$ fixed (say $\delta_0=\tfrac18$), Theorem 14.2
gives $1/\zeta\ll_{\delta_0}U^\varepsilon$; its implied constant depends on $\delta_0$ and must
not be asked to survive $\delta_0\to0$. On the remaining stretch $1<\sigma\le\tfrac98$ the Euler
product gives $1/\zeta\ll 1$ with no $t$. The
remaining short strip down to $\tfrac12+\eta$ is controlled by the upper bounds of 14.14,
Cauchy's formula, and the lower bound of 14.16 at the selected height. The $1/t$ floor
dominates $|\chi|/t$ throughout this wing, so
$|\Phi(\sigma+iU)|\ll U^{-2+\varepsilon}$. The product $Y^\sigma U^{-2+\varepsilon}$ is maximised at
$\sigma=\tfrac98$ since $Y>U$, hence $\ll Y^{9/8}U^{-2+\varepsilon}$. With $U=Y^{2/3}$ this is
$Y^{9/8-4/3+o(1)}=Y^{-5/24+o(1)}$, which is $o(Y^{1/2})$.

(ii) *Left wing,* $\sigma\in[\sigma_0,\tfrac12-\eta]=[\tfrac14,\tfrac12-\eta]$. The functional
equation reduces $1/\zeta(s)$ to $1/(\chi(s)\zeta(1-s))$ with $\Re(1-s)\in[\tfrac12+\eta,\tfrac34]$,
where the right-half estimates apply — Theorem 14.2 on the compact piece and the
selected-height bound on the short strip — *gaining* a factor $|\chi(s)|^{-1}$. The remainder in $H_A$
still carries $|\chi|/|t|$, which dominates, and the table gives
$|\Phi|\ll U^{-3/2-\sigma+\varepsilon}$. The product
$Y^\sigma U^{-3/2-\sigma+\varepsilon}=U^{-3/2+\varepsilon}(Y/U)^\sigma$ has $Y/U=Y^{1/3}>1$,
so is maximised at the right endpoint $\sigma=\tfrac12-\eta$:
$\ll Y^{1/2}U^{-2+\varepsilon}=Y^{1/2-4/3+o(1)}=Y^{-5/6+o(1)}$. (At the left
endpoint the same product is $\ll Y^{1/4}U^{-7/4+\varepsilon}=Y^{1/4-(7/4)(2/3)+o(1)}=Y^{-11/12+o(1)}$,
smaller.)

(iii) *The band,* $|\sigma-\tfrac12|\le\eta$. Theorem 14.16 already includes $\sigma=\tfrac12$:
$1/\zeta(\sigma+iU)\ll U^\varepsilon$ uniformly for $\tfrac12\le\sigma\le\tfrac12+\eta$. To the
left of the line, $1/\zeta(s)=1/(\chi(s)\zeta(1-s))$ and $|\chi|\ll U^\eta=U^{o(1)}$, so
$1/\zeta\ll U^{o(1)}$ as well; equivalently the leading term of $H_A$ is $1/((z-1)\zeta(1-z))$.
Thus $|\Phi(\sigma+iU)|\ll U^{-2+\varepsilon}$ on the band, $Y^\sigma=Y^{1/2}U^{o(1)}$, and the
contribution is $\ll Y^{1/2}U^{-2+\varepsilon}/\log\log U=Y^{-5/6+o(1)}$ at $U=Y^{2/3}$.

All three pieces are $o(Y^{1/2}\log Y)$. Taking $U\to\infty$ independently of $Y$, as in [G]
Lemma 2.5, is for a different integral (his $H_{a,c}(T)$, length $O(1)$, integrand
$\mathrm{ERi}-\zeta'/\zeta$) and does not apply here: the right-wing term $Y^{\sigma_1}U^{-2}$ grows
with $Y$ unless $U$ grows with $Y$. The choice $U=Y^{2/3}$ is large enough for the tails and
the right wing at $\sigma_1=\tfrac98$, and small enough that $U<Y$.

Collecting: for every $Y>2$,
$$
\tag{4.3}
P_A(Y)
  =\sum_{|\gamma|<Y^{2/3}}Y^\rho\bigl(c_2\log Y+c_1\bigr)
   +o\bigl(Y^{1/2}\log Y\bigr),
$$
the $o$ as $Y\to\infty$, with no restriction on $Y$ beyond largeness. Dividing by
$Y^{1/2}\log Y$ and using $c_2=r_\rho/\rho$, $Y^\rho=Y^{1/2}Y^{i\gamma}$ under RH, gives
$$
\frac{P_A(Y)}{Y^{1/2}\log Y}
  =\sum_{\gamma>0}\frac{r_\rho}{\rho}\,Y^{i\gamma}+\frac{\overline{r_\rho}}{\overline{\rho}}\,Y^{-i\gamma}
   +o(1).
$$

### 4.5 From good ordinates to all $T$

Two contours, two notions of good ordinate, not to be conflated.

*[G]'s rectangle,* identity $(\dagger)$. This holds along good $T'\in[T,T+1]$, with $o(1)$ as
$T'\to\infty$. For a general $T$, two errors appear. First, $P_A(T\log x)-P_A(T'\log x)$: the
interval of integration has length $O(1)$ and the integrand is $O(\log T)$, so the difference
is $O(\log T)$. Second, $\Sigma R_T-\Sigma R_{T'}$ is a sum over the $O(\log T)$ zeros in
$[T,T+1]$, each $\mathrm{ERi}(\rho\log x)=O_x(1)$ ($x$ fixed). Both are $O(\log T)$ unnormalised,
hence $O(T^{-1/2})$ after dividing by $\sqrt T\log T$. The $o(1)$ in $(\dagger)$ itself, divided
by $\sqrt T\log T$, is smaller still. Thus $(\dagger)$ transfers to every large $T$.

*The Mellin rectangle* of §4.4. Identity (4.3) already holds for every large $Y$, because $P_A$
is $C^1$ and the height $U$ was chosen by Theorem 14.16 in $[Y^{2/3},Y^{2/3}+1]$, independently
of whether $Y/\log x$ is a good ordinate for $(\dagger)$.

Set $Y=T\log x$. Then $Y^{1/2}=\sqrt T\sqrt{\log x}$ and
$\log Y=\log(T/2\pi)+\log C$, so $\log Y/\log(T/2\pi)\to1$. The constant $\log C$ produces a
further $Y^\rho$ term, $O(Y^{1/2})$ after summing, which dies on dividing by
$\sqrt T\log(T/2\pi)$. Combining with $(\dagger)$,
$$
\tag{4.4}
\frac{\Sigma R_T-L}{\sqrt T\,\log(T/2\pi)}
  =\sum_{\gamma>0}a_\rho e^{i\gamma\log T}+\overline{a_\rho}\,e^{-i\gamma\log T}+\varepsilon(\log T),
$$
where
$$
a_\rho=\frac{r_\rho}{\rho}\cdot\frac{(\log x)^{i\gamma}}{\pi\sqrt{\log x}},
\qquad
|a_\rho|=\frac{|r_\rho|}{|\rho|\,\pi\sqrt{\log x}},
$$
and $\varepsilon(u)\to0$ as $u\to\infty$, for all real $u$, not merely along a sequence. This
is the expansion that §4.1 feeds to [ANS] (1.6). Truncating the series at any $X$ costs
$O(\sum_{\gamma>X}|a_\rho|)$ uniformly in $T$, by absolute convergence.

### 4.6 Consequences

Along every large $T$,
$$
\frac{\Sigma R_T-L}{\sqrt T\,\log(T/2\pi)}
  =2\Re\sum_{0<\gamma\le X}a_\rho e^{i\gamma\log T}
   +O\Bigl(\sum_{\gamma>X}2|a_\rho|\Bigr)+\varepsilon(\log T),
$$
the first error $o(1)$ as $X\to\infty$ uniformly in $T$, and $\varepsilon(u)\to0$ as
$u\to\infty$ — both under WMC. (The contour bounds of §4.2–4.4 use RH; the residue tail (4.1)
uses WMC; WMC is in force throughout §4 and implies RH.) This is [ANS] (1.6); (1.7) with
$X=e^Y$ follows as in §4.1. Their density condition (1.11),
$\sum_{T<\lambda_n\le T+1}1\ll\log T$, is Riemann–von Mangoldt, unconditionally. The
amplitude condition (1.12) holds with $\theta=0$, since
$\sum_{\gamma\le T}\gamma^2|a_\rho|^2\asymp\sum_{\gamma\le T}1/|\rho\zeta'(\rho)|^2=O(1)$,
comfortably inside the admissible range $\theta<3-\sqrt3$. [ANS] Corollary 1.3(b) gives that
$F$ is $B^2$-almost periodic (hence has a limiting distribution); Theorem 1.9, under LI, gives
the Bessel product. Absolute convergence is stronger than $B^2$: the series
$\phi(u)=\sum_{\gamma>0}2\Re\bigl(a_\rho e^{i\gamma u}\bigr)$ is Bohr uniformly almost periodic,
so $F(u)=\phi(u)+\varepsilon(u)$ with $\varepsilon(u)\to0$. A $o(1)$ perturbation of a Bohr
function need not itself be Bohr, but it has the same limiting distribution in the mean
$\lim_{U\to\infty}U^{-1}\int_0^U g(F(u))\,du$ for bounded continuous $g$; its support is
contained in $[-S,S]$, where $S=\sum_{\gamma>0}2|a_\rho|$. Under LI, finite-dimensional
Kronecker equidistribution together with the uniform tail makes the support the full interval
$[-S,S]$, gives the Bessel product, and makes $\limsup|F|=S$:

> **Theorem 2.** Assume WMC and LI. Then
> $F(u)=(\Sigma R_{e^u}-L)/(e^{u/2}\log(e^u/2\pi))$ equals a Bohr uniformly almost periodic
> function of $u$ plus $o(1)$; it possesses a limiting distribution $\mu$ with compact support
> and characteristic function $\hat\mu(\xi)=\prod_{\gamma>0}J_0(2|a_\rho|\xi)$; and
> $$
> \limsup_{T\to\infty}\frac{|\Sigma R_T(x)-L|}{\sqrt T\,\log(T/2\pi)}
>   =\sum_{\gamma>0}2|a_\rho|
>   =\frac1{\pi\sqrt{\log x}}\sum_{\gamma>0}\frac{2|r_\rho|}{|\rho|},
> \]
> an equality. At $x=1.1$ the first sixty zeros contribute $0.039108$ to the right-hand side;
> the tail is positive and no numerical upper bound for it is claimed.

**What this does not give.** No Lindelöf-free version: RH (via WMC) is used for $1/\zeta$ on
$\sigma>\tfrac12$. Off RH the contour would stop at $\Theta-\delta$ and the remainder would be
$Y^{\Theta-\delta}$, which is the trichotomy of §3 rather than a limiting distribution. The tail
(4.1) does carry a rate, $\ll\sqrt{\log X/X}$, and that is what turns the dropped residues of
§4.4 into $O(Y^{1/6}(\log Y)^{3/2})$.

**Numerics.** Four hundred residues against the $100\,000$-zero sum: absolute error $O(1)$ and
flat while $\Sigma R_T-L$ grows to $30$; after normalisation, $10^{-4}$ to $10^{-5}$, three
orders below the signal. The pointwise tail the argument predicts, with no averaging
(Appendix B).

---

## 5. Riesz means of the zero-sum

The obvious repair for a divergent partial sum is to replace the sharp cutoff by a Riesz weight
of order $k>0$,
$$
\Sigma R_T^{(k)}(x):=\sum_{0<|\Im\rho|\le T}\Bigl(1-\frac{|\gamma|}T\Bigr)^k m_\rho\,
\mathrm{ERi}(\rho\log x)
$$
($k=0$ is the original sum, $k=1$ is Fejér). One expects, by a careless reading of Hardy–Riesz,
that these converge as soon as $k>\Theta$, with error $\asymp T^{\Theta-k}\log T$. They do not.
The oscillation is logarithmic, its instantaneous frequency at height $T$ is $\gamma/T$, and a
kernel of width $T$ therefore suppresses each mode by a factor independent of $T$. The weight
reduces each residue amplitude by a computable factor and changes nothing else at the
$\sqrt T\log T$ scale.

The contour centre is $D_{a,x}$; identifying it with the numerical candidate $L$ remains
unproved. Every normalised result below is unchanged by replacing either by another fixed
constant. No Riesz sum of fixed order converges.

### 5.1 Why Hardy–Riesz does not apply

A Riesz kernel of width $T$ suppresses a frequency $\nu$ by $\asymp(\nu T)^{-k}$. By §4 the
model is almost periodic in $u=\log T$, and at $u=\log T$ the instantaneous frequency of the mode
$\gamma$ is $\gamma/T$, so $\nu T=\gamma$ is **constant in $T$**. The kernel widens exactly as
fast as the frequency shrinks.

The textbook contrast is $\sum_{n\ge1}(-1)^{n-1}n^a$: its partial sums grow like $n^a$ and *are*
$(C,k)$-summable for $k>a$, because the alternation sits at frequency $\pi$ in the type variable
— fixed — so the suppression $\omega^{-k}$ genuinely grows. Logarithmic oscillation is a
different regime.

In Hardy–Riesz's own terms, with $A^K(\omega)=\sum_{\gamma<\omega}c_\gamma(\omega-\gamma)^K$, the
model suggests $A^K(\omega)\asymp\omega^K(L+\omega^{1/2}\cdot\mathrm{osc})$, and Theorem 31 of the
tract (p. 45) would give
$$
\sigma_K=\limsup_{\omega\to\infty}\frac{\log|A^K(\omega)|}{\log\omega}-K=\tfrac12
\qquad\text{for every }K .
$$
This is a model comparison, not a second proof. Theorem 3 below is the statement for Riesz means
of each *fixed* finite order: the double pole
survives every such weight. Logarithmic-type means fail for the same reason (suppression
$(\gamma\log T)^{-k}$, still beaten by $\sqrt T$), on the residue model; Abel summation of a
term $T^\rho$ with $\Re\rho>0$ blows up like $\varepsilon^{-\rho}$, likewise on the model, not
by a Laplace-transform twin of Lemma A. (The same shape $\sum_\rho x^\rho/(\rho\zeta'(\rho))$
for $\sum_{n\le x}\mu(n)$ is why the literature calls $\sum\mu(n)=1/\zeta(0)=-2$
zeta-regularized rather than Cesàro-summed.)

A Fejér-weighted $\sum\mathrm{ERi}(\rho\log x)$ is also not a Weil–Guinand identity. The
triangular function $h(t)=(1-|t|/T)_+$ is a standard admissible test function, but it weights
$\sum_\rho\Phi(\rho)$, not the Riemann–von Mangoldt summands $\mathrm{ERi}(\rho\log x)$.
Different object.

### 5.2 The weighted explicit formula

The unconditional argument uses the primitive $P_A$, not the conditional expansion of §4.
For fixed $k>0$, Stieltjes integration by parts gives
$$
\Sigma R_T^{(k)}=\frac kT\int_0^T\Bigl(1-\frac uT\Bigr)^{k-1}\Sigma R_u\,du.
$$
The kernel is integrable even when $0<k<1$. The bridge $(\dagger)$ extends from good ordinates
to every large $u$ with error $O_x(\log u)$: $P_A$ changes by $O(\log u)$ across a unit interval,
and the $O(\log u)$ zeros there, counted with multiplicity, each contribute $O_x(1)$ by [G]
Prop. 1.4 and the Riemann--von Mangoldt zero count (Titchmarsh Thm. 9.2). Put
$$
(\mathcal R_kP_A)(Y):=\frac{k}{Y}\int_1^Y(1-v/Y)^{k-1}P_A(v)\,dv .
$$
The bounded initial interval contributes $O(1/T)$. Thus

> **Proposition.** For each fixed $x>1$ and real $k>0$, unconditionally,
> $$
> \Sigma R_T^{(k)}-D_{a,x}
> =\frac1{\pi\log x}(\mathcal R_kP_A)(T\log x)+O_x(\log T).
> \]
> In the half-plane of absolute convergence, the Mellin transform of $\mathcal R_kP_A$ is
> $m_k(s)\Phi(s)$, where
> $m_k(s)=kB(s+1,k)=\Gamma(k+1)\Gamma(s+1)/\Gamma(s+1+k)$.

The identity follows by Fubini from $P_A(Y)=O(Y\log Y)$ and the Beta integral, and continues
meromorphically. At $\rho_1$, $m_k$ is holomorphic and nonzero, so the double pole of $\Phi$
retains leading coefficient $(r_{\rho_1}/\rho_1)m_k(\rho_1)$. The multiplier introduces no real
pole with real part at least $1/2$. Lemma A therefore applies to $\mathcal R_kP_A$, and the
$O(\log T)$ transfer error is negligible at the $\sqrt T\log T$ scale. The case $k=0$ is
Theorem 1; the literal beta integral is only for $k>0$, while its gamma-ratio multiplier extends
to $m_0(s)=1$.

The §4 residue expansion gives a *conditional* weighted model under WMC, but it is not used in
this proof. A numerical fit suggests a smaller transient $k\tau(x)/T$ with
$\tau(1.1)\approx510$; its coefficient has not been derived. An error $o(\sqrt T)$ cannot
identify any fixed $T^{-1}$ coefficient. The fitted $\tau$ is distinct from the regularised
trivial-zero value $I(x)$.

With $\rho=\tfrac12+i\gamma$, the amplitude multiplier is the gamma ratio $m_k(\rho)$; for
$k>0$ it equals $kB(\rho+1,k)$. Thus
$$
W_k(\gamma):=\Bigl|\frac{\Gamma(k+1)\Gamma(\tfrac32+i\gamma)}{\Gamma(\tfrac32+k+i\gamma)}\Bigr|
  \ \sim\ \frac{k!}{\gamma^k}\qquad(\gamma\to\infty),
$$
with $W_0\equiv1$ and $W_k(\gamma)>0$ for every finite $k$ and every zero — $\Gamma$ has no
zeros, and $\tfrac32+k+i\gamma$ is never a non-positive integer. So $k\,B(\rho_1+1,k)$ is a
non-zero holomorphic multiplier at the pole $\rho_1$ of the Mellin transform of $P_A$, the double
pole of §3.2 survives the weighting, and Lemma A applies verbatim to the weighted primitive:

> **Theorem 3.** For every fixed $k\ge0$ and every fixed $x>1$,
> $\Sigma R_T^{(k)}(x)=\Omega_\pm\bigl(\sqrt T\log T\bigr)$, unconditionally, with constant
> $0.008144\cdot W_k(\gamma_1)$ at $x=1.1$.

Under WMC, multiply the inverse Mellin integrand of §4.4 by $m_k(s)$ for $k>0$. Stirling gives
$|m_k(\sigma+it)|\ll_k(1+|t|)^{-k}$ in the contour strip, so the line and horizontal bounds
remain integrable or improve. The residues have leading coefficients
$(r_\rho/\rho)m_k(\rho)$ and an absolutely summable lower-order series. The all-$T$ bridge in
the Proposition then yields an absolutely convergent Bohr series with amplitudes multiplied by
$W_k(\gamma)$. For $k=0$, use Theorem 2. With LI, Kronecker applies and the $\limsup$ is an
equality:

> **Theorem 4.** Assume WMC and LI. For every fixed $k\ge0$,
> $$
> \limsup_{T\to\infty}\frac{\bigl|\Sigma R_T^{(k)}(x)-L\bigr|}{\sqrt T\,\log(T/2\pi)}
>   =\frac1{\pi\sqrt{\log x}}\sum_{\gamma>0}\frac{2|r_\rho|}{|\rho|}\,W_k(\gamma)\ >\ 0 ,
> \]
> independent of $T$.

At $x=1.1$, from the first 60 residues, and with the unconditional constants of Theorem 3
alongside:

| $k$ | first-sixty weighted sum | ratio of sixty-term sums | $W_k(\gamma_1)$ | uncond. constant |
|---|---|---|---|---|
| 0 | 0.0379306 | 1 | 1 | $8.144\cdot10^{-3}$ |
| 1 | 0.00175808 | 21.6 | 0.07035271 | $5.730\cdot10^{-4}$ |
| 2 | 0.00020064 | 189.0 | 0.00980245 | $7.984\cdot10^{-5}$ |
| 3 | 0.00003728 | 1017.5 | 0.00201951 | $1.645\cdot10^{-5}$ |

A factor of 189 at $k=2$ is a genuine numerical improvement. It is a bounded improvement, not
convergence.

### 5.3 Why $T=10^4$ looked like convergence

The transient $k\tau/T$ dominates the growing oscillation until the two are comparable; at $x=1.1$
the crossover is $T\approx1.4\cdot10^3$ ($k=1$), $7.8\cdot10^3$ ($k=2$), $2.8\cdot10^4$ ($k=3$).
At $T=10^4$ the $k=3$ discrepancy $0.14$ is the transient dying — its oscillatory part is still
only $0.03$. Three checkpoints at $T=2048$, $4096$, $9999$ hide this: $T=2048$ for $k=1$ lands
near a zero of the oscillation ($|\cdot-L|=0.058$), while the per-octave maxima of
$|\Sigma R_T^{(1)}-L|$ are already flat at $0.58$–$0.87$ across $T\in[10^3,10^4]$.

### 5.4 The falsifier

Möbius–Ei evaluation of $\mathrm{ERi}(\rho\log1.1)$ for the first $10^5$ Odlyzko zeros
($\gamma\le74920.8$) in double precision, with maximum absolute error $4.4\cdot10^{-13}$ against
Gram's series at 460 digits on the first 10142 zeros (Appendix B). Per-octave maxima of
$|\Sigma R_T^{(2)}-L|$:

| octave | observed | aligned envelope (§5.2) | refuted $T^{-1/2}\log T$ |
|---|---|---|---|
| $9999$–$19998$ | 0.2333 | 0.2360 | 0.1180 |
| $19998$–$39996$ | 0.3504 | 0.3624 | 0.0906 |
| $39996$–$74920$ | 0.3511 | 0.5315 | 0.0709 |

Observed **grows** by $50\%$; the refuted law predicts a **$40\%$ decrease**. The two are a
factor of five apart by the last octave. Ratios to the aligned bound are $0.99$, $0.97$, $0.66$,
all below $1$, with the finite-window undershoot of §4 opening as the envelope grows.
Corroboration: $k=1$ reaches $|\Sigma R^{(1)}-L|=2.18$ at $T=74920$ against $0.65$ at $T=10^4$;
the raw sum reaches $+21.17$.

### 5.5 What the weight is for

To recover $L$ from a truncated zero-sum one subtracts the residue main term of §4 — sixty
residues track the 10142-zero sum to $0.3$–$3\%$ with no drift. No weighting will do it. The
Riesz weight remains useful as a **preconditioner**: $k=2$ cuts the envelope by 189, which is the
difference between an $O(1)$ and an $O(10^{-2})$ error at heights where zeros are routinely used
— provided one remembers that the error then grows like $\sqrt T\log T$ from that smaller
baseline.

### 5.6 The trivial-zero series, at the same resolution

[G] Prop. 4.1 shows that $\Sigma R^{\mathrm{triv}}_T(x):=\sum_{k\le T}\mathrm{ERi}(-2k\log x)$ is
not $O(T^\theta)$ for any $\theta<\Theta$, by the mechanism of his main theorem one level down. Its
Mellin transform is $D_x(s)/s$ with $D_x(s)=\sum_k\mathrm{ERi}(-2k\log x)\,k^{-s}$ ([G] (4.4)), and
his (4.3),
$$
D_x(s)=\frac{\Gamma(1-s)\,(2\log x)^{s-1}}{(s-1)\,\zeta(s)}+H_{x,\delta}(s),
$$
$H_{x,\delta}$ holomorphic on $\mathbb C\setminus\{1,2,3,\dots\}$, gives $D_x$ a pole of order
$m_\rho$ at every non-trivial zero $\rho$ — at $s=\rho$ itself, from $1/\zeta$, with a bare
$\Gamma(1-\rho)$ in the residue and no $P_A$ behind it. Lemma A applies as it stands:
$\Sigma R^{\mathrm{triv}}_T(x)$ is real and locally bounded; $\mathrm{ERi}(-2k\log x)\to0$, so the
partial sums are $o(T)$ and $(s-1)D_x(s)\to0$ as $s\downarrow1$ on the real axis, hence $D_x$ has
residue zero at $s=1$ and there is no pole to subtract. On $[\tfrac12,1)$, $\zeta$ has no real
zero, and beyond $1$ the Dirichlet integral is holomorphic. Division by $s$ introduces no pole
with real part at least $\tfrac12$, so the real main term in Lemma A is zero;
and at $\rho_1$ the pole is simple, with leading coefficient
$c_1=\Gamma(1-\rho_1)(2\log x)^{\rho_1-1}/\bigl((\rho_1-1)\zeta'(\rho_1)\rho_1\bigr)$. So
$$
\Sigma R^{\mathrm{triv}}_T(x)=\Omega_\pm(\sqrt T),\qquad
\limsup_{T\to\infty}\frac{\Sigma R^{\mathrm{triv}}_T(x)}{\sqrt T}\ge|c_1|,
$$
with no logarithm, the pole being simple. At $x=1.1$, $|c_1|=8.241\cdot10^{-12}$;
$2|c_1|=1.65\cdot10^{-11}$ is a heuristic conjugate-pair amplitude, not the one-pole lower
bound (numerical evaluation of the coefficient above). The estimate $T\approx4\cdot10^{21}$ for a
unit excursion from $I(x)$ assumes an unproved centering identity and is only illustrative.
The Beta factor of
§5.2 multiplies $c_1$ by $k\,B(\rho_1+1,k)\ne0$ at every Riesz order, so no Cesàro or Riesz mean
sums $\sum_kR(x^{-2k})$ either — which answers the recollection in [MO] that Cesàro summation
succeeds for the trivial zeros: it does not, but the oscillation it fails on is eleven orders of
magnitude below the one from the non-trivial zeros, and no computation would see it. The two
halves of the folklore identity diverge for one reason at two scales.

**The value $I(x)$ is $D_x(0)$.** At $s=0$ the first term of (4.3) is $1/\log x$, and [G]'s
remainder,
$$
H_{x,\delta}(s)=\Gamma(1-s)(2\log x)^{s-1}\int_\delta^1Q(u)u^{s-1}\,du-\int_0^\delta Q(u)\,G(s,2u\log x)\,du
-\int_\delta^1Q(u)\,\mathrm{Li}_s(e^{-2u\log x})\,du,
$$
$G(s,v)=\sum_{j\ge0}(-1)^j\zeta(s-j)v^j/j!$ and $2\delta\log x<2\pi$, collapses there:
$\mathrm{Li}_0(e^{-v})=1/(e^v-1)$ and $G(0,v)=1/(e^v-1)-1/v$, so the three terms are one,
$$
H_{x,\delta}(0)=-\int_0^1Q(u)\Bigl(\frac1{e^{2u\log x}-1}-\frac1{2u\log x}\Bigr)du,
$$
independent of $\delta$, integrable because the bracket tends to $-\tfrac12$ at $u=0$ and
$Q\in L^1$, and valid for every $x>1$. (The polylogarithm integral on its own is not integrable at
$u=0$ when $s=0$; the split is what makes the value exist.) For $1<x<e^\pi$ the series for
$G(0,2u\log x)$ converges uniformly on $[0,1]$ — [G]'s argument with $v_0=2\log x<2\pi$ — and
termwise integration against $Q$, with [G] Lemma 1.2, $\int_0^1Q(u)u^{s-1}du=-1/((s-1)\zeta(s))$
at $s=j+1$, and $\int_0^1Q=-\mathrm{ERi}(0)=-1$, gives
$$
H_{x,\delta}(0)=-\frac12+\sum_{j\ge1}\frac{(-1)^j\,\zeta(-j)\,(2\log x)^j}{j\cdot j!\,\zeta(j+1)} .
$$
Only odd $j=2k-1$ survive, since $\zeta(-2k)=0$; with $\zeta(1-2k)=-B_{2k}/2k$ and Euler's
$\zeta(2k)=(-1)^{k+1}B_{2k}(2\pi)^{2k}/(2\,(2k)!)$ the $k$-th term is
$\frac1\pi\frac{(-1)^{k+1}}{2k-1}\bigl(\frac{\log x}\pi\bigr)^{2k-1}$, so
$$
H_{x,\delta}(0)=-\frac12+\frac1\pi\arctan\frac{\log x}\pi=-\frac1\pi\arctan\frac\pi{\log x},
\qquad
D_x(0)=\frac1{\log x}-\frac1\pi\arctan\frac\pi{\log x}=I(x),
$$
for $1<x<e^\pi$, and for all $x>1$ because both sides are real-analytic in $x$ there. So the
value $I(x)$ is the value at $s=0$ of the continuation of their Dirichlet series. This is a
regularisation, not a proved centre of the partial sums. The series itself still diverges, as
above, with a one-pole $\Omega_\pm(\sqrt T)$ constant $8.241\cdot10^{-12}$ at $x=1.1$.
The evaluation is Euler on
the even zeta values against the Mellin transform of $Q$; it uses [G] (4.3) and Lemma 1.2 and
nothing from §§2–5, and it is not the residue evaluation of Paper A §2.3.

---

## Appendix A. The Pintz–Révész $\pi/2$ does not transfer

Can the Landau constant $0.008144$ of Theorem 1 be replaced by
$(\pi/2)|c_2|/(\pi\sqrt{\log x})=0.012793$, still unconditional from $\rho_1$? **No.** The extra
$\log Y$ is not the obstruction. Two other things are, and a transplanted kernel argument gives at
most $21\%$, to $0.00989$, not $57\%$.

**Where $\pi/2$ comes from.** Révész (arXiv:2202.01837v3, Theorem 1 = Pintz, Acta Arith. 36
(1980), Thm 1) evaluates a Gaussian-weighted integral $S=U(\rho_0)+U(\overline{\rho_0})$ two
ways, with $D(s)=-\zeta'/\zeta(s)-s/(s-1)$ and
$U(w)=\frac1{2\pi i}\int_{(2)}D(s+w)e^{ms^2+Ms}\,ds$. Residues of $D$ at a simple zero are $-1$,
real. The upper bound pulls out $K$ from $|\Delta(x)|\le Kx^{\beta_0}$ and uses
$\frac1{\sqrt\pi}\int|\cos(Py+R)|e^{-y^2}dy\to2/\pi$, the mean of $|\cos|$ — this is the Gibbs
factor, and it has nothing to do with Mellin pole order. The lower bound shifts the contour, so
$U(w)$ becomes a sum of pure exponentials $\exp(m(\rho-w)^2+M(\rho-w))$, and Cassels (his
Lemma 3) on those power sums, using that $\rho_0$ contributes two terms equal to $1$, gives
$|S|\ge2-o(1)$. Comparing, $K\ge\pi/(2|\rho_0|)$. His Theorem 3 exhibits Beurling systems
attaining $\pi/2+\varepsilon$, so $\pi/2$ is optimal **in that simple-pole, residue-$(-1)$
world**. Note also that the conclusion is $|\Delta(x)|\ge\cdots$ for arbitrarily large $x$: an
$\Omega$ for the absolute value, where Lemma A already gives both signs.

**What changes at a double pole.** The extra $\log Y$ is slowly varying on a Gaussian window of
width $O(\sqrt m)$ in the log-scale, with relative variation $O(m^{-1/2})\to0$; the cosine, and
therefore the $2/\pi$, survives, and the factors $M\sim\log Y$ cancel between the two sides. The
residue side does not survive. A double pole of $F(s+w)e^{ms^2+Ms}$ at $s=\rho-w$ contributes
$\bigl(c_1+c_2(2m(\rho-w)+M)\bigr)\exp(m(\rho-w)^2+M(\rho-w))$, a linear polynomial in $m$ times
an exponential; Cassels' lemma is about $\sum r_\ell^Le^{i\alpha_\ell L}$ with some $r_\ell=1$ and
does not apply to $(a+bL)z^L$. For $\rho_1$ this is nearly moot: the window $|\gamma-\gamma_1|<5$
contains only $\rho_1$, the next zero being at $21.02$, a gap of $6.89$. Isolation substitutes for
Cassels — and leaves nothing to rotate.

**The actual obstruction: $c_2$ is not real.** At $\rho_1$,
$$
c_2=\frac{r_{\rho_1}}{\rho_1}=0.00610589-0.00501147\,i,
\qquad
\frac{|\Re c_2|}{|c_2|}=0.773,
\qquad
\arg c_2=-0.219\pi .
$$
Put $w=\rho_1$; the $\rho_1$-term in $U(w)$ is $c_1+Mc_2$ with $M=16m>0$ frozen in sign, so
$$
S=U(\rho_1)+U(\overline{\rho_1})=2\Re c_1+2M\Re c_2\ \sim\ 2M\,\Re c_2,
$$
**not** $2M|c_2|$. Varying $m$ scales $M$ but cannot rotate $\arg c_2$: the exponential at
$\rho-w=0$ is identically $1$. The same inner product appears on the $x$-side. Comparing
$|S|\le(4KM/\pi)(1+o(1))$ with $|S|\ge2M|\Re c_2|(1+o(1))$ gives
$$
K\ \ge\ \frac\pi2|\Re c_2|=1.214\,|c_2|,
\qquad\text{i.e.}\qquad 0.008144\times1.214=0.00989
$$
at $x=1.1$: $21\%$, not $57\%$. Bounding $|U(\rho_1)|$ instead of $|S|$ sees the full $|c_2|$ on
the residue side, but then there is no cosine and no $2/\pi$, and one recovers exactly Landau's
$|c_2|$. Adding the conjugate introduces both the factor $2$ and the factor $2/\pi$, net $\pi/2$
— and, when the residue is complex, $\pi/2$ times the real part.

**The ladder.**

| | constant at $x=1.1$ | what it is |
|---|---|---|
| 1 | **0.008144** | Lemma A at $\rho_1$. $\Omega_\pm$, unconditional. **Theorem 1.** |
| 1b | 0.00989 | transplanted kernel at the isolated double pole, $(\pi/2)|\Re c_2|$. $\Omega$ for the absolute value only; unconditional, but a separate write-up. |
| 2 | 0.012793 | $(\pi/2)|c_2|$. **Does not transfer** — would need $\arg c_2=0$, the PNT situation. |
| 3 | 0.016289 | naive aligned pair $2|c_2|$. Not a theorem from one named zero. |
| 4 | 0.039108 | First sixty positive terms of Theorem 2's conditional infinite envelope; not its full value. |

Row 1b is not worth the page: Lemma A already gives both signs at $0.008144$, and the kernel
argument would rewrite Révész §§4–5 with a double-pole residue and no Cassels to win $21\%$ on
the absolute $\limsup$ only. Révész's Theorem 3 is sharpness for residue $-1$ and does not say
that $0.012793$ is a ceiling here. In one sentence: *the PNT residue is $-1$; ours is
$r_{\rho_1}/\rho_1$, with argument $-0.219\pi$, and Cassels cannot rotate it.*

---

## Appendix B. Computation

**Evaluating $\mathrm{ERi}$ at large complex argument.** Gram's series suffers catastrophic
cancellation: intermediate terms reach $e^{|z|}$, so it needs about $|z|/\ln10$ digits — 460 at
$\gamma=10^4$, and roughly 3100 at $\gamma=74920$. Instead use the defining Möbius–$\mathrm{Ei}$
sum with an explicit tail. With $N$ chosen so that $|z|/N$ is small,
$$
\mathrm{ERi}(z)=\sum_{n\le N}\frac{\mu(n)}n\mathrm{Ei}\Bigl(\frac zn\Bigr)
  +(\gamma_E+\mathrm{Log}(-z))A_0-A_1+\sum_{k\ge1}\frac{z^k}{k\,k!}G_k ,
$$
$$
A_0=\sum_{n>N}\frac{\mu(n)}n=-m(N),\quad
A_1=\sum_{n>N}\frac{\mu(n)\log n}n=-1-\sum_{n\le N}\frac{\mu(n)\log n}n,\quad
G_k=\sum_{n>N}\frac{\mu(n)}{n^{k+1}} ,
$$
using $\sum\mu(n)/n=0$ and $\sum\mu(n)\log n/n=-1$. The branch is Elliott's,
$\mathrm{Ei}(w)=\gamma_E+\mathrm{Log}(-w)+\sum_kw^k/(k\,k!)$, equivalently
$\mathrm{Ei}(w)=-E_1(-w)$ for $\Im w>0$. Mathematica's `ExpIntegralEi` is larger by
$i\pi\,\mathrm{sgn}(\Im w)$ (Elliott §4.4, pp. 125–127). Using that alternative branch
consistently in both head and tail gives the same entire $\mathrm{ERi}$, because
the shifts cancel by $\sum\mu(n)/n=0$. This is the route of
Riesel–Göhl, Math. Comp. 24 (1970), 969–983, pp. 978–979, at a larger truncation.

Two implementation points. $E_1$ is evaluated by Gauss–Laguerre,
$E_1(v)=e^{-v}\int_0^\infty e^{-t}/(t+v)\,dt$, for $|w|$ above a threshold and by the power
series below it. And the tail coefficients $G_k=1/\zeta(k+1)-\sum_{n\le N}\mu(n)n^{-(k+1)}$ must
themselves be formed in high precision: both terms are $\approx1$ and their difference is
$\approx N^{-k}$ ($10^{-139}$ at $k=60$, $N=191$), so in double they return $10^{-16}$ of
cancellation noise, which the factor $z^k/k!\sim10^{96}$ amplifies to $10^{80}$. The same disease
as Gram's series, one level down. Quantising $N$ to powers of two lets the tail data be cached.

**The working precision must scale with both $k$ and $N$**, not merely be "high". Resolving
$G_k\approx N^{-k}$ takes about $k\log_{10}N$ significant digits, so a fixed precision is
silently adequate for small $|z|$ and silently wrong past a threshold, the factor $z^k/k!$
amplifying whatever noise remains. Take $\mathrm{dps}=\lceil K\log_{10}N\rceil+60$ with the tail
series truncated at $K=48$ — ample, since $(|z|/N)^k/(k^2k!)<10^{-23}$ by $k=40$ once
$|z|/N\le5$. A fixed $\mathrm{dps}=60+2K$ is not enough: at $N=4096$, $K=80$ it resolves $220$ of
the $289$ digits required, and the computed $\mathrm{ERi}$ then leaves the bound
$|\mathrm{ERi}(z)|\le\|Q\|_1e^{\Re z}$ near $|z|\approx1.9\cdot10^{4}$, reaching $10^{15}$ by
$|z|=3\cdot10^{4}$. That bound is the cheapest check available on any implementation: along a
vertical ray $\mathrm{ERi}$ is bounded, so growth of any kind is a precision failure and not a
feature of the function.

**Validation.** Against the 10142 values computed from Gram's series at 460 digits: maximum
absolute error $4.4\cdot10^{-13}$, median $6.9\cdot10^{-15}$, maximum relative error
$2.6\cdot10^{-11}$; all four Riesz sums at $T=9999$ reproduced. $58$ s for $100\,000$ zeros in
double precision, $0.6$ ms per zero.

**Residues.** $|r_\rho|$ from $(\ast)$ at 40 digits reproduces the eight values quoted in §2 to
all printed digits, with $\Gamma$-form and $\chi$-form agreeing to $10^{-40}$;
$\sum_{n\le60}2|r_\rho|/|\rho|=0.037931$, and dividing by $\pi\sqrt{\log1.1}$ gives the envelope
$R_{60}=0.039108$, the sixty-term partial envelope associated with Theorem 2.

**Pointwise convergence of the residue expansion** (§4.5), 400 residues against the
$100\,000$-zero sum:

| $T$ | $\Sigma R_T-L$ | err $X=60$ | $X=100$ | $X=200$ | $X=400$ | normalised |
|---|---|---|---|---|---|---|
| 2048 | $-2.11739$ | 1.4e−1 | 2.2e−2 | 2.7e−2 | 2.7e−2 | 1.0e−4 |
| 9999 | $-3.75652$ | 1.2e−1 | 3.9e−2 | 4.2e−2 | 5.7e−3 | 7.8e−6 |
| 30000 | $-8.97409$ | 4.6e−1 | 3.2e−2 | 1.0e−2 | 3.9e−2 | 2.7e−5 |
| 74920 | $+30.11368$ | 2.4e−1 | 6.2e−2 | 3.0e−2 | 3.7e−2 | 1.4e−5 |

**Files.** `compute/mobius_ei.py`, `compute/run_k2.py`; data `data/eri_list.csv`,
`data/residues_60.csv`, `data/residues_400.npy`, `data/pair_100k.npy`, [data provenance](../docs/data-provenance.md); figures `data/figures/`.

---

## References

- [G] H. Grobner, *On divergence related to Riemann–von Mangoldt's explicit formula of the
  prime-counting function*, arXiv:2609.02713v1 (2026).
- [La] E. Landau, *Nouvelle démonstration pour la formule de Riemann sur le nombre des nombres
  premiers inférieurs à une limite donnée et démonstration d'une formule plus générale pour le cas
  des nombres premiers d'une progression arithmétique*, Ann. Sci. École Norm. Sup. (3) 25 (1908),
  399–442; Numdam, DOI 10.24033/asens.595. Part I, (25)–(26), p. 425 (Riemann's formula, the
  zero-sum in conjugate pairs by increasing ordinate); Part II §XV, (38), pp. 441–442 (the formula
  for progressions).
- [MO] *Is $\pi(x)=\mathrm{R}(x)-\sum_\rho\mathrm{R}(x^\rho)$ correct at all?*, MathOverflow
  question 386213 (2021), with the accepted answer (2021), an answer by I. Zakharevich (2024) and
  one by H. Grobner (2026-09-05).
- [ANS] A. Akbary, N. Ng, M. Shahabi, *Limiting distributions of the classical error terms of
  prime number theory*, Quart. J. Math. 65 (2014), 743–780; arXiv:1306.1657.
- N. Ng, *The distribution of the summatory function of the Möbius function*, Proc. LMS (3) 89
  (2004), 361–389; arXiv:math/0310381.
- D. Lowry-Duda, *Non-real poles and irregularity of distribution*, J. Number Theory 217 (2020),
  23–35; arXiv:1910.09969.
- A. E. Ingham, *The Distribution of Prime Numbers*, Cambridge Tracts 30 (1932), Theorem H,
  pp. 88–89.
- G. H. Hardy, M. Riesz, *The General Theory of Dirichlet's Series*, Cambridge Tract 18 (1915),
  Theorem 31, p. 45.
- Sz. Gy. Révész, *Oscillation of the remainder term in the prime number theorem of Beurling,
  "caused by a given zeta-zero"*, arXiv:2202.01837v3 (2022);
  J. Pintz, Acta Arith. 36 (1980), 341–365.
- D. J. Platt, *Isolating some non-trivial zeros of zeta*, Math. Comp. 86 (2017), 2449–2467.
- E. C. Titchmarsh, *The Theory of the Riemann Zeta-Function*, 2nd ed. (Heath-Brown), 1986.
- H. L. Montgomery, R. C. Vaughan, *Multiplicative Number Theory I*, CUP (2007).
- J. Elliott, *Analytic Number Theory and Algebraic Asymptotic Analysis*, World Scientific (2025);
  arXiv:2407.17820.
- H. Riesel, G. Göhl, *Some calculations related to Riemann's prime number formula*, Math. Comp.
  24 (1970), 969–983.
- X. Cao, Y. Tanigawa and W. Zhai, *Some mean value results related to Hardy's function*,
  arXiv:2003.11349 (2020).
- D. Zagier, *The first 50 million prime numbers*, Math. Intelligencer 0 (1977), 7–19.
- A. M. Odlyzko, tables of zeros of the Riemann zeta function,
  `www.dtc.umn.edu/~odlyzko/zeta_tables/`.
