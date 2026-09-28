# The sum $\sum_\rho m_\rho\,\mathrm{ERi}(\rho\log x)$ over the zeros of a Dirichlet $L$-function diverges

**Open**

## Abstract

The untwisted entire ERi function, summed over the zeros of a primitive Dirichlet
$L$-function, retains the zeta-driven oscillation of its companion zeta series. For real
characters the partial sums are unconditionally $\Omega_\pm(\sqrt T\log T)$; for complex
characters the real part has the same oscillation and the imaginary part converges to
an explicitly identified finite prime-counting expression. No fixed nonnegative real
Riesz order removes the oscillation. Exact envelopes transfer under the Weak Mertens
Conjecture and linear independence for the zeta ordinates. No simplicity assumption
on the Dirichlet zeros is used.

## 0. Notation

$\chi$ is a primitive character mod $q>1$. $\mathfrak a\in\{0,1\}$ is its parity,
$\chi(-1)=(-1)^{\mathfrak a}$; $\mathfrak a=0$ is *even*, $\mathfrak a=1$ is *odd*. (The letter $a$ is
reserved for the contour abscissa below; the two are unrelated.) $\rho=\beta+i\gamma$ runs over the
non-trivial zeros of $L(s,\chi)$, $m_\rho$ is the multiplicity of $\rho$, and for $T>0$

$$
\Sigma R_T^\chi(x):=\sum_{|\gamma|<T} m_\rho\,\mathrm{ERi}(\rho\log x).
$$

All objects below are [G]'s, with $B(y)=-\sum_{n\le y}\mu(n)/n$, $Q(u)=u^{-1}B(u^{-1})\in L^1(0,1)$,
$q_A(u)=e^{-Au}Q(u)$ for fixed $A>0$:

| | | |
|---|---|---|
| $\mathrm{ERi}(z)$ | $=-\int_0^1 Q(u)e^{zu}\,du$; entire; $\mathrm{ERi}(0)=1$ | [G] Prop. 1.4 |
| $f_A(Y)$ | $=\Re\,\mathrm{ERi}(-A+iY)=-\int_0^1 q_A(u)\cos(Yu)\,du$; real, bounded | [G] §1.3 |
| $H_A(z)$ | $=\int_1^\infty f_A(Y)Y^{-z}\,dY$ | [G] §1.3 |
| $P_A(Y)$ | $=\int_1^Y \log(y/C)\,f_A(y)\,dy$, with $C$ fixed in §5 | [G] §1.4 |

**None of these carries $\chi$.** That is the content of §1, and it is why the proofs of [G] §1 are
cited rather than repeated.

## 1. The series, and why its summand is untwisted

Let $\chi$ be real primitive mod $q$, write $\Pi^*(x,\chi)=\sum_{p^k\le x}\chi(p^k)/k$ with the
average taken at jumps, and $\pi^*(x,\chi)$ for the prime count alone. Paper A, Proposition 3, is
the finite identity
$$
\pi^*(x,\chi)+\pi_-^*(\sqrt x)=\sum_{n\le N}\frac{\mu(n)}n\,\Pi^*(x^{1/n},\chi),
\qquad N=\lfloor\log x/\log2\rfloor,
$$
obtained by inverting $\Pi^*(x,\chi)=\sum_k\pi^*(x^{1/k},\chi^k)/k$ in the exponent $k$. The
correction $\pi_-^*(\sqrt x)$ is the trace of the even $k$, for which $\chi^k=\chi_0$; its form is
Paper A, Lemmas 1–2, and is not re-derived here. What matters for this paper is the weight: it is
$\mu(n)/n$, it counts how prime powers of $x$ are stripped, and there is nothing for it to twist
against. Substituting the explicit formula for $\Pi^*(y,\chi)$ term by term produces, at each zero
$\rho$, the term $m_\rho\sum_{n\le N}(\mu(n)/n)\,\mathrm{li}(x^{\rho/n})=m_\rho R_N(x^\rho)$, with
$R_N$ the $N$-th partial sum of
$\mathrm{ERi}(\rho\log x)=\sum_{n\ge1}(\mu(n)/n)\,\mathrm{Ei}(\rho\log x/n)$, Riemann's $R(x^\rho)$,
untwisted. (The explicit formula at this level — $\Pi^*(y,\chi)$ as $-\sum_\rho\mathrm{li}(y^\rho)$
plus the trivial-zero integral plus a constant, $\log L(0,\chi)$ for odd $\chi$ and
$\log L'(0,\chi)-C_0$ for even, $C_0$ Euler's constant — is Landau's (38) of 1908 [La], for every
character of either parity; Paper A's (†) is its truncation, with the remainder from [MV] 12.10.)

Paper A stops there, with a finite sum in $n$ against a finite sum over zeros, and proves the
resulting identity with an explicit remainder (its Theorem A). What it does not do, and what the
sources it corrects do, is replace the partial sum by $\mathrm{ERi}$ itself and let $T\to\infty$.
Paper A §3 shows that this replacement is licensed if and only if
$$
\sum_\rho m_\rho\,\mathrm{ERi}(\rho\log x),\qquad\rho\ \text{over the zeros of }L(s,\chi),
\ \text{ordered by }|\gamma|,
$$
converges. That is the series of this paper. Its summand is [G]'s function; only the index set
has changed.

Two consequences shape what follows. Every object of [G] §1 is available unchanged, since none of
them sees $\chi$; in particular $\sum_{n\ge1}\mu(n)/n=0$, which [G] uses at (2.11)–(2.12), is still
the relevant sum. And there is no twisted analogue to construct. The natural guess — replace
$\mu(n)$ by $\mu(n)\chi(n)$ throughout, so that $1/\zeta$ becomes $1/L(\cdot,\chi)$ — describes a
different function, one for which $\sum\mu(n)\chi(n)/n=1/L(1,\chi)\ne0$ and the argument of §4
fails. The character enters this paper in exactly one place: the second factor of the integrand,
$-L'/L(s,\chi)$ in place of $-\zeta'/\zeta(s)$.

## 2. The contour

Fix $x>1$, $1<c\le2$, $0<a<1$, and let $\Delta_T(a,c)$ be the rectangle with vertices $-a\pm iT$,
$c\pm iT$. Write $I_c^\chi(T)$, $I_{-a}^\chi(T)$ for the right and left vertical contributions and
$H_{a,c}^\chi(T)$ for the two horizontals, exactly as in [G] §2.

> **Proposition 1 (residues; [G] Prop. 2.1).** Let $T$ not be an ordinate of a zero of $L(s,\chi)$.
> Then
> $$
> -\Sigma R_T^\chi(x)-\mathfrak e_\chi
> =\frac1{2\pi i}\oint_{\Delta_T(a,c)}\mathrm{ERi}(s\log x)\Bigl(-\frac{L'}{L}(s,\chi)\Bigr)ds
> = I_c^\chi(T)-I_{-a}^\chi(T)+H_{a,c}^\chi(T),
> $$
> where $\mathfrak e_\chi=0$ if $\chi$ is odd and $\mathfrak e_\chi=1$ if $\chi$ is even.

*Proof.* $\mathrm{ERi}(s\log x)$ is entire, so the integrand is singular only at the poles of
$-L'/L(\cdot,\chi)$, which for primitive non-principal $\chi$ are exactly the zeros of $L(s,\chi)$
($L$ is entire and $L(1,\chi)\ne0$). Inside $\Delta_T(a,c)$ these are the non-trivial zeros with
$|\gamma|<T$, since $0<\beta<1\subset(-a,c)$, together with $s=0$ when $\chi$ is even. No other
trivial zero is enclosed: for odd $\chi$ they lie at $-1,-3,\dots$ and for even $\chi$ at
$-2,-4,\dots$, all with real part $\le-1<-a$. The sides meet no pole: $\Re s=c>1$ is zero-free,
$\Re s=-a\in(-1,0)$ passes between $0$ and $-1$, and $T$ is off the ordinates. At a zero $\rho$ of
multiplicity $m_\rho$ the residue is $-m_\rho\,\mathrm{ERi}(\rho\log x)$. For even $\chi$ the zero at
$0$ is simple, $-L'/L$ has residue $-1$ there, and $\mathrm{ERi}(0)=1$; the residue is $-1$. $\square$

*Remarks.* (i) The interval for $a$ is open because $a=1$ would put the left side through the odd
trivial zeros; for $\zeta$ nothing lives on $\Re s=-1$ and [G] does not need to say this.
(ii) There is no $s=1$ term; [G]'s $\mathrm{ERi}(\log x)$ is the residue at the pole of $\zeta$.
(iii) Two conventions to keep apart. [G] writes $|\Im\rho|<T$ and Montgomery–Vaughan write
$|\gamma|\le T$; with $T$ off the ordinates these agree. More importantly, the $\zeta$ paper's
$\Sigma R_T(x)$ is over $0<|\Im\rho|\le T$, and the $0<$ excludes $\gamma=0$ — vacuously, since
$\zeta$ has no real zero in the strip. Here $\Sigma R_T^\chi(x)$ is over $|\gamma|<T$ and real zeros
of $L(s,\chi)$ are **included**: they lie inside $\Delta_T(a,c)$, and the proof above has counted
them. A reader with both papers open should not take the two sums to have the same range.

## 3. The horizontals

> **Proposition 2 (horizontals; [G] Lem. 2.5).** For every sufficiently large $T$ there is
> $T'\in[T,T+1]$, not an ordinate of a zero of $L(s,\chi)$, such that for every $B>3$
> $$
> H_{a,c}^\chi(T')\ll_{B,a,c,x,q}(\log T')^{2-B},
> $$
> in particular $H_{a,c}^\chi(T')=o(1)$.

*Proof.* By [MV] Lemma 12.7 there is $T_1\in[T,T+1]$ with
$(L'/L)(\sigma\pm iT_1,\chi)\ll(\log qT)^2$ uniformly for $-1\le\sigma\le2$; its proof chooses $T_1$
with $|T_1-\gamma|\gg1/\log qT$ and $|T_1+\gamma|\gg1/\log qT$ for every zero, so $T_1$ is not an
ordinate; take $T':=T_1$. [G] Lemma 2.3 is a bound on $\mathrm{ERi}$ alone and carries no character:
$\mathrm{ERi}(\log x\,(\sigma\pm iT'))\ll_{B,[-a,c],\log x}(\log T')^{-B}$ uniformly for
$\sigma\in[-a,c]$. For fixed $q$, $(\log qT)^2=O_q(\log^2T)$. Multiply and integrate over two
segments of length $a+c$. $\square$

*Remark.* [G] obtains the lower horizontal from the upper "by complex conjugation". That uses the
conjugate symmetry of the zeros of $\zeta$, which the zeros of $L(s,\chi)$ do not have in general
($\rho\mapsto1-\bar\rho$ preserves them; $\rho\mapsto\bar\rho$ sends them to zeros of
$L(s,\bar\chi)$). [MV] 12.7 supplies both sides directly, which is why it is quoted rather than
[MV] 12.2 with conjugation. For real $\chi$ conjugation would in fact suffice.

## 4. The right vertical

> **Proposition 3 (right vertical; [G] Prop. 2.6).** For fixed $1<c\le2$,
> $$
> \lim_{T\to\infty}I_c^\chi(T)=\sum_{n\le N}\frac{\mu(n)}n\,\Pi^*(x^{1/n},\chi),
> \qquad N=\lfloor\log x/\log2\rfloor .
> $$
> If $\chi$ is real, this equals $\pi^*(x,\chi)+\pi_-^*(\sqrt x)$ (Paper A, Proposition 3).

*Proof.* On $\Re s=c$ write $s=c+it$. Since $c>1$,
$$
-\frac{L'}{L}(c+it,\chi)=\sum_{m\ge2}\frac{\Lambda(m)\chi(m)}{m^{c}}\,e^{-it\log m},
$$
absolutely convergent, since $|\chi(m)|\le1$ and $\sum\Lambda(m)m^{-c}<\infty$. By [G] Prop. 1.4,
$$
\mathrm{ERi}(\log x\,(c+it))=-\int_0^1Q_c(u)\,e^{i\log x\,tu}\,du,\qquad Q_c(u):=Q(u)\,e^{cu\log x},
$$
and $Q_c\in L^1(0,1)$ with $\|Q_c\|_1\le e^{c\log x}\|Q\|_1$; on every compact subinterval of
$(0,1]$ it is of bounded variation, its only jumps being those of $Q$ at $u=1/n$. Both are facts
about $Q$ alone. For fixed $T$ the double integral
$\frac1{2\pi}\int_{-T}^T\int_0^1\sum_m|Q_c(u)|\Lambda(m)|\chi(m)|m^{-c}\,du\,dt$ is at most
$\frac T\pi\|Q_c\|_1\sum_m\Lambda(m)m^{-c}<\infty$, so Fubini gives
$$
I_c^\chi(T)=-\sum_{m\ge2}\frac{\Lambda(m)\chi(m)}{m^{c}}\int_0^1Q_c(u)\,
\frac{\sin\bigl(T(u\log x-\log m)\bigr)}{\pi\,(u\log x-\log m)}\,du,
\tag{4.1}
$$
which is [G] (2.7) with $\Lambda(m)$ replaced by $\Lambda(m)\chi(m)$. We evaluate the limit of each
term and justify the interchange of limit and sum.

*Terms with $m>x$.* Here $\log m-\log x>0$, so
$h_m(u):=Q_c(u)/(\pi(u\log x-\log m))$ lies in $L^1(0,1)$ with
$\|h_m\|_1\le\|Q_c\|_1/(\pi(\log m-\log x))$, and the integral in (4.1) is
$\Im\bigl(e^{-iT\log m}\int_0^1h_m(u)e^{iT\log x\,u}\,du\bigr)$, which tends to $0$ by
Riemann–Lebesgue since $\log x>0$. The bound $\|h_m\|_1$ is independent of $T$; multiplied by
$\Lambda(m)|\chi(m)|m^{-c}\le\Lambda(m)m^{-c}$ and summed over $m>x$ it converges, because
$\Lambda(m)\le\log m$ and $\log m-\log x\ge\tfrac12\log m$ for $m\ge x^2$. Dominated convergence
therefore shows the $m>x$ part of (4.1) tends to $0$.

*Terms with $m<x$.* Put $u_m=\log m/\log x\in(0,1)$ and extend $Q_c$ by zero outside $(0,1)$.
Fourier's single-integral theorem ([Tit48], Thm 12 of §1.14) at the point $u_m$, where $Q_c$ is of
bounded variation, gives
$$
\lim_{T\to\infty}\int_0^1Q_c(u)\frac{\sin(T(u\log x-\log m))}{\pi(u\log x-\log m)}\,du
=\frac1{\log x}\cdot\frac{Q_c(u_m^+)+Q_c(u_m^-)}2 .
$$
*The term $m=x$*, present only if $x$ is a prime power: $u_m=1$, and the same theorem at the jump
$u=1$ of the extended $Q_c$ gives the limit $Q_c(1^-)/(2\log x)$.

The finitely many terms with $m\le x$ may be passed to the limit individually. Hence
$$
\lim_{T\to\infty}I_c^\chi(T)
=-\sum_{m<x}\frac{\Lambda(m)\chi(m)}{m^c}\cdot\frac{Q_c(u_m^+)+Q_c(u_m^-)}{2\log x}
-\bigl[x=m\bigr]\frac{\Lambda(m)\chi(m)}{m^c}\cdot\frac{Q_c(1^-)}{2\log x}.
$$
Now $Q_c(u)=e^{cu\log x}\,u^{-1}B(u^{-1})$, and at $u_m$ one has $e^{cu_m\log x}=m^c$ and
$u_m^{-1}=\log x/\log m$. [G] (2.11)–(2.12) defines
$A_0(y):=\sum_{n<y}\mu(n)/n+\tfrac12\mathbf 1_{\mathbb N}(y)\,\mu(y)/y$ and proves
$-\tfrac12\bigl(B(y^-)+B(y^+)\bigr)=A_0(y)$ for all $y\ge1$; the proof uses
$\sum_{n\ge1}\mu(n)/n=0$, and this sum is untwisted and still vanishes. Each $m<x$ therefore
contributes $\frac{\Lambda(m)\chi(m)}{\log m}A_0\bigl(\frac{\log x}{\log m}\bigr)$. For $m=x$,
$Q_c(1^-)=e^{c\log x}\cdot B(1^+)=-m^c$ (here $m=x$, and $B(1^+)=-\mu(1)/1=-1$), so the contribution is
$\frac{\Lambda(m)\chi(m)}{2\log m}=\frac{\Lambda(m)\chi(m)}{\log m}A_0(1)$, as $A_0(1)=\tfrac12$.
Altogether,
$$
\lim_{T\to\infty}I_c^\chi(T)=\sum_{m\le x}\frac{\Lambda(m)\chi(m)}{\log m}\,A_0\Bigl(\frac{\log x}{\log m}\Bigr)
=\sum_{\substack{m\ge2,\ n\ge1\\ m^n<x}}\frac{\Lambda(m)\chi(m)}{\log m}\frac{\mu(n)}n
+\frac12\sum_{\substack{m\ge2,\ n\ge1\\ m^n=x}}\frac{\Lambda(m)\chi(m)}{\log m}\frac{\mu(n)}n,
\tag{4.2}
$$
which is [G] (2.13) with the character inserted.

Group (4.2) by $n$. For fixed $n$, the condition $m^n<x$ is $m<x^{1/n}$, and with $m=p^k$ one has
$\Lambda(m)/\log m=1/k$ and $\chi(m)=\chi(p^k)$, so
$$
\sum_{m<x^{1/n}}\frac{\Lambda(m)\chi(m)}{\log m}+\frac12\sum_{m=x^{1/n}}\frac{\Lambda(m)\chi(m)}{\log m}
=\sum_{p^k<x^{1/n}}\frac{\chi(p^k)}k+\frac12\sum_{p^k=x^{1/n}}\frac{\chi(p^k)}k
=\Pi^*(x^{1/n},\chi),
$$
the $\tfrac12$ being exactly the averaging at a jump. Terms with $m\ge2$ require $x^{1/n}\ge2$, so
only $n\le N$ occur. This is the first display. If $\chi$ is real, Paper A, Proposition 3, gives
the second. $\square$

*Remarks.* (i) This is where a twisted summand would fail: with $\mathrm{ERi}_\chi$ in place of
$\mathrm{ERi}$, the step (2.11)–(2.12) would need $\sum\mu(n)\chi(n)/n$, which is
$1/L(1,\chi)\ne0$. (ii) The proof stops one step short of [G]'s. His last step puts $m=p^k$,
$r=kn$, and collapses (2.13) to $\pi_0(x)$ by $\sum_{n\mid r}\mu(n)=\delta_{r,1}$. With $\chi$
present the summand of (4.2) is $\chi(p)^k\mu(n)/r$, which depends on $k$ and not only on $r$, and
the collapse fails. That failure is the $\chi^k$ obstruction of Paper A §1, and Paper A's
Proposition 3 is its resolution for real $\chi$ — which is why the second display is a citation
and not a continuation of the computation.

## 5. The left vertical, and the constant $C$

> **Lemma (the [G] (2.14) analogue).** For fixed $0<a<1$ and $t\to\infty$,
> $$
> -\frac{L'}{L}(-a+it,\chi)=\log\frac{qt}{2\pi}+\frac{i(a+\frac12)}{t}+\frac{L'}{L}(1+a-it,\bar\chi)+O_{a,q}(t^{-2}),
> $$
> for $\chi$ of either parity.

*Proof.* From $\Lambda(s,\chi)=(q/\pi)^{(s+\mathfrak a)/2}\Gamma((s+\mathfrak a)/2)L(s,\chi)$ and
$\Lambda(s,\chi)=\varepsilon(\chi)\Lambda(1-s,\bar\chi)$, logarithmic differentiation gives
$$
-\frac{L'}{L}(s,\chi)=\log\frac q\pi+\tfrac12\psi\Bigl(\frac{s+\mathfrak a}2\Bigr)+\tfrac12\psi\Bigl(\frac{1-s+\mathfrak a}2\Bigr)+\frac{L'}{L}(1-s,\bar\chi).
$$
Put $s=-a+it$ and use $\psi(z)=\log z-\frac1{2z}+O(|z|^{-2})$ in $|\arg z|\le\pi-\delta$. With
$z_1=(-a+\mathfrak a+it)/2$ and $z_2=(1+a+\mathfrak a-it)/2$,
$$
\psi(z_1)=\log\tfrac t2+\tfrac{i\pi}2-\tfrac{i(\mathfrak a-a)}t+\tfrac it+O(t^{-2}),\qquad
\psi(z_2)=\log\tfrac t2-\tfrac{i\pi}2+\tfrac{i(1+a+\mathfrak a)}t-\tfrac it+O(t^{-2}).
$$
The $\pm i\pi/2$ cancel, the $\pm i/t$ cancel, and
$\tfrac12\psi(z_1)+\tfrac12\psi(z_2)=\log\tfrac t2+\tfrac i{2t}\bigl[-(\mathfrak a-a)+(1+a+\mathfrak a)\bigr]+O(t^{-2})
=\log\tfrac t2+\tfrac{i(a+\frac12)}t+O(t^{-2})$; the parity cancels. Add $\log(q/\pi)$. $\square$

*Remarks.* (i) This is [G] (2.14) with $\log(t/2\pi)\mapsto\log(qt/2\pi)$ and nothing else changed;
there is no cotangent because the two gamma factors pair directly. (ii) That the $1/t$ term is
purely imaginary is used in the next proof and is not automatic; the computation shows its real part
is exactly $0$.

Now specify $A:=a\log x$ and
$$
C:=\frac{2\pi}{q}\log x
$$
in the definition of $P_A$.

> **Proposition 4 (left vertical; [G] Prop. 2.15).** Let $\chi$ be real. There is a constant
> $C_{a,x,\chi}$, independent of $T$, such that
> $$
> I_{-a}^\chi(T)=\frac1{\pi\log x}\,P_A(T\log x)+C_{a,x,\chi}+o(1)\qquad(T\to\infty).
> $$

*Proof.* Write $g(t):=\mathrm{ERi}(\log x\,(-a+it))\,(-L'/L)(-a+it,\chi)$, so that
$I_{-a}^\chi(T)=\frac1{2\pi}\int_{-T}^Tg(t)\,dt$. Since $\chi$ is real, $L(\bar s,\chi)=\overline{L(s,\chi)}$
and $\mathrm{ERi}(\bar z)=\overline{\mathrm{ERi}(z)}$, whence $g(-t)=\overline{g(t)}$ and
$$
I_{-a}^\chi(T)=\frac1\pi\,\Re\int_0^Tg(t)\,dt .
$$
(For complex $\chi$, Proposition 4′ below.) Note
$\mathrm{ERi}(\log x\,(-a+it))=\mathrm{ERi}(-A+i\log x\,t)$ with $A=a\log x$.

Choose $t_0\ge2$ with $t_0\log x\ge1$ and such that the Lemma's expansion holds for $t\ge t_0$.
The integral over $[0,t_0]$ is a constant depending only on $a$, $x$, $\chi$ and is absorbed into
$C_{a,x,\chi}$. Inserting the Lemma, it remains to treat
$$
\frac1\pi\,\Re\int_{t_0}^T\mathrm{ERi}(-A+i\log x\,t)\Bigl(\log\frac{qt}{2\pi}+\frac{i(a+\frac12)}t
+\frac{L'}{L}(1+a-it,\chi)+O_{a,q}(t^{-2})\Bigr)dt .
\tag{5.1}
$$

*The $O(t^{-2})$ term.* $\mathrm{ERi}(-A+i\log x\,t)$ is bounded in $t$ ([G] Prop. 1.4), so the
integrand is absolutely integrable on $[t_0,\infty)$; the term contributes a constant plus $o(1)$.

*The $i(a+\tfrac12)/t$ term.* Since $\Re(iz)=-\Im z$, its contribution is
$-\frac{a+1/2}\pi\int_{t_0}^T\Im\,\mathrm{ERi}(-A+i\log x\,t)\,\frac{dt}t$. By [G] Prop. 1.4,
$\Im\,\mathrm{ERi}(-A+i\log x\,t)=-\int_0^1q_A(u)\sin(\log x\,tu)\,du$ with $q_A\in L^1(0,1)$, and
Fubini for fixed $T$ gives
$$
\int_{t_0}^T\Im\,\mathrm{ERi}(-A+i\log x\,t)\,\frac{dt}t
=-\int_0^1q_A(u)\bigl(\mathrm{Si}(T\log x\,u)-\mathrm{Si}(t_0\log x\,u)\bigr)du,
$$
$\mathrm{Si}(y)=\int_0^y\sin v/v\,dv$. The bracket is bounded uniformly in $T\ge t_0$ and
$u\in(0,1)$, and for each fixed $u>0$ tends to $\frac\pi2-\mathrm{Si}(t_0\log x\,u)$; dominated
convergence gives a finite limit. Constant plus $o(1)$. (This is where the purely imaginary
coefficient in the Lemma is used: a real part $\alpha/t$ would replace $\mathrm{Si}$ by the cosine
integral, whose logarithmic singularity at $u=0$ is not controlled by $q_A\in L^1$.)

*The Dirichlet-series term.* For $a>0$,
$\frac{L'}{L}(1+a-it,\chi)=-\sum_{m\ge2}\Lambda(m)\chi(m)m^{-1-a}e^{it\log m}$, absolutely
convergent. With $\mathrm{ERi}(-A+i\log x\,t)=-\int_0^1q_A(u)e^{i\log x\,tu}\,du$, the contribution
of this term to (5.1) is
$$
\Re\sum_{m\ge2}b_m\,\frac1\pi\int_0^1q_A(u)\,\frac{e^{iT\lambda_m(u)}-e^{it_0\lambda_m(u)}}{i\lambda_m(u)}\,du,
\qquad b_m:=\Lambda(m)\chi(m)m^{-1-a},\ \ \lambda_m(u):=u\log x+\log m,
$$
the interchange being justified by $\sum|b_m|\le\sum\Lambda(m)m^{-1-a}<\infty$ and $q_A\in L^1$.
The $t_0$-exponential is independent of $T$ and gives a constant. For the $T$-exponential,
$\lambda_m(u)\ge\log m\ge\log2>0$ uniformly in $u$, so for each $m$ the integral tends to $0$ by
Riemann–Lebesgue, and it is bounded by $\|q_A\|_1/\log m$ uniformly in $T$. Since
$$
\sum_m|b_m|\cdot\frac{\|q_A\|_1}{\log m}
\le\|q_A\|_1\sum_m\frac{\Lambda(m)}{m^{1+a}\log m}
\le\|q_A\|_1\sum_m m^{-1-a}<\infty,
$$
dominated convergence gives $0$. Constant plus $o(1)$.

*The logarithmic term.* $\Re\,\mathrm{ERi}(-A+i\log x\,t)=f_A(\log x\,t)$ by definition, so the
contribution is $\frac1\pi\int_{t_0}^Tf_A(\log x\,t)\log\frac{qt}{2\pi}\,dt$. Substitute
$y=t\log x$:
$$
\log\frac{qt}{2\pi}=\log\frac{qy}{2\pi\log x}=\log\frac yC,\qquad C=\frac{2\pi}q\log x,
$$
and the contribution becomes
$$
\frac1{\pi\log x}\int_{t_0\log x}^{T\log x}\log\frac yC\,f_A(y)\,dy
=\frac1{\pi\log x}\Bigl(P_A(T\log x)-P_A(t_0\log x)\Bigr),
$$
the second term a constant. Collecting the constants into $C_{a,x,\chi}$ proves the proposition.
$\square$

*Remarks.* (i) $C$ is forced by the substitution, exactly as $C=2\pi\log x$ is forced in [G]; it is
not a parameter. (ii) Every step other than the Lemma is [G]'s with $\Lambda(m)\mapsto\Lambda(m)\chi(m)$
in one Dirichlet series; the Lemma is the only place the functional equation of $L(s,\chi)$ is used,
and its parity cancels there.

> **Proposition 4′ (left vertical, complex $\chi$).** Let $\chi$ be primitive mod $q>1$. There is a
> constant $C_{a,x,\chi}\in\mathbb C$, independent of $T$, such that
> $$
> I_{-a}^\chi(T)=\frac1{\pi\log x}\,P_A(T\log x)+C_{a,x,\chi}+o(1)\qquad(T\to\infty).
> $$

*Proof.* Write $g_\chi(t):=\mathrm{ERi}(\log x\,(-a+it))\,(-L'/L)(-a+it,\chi)$. Now
$L(\bar s,\chi)=\overline{L(s,\bar\chi)}$, so $g_\chi(-t)=\overline{g_{\bar\chi}(t)}$ and
$$
I_{-a}^\chi(T)=\frac1{2\pi}\int_0^Tg_\chi(t)\,dt+\frac1{2\pi}\,\overline{\int_0^Tg_{\bar\chi}(t)\,dt},
$$
which for real $\chi$ is the fold of Proposition 4. In each half the integral over $[0,t_0]$ is a
constant. Insert the Lemma — for $\chi$ in the first half, for $\bar\chi$ in the second. The two
characters have the same parity, so the terms $\log(qt/2\pi)$, $i(a+\tfrac12)/t$ and
$O_{a,q}(t^{-2})$ of the two expansions are identical, and their contributions to the two halves are
a quantity and its conjugate: together they are $\frac1\pi\Re$ of
$\int_{t_0}^T\mathrm{ERi}(-A+i\log x\,t)\bigl(\log\frac{qt}{2\pi}+\frac{i(a+\frac12)}t+O(t^{-2})\bigr)dt$,
which is (5.1) without its Dirichlet-series term, and is treated as there — $P_A$, the
$\mathrm{Si}$-term, a constant plus $o(1)$. The halves differ only in the Dirichlet-series term,
$\frac{L'}L(1+a-it,\bar\chi)$ in the first and $\frac{L'}L(1+a-it,\chi)$ in the second, i.e.
$b_m=\Lambda(m)\bar\chi(m)m^{-1-a}$ against $\Lambda(m)\chi(m)m^{-1-a}$. The treatment of that
term in the proof of Proposition 4 used $\sum_m|b_m|<\infty$, $q_A\in L^1$ and Riemann–Lebesgue,
and at no point the real part; applied to each half without the $\Re$ it gives a complex constant
plus $o(1)$. Collecting the constants proves the proposition. $\square$

*Remark.* $P_A$ is untouched: it comes from the gamma-factor term, which $\chi$ and $\bar\chi$
share. $C_{a,x,\chi}$ is the only place where $\chi\ne\bar\chi$ can show; the Corollary of §6
evaluates its imaginary part, which is $0$.

## 6. Assembly, and the divergence

Along the ordinates $T'$ of Proposition 2, Propositions 1–4 give, for real primitive $\chi$,
$$
\Sigma R_{T'}^\chi(x)=\frac1{\pi\log x}P_A(T'\log x)+\ell_\chi(x)+o(1),
\qquad
\ell_\chi(x)=C_{a,x,\chi}-\pi^*(x,\chi)-\pi_-^*(\sqrt x)-\mathfrak e_\chi .
$$
Here $\ell_\chi(x)$ plays the role of the constant written $L$ in the $\zeta$ paper (its §1.1); the
letter is changed because $L$ is the $L$-function throughout this paper.
The left side does not depend on $a$; the $a$-dependence of $P_A$ (through $A$) and of
$C_{a,x,\chi}$ cancels, as in [G]. For complex primitive $\chi$, Proposition 4′ in place of
Proposition 4 and the first display of Proposition 3 give the same identity with
$\ell_\chi(x)=C_{a,x,\chi}-\sum_{n\le N}\frac{\mu(n)}n\Pi^*(x^{1/n},\chi)-\mathfrak e_\chi$, now
complex; the step from $T'$ to all $T$ below bounds $\Sigma R_{T'}^\chi-\Sigma R_T^\chi$ in modulus
and is unchanged.

**From $T'$ to all $T$.** $P_A'(Y)=\log(Y/C)f_A(Y)=O(\log Y)$ with $f_A$ bounded, and changing $C$
shifts $\log(Y/C)$ by a constant, so $T'=T+O(1)$ moves $P_A(T\log x)$ by $O(\log T)$. On the
other side, $\Sigma R_{T'}^\chi-\Sigma R_T^\chi$ is the sum over the $\ll\log qT$ zeros (with
multiplicity; the local count recorded in [MV], proof of Lemma 12.6, p. 402) with ordinates in
$[T,T']$, each $\ll_B(\log T)^{-B}$ by [G] Lemma 2.3
at $\sigma=\beta\in(0,1)$; the difference is $o(1)$. Both are $o(\sqrt T\log T)$.

> **Theorem.** Let $\chi$ be a real primitive character mod $q>1$ and $x>1$. Then, unconditionally,
> $$
> \Sigma R_T^\chi(x)=\sum_{|\gamma|<T}m_\rho\,\mathrm{ERi}(\rho\log x)=\Omega_\pm\bigl(\sqrt T\log T\bigr)
> \qquad(T\to\infty),
> $$
> and in particular the series $\sum_\rho m_\rho\,\mathrm{ERi}(\rho\log x)$ over the zeros of
> $L(s,\chi)$, ordered by $|\gamma|$, diverges. The $\Omega$ comes from the first zero $\rho_1$ of
> $\zeta$, with the constant of the $\zeta$ paper's Theorem 1, $|r_{\rho_1}/\rho_1|/(\pi\sqrt{\log x})$;
> at $x=1.1$ this is $0.008144$.

*Proof.* By the display and the transfer, $\Sigma R_T^\chi(x)-\ell_\chi(x)=\frac1{\pi\log x}P_A(T\log x)
+O(\log T)$ for all large $T$. The $\zeta$ paper, §3.2–3.5, proves $P_A(Y)=\Omega_\pm(Y^{1/2}\log Y)$
from the double pole of the Mellin transform of $P_A$ at $\rho_1$, with leading coefficient
$r_{\rho_1}/\rho_1$; the argument uses $C$ only through the term $(\log C)H_A$, whose poles are
simple and which cannot cancel the double pole, so it holds for the present $C$ verbatim. The
constant of that $\Omega$ is $|r_{\rho_1}/\rho_1|$; dividing by $\pi\log x$ and putting $Y=T\log x$,
with $\log Y\sim\log T$, turns it into $|r_{\rho_1}/\rho_1|/(\pi\sqrt{\log x})$. $\square$

*Remarks.* (i) The oscillation is produced by the summand, which is a $\zeta$-object; the index set
is the zeros of $L(s,\chi)$. No zero of $L(s,\chi)$ is named in the argument. (ii) By Paper A §3,
the divergence is equivalent to the failure of the identity written with "$=$" in the sources cited
there. (iii) Where the hypotheses enter. *Primitive* is used throughout: Proposition 1 ($L$ entire,
$L(1,\chi)\ne0$, the trivial zeros), Proposition 2 ([MV] 12.7), the Lemma (the completed
$\Lambda$), and Proposition 4. *Real* is used in exactly two places: the second sentence of
Proposition 3, which rewrites $\lim I_c^\chi$ as $\pi^*(x,\chi)+\pi_-^*(\sqrt x)$ via Paper A, and
the fold in the proof of Proposition 4. The first is cosmetic — for any primitive non-principal
$\chi$ the assembly reads
$$
\Sigma R_T^\chi(x)=\frac1{\pi\log x}P_A(T\log x)+C_{a,x,\chi}-\sum_{n\le N}\frac{\mu(n)}n\Pi^*(x^{1/n},\chi)-\mathfrak e_\chi+o(1),
$$
with the same $P_A$; only the naming of the constant needs $\chi$ real. The second is removed by
Proposition 4′, and the Corollary below says what survives for complex $\chi$: the real part carries
the $\Omega$, the imaginary part converges. The Theorem is kept as stated, for real $\chi$, where
$\Sigma R_T^\chi$ is real and the statement is about the sum itself. (iv) For $\zeta$ the
corresponding identity is the one whose convergence was asked on MathOverflow in 2021 [MO], and
whose accepted answer passes to the infinite form without proof; the $\zeta$ paper's §1 records the
thread. The untruncated $\mathrm{li}$-level identity for $\Pi^*(x,\chi)$, from which Paper A's (†)
is the truncation, is Landau's (38) of 1908 [La].

> **Corollary (complex $\chi$).** Let $\chi$ be a primitive character mod $q>1$ and $x>1$. Then,
> unconditionally,
> $$
> \Re\,\Sigma R_T^\chi(x)=\Omega_\pm\bigl(\sqrt T\log T\bigr),\qquad
> \Im\,\Sigma R_T^\chi(x)\to-\Im\sum_{n\le N}\frac{\mu(n)}n\,\Pi^*(x^{1/n},\chi)\qquad(T\to\infty),
> $$
> the first with the constant of the Theorem, the second because $\Im\,\ell_\chi(x)$ is that sum:
> $\Im\,C_{a,x,\chi}=0$. In particular the series
> $\sum_\rho m_\rho\,\mathrm{ERi}(\rho\log x)$ over the zeros of $L(s,\chi)$, ordered by
> $|\gamma|$, diverges for every primitive $\chi$ mod $q>1$.

*Proof.* By §6 with Proposition 4′, along the ordinates $T'$ of Proposition 2,
$\Sigma R_{T'}^\chi(x)-\ell_\chi(x)=\frac1{\pi\log x}P_A(T'\log x)+o(1)$, and $P_A$ is real. Real
parts: the proof of the Theorem applies verbatim to $\Re\,\Sigma R_T^\chi$. Imaginary parts:
$\Im\,\Sigma R_{T'}^\chi(x)-\Im\,\ell_\chi(x)=o(1)$ along $T'$, and
$\Sigma R_{T'}^\chi-\Sigma R_T^\chi=o(1)$ for $T'=T+O(1)$ by the step from $T'$ to all $T$, so
$\Im\,\Sigma R_T^\chi(x)\to\Im\,\ell_\chi(x)$ through all $T$. It remains to show
$\Im\,C_{a,x,\chi}=0$, i.e. $\Im\,I_{-a}^\chi(T)\to0$, $P_A$ being real. By the fold,
$\Im\,I_{-a}^\chi(T)=\frac1{2\pi}\Im\int_0^T(g_\chi-g_{\bar\chi})(t)\,dt$, and the difference is
exact, not asymptotic: the gamma terms of the two functional equations coincide, so for every real $t$
$$
(g_\chi-g_{\bar\chi})(t)=\mathrm{ERi}(-A+i\log x\,t)\Bigl(\frac{L'}L(1+a-it,\bar\chi)-\frac{L'}L(1+a-it,\chi)\Bigr)
=\mathrm{ERi}(-A+i\log x\,t)\sum_{m\ge2}2i\,\Lambda(m)\,\Im\chi(m)\,m^{-1-a}e^{it\log m},
$$
the series absolutely convergent. With $\mathrm{ERi}(-A+iY)=-\int_0^1q_A(u)e^{iYu}\,du$ and Fubini
as in the proof of Proposition 4,
$$
\int_0^T(g_\chi-g_{\bar\chi})\,dt
=2\sum_{m\ge2}\Lambda(m)\,\Im\chi(m)\,m^{-1-a}\int_0^1\frac{q_A(u)}{\lambda_m(u)}\,du
-2\sum_{m\ge2}\Lambda(m)\,\Im\chi(m)\,m^{-1-a}\int_0^1\frac{q_A(u)\,e^{iT\lambda_m(u)}}{\lambda_m(u)}\,du .
$$
The first sum is real — $q_A$, $\lambda_m$ and $\Im\chi(m)$ are — and does not depend on $T$; the
second tends to $0$ by Riemann–Lebesgue and dominated convergence, exactly as the $T$-exponential
in Proposition 4. Hence $\Im\int_0^T(g_\chi-g_{\bar\chi})\,dt\to0$. $\square$

*Remarks.* (i) The zeros of $L(s,\bar\chi)$ are the conjugates of those of $L(s,\chi)$, so
$\Sigma R_T^\chi+\Sigma R_T^{\bar\chi}=2\Re\,\Sigma R_T^\chi$ is the ERi series over the zeros of
$L(s,\chi)L(s,\bar\chi)$, and the divergence of the real part is the Theorem for that real
$L$-function in all but name; the imaginary part is the antisymmetric combination
$(\Sigma R_T^\chi-\Sigma R_T^{\bar\chi})/2i$, from which the $\zeta$-driven oscillation — symmetric,
because $P_A$ is — cancels. (ii) Numerically, for the character of order $4$ mod $5$ and its
conjugate on PARI lists to $T=3000$ (3237 and 3236 positive ordinates; `compute/eri_series_chi5.py`, output
`data/eri_series_chi5_pari3000.txt`): at $x=1.1$ the per-octave width of $\Re\,\Sigma R_T^\chi$
over the octaves from $[100,200)$ to $[1600,3000)$ is $2.8,\ 4.2,\ 5.4,\ 10.0,\ 15.9$, that of
$\Im\,\Sigma R_T^\chi$ is $0.19,\ 0.23,\ 0.18,\ 0.14,\ 0.11$; at $x=2$ the imaginary part reads
$-0.498,\ -0.495,\ -0.498$ at $T=300,\ 1000,\ 3000$, which is the Corollary's
$-\Im\,\Pi^*(2,\chi)=-\Im\,\chi(2)/2=-\tfrac12$ ($\chi(2)=i$, $N=1$); at $x=1.1$, $N=0$ and the
predicted limit is $0$: $-0.002$ at $T=3000$, the last octave's width being $0.11$.

The positive spectra of $\chi$ and $\bar\chi$ need not coincide: conjugation relates
positive ordinates of one to negative ordinates of the other. Their one-zero count
difference is consistent with this asymmetry and provides no completeness certificate.

### 6.1 Riesz means

The same identity settles the Riesz-mean question, and settles it through $P_A$ rather than through
any abscissa of the $L$-zeros.

> **Corollary (Riesz refutation).** Let $\chi$ be real primitive mod $q>1$. For every fixed $k\ge0$
> and every fixed $x>1$, writing
> $$
> \Sigma R_T^{(k),\chi}(x):=\sum_{|\gamma|<T}\Bigl(1-\frac{|\gamma|}T\Bigr)^k m_\rho\,\mathrm{ERi}(\rho\log x),
> $$
> one has $\Sigma R_T^{(k),\chi}(x)=\Omega_\pm(\sqrt T\log T)$ unconditionally, with the same
> constant as the $\zeta$ paper's Theorem 3, $|r_{\rho_1}/\rho_1|\,W_k(\gamma_1)/(\pi\sqrt{\log x})$ —
> at $x=1.1$, $0.008144\cdot W_k(\gamma_1)$ — where
> $\gamma_1$ is the ordinate of the first zero of $\zeta$ and
> $W_k(\gamma)=\bigl|\Gamma(k+1)\Gamma(\tfrac32+i\gamma)/\Gamma(\tfrac32+k+i\gamma)\bigr|$. No Riesz
> order sums the series.

*Proof.* For $k=0$ this is the Theorem. For $k>0$ the Riesz mean is the weighted integral of the
partial sum: writing $S(u)=\Sigma R_u^\chi(x)$, one has
$\frac kT\int_0^T(1-u/T)^{k-1}\,du=1$, and
$$
\Sigma R_T^{(k),\chi}(x)-\ell_\chi(x)=\frac kT\int_0^T\Bigl(1-\frac uT\Bigr)^{k-1}\bigl(S(u)-\ell_\chi(x)\bigr)\,du.
$$
For $0<k<1$ the kernel has an integrable singularity at $u=T$; the identity remains valid.
By §6, $S(u)-\ell_\chi=\frac1{\pi\log x}P_A(u\log x)+O(\log u)$ for large $u$. Split the integral
at a fixed $T_0$. On $[0,T_0]$ only finitely many zeros enter, $S(u)-\ell_\chi$ is $O_{x,T_0}(1)$,
and the contribution is $O(1/T)$. On $[T_0,T]$ the $O(\log u)$ integrates to $O(\log T)$, since
the kernel has mass $1$. Hence
$$
\Sigma R_T^{(k),\chi}(x)-\ell_\chi(x)=\frac1{\pi\log x}\cdot\frac kT\int_0^T\Bigl(1-\frac uT\Bigr)^{k-1}P_A(u\log x)\,du+O(\log T).
$$
The main term is the Riesz mean of $P_A$ in the variable $Y=T\log x$ (substitute $y=u\log x$).
The $\zeta$ paper's §5.2 records that the Beta integral
$\frac kT\int_0^T(1-u/T)^{k-1}u^\rho\,du=T^\rho\,k\,B(\rho+1,k)$ multiplies the double pole of the
Mellin transform of $P_A$ at $\rho_1$ by $k\,B(\rho_1+1,k)=\Gamma(k+1)\Gamma(\rho_1+1)/\Gamma(\rho_1+1+k)$,
which is non-zero ($W_k(\gamma_1)=|k\,B(\rho_1+1,k)|$), so Lemma A applies to that weighted
primitive. Its Theorem 3 is the resulting $\Omega_\pm(\sqrt T\log T)$ with the stated constant.
The $O(\log T)$ is $o(\sqrt T\log T)$ and does not affect it. $\square$

*Remarks.* (i) It is the same $P_A$ and the same $\rho_1$ of $\zeta$, so the constant is not merely
of the same form as the $\zeta$ paper's — it is identical, $W_k$ evaluated at $\zeta$'s first
ordinate. The zeros of $L(s,\chi)$ do not enter the leading oscillation at any Riesz order.
(ii) The naive expectation — Riesz summability once $k$ exceeds the abscissa
$\Theta_\chi=\sup\{\Re\rho:L(\rho,\chi)=0\}$ — is false, exactly as it is for $\zeta$ (the $\zeta$
paper's §5 opens by refuting it). The oscillation is logarithmic in $\log T$, produced by $P_A$,
so its instantaneous frequency is $\gamma_1/T$ from $\zeta$, not from $L$; a width-$T$ kernel
widens as fast as that frequency shrinks, and the suppression $W_k(\gamma_1)$ is independent of
$T$. The refutation goes through $P_A$, not through $\Theta_\chi$.
(iii) The small-$u$ contribution is $O(1/T)$. It is not identified with the $\zeta$ paper's
numerically fitted transient $k\tau(x)/T$. Its coefficient has not been derived.
(iv) For complex primitive $\chi$ the same proof, with Proposition 4′, gives the statement for
$\Re\,\Sigma R_T^{(k),\chi}(x)$ with the same constant; the imaginary part of the Riesz mean is
the weighted integral of a convergent function and converges.

### 6.2 The Theorem at $q=4$, numerically

On the 12 349 zeros of $L(s,\chi_{-4})$ to $T=10^4$ of Paper A §4 (PARI, count-checked, and
cross-checked zero for zero against an independent finder on $[0,3000]\cup[6000,8000]$),
$\Sigma R_T^\chi(1.1)=\sum_{|\gamma|<T}\mathrm{ERi}(\rho\log1.1)$ was evaluated with the $\zeta$
paper's evaluator (its Appendix B; `compute/paperB_q4.py`, output `data/paperB_q4.txt`, the
pair sums in `data/eri_list_chi4.csv`) and compared with the $\zeta$ paper's residue model of its
§4 — the sixty $r_\rho$ at the first sixty zeros of $\zeta$, nothing of $L$ — with
$C=\frac\pi2\log x$ in place of $2\pi\log x$:
$$
\Sigma R_T^\chi(x)-\ell_\chi(x)\approx\frac1{\pi\log x}\,2\Re\sum_{\rho}\frac{r_\rho\,Y^\rho}\rho
\Bigl(\log\frac YC-\frac1\rho\Bigr),\qquad Y=T\log x .
$$
The centre $\ell_\chi(1.1)=-10.487$ is fitted, being the absorbed constant; the same fit on the
$\zeta$ data returns $-8.9414$ against that paper's $L=-8.9418$. Over $T\in[64,10^4]$ the model's
rms residual is $0.11$ and its maximum $0.41$, against swings of the sum up to $18.6$; at the ten
checkpoints of the $\zeta$ paper the relative error is $0.4$–$3.2\%$, as $0.3$–$3\%$ there. With
$C=2\pi\log x$ kept instead, the rms residual is $1.26$ and the maximum $3.4$: the constant of §5
is what the data see, at the twelve-fold level. Normalised by the logarithm $P_A$ carries,
$F_q(T):=(\Sigma R_T^\chi-\ell_\chi)/(\sqrt T\log(qT/2\pi))$, the per-octave maxima of $|F_q|$ on
the octaves from $[128,256)$ to $[8192,10^4)$ are
$$
0.0225,\ 0.0237,\ 0.0227,\ 0.0234,\ 0.0223,\ 0.0263,\ 0.0128,
$$
and the $\zeta$ paper's $|F|$ on the same octaves are $0.0224,\ 0.0237,\ 0.0227,\ 0.0235,\ 0.0223,\
0.0263,\ 0.0128$. The two normalised functions have correlation $0.998$ on $[64,10^4]$ and rms
difference $0.0008$ against rms $0.0130$ each; the unnormalised difference
$(\Sigma R_T^\chi-\ell_\chi)-(\Sigma R_T-L)$ has per-octave maxima $0.031$–$0.037$ times $\sqrt T$,
the $O(\sqrt T)$ that $\log q$ times the integral of $f_A$ predicts under WMC. The first-sixty
partial amplitude guide is
$\sum_{0<\gamma\le\gamma_{60}}2|r_\rho/\rho|\sqrt{\log x}\,\log q/(\pi\log x)\approx0.054$;
it is neither the infinite envelope nor a rigorous upper bound. This is remark (i) after the
Theorem in numbers: the zeros of $L(s,\chi)$ index the sum and do not shape it.

## 7. What remains conditional

The Theorem and the Corollaries are unconditional. All use nothing about the zeros of $L(s,\chi)$
beyond the standard local count ([MV], proof of Lemma 12.6, p. 402) and the height choice in
[MV] Lemma 12.7, with multiplicity: the oscillation is
produced by $\rho_1$ of $\zeta$, through $P_A$.

The exact $\limsup$ — the value of $\limsup|\Sigma R_T^\chi-\ell_\chi|/(\sqrt T\log T)$ rather than
its positivity — transfers the same way. Section 6 gives
$\Sigma R_T^\chi-\ell_\chi=\frac1{\pi\log x}P_A(T\log x)+O(\log T)$, and $O(\log T)$ dies on
dividing by $\sqrt T\log T$, so the normalised $\limsup$ **equals** that of the $\zeta$ paper's
$|\Sigma R_T-L|$. It is therefore the $\zeta$ paper's Theorem 2, conditional on the Weak Mertens
Conjecture and linear independence for $\zeta$ (Titchmarsh 14.29 for $M(x)=\sum_{n\le x}\mu(n)$,
and LI of the positive ordinates of $\zeta$). No twin for $M(x,\chi)$ is required. The weighted
exact $\limsup$ is likewise the $\zeta$ paper's Theorem 4, the same hypotheses, with $W_k(\gamma)$
multiplying the sum over **zeros of $\zeta$**, not of $L$.

These envelope formulas use the infinite residue sum. The sixty-term numerical comparisons
in §6.2 are partial sums and do not determine its exact value.

A Bohr series for $\Sigma R_T^\chi$ over the zeros of $L(s,\chi)$ would be a different description,
and that is what would need a Weak Mertens twin for $M(x,\chi)$ (simplicity of those zeros, and
$\sum|\rho\,L'(\rho,\chi)|^{-2}<\infty$) together with Pearce–Crump Lemma 19 in the role of
Titchmarsh 14.16. Such an alternative expansion is not established here; those
character-specific conditions are not used above. The $m_\rho$ throughout this paper
is present because simplicity of the $L$-zeros is not assumed, not because a 14.29 twin is missing.

---

## References

- [G] H. Grobner, *On divergence related to Riemann–von Mangoldt's explicit formula of the
prime-counting function*, arXiv:2609.02713.

- [MV] H. L. Montgomery and R. C. Vaughan,
*Multiplicative Number Theory I: Classical Theory*, Cambridge Studies in Advanced Mathematics 97,
Cambridge University Press (2007), Ch. 12: proof of Lemma 12.6 and Lemma 12.7, p. 402;
Theorem 12.10, pp. 403–405. The local zero-count estimate is cited as recorded in the proof of
Lemma 12.6; that passage attributes it to Theorem 10.17.

- [Tit48] E. C. Titchmarsh, *Introduction to the Theory of Fourier
Integrals*, 2nd ed., Oxford 1948 (Theorem 12 of §1.14, as cited in [G]).

- [Tit] E. C. Titchmarsh, *The
Theory of the Riemann Zeta-Function*, 2nd ed., revised by D. R. Heath-Brown,
Oxford University Press (1986), §14.16 and §14.29.

- [PC] A. Pearce-Crump, *Negative discrete second moments of Dirichlet $L$-functions*, arXiv:2606.25094
(Lemma 19, the Dirichlet twin of [Tit] 14.16; §7 only).

- [La] E. Landau, *Nouvelle démonstration
pour la formule de Riemann sur le nombre des nombres premiers inférieurs à une limite donnée et
démonstration d'une formule plus générale pour le cas des nombres premiers d'une progression
arithmétique*, Ann. Sci. École Norm. Sup. (3) 25 (1908), 399–442, Part II §XV (38), pp. 441–442
(Numdam, DOI 10.24033/asens.595).

- [MO] *Is $\pi(x)=\mathrm{R}(x)-\sum_\rho\mathrm{R}(x^\rho)$
correct at all?*, MathOverflow question 386213 (2021), with its accepted answer (2021) and the
answers of I. Zakharevich (2024) and H. Grobner (2026-09-05). Paper A: [the real-character manuscript](real-character-prime-counts.md). The
$\zeta$ paper: [the zeta manuscript](zeta-zero-sums.md).
