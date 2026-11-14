# Piecewise-Quadratic Approximations of Group Likelihoods: Construction and Accuracy

*Draft. Companion to "Convergence of Two-Block Gibbs Samplers with Piecewise-Quadratic Group Likelihoods" (`Convergence_Piecewise_Quadratic_Gibbs.md`), which takes the approximations constructed here as given and studies the sampler.*

---

## Abstract

We construct concave, continuously differentiable, piecewise-quadratic approximations to the group log-likelihoods of a hierarchical generalized linear model, and we argue that they are close to the original likelihoods.

- **Scalar group effects.** Three tangent quadratics, joined at switch points where the tangent lines of the score cross, with level shifts that make the joined function continuous. Log-concavity is automatic.
- **Placement.** The tangency points are placed by minimizing the total-variation distance between the approximate and exact conditional posteriors. The first-order conditions are exact. An explicit positioning rule, in the spirit of Nygren and Nygren (2006), is organized around two standardized shape numbers: the skewness $\varepsilon$ and the kurtosis $\delta$.
- **Hierarchical use.** For use inside the hierarchical sampler, the placement is made robust across the posterior support of the population parameter by a deterministic cavity construction.
- **Accuracy.** We bound the distance between the marginal posteriors under the approximation and under the original likelihood.
- **Vector group effects.** A direct multivariate generalization is impossible. The approximation is instead built from one-dimensional pieces along the linear predictors.

Results are labeled **Theorem/Proposition/Lemma** (proved), **Formal result** (asymptotic expansion without remainder bounds) or **Computation** (a numerically evaluated property of an explicitly defined problem).

---

## 1. Purpose and scope

The two-block Gibbs sampler for a hierarchical generalized linear model needs exact draws of each group effect given the population parameter. When the group likelihood is replaced by a piecewise-quadratic function, that conditional is a finite mixture of truncated normals and can be drawn exactly. This document is concerned with **which** piecewise-quadratic function to use, and **how close** it is. Section 2 fixes notation. Sections 3–4 treat a scalar group effect. Section 5 bounds the error in the marginal posteriors and makes the placement robust across $\gamma$. Section 6 treats a vector population parameter, and §7 a vector group effect.

## 2. Model and notation

### 2.1 Model

Let $\gamma\in\mathbb R^p$ and, for groups $i=1,\dots,J$, let $\beta_i\in\mathbb R^q$:
$$\gamma\sim\mathcal N_p(\mu_\gamma,\Lambda^{-1}),\qquad\beta_i\mid\gamma\sim\mathcal N_q(H_i\gamma,\Psi_i),\qquad p(y_i\mid\beta_i)=\exp\{\ell_i(\beta_i)\}.$$
Here $H_i\in\mathbb R^{q\times p}$, $\Psi_i\succ0$ and $\Lambda\succ0$. For a canonical exponential family with cumulant function $b$,
$$\ell_i(\beta_i)=\sum_{j=1}^{m_i}\bigl[y_{ij}\,x_{ij}^\top\beta_i-n_{ij}\,b(x_{ij}^\top\beta_i)\bigr],$$
where $j$ runs over the $m_i$ distinct covariate patterns of group $i$ (counts aggregated). When $q=1$ we write $\ell_i(\beta)=y_i\beta-n_ib(\beta)$.

**Logit.** $b(u)=\log(1+e^u)$, with $p(u)=b'(u)=\operatorname{expit}(u)$ and $w(u)=b''(u)=p(u)(1-p(u))$.

### 2.2 Two-block Gibbs sampler

One iteration from $\gamma$:

1. **(Block 1)** draw $\beta_i\sim f_i(\cdot\mid H_i\gamma)$ independently for $i=1,\dots,J$, where
$$f_i(\beta\mid\eta)\propto e^{\ell_i(\beta)}\,\phi_q(\beta;\eta,\Psi_i).$$
2. **(Block 2)** draw $\gamma'\sim\mathcal N_p\bigl(\bar\gamma(\beta),P_0^{-1}\bigr)$, where
$$P_0=\Lambda+\sum_iH_i^\top\Psi_i^{-1}H_i,\qquad\bar\gamma(\beta)=P_0^{-1}\Bigl(\Lambda\mu_\gamma+\sum_iH_i^\top\Psi_i^{-1}\beta_i\Bigr).$$

The $\gamma$-marginal chain has transition kernel $K$, is reversible with respect to the marginal posterior $\pi_\gamma$, and has a non-negative spectrum (Liu, Wong and Kong 1994). $\mathcal L(\gamma^{(k)})$ denotes the law after $k$ iterations from a fixed $\gamma^{(0)}$.

### 2.3 General notation

- $\|\mu-\nu\|_{TV}=\sup_A|\mu(A)-\nu(A)|$.
- $\|v\|_A=(v^\top Av)^{1/2}$ for $A\succ0$.
- $\lambda_{\max}$ and $\lambda_{\min}$ denote extreme eigenvalues.
- $\varphi$, $\Phi$ and $\bar\Phi=1-\Phi$ denote the standard normal density, distribution and survival functions.
- $\preceq$ is the Loewner order.
- A density $e^{-U}$ is **$K$-strongly log-concave** if $\nabla^2U\succeq K$ (in the sense of distributions where $U$ is only $C^1$).

---

---

## 3. Three-piece approximations

### 3.1 Tangent quadratics

**Definition 3.1 (tangent quadratic).** For $t\in\mathbb R$,
$$T_t(\beta)=\ell(t)+\ell'(t)(\beta-t)+\tfrac12\ell''(t)(\beta-t)^2 .$$

**Lemma 3.2 (moving the tangency point).** For every $\beta$,
$$\partial_tT_t(\beta)=\tfrac12\ell'''(t)(\beta-t)^2 .$$
More generally, for the order-$k$ Taylor polynomial, $\partial_tT_t^{(k)}(\beta)=\ell^{(k+1)}(t)(\beta-t)^k/k!$.

*Proof.* Differentiate term by term; the sum telescopes. $\square$

**Proposition 3.3 (tangent quadratics generally do not intersect).** For $t_0<t_1$ and all $\beta$,
$$T_{t_1}(\beta)-T_{t_0}(\beta)=\tfrac12\int_{t_0}^{t_1}\ell'''(v)(\beta-v)^2\,dv .$$
If $\ell'''$ has constant sign on $[t_0,t_1]$, the difference has a strict sign for every $\beta$, and the two quadratics do not intersect.

*Proof.* Integrate Lemma 3.2 in $t$. $\square$

For the logit, $\ell'''(u)=n\,w(u)(2p(u)-1)$ has the sign of $u$. So two tangency points on the same side of $p=\tfrac12$ give quadratics that never meet: one lies above $\ell$ and the other below between them (Proposition 3.7). Joining them therefore requires level shifts.

### 3.2 The C¹ construction

**Definition 3.4 (three-piece approximation).** Given tangency points $t_-<t_0<t_+$ and switch points $\tau^-<\tau^+$ with $t_-<\tau^-<t_0<\tau^+<t_+$, define
$$q(\beta)=\begin{cases}T_{t_-}(\beta)+c_-,&\beta\le\tau^-,\\T_{t_0}(\beta),&\tau^-<\beta\le\tau^+,\\T_{t_+}(\beta)+c_+,&\beta>\tau^+,\end{cases}\qquad c_\pm=T_{t_0}(\tau^\pm)-T_{t_\pm}(\tau^\pm).$$
$q$ is continuous by construction. A switch point is **C¹** if the one-sided derivatives of $q$ agree there.

**Proposition 3.5 (C¹ switch points; log-concavity).** The C¹ switch point between tangency points $a<c$ is the intersection of the tangent lines to $\ell'$ at $a$ and $c$:
$$\tau(a,c)=\frac{b'(c)-b'(a)+a\,b''(a)-c\,b''(c)}{b''(a)-b''(c)}\qquad(b''(a)\ne b''(c)).$$
With C¹ switch points, $q'$ is continuous and piecewise linear with slopes $-n\,b''(t_k)<0$. Hence $q$ is concave, and $q-\beta^2/(2\Psi)$ has curvature at least $n\min_kb''(t_k)+1/\Psi$ everywhere.

*Proof.* $T_t'(\beta)=\ell'(t)+\ell''(t)(\beta-t)$ is the tangent line to $\ell'=y-nb'$ at $t$. Equating two such lines gives the formula. A continuous function that is piecewise linear with negative slopes is decreasing. $\square$

**Lemma 3.6 (switch points move with the tangency points).** With $\tau=\tau(a,c)$,
$$\frac{\partial\tau}{\partial a}=\frac{\ell'''(a)(\tau-a)}{n\,(b''(a)-b''(c))},\qquad\frac{\partial\tau}{\partial c}=\frac{\ell'''(c)(c-\tau)}{n\,(b''(a)-b''(c))}.$$

*Proof.* Implicit differentiation of $T'_a(\tau)=T'_c(\tau)$. $\square$

For the logit, both derivatives are positive when $a$ and $c$ lie on the same side of $p=\tfrac12$.

**Proposition 3.7 (position relative to $\ell$; logit).** Let $D_t=\ell-T_t$. Then $D_t''(\beta)=n\bigl(w(t)-w(\beta)\bigr)$. For $t>0$:

- $T_t<\ell$ for $\beta>t$;
- $T_t>\ell$ on $(\beta_c(t),t)$, for a unique $\beta_c(t)<-t$;
- $T_t<\ell$ again for $\beta<\beta_c(t)$.

For $t<0$, mirror. For $t=0$, $T_0\le\ell$ everywhere (Böhning–Lindsay). With all tangency points on one side of $p=\tfrac12$, the joined approximation lies **above** $\ell$ on the side toward $p=\tfrac12$ and **below** it on the side away from $p=\tfrac12$. It crosses $\ell$ at $t_0$ and once more, beyond the mirror point of $t_-$.

*Proof.* $D_t$ is convex where $|\beta|>|t|$ and concave where $|\beta|<|t|$, with $D_t(t)=D_t'(t)=D_t''(t)=0$. The sign pattern follows from the shape on each region and from $D_t\to+\infty$ as $\beta\to-\infty$. The joined statement follows from the signs of $c_\pm$ obtained from Proposition 3.3. $\square$

**Proposition 3.8 (when a C¹ switch point lies between its tangency points; logit).** Let $E(t,u)=p(u)-p(t)-w(t)(u-t)$. Then $\tau(a,c)\in(a,c)$ if and only if $E(a,c)\,E(c,a)>0$.

1. If $a$ and $c$ lie on the same side of $0$, this always holds.
2. If $a<0<c$, it holds if and only if
$$a\in\bigl(a^{**}(c),0\bigr)\ \cup\ \bigl(-\infty,a^*(c)\bigr),$$
where $a^{**}(c)\in(-c,0)$ is the unique root of $E(a,c)=0$ and $a^*(c)<-c$ is the unique root of $E(c,a)=0$.

Both boundaries decrease in $c$. As $c\downarrow0$, $a^{**}(c)=-\tfrac c2+o(c)$ and $a^*(c)=-2c+o(c)$. As $c\to\infty$, $a^{**}(c)\sim-\log c$ and $a^*(c)\sim-e^{c}$. The interval $[a^*(c),a^{**}(c)]$ (the **forbidden band**) contains the mirror point $-c$.

*Proof.* The difference of the two score tangent lines is linear, with endpoint values $-nE(c,a)$ at $a$ and $nE(a,c)$ at $c$. The case analysis uses concavity of $p$ on $[0,\infty)$ and convexity on $(-\infty,0]$, the monotonicity of $E$ in each argument, and the identity $E(-c,c)\propto\sinh c-c>0$. The asymptotics come from expanding $p$ and $w$ at $0$ and from their exponential tails. $\square$

**Corollary 3.9 (which piece is most precise).** Piece precisions are $c_k=n\,w(t_k)$, ordered by $|t_k|$. Take $t_0\ge0$. Then exactly one of the following holds:

- **(B)** $|t_-|<t_0$: the precisions are monotone across the pieces, $c_->c_0>c_+$;
- **(A)** $t_-<a^*(t_0)$: the center piece is the most precise.

No feasible configuration lies between them.

---

---

## 4. Optimal positioning

### 4.1 Criterion and first-order conditions

Fix a Gaussian prior $\mathcal N(\eta,\Psi)$. Let
$$f\propto e^{\ell}\phi(\cdot;\eta,\Psi),\qquad g_\theta\propto e^{q_\theta}\phi(\cdot;\eta,\Psi),\qquad\theta=(t_-,t_0,t_+),$$
with C¹ switch points.

**Lemma 4.1 (derivative of TV).** If $\{f=g_\theta\}$ has Lebesgue measure zero, then
$$\partial_\theta\|f-g_\theta\|_{TV}=-\operatorname{Cov}_{g_\theta}\bigl(\partial_\theta q_\theta,\ \mathbf 1\{f>g_\theta\}\bigr).$$

*Proof.* Write $\|f-g\|_{TV}=\int_{\{f>g\}}(f-g)$ and apply Leibniz' rule; the boundary term vanishes because $f=g$ on the boundary. Then use $\partial_\theta g=g\,(\partial_\theta q-E_g\partial_\theta q)$. $\square$

The same identity holds for forward KL and squared Hellinger distance, with $\mathbf 1\{f>g\}$ replaced by $f/g$ and $\tfrac12\sqrt{f/g}$ respectively.

**Theorem 4.2 (first-order conditions).** Let $A=\{f>g\}$ and $\operatorname{cl}(\beta)=\min(\max(\beta,\tau^-),\tau^+)$. Define
$$\psi_0=(\operatorname{cl}(\beta)-t_0)^2,\qquad\psi_+=(\beta-\tau^+)_+^2-2(t_+-\tau^+)(\beta-\tau^+)_+,\qquad\psi_-=(\tau^--\beta)_+^2-2(\tau^--t_-)(\tau^--\beta)_+ .$$
At an interior minimizer of $\|f-g_\theta\|_{TV}$, for each $k\in\{-,0,+\}$,
$$\ell'''(t_k)\cdot\Bigl(E_g[\psi_k\mid A]-E_g[\psi_k]\Bigr)=0 .$$

*Proof.* By Lemma 3.2, $\partial_{t_k}q=\tfrac12\ell'''(t_k)\psi_k$; the constants come from the level shifts. Moving a C¹ switch point changes $q$ by the slope jump, which is zero there, so the chain rule through Lemma 3.6 contributes nothing. Apply Lemma 4.1. $\square$

**Proposition 4.3 (the root at $p=\tfrac12$; logit).** Let $G(t)=\tfrac12\ell'''(t)\Delta(t)$ be the gradient in one tangency coordinate, with $\Delta(t)=-\operatorname{Cov}_g(\psi,\mathbf 1_A)$. Then $G(0)=0$ and $G'(0)=\tfrac{n}{16}\Delta(0)$, since $\ell''''(0)=n/8$. So $t=0$ is always stationary in that coordinate: a local minimum if $\Delta(0)>0$, a local maximum if $\Delta(0)<0$. $\square$

### 4.2 Standardization and shape numbers

Let $m$ be the mode of $f$, $P^*=-\ell''(m)+1/\Psi$, $s=1/\sqrt{P^*}$ and $\lambda_0=1/(\Psi P^*)$. Put $z=(\beta-m)/s$. Then
$$\log f(m+sz)=\text{const}-\tfrac12z^2+\tfrac\varepsilon6z^3+\tfrac\delta{24}z^4+\tfrac{\mu\varepsilon^3}{120}z^5+\cdots,$$
$$\varepsilon=\ell'''(m)s^3,\qquad\delta=\ell''''(m)s^4,\qquad\kappa=\delta/\varepsilon^2,\qquad\mu=\ell^{(5)}(m)s^5/\varepsilon^3,\qquad r=\delta/\varepsilon .$$
For the logit, with $\vartheta=\tanh(m/2)$ and $w=w(m)$:
$$\varepsilon=(1-\lambda_0)\vartheta s,\qquad\delta=-(1-\lambda_0)(1-6w)s^2,\qquad r=-\frac{(1-6w)\,s}{(1-\lambda_0)\vartheta}.$$
We call $\varepsilon$ the **standardized skewness** and $\delta$ the **standardized kurtosis**. Their ratio $r$ is near $0$ when the mode is far from $p=\tfrac12$ and tends to $+\infty$ as the mode approaches $p=\tfrac12$.

### 4.3 Leading order: the skew-led layout

**Theorem 4.4 (leading-order positions; formal in the reduction, rigorous in the analysis).** Let $\varepsilon\to0$ with $\delta=O(\varepsilon^2)$. To leading order,
$$\|f-g\|_{TV}=\tfrac\varepsilon2\,E_\varphi|R-E_\varphi R|+O(\varepsilon^2),$$
with $R$ the piecewise cubic residual. The layout $z=(-2b,0,2b)$ with switch points $\pm b$ satisfies all three first-order conditions if and only if
$$(1+3b^2)\,\bar\Phi(b)=3b\,\varphi(b).$$
This equation has a unique root $b^*\in(0,\sqrt{2/3})$; numerically $b^*=0.52553$.

*Proof.* For a cubic log-density, C¹ switch points are midpoints of adjacent tangency points, since $T'_t(u)=\ell'(u)-\tfrac\varepsilon2(u-t)^2$. With the symmetric layout, $R$ is odd and $\{R>0\}=\{z>0\}$. So the center condition holds by symmetry, and the side condition reduces to $E_\varphi[\psi_+]=0$, which is the displayed equation.

For uniqueness, let $h(b)=(1+3b^2)\bar\Phi(b)-3b\varphi(b)$. Then $h(0)=\tfrac12>0$ and $h(b)\sim-2\varphi(b)/b<0$ as $b\to\infty$. At any root, $h'(b)=2\varphi(b)(3b^2-2)/(1+3b^2)$. So every root below $\sqrt{2/3}$ is a down-crossing and every root above it an up-crossing. Since $h$ ends negative, no up-crossing can occur, and two down-crossings cannot occur without an up-crossing between them. $\square$

This is the quadratic-tangent analogue of the layout of Nygren and Nygren (2006, §3.2): center tangency at the mode, symmetric outer tangency points, and switch points halfway between.

### 4.4 First order: the skew correction

**Lemma 4.5 (normalization condition; formal).** At the optimum, the approximation built with its center piece unshifted has the same normalizing constant as $f$, to second order in $\varepsilon$:
$$\int e^{q}\phi(\cdot;\eta,\Psi)=\int e^{\ell}\phi(\cdot;\eta,\Psi)\,\bigl(1+o(\varepsilon^2)\bigr).$$

*Derivation.* On the center piece, $\log f-\log g=\tfrac\varepsilon6(z-z_0)^3-C$ with $C=\log(Z_g/Z_f)$. If $C\ne0$ at order $\varepsilon^2$, the crossing of $f$ and $g$ near $z_0$ lies at distance of order $\varepsilon^{1/3}$, and the center condition of Theorem 4.2 is violated at that order. $\square$

**Formal result 4.6 (first-order positions; closed form).** Write $z=(-2b^*,0,2b^*)+\varepsilon\,(c_-,c_0,c_+)$, and let $b=b^*$, $p=\varphi(b^*)$. Then
$$c_0=\frac{5+15b^2-96b^3p-48b^5p+\kappa\,(3+9b^2-64b^3p)}{12\,(1+3b^2-8bp)}=0.8234+0.5118\,\kappa,$$
$$c_\pm=\frac{2b^4(3+\kappa)+3b^2c_0}{2-3b^2}=0.9731+0.4922\,\kappa .$$

*Derivation.* Expand the three conditions of Theorem 4.2 and Lemma 4.5 to the next order. This gives four linear equations, whose coefficients are Gaussian integrals over intervals with endpoints $0,\pm b^*,\pm2b^*$. Using the Mills identity $\bar\Phi(b^*)=3b^*p/(1+3b^{*2})$, the normalization equation involves $c_0$ alone, and each tail equation reads $(3b^2-2)c_\pm+3b^2c_0=-2b^4(3+\kappa)$. $\square$

To first order the layout shifts rigidly toward the flat side; its half-width is unchanged.

**Computation 4.7 (second order).** The second-order term widens the layout symmetrically:
$$z_\pm=\pm\bigl(2b^*+\varepsilon^2h_2\bigr)+\varepsilon c_\pm,\qquad h_2=0.5607+0.4014\,\kappa-0.0532\,\kappa^2+0.0659\,\mu,$$
with no second-order shift (to within $0.01$ for $|\kappa|\le1$). The coefficients are Taylor coefficients of the exact solution of the conditions of Theorem 4.2 in the model of §4.2. They were computed by high-accuracy quadrature; the same procedure reproduces Formal result 4.6 to four decimals.

### 4.5 Kurtosis-led regime and the unified problem

When the mode approaches $p=\tfrac12$, $\varepsilon\to0$ but $\delta$ does not, so the assumption $\delta=O(\varepsilon^2)$ fails.

**Definition 4.8 (unified problem).** For $r\in\mathbb R$, let $Z(r)$ minimize $E_\varphi|R_r-E_\varphi R_r|$ over layouts $z$, where $R_r$ is the piecewise residual of $\tfrac16z^3+\tfrac r{24}z^4$ with C¹ switch points. With $r=\delta/\varepsilon$ held fixed, $\|f-g\|_{TV}=\tfrac\varepsilon2\min_zE|R_r-ER_r|+O(\varepsilon^2)$.

**Proposition 4.9 (mirror symmetry).** $Z(-r)=-\operatorname{rev}Z(r)$, where $\operatorname{rev}$ reverses the order of the three points.

*Proof.* $R_{-r}(-z)=-R_r(z)$ under reflection of the layout, and $|\cdot|$ is even. $\square$

**Computation 4.10 (structure of $Z(r)$).**

1. Near $r=0$, $Z(r)=(-2b^*,0,2b^*)+r\,(0.492,0.512,0.492)+O(r^2)$, in agreement with the $\kappa$-coefficients of Formal result 4.6.
2. There are three branches: X (continuing the skew-led layout), Y (a crossing type, with the steep-side tail far across $p=\tfrac12$) and K (continuing the kurtosis-led layout). They are optimal on $r<0.37$, $0.37<r<2.14$ and $r>2.14$ respectively.
3. As $r\to\infty$, $Z_K(r)=Z_\infty+a/r+b_K/r^2+O(r^{-3})$, with
$$Z_\infty=(-1.5392,-0.4749,1.7925),\qquad a=(-0.6785,-0.9485,-1.1912),\qquad b_K=(-0.1150,-0.1247,0.5671).$$
$Z_\infty$ is **asymmetric**: for a symmetric kurtotic target the optimal three-piece layout breaks the symmetry. The best symmetric layout is 35% worse.

### 4.6 The positioning rule

**Definition 4.11 (three-case rule; logit).** After mirroring so that $m\ge0$:

| case | condition | positions |
|---|---|---|
| **1** | otherwise (including all $r<0$) | $z_0=\varepsilon c_0$, $z_\pm=\pm(2b^*+\varepsilon^2h_2)+\varepsilon c_\pm$ |
| **2** | $r\ge2.14$ | $z=Z_\infty+a/r+b_K/r^2$ |
| **3** | $0.37\le r<2.14$ and $\varepsilon\lvert C_{Y,-}(r)\rvert\le0.3$ | $z=Z_Y(r)+\varepsilon C_Y(r)$ |

In every case, $t_k=m+sz_k$ and the switch points come from Proposition 3.5. $Z_Y$ and $C_Y$ are cubic approximations in $1/r$ to the branch-Y solution and its first-order correction. Their coefficients, and the provenance of all constants, are in `Piecewise_Gaussian_Positioning_Rule.md`.

---

### 4.7 How the accuracy depends on the number of pieces

Throughout this subsection: skew-led regime ($\delta=O(\varepsilon^2)$), tangent pieces joined with C¹ switch points, standardized coordinates $z=(\beta-m)/s$ of §4.2, and $\varphi_0=\varphi(0)=1/\sqrt{2\pi}$.

**Reduction.** For any number of pieces $K$, the minimal TV over placements is
$$\min_{\text{placements}}\|f-g_K\|_{TV}=\varepsilon\,J_K+O(\varepsilon^2),\qquad J_K=\min_{t_1<\dots<t_K}\tfrac12E_\varphi\bigl|R_t-E_\varphi R_t\bigr|,$$
where $R_t$ is the piecewise residual of $z^3/6$. Since $\log f-T_t=\tfrac\varepsilon6(z-t)^3$ for a cubic log-density, the C¹ switch points are the midpoints of adjacent tangency points. $J_K$ is a pure number. **All dependence on the data is through $\varepsilon$.**

**Proposition 4.12 ($K=1,2,3$ in closed form).**
$$J_1=\frac{\varphi_0}3,\qquad J_2=\frac{\varphi_0}3-\frac1{32\,\varphi_0},\qquad J_3=\frac{\varphi_0}3-\frac{2b^{*2}\varphi(b^*)}{1+3b^{*2}},$$
attained respectively at
- $t=0$;
- $t=\pm a^*$ with $a^*=\dfrac1{4\varphi_0}=\dfrac{\sqrt{2\pi}}4$, switch point at $0$;
- $t=(-2b^*,0,2b^*)$ with $b^*$ the root of $(1+3b^2)\bar\Phi(b)=3b\,\varphi(b)$ (Theorem 4.4).

Numerically, $J_1=0.1330$, $J_2=0.0546$ and $J_3=0.0280$.

*Proof.* For symmetric layouts the residual is odd and positive exactly for $z>0$, so $\tfrac12E|R-ER|=E[R;z>0]$. We use the partial moments $\int_0^\infty z\varphi=\varphi_0$, $\int_0^\infty z^2\varphi=\tfrac12$ and $\int_0^\infty z^3\varphi=2\varphi_0$.

- *$K=1$:* $E[z^3/6;\,z>0]=\varphi_0/3$.
- *$K=2$:* for $z>0$, $R=\bigl((z-a)^3+a^3\bigr)/6$, so
$$E[R;z>0]=\tfrac16\bigl(2\varphi_0-\tfrac32a+3a^2\varphi_0\bigr).$$
Its derivative in $a$ is $\tfrac12(2a\varphi_0-\tfrac12)$, which vanishes at $a^*=1/(4\varphi_0)$. Substituting gives $J_2$.
- *$K=3$:* $R=z^3/6$ on $[0,b]$ and $R=\bigl((z-2b)^3+2b^3\bigr)/6$ on $(b,\infty)$. Integrating,
$$E[R;z>0]=\frac{\varphi_0}3-b(1+b^2)\bar\Phi(b)+b^2\varphi(b).$$
At the Mills root, $\bar\Phi(b^*)=3b^*\varphi(b^*)/(1+3b^{*2})$, which gives the stated form. $\square$

**Formal result 4.13 (large-$K$ law and optimal piece density).** For large $K$,
$$J_K=\frac{I^3}{48\,K^2}\bigl(1+o(1)\bigr),\qquad I=\int_{-\infty}^\infty\bar\Phi(|u|)^{1/3}du\ (=2.7075),$$
so $J_K\approx0.413/K^2$. The asymptotically optimal layout gives a piece centered at $z$ the width
$$h(z)=\frac IK\,\bar\Phi(|z|)^{-1/3},$$
that is, it places tangency points at the quantiles of the density proportional to $\bar\Phi(|z|)^{1/3}$.

*Derivation.* On piece $k$, the score of the approximation falls short of the truth by $\tfrac\varepsilon2(z-t_k)^2\ge0$ (Taylor's theorem for a cubic log-density; the C¹ joins prevent any correction across pieces). So the residual $r=\log f-q$ is non-decreasing, with average slope $\varepsilon h^2/24$ over a piece of width $h$. This systematic tilt dominates the within-piece oscillation, which is $O(h^3)$. With a smoothly varying width $h(u)$, $r(z)\approx\tfrac\varepsilon{24}\int_0^zh(u)^2du$. Then
$$\tfrac12E|r-\operatorname{med}r|=\frac\varepsilon{48}\int h(u)^2\,\bar\Phi(|u|)\,du,$$
to be minimized subject to $\int h(u)^{-1}du=K$. The Lagrange condition gives $h^3\propto\bar\Phi^{-1}$, which yields the stated $h$ and value. $\square$

*Consequences.*

1. **Diminishing returns are structural.** Each extra piece improves the TV by about $\varepsilon\,I^3/(24K^3)$, and the improvement shrinks like $K^{-3}$.
2. **Most of the gain comes early.** The closed forms show that the first extra piece removes $(J_1-J_2)/J_1=3/(32\varphi_0^2)=59\%$ of the single-Gaussian error. The second removes a further $20\%$, so two and three pieces capture about $80\%$.
3. **The tilt explains why the rate is $K^{-2}$ rather than $K^{-3}$.** Tangent pieces with C¹ joins always underestimate the score on the flat side, so the error is a tilt of order $h^2$ rather than an oscillation of order $h^3$.

*Check of Formal result 4.13* (`k_pieces/asym_check.py`, quadrature): $J_K$ evaluated at the asymptotic layout, divided by $I^3/(48K^2)$, is $0.79$, $0.86$, $0.90$, $0.94$ and $0.97$ for $K=3,4,6,10,20$, so the formula is approached from below. The exactly optimal values for $K\le3$ are those of Proposition 4.12.

**How the data enter (analytic).** The minimal TV is $\varepsilon J_K$, with
$$\varepsilon=(1-\lambda_0)\tanh(m/2)\,s,\qquad s^2=\frac1{n\,w(m)+1/\Psi_{\rm eff}},\qquad\lambda_0=\frac{s^2}{\Psi_{\rm eff}},\qquad\Psi_{\rm eff}=\Psi_i+H_iVH_i^\top .$$
For a covariate pattern in Part III, $n,y$ are $n_{ij},y_{ij}$ and the cavity variance replaces $\Psi_{\rm eff}$. So the **absolute** benefit of extra pieces, $\varepsilon(J_K-J_{K+1})$, increases:

- as $n_{ij}$ decreases (for a weak prior, $s\approx(n\,w)^{-1/2}$, so the benefit scales like $n^{-1/2}$);
- as $y_{ij}/n_{ij}$ becomes extreme ($\tanh(m/2)\to1$);
- as the group prior $\Psi_i$ weakens;
- as the posterior variance $V$ of $\gamma$ grows, through $\Psi_{\rm eff}$, which raises $s$ and lowers $\lambda_0$.

Near $p=\tfrac12$, where $r$ is large, the same structure holds with the kurtosis-led constants (§4.5). Those are not derived here for general $K$.

**The cost side.** More pieces place the outermost tangency point farther out: $z_{\max}=a^*=0.63$ for $K=2$, $2b^*=1.05$ for $K=3$, and growing for larger $K$ (for the asymptotic layout, the outermost quantile of $\bar\Phi(|z|)^{1/3}$). The flattest curvature falls accordingly, to first order $c_{\min}s^2\approx(1-\lambda_0)-|\varepsilon|z_{\max}-\tfrac12\delta z_{\max}^2$, and the worst-case rate $\rho_{\rm hi}$ rises. Accuracy improves like $K^{-2}$ while $c_{\min}$ keeps decreasing. That favors small $K$.

### 4.8 Adaptive refinement: start with three pieces, add pieces where needed

The idea is to give each group, or each covariate pattern in Part III, as many pieces as its own shape requires. Start with $K_i=3$, measure each group's approximation error deterministically, and refine only the groups whose error exceeds a threshold. Everything below is deterministic.

**Step 0: targets.** Fix an approximation budget $\tau_{\rm apx}$ for the marginal posteriors, and split it as follows.

- *Per-group budget $\tau_\beta$,* for the first term of Theorem 5.2: the averaged conditional TV of each group.
- *Budget $\tau_\gamma$ for the $\gamma$-marginal.* By Theorem 5.1 this requires the weighted conditional-mean errors to satisfy
$$\sum_i\frac{|H_i|}{\Psi_i}\,\bigl|\tilde\mu_i-\mu_i\bigr|\le2\sqrt{\kappa_\gamma}\,\tau_\gamma .$$
For vector $\gamma$, use the norm of Theorem 7.7. A simple equal allocation per group is $|\tilde\mu_i-\mu_i|\le2\sqrt{\kappa_\gamma}\,\tau_\gamma\,\Psi_i/(J|H_i|)$.

**Step 1: predicted number of pieces (closed form, no integration).** By §4.7 the error with $K$ well-placed pieces is about $\varepsilon_iJ_K$. With $J_1=\varphi_0/3$, $J_2$ and $J_3$ from Proposition 4.12, and $J_K\approx I^3/(48K^2)$ for $K\ge3$ (Formal result 4.13),
$$\hat K_i=\min\{K:\ \varepsilon_iJ_K\le\tau_\beta\}\ \approx\ \max\Bigl\{3,\ \Bigl\lceil\sqrt{\tfrac{I^3}{48}}\sqrt{\tfrac{\varepsilon_i}{\tau_\beta}}\Bigr\rceil\Bigr\}=\max\Bigl\{3,\ \Bigl\lceil0.643\sqrt{\varepsilon_i/\tau_\beta}\Bigr\rceil\Bigr\}.$$
Here $\varepsilon_i$ is computed at the group's effective prior $\mathcal N(H_i\gamma^\star,\Psi_i+H_iVH_i^\top)$ (Procedure 5.4). For $K\le3$, use the exact $J_K$. For example, $\tau_\beta=0.005$ gives $\hat K_i\le3$ whenever $\varepsilon_i\le\tau_\beta/J_3=0.18$. For $K\ge4$ the asymptotic formula overstates $J_K$ (by about $30\%$ at $K=4$; §4.7), so it errs toward more pieces. Groups with many trials, or with $y_{ij}/n_{ij}$ near $\tfrac12$, need no refinement.

**Step 2: place $K_i$ pieces.**

- $K_i=3$: the positioning rule of §4.6.
- $K_i\ge4$: tangency points at the quantiles $(k-\tfrac12)/K_i$, $k=1,\dots,K_i$, of the density proportional to $\bar\Phi(|z|)^{1/3}$ (Formal result 4.13), shifted toward the flat side by the first-order correction of Formal result 4.6 as a common shift $\varepsilon c$. Optionally polish by Newton steps on the first-order conditions of Theorem 4.2, which have the same per-piece form for any $K$: one condition $\ell'''(t_k)\,(E_g[\psi_k\mid A]-E_g[\psi_k])=0$ per tangency point.
- Switch points by Proposition 3.5. They are C¹, so log-concavity is automatic for every $K$. For the logit, check feasibility across $p=\tfrac12$ (Proposition 3.8).

**Step 3: verify deterministically.** For each group, compute by one-dimensional quadrature:

- the conditional TV $\|f_i-g_{i,K_i}\|_{TV}$ at the effective prior (or averaged over $\mathcal N(H_i\gamma^\star,H_iVH_i^\top)$ by Gauss–Hermite nodes, which is the averaged criterion);
- the conditional-mean error $|\tilde\mu_i-\mu_i|$ at the same prior.

**Step 4: refine.** For each group failing either target, set $K_i\leftarrow K_i+1$ and return to Step 2 for that group only. Stop when all groups pass, or when the cost criterion below says further pieces do not pay.

**Step 5: update the convergence constants.** Recompute $c_{i,\min}$ (the curvature of group $i$'s flattest piece), $\rho_{\rm hi}$, and the run length $k$ of the convergence document (Corollary 9.1).

**When further pieces stop paying.** Adding a piece to group $i$ reduces its error by about $\varepsilon_i(J_K-J_{K+1})\approx\varepsilon_iI^3/(24K^3)$. It also moves the outermost tangency point outward: $z_{\max}=1.05$ s.d. for $K=3$ (optimal layout), and about $1.79$, $2.18$ and $2.63$ s.d. for $K=4,6,10$ (asymptotic layout). That lowers $c_{i,\min}$ and raises $\rho_{\rm hi}$ by some $\Delta\rho_i$. At a fixed run length $k$, the mixing part of the bound, $B=\rho_{\rm hi}^kD/\sqrt{2\pi}$ (convergence document, Theorem 7.2), increases by about $B\,k\,\Delta\rho_i/\rho_{\rm hi}$. So refine group $i$ only while
$$\varepsilon_i\bigl(J_K-J_{K+1}\bigr)\ >\ B\,k\,\frac{\Delta\rho_i}{\rho_{\rm hi}} .$$
Refinement pays mainly for groups with large $\varepsilon_i$ and small loading on $\gamma$, whose flatter pieces barely move $\rho_{\rm hi}$.

**Properties.**

- Every step is deterministic: closed-form predictions, quadrature checks, Newton polishing.
- The convergence theorems apply for any per-group $K_i$, since they use only $c_{i,\min}$ and $c_{i,\max}$.
- Step 1 alone, without the Step 3 checks, gives a fully closed-form allocation of pieces. Steps 3–4 replace the leading-order predictions by exact values where $\varepsilon_i$ is large and the expansion is least reliable.
- *Not covered:* closed-form $J_K(r)$ for the kurtosis-led regime ($r$ large). For those groups, Step 1 can use the $K\le3$ rule and rely on Steps 3–4.

---

## 5. Accuracy of the approximation: marginal posteriors

Write $\mu_i(\eta)=E[\beta_i\mid\eta]$ under $e^{\ell_i}$ and $\tilde\mu_i(\eta)$ under $e^{q_i}$.

**Theorem 5.1 (marginal posterior of $\gamma$).** Let $D(\gamma)=\sum_i(H_i/\Psi_i)\bigl(\tilde\mu_i(H_i\gamma)-\mu_i(H_i\gamma)\bigr)$ and suppose $|D|\le L$. Then
$$\|\pi_\gamma-\tilde\pi_\gamma\|_{TV}\le\frac{L}{2\sqrt{\kappa_\gamma}} .$$

*Proof.* $\tilde\pi_\gamma\propto\pi_\gamma e^{R}$ with $R'=D$, by Tweedie. $\tilde\pi_\gamma$ is $\kappa_\gamma$-strongly log-concave (Convergence document, Lemma 3.3), so Herbst's argument gives $\mathrm{KL}(\tilde\pi_\gamma\|\pi_\gamma)=\log E_{\tilde\pi}e^{-(R-ER)}\le L^2/(2\kappa_\gamma)$. Then apply Pinsker. $\square$

Thus the **conditional means**, not the per-group TV fit, control the $\gamma$-marginal. Errors accumulate over groups with weights $H_i/\Psi_i$.

**Theorem 5.2 (marginal posterior of $\beta_i$).**
$$\|\pi_{\beta_i}-\tilde\pi_{\beta_i}\|_{TV}\le E_{\tilde\pi_\gamma}\bigl\|f_i(\cdot\mid H_i\gamma)-\tilde f_i(\cdot\mid H_i\gamma)\bigr\|_{TV}+\|\pi_\gamma-\tilde\pi_\gamma\|_{TV}.$$

*Proof.* Insert $\int f_i\,\tilde\pi_\gamma$ and use the triangle inequality; mixing through a common kernel does not increase TV. $\square$

### 5.1 Positioning across the support of $\gamma$

Theorem 5.2 motivates the averaged criterion $A_i(t)=E_{\gamma}\|f_i(\cdot\mid H_i\gamma)-\tilde f_i(\cdot\mid H_i\gamma;t)\|_{TV}$.

**Proposition 5.3 (cavity equivalence; first order).** Suppose the posterior of $\gamma$ given $y_{-i}$ is $\mathcal N(\bar\gamma,V)$. The marginal of $\beta_i$ is then the single-group posterior with prior $\mathcal N(H_i\bar\gamma,\Psi_i+H_i^2V)$. Averaging the conditionals over $\gamma$ widens the conditional variance by $\rho_i^2H_i^2V$ with $\rho_i=(1+c_i\Psi_i)^{-1}$, and
$$\frac{\partial}{\partial(H_i^2V)}\Bigl(c_i+\frac1{\Psi_i+H_i^2V}\Bigr)^{-1}\Big|_{V=0}=\rho_i^2 .$$
So both reductions agree to first order in $H_i^2V$.

*Proof.* Direct computation. $\square$

**Procedure 5.4 (deterministic positioning).**

1. Position every group at a fixed $\gamma^{(0)}$.
2. Compute the mode $\gamma^\star$ of $\tilde\pi_\gamma$ by Newton's method, using the closed-form $\tilde\mu_i$ and $\widetilde{\operatorname{Var}}_i$, and set $V=\bigl(P_0-\sum_iH_i^2\widetilde{\operatorname{Var}}_i/\Psi_i^2\bigr)^{-1}$.
3. Reposition each group with prior $\mathcal N(H_i\gamma^\star,\Psi_i+H_i^2V)$, and iterate.
4. *Refinement:* replace $\gamma^\star$ and $V$ by the leave-one-out mean and variance of $\gamma$ under the other groups' three-piece likelihoods. The cavity then matches the true cavity's first two moments exactly.

The sampler always uses the true $\Psi_i$ and the current $\gamma$.

---

---

## 6. Vector population parameter $\gamma$ (scalar $\beta_i$)

### 6.1 Positioning is unchanged

**Proposition 6.1.** Sections 3–4 apply verbatim to each group. In Procedure 5.4 the effective prior is $\mathcal N(H_i\gamma^\star,\ \Psi_i+H_iVH_i^\top)$, with $V\in\mathbb R^{p\times p}$, and the averaged criterion reduces as in Proposition 5.3 with $H_i^2V$ replaced by $H_iVH_i^\top$.

*Proof.* Group $i$'s conditional depends on $\gamma$ only through the scalar $\eta_i$, whose distribution under $\mathcal N(\gamma^\star,V)$ is $\mathcal N(H_i\gamma^\star,H_iVH_i^\top)$. $\square$

**Lemma 6.2 (Laplace step and leave-one-out downdate).**
$$\nabla\log\tilde\pi_\gamma=-\Lambda(\gamma-\mu_\gamma)+\sum_iH_i^\top\frac{\tilde\mu_i-\eta_i}{\Psi_i},\qquad-\nabla^2\log\tilde\pi_\gamma=P_0-\sum_iH_i^\top H_i\,\frac{\widetilde{\operatorname{Var}}_i}{\Psi_i^2}.$$
Removing group $i$ is a rank-one change to the Hessian, with $h_i=1/\Psi_i-\widetilde{\operatorname{Var}}_i/\Psi_i^2$. So
$$V_{-i}=V+\frac{h_i\,VH_i^\top H_iV}{1-h_iH_iVH_i^\top}\qquad(\text{Sherman–Morrison}).$$
$\square$

---

## 7. Vector group effects $\beta_i$

### 7.1 The obstruction to a direct generalization

**Lemma 7.1 (C¹ joins in $q$ dimensions).** Let $Q_a$ and $Q_b$ be quadratics on $\mathbb R^q$ that agree in value and gradient on the hyperplane $\{x:n^\top x=c\}$. Then
$$Q_a(x)-Q_b(x)=\alpha\,(n^\top x-c)^2$$
for some $\alpha\in\mathbb R$. In particular, $\nabla^2Q_a-\nabla^2Q_b=2\alpha\,nn^\top$ has rank at most one.

*Proof.* $Q_a-Q_b$ is a quadratic that vanishes together with its gradient on the hyperplane. Restricted to any line crossing it, the difference is a one-variable quadratic with a double root at the crossing, so it is a multiple of $(n^\top x-c)^2$. $\square$

So tangent quadratics of a general $\ell_i$ at points of a $q$-dimensional grid cannot be joined smoothly: their Hessians differ in full rank. A C¹, log-concave construction must change curvature in one direction at a time.

### 7.2 Construction by projections

Using the GLM structure $\ell_i(\beta)=\sum_{j=1}^{m_i}\ell_{ij}(x_{ij}^\top\beta)$:

**Definition 7.2.** For each pattern $j$, let $q_{ij}$ be a one-dimensional three-piece C¹ approximation to $u\mapsto\ell_{ij}(u)$ (Definition 3.4), and set
$$q_i(\beta)=\sum_{j=1}^{m_i}q_{ij}(x_{ij}^\top\beta).$$

**Proposition 7.3.** $q_i$ is concave and $C^1$. Its Hessian exists off the hyperplanes $\{x_{ij}^\top\beta=\tau^\pm_{ij}\}$ and satisfies
$$-\nabla^2q_i(\beta)=\sum_jc_{ij,k(j,\beta)}\,x_{ij}x_{ij}^\top\ \succeq\ C_{i,\min}:=\sum_jc_{ij,\min}\,x_{ij}x_{ij}^\top .$$
Across each hyperplane, the Hessian changes by a multiple of $x_{ij}x_{ij}^\top$, in accordance with Lemma 7.1.

*Proof.* A sum of concave $C^1$ functions of linear forms is concave and $C^1$. Differentiate twice. $\square$

**Positioning (cavity per projection).** Pattern $j$'s target is the distribution of $u_{ij}=x_{ij}^\top\beta_i$ under the effective prior $\mathcal N(H_i\gamma^\star,\ \Psi_i+H_iVH_i^\top)$ and the other patterns of the same group. Gaussian-approximating the rest gives a one-dimensional cavity $\mathcal N(m_{ij},v_{ij})$ for $u_{ij}$. Section 4 then applies verbatim to $\ell_{ij}$, with $\mathcal N(m_{ij},v_{ij})$ as its prior, and with shape numbers $\varepsilon_{ij}$ and $\delta_{ij}$ along the direction $x_{ij}$. The procedure is deterministic and has the structure of expectation propagation with one-dimensional sites.

### 7.3 Sampling the conditional

The conditional $\tilde f_i(\beta\mid\eta)\propto e^{q_i(\beta)}\phi_q(\beta;\eta,\Psi_i)$ is a normal truncated to the cells of a hyperplane arrangement.

**Proposition 7.4.** Along any line $\beta=\beta_0+sv$, $\log\tilde f_i$ is a concave, $C^1$, piecewise-quadratic function of $s$ with at most $2m_i+1$ pieces.

*Proof.* Each term $q_{ij}(x_{ij}^\top\beta_0+s\,x_{ij}^\top v)$ contributes at most two breakpoints in $s$. $\square$

Consequently, coordinate-wise or direction-wise Gibbs sweeps within block 1 use exact one-dimensional draws by inversion, generalized from 3 to $2m_i+1$ pieces. Exact joint draws are possible when $m_i$ is small, by enumerating the cells.

### 7.4 Positioning each pattern, and an example

**Per-pattern cavity.** Each pattern's approximation is positioned by the one-dimensional rule of §4. The prior is replaced by the distribution of its own linear predictor $u_{ij}=x_{ij}^\top\beta_i$ under everything else in the model:

1. Start from the group's effective prior $\mathcal N(H_i\gamma^\star,\ \Psi_i+H_iVH_i^\top)$ (Procedure 5.4).
2. Approximate the other patterns' contributions by Gaussians (Laplace at the group's conditional mode, then, once all pieces exist, by their exact first two moments under the pieces).
3. Project onto $x_{ij}$. This gives a one-dimensional cavity $\mathcal N(m_{ij},v_{ij})$ for $u_{ij}$.
4. Apply §4 to $\ell_{ij}$ with prior $\mathcal N(m_{ij},v_{ij})$, giving shape numbers $\varepsilon_{ij}$ and $\delta_{ij}$ along $x_{ij}$.

The procedure is deterministic. It has the structure of expectation propagation with one-dimensional sites.

**Example: random intercept and slope, binary covariate.** The patterns are $x_{i1}=(1,0)^\top$ and $x_{i2}=(1,1)^\top$, so
$$q_i(\beta_0,\beta_1)=q_{i1}(\beta_0)+q_{i2}(\beta_0+\beta_1).$$
- The four switch lines $\beta_0=\tau_{i1}^\pm$ and $\beta_0+\beta_1=\tau_{i2}^\pm$ cut $\mathbb R^2$ into nine cells.
- On each cell, $-\nabla^2q_i=c_{i1,k}\,e_1e_1^\top+c_{i2,k'}\,\mathbf 1\mathbf 1^\top$, for the pieces $k$ and $k'$ active there.
- Crossing a line $\beta_0=\tau_{i1}^\pm$ changes only the $e_1e_1^\top$ term, and crossing $\beta_0+\beta_1=\tau_{i2}^\pm$ only the $\mathbf 1\mathbf 1^\top$ term, as Lemma 7.1 requires.
- The curvature floor $C_{i,\min}=c_{i1,\min}e_1e_1^\top+c_{i2,\min}\mathbf 1\mathbf 1^\top$ is positive definite, because the patterns span $\mathbb R^2$.

**Patterns that do not span $\mathbb R^q$.** The group prior supplies curvature in the remaining directions: $C_{i,\min}+\Psi_i^{-1}\succ0$ always.

**Many patterns.** With a continuous covariate, $m_i$ may be large. The arrangement then has many cells, but neither positioning nor line-wise sampling enumerates them: each pattern is positioned separately, and each line meets at most $2m_i$ switch hyperplanes. Binning the covariate trades accuracy for fewer patterns.

### 7.5 Accuracy for vector group effects

Theorems 5.1 and 5.2 carry over with $\mu_i,\tilde\mu_i\in\mathbb R^q$:
$$D(\gamma)=\sum_iH_i^\top\Psi_i^{-1}\bigl(\tilde\mu_i(H_i\gamma)-\mu_i(H_i\gamma)\bigr)\in\mathbb R^p,\qquad\|\pi_\gamma-\tilde\pi_\gamma\|_{TV}\le\tfrac12\sup_\gamma\|D(\gamma)\|_{K_\gamma^{-1}},$$
with $K_\gamma$ the curvature floor of $\tilde\pi_\gamma$ (Convergence document, Lemma 3.3). The proof is the change of variables $K_\gamma^{1/2}\gamma$, then Herbst and Pinsker. Per-pattern errors in the conditional means accumulate along the directions $H_i^\top\Psi_i^{-1}x_{ij}$.

---

## 8. Why the approximations are close: summary of the arguments

1. **Exact where the likelihood is Gaussian.** Each piece matches value, slope and curvature at its tangency point. The error grows like $\ell'''(t)(\beta-t)^3/6$, so it is third order in the distance from the tangency point.
2. **Placement is optimal in total variation** (Theorem 4.2). The positioning rule approximates the optimum: within a few percent at moderate and large $n$, and within a factor of about $1.5$ at small $n$, for the logit with the glmbayes default prior (computations in `Piecewise_Gaussian_Positioning_Rule.md`).
3. **Marginals.**
   - The $\gamma$-marginal error is controlled by the errors in the groups' conditional **means** (Theorem 5.1).
   - Each $\beta_i$-marginal error is controlled by the group's conditional fit averaged over the posterior of $\gamma$, plus the $\gamma$-marginal error (Theorem 5.2).
   - The cavity placement of Procedure 5.4 targets exactly that average (Proposition 5.3).
4. **Vector group effects.** The projection construction reduces every approximation to a one-dimensional one, to which items 1–3 apply pattern by pattern.

## 9. Open problems

1. Remainder bounds for the formal results of §4.
2. Closed forms for the computed constants, and a derived replacement for the Case 3 tolerance.
3. Tail control in Theorem 5.1: $\mu_i-\tilde\mu_i$ grows linearly in $|\eta|$ for the logit.
4. A mean-matching constraint in the placement, targeting Theorem 5.1 directly.
5. Positioning constants for the Poisson, where $b''=e^u$ is monotone: no mirror point and no forbidden band.

## References

- Nygren, K. N. and Nygren, L. M. (2006). Likelihood subgradient densities. *Journal of the American Statistical Association* 101. *(Page range to be checked against the published version.)*
- Supporting documents: `Piecewise_Gaussian_Positioning_Rule.md` (rule, constants, validation), `Piecewise_Gaussian_GLMM_Specification.md` (derivations), scripts in `positioning_rule/`.
