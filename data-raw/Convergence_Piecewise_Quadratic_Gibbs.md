# Convergence of Two-Block Gibbs Samplers with Piecewise-Quadratic Group Likelihoods

*Draft. Companion to "Piecewise-Quadratic Approximations of Group Likelihoods: Construction and Accuracy" (`Construction_Piecewise_Quadratic_Approximations.md`).*

---

## Abstract

We give explicit convergence bounds for the two-block Gibbs sampler in a hierarchical model with a Gaussian population layer, when each group's log-likelihood is concave, continuously differentiable and piecewise quadratic, with Hessian bounded between two positive definite matrices. Such likelihoods arise as approximations to generalized-linear-model likelihoods (see the companion document), but here they are taken as the model.

All bounds are expressed through two **rate matrices**, built from the flattest and the steepest curvatures of the pieces. We treat four cases:

1. scalar $\gamma$, scalar $\beta_i$;
2. vector $\gamma$, scalar $\beta_i$;
3. scalar $\gamma$, vector $\beta_i$;
4. vector $\gamma$, vector $\beta_i$.

In Cases 1–3 the contraction is pathwise, obtained by quantile coupling of a one-dimensional statistic. In Case 4 it is pathwise for groups whose loading matrix has rank one and holds in $L^2$ for the others. In all four cases the rate is the largest eigenvalue $\rho_{\rm hi}<1$ of the same rate matrix. In every case we give a total-variation bound with explicit constants and a sandwich for the exact geometric rate. All results reduce to the Gaussian case when the pieces have a single curvature.

---

## 1. Setting

### 1.1 Model

For $\gamma\in\mathbb R^p$ and groups $i=1,\dots,J$ with $\beta_i\in\mathbb R^q$:
$$\gamma\sim\mathcal N_p(\mu_\gamma,\Lambda^{-1}),\qquad\beta_i\mid\gamma\sim\mathcal N_q(H_i\gamma,\Psi_i),\qquad p(y_i\mid\beta_i)=e^{q_i(\beta_i)},$$
with $H_i\in\mathbb R^{q\times p}$, $\Psi_i\succ0$ and $\Lambda\succ0$. Write
$$A_i=\Psi_i^{-1}H_i\in\mathbb R^{q\times p},\qquad P_0=\Lambda+\sum_iH_i^\top\Psi_i^{-1}H_i .$$

**What $P_0$ is.** $P_0=\Lambda+\sum_iH_i^\top\Psi_i^{-1}H_i$ is the precision of $\gamma$ **given all the group effects** $\beta_1,\dots,\beta_J$: the conditional distribution used in block 2. It combines the population prior precision $\Lambda$ with the information the $\beta_i$ carry about $\gamma$ through the hierarchy, as if the $\beta_i$ were observed. It is:

- **not** the prior precision of $\gamma$ ($\Lambda$ alone), which is usually much smaller;
- **not** the posterior precision of $\gamma$. The marginal posterior has curvature $P_0-R(\gamma)$ (Lemma 3.3), between $(1-\rho_{\rm hi})P_0$ and $(1-\rho_{\rm lo})P_0$, smaller than $P_0$ because the $\beta_i$ are themselves uncertain.

$P_0$ does not depend on the data $y$. It appears as the precision of the block-2 noise, as the metric in which distances and contractions are measured, and in $\kappa_\gamma=P_0(1-\rho_{\rm hi})$.

*Example* (the 20-group model of §9.2: $\Lambda=0.01$, $\Psi_i=0.5$, $H_i$ evenly spaced on $[-1,1]$): $P_0=0.01+\sum_iH_i^2/0.5=14.75$. The prior precision is $0.01$. The actual posterior precision of $\gamma$, by quadrature, is $1/0.1303=7.67$, which lies between $(1-\rho_{\rm hi})P_0=6.72$ and $P_0$.

### 1.2 Assumption on the group likelihoods

**Assumption PQ.** Each $q_i:\mathbb R^q\to\mathbb R$ is concave, $C^1$, and piecewise quadratic on finitely many polyhedral cells, with
$$C_{i,\min}\preceq-\nabla^2q_i(\beta)\preceq C_{i,\max}\qquad\text{for almost every }\beta,$$
for symmetric matrices $0\preceq C_{i,\min}\preceq C_{i,\max}$. When $q=1$ these are scalars $0\le c_{i,\min}\le c_{i,\max}$.

The three-piece and projection constructions of the companion document satisfy Assumption PQ:

- for $q=1$, $c_{i,\min}=\min_kc_{ik}$ and $c_{i,\max}=\max_kc_{ik}$;
- for $q\ge2$, $C_{i,\min}=\sum_jc_{ij,\min}x_{ij}x_{ij}^\top$ and $C_{i,\max}=\sum_jc_{ij,\max}x_{ij}x_{ij}^\top$.

**Remark (number of pieces).** Nothing in §3–§9 uses the number of pieces. The proofs use only Assumption PQ:

- concavity;
- $C^1$ junctions, needed for the Cramér–Rao lower bound in Lemma 3.2 (Lemma 3.4 also tolerates concave kinks);
- the two curvature bounds.

A $K$-piece construction, for any $K\ge1$, gives the same theorems with $c_{i,\min}=\min_kc_{ik}$ and $c_{i,\max}=\max_kc_{ik}$ over its $K$ pieces. Sampling stays exact: a $K$-component truncated-normal mixture, drawn by inversion. The number of pieces affects the bounds **only through these two constants**:

- $K=1$ (a single Gaussian): $c_{\min}=c_{\max}$, so $\rho_{\rm lo}=\rho_{\rm hi}$, the Laplace-type rate, which is typically optimistic for the true posterior.
- **Larger $K$:** the extra pieces sit farther into the tails, where the true curvature is smaller. So $c_{\min}$ decreases, and $\rho_{\rm hi}$ increases, while the approximation improves. As $K\to\infty$ with pieces covering the whole line, $c_{\min}\to\inf_\beta n\,b''(\beta)=0$ for the logit, and $\rho_{\rm hi}$ tends to the prior-only rate. That is, the bound becomes vacuous exactly as the approximation becomes exact.
- **Smaller $K$** (for example $K=2$: one junction, so only one tail can be represented) keeps $c_{\min}$ larger but fits the likelihood less well.

The worst-case rate is driven by the flattest piece **regardless of how much posterior mass that piece carries**. So the choice of $K$ and of the outer tangency points trades accuracy of the approximation against sharpness of the bound. Three pieces, with outer tangency points about one to two conditional s.d. from the mode, is a deliberate point on that trade-off. The localized rate (Specification §9A.4) is one way to discount flat pieces that carry little mass.

### 1.3 Sampler

One iteration from $\gamma$:

1. **(Block 1)** draw $\beta_i\sim f_i(\cdot\mid\gamma)\propto e^{q_i(\beta)}\phi_q(\beta;H_i\gamma,\Psi_i)$, independently over $i$. These draws are exact; under Assumption PQ, $f_i$ is a normal truncated to the cells of a polyhedral partition.
2. **(Block 2)** draw $\gamma'\sim\mathcal N_p\bigl(P_0^{-1}(\Lambda\mu_\gamma+\sum_iA_i^\top\beta_i),\ P_0^{-1}\bigr)$.

$K$ is the transition kernel of the $\gamma$-chain, and $\pi_\gamma$ is the marginal posterior of $\gamma$. $K$ is reversible with respect to $\pi_\gamma$ and has non-negative spectrum (Liu, Wong and Kong 1994).

### 1.4 Notation

- $\|v\|_M=(v^\top Mv)^{1/2}$.
- $\preceq$ is the Loewner order.
- $\lambda_{\min}$ and $\lambda_{\max}$ denote extreme eigenvalues.
- $\|\cdot\|_{TV}$ is the total-variation distance.
- $\mathcal L(\gamma^{(k)})$ is the law after $k$ iterations from a fixed $\gamma^{(0)}$.

---

## 2. Rate quantities

**Definition 2.1.**
$$K_i^{\rm lo}=C_{i,\min}+\Psi_i^{-1},\qquad K_i^{\rm hi}=C_{i,\max}+\Psi_i^{-1},$$
$$R_{\rm hi}=\sum_iA_i^\top(K_i^{\rm lo})^{-1}A_i,\qquad R_{\rm lo}=\sum_iA_i^\top(K_i^{\rm hi})^{-1}A_i,$$
$$\tilde A_{\rm hi}=P_0^{-1/2}R_{\rm hi}P_0^{-1/2},\qquad\tilde A_{\rm lo}=P_0^{-1/2}R_{\rm lo}P_0^{-1/2},$$
$$\rho_{\rm hi}=\lambda_{\max}(\tilde A_{\rm hi}),\qquad\rho_{\rm lo}=\lambda_{\min}(\tilde A_{\rm lo}),\qquad K_\gamma=P_0-R_{\rm hi}\succeq(1-\rho_{\rm hi})P_0 .$$


*Interpretation.* $\tilde A_{\rm hi}$ is the rate matrix of the Gaussian model in which every group's log-likelihood has curvature $C_{i,\min}$ (all groups on their flattest pieces), and $\tilde A_{\rm lo}$ the same with $C_{i,\max}$. For $p=q=1$,
$$\rho_{\rm hi}=a_{\rm hi}=\frac1{P_0}\sum_i\frac{H_i^2}{\Psi_i(1+c_{i,\min}\Psi_i)}.$$
Always $0<\rho_{\rm hi}<1-\lambda_{\min}(P_0^{-1/2}\Lambda P_0^{-1/2})<1$.

**Definition 2.2 (admissible reference set).** Let $\bar\gamma$ and $\gamma^\star$ be the mean and the mode of $\pi_\gamma$. Define
$$\mathcal C=\Bigl\{c\in\mathbb R^p:\ E_{\pi_\gamma}\|\tilde\gamma-c\|_{P_0}^2\le\frac p{1-\rho_{\rm hi}}\Bigr\}=\Bigl\{c:\ \|c-\bar\gamma\|_{P_0}^2\le\frac p{1-\rho_{\rm hi}}-\operatorname{tr}\bigl(P_0\operatorname{Cov}_{\pi_\gamma}\bigr)\Bigr\},$$
using $E\|\tilde\gamma-c\|^2_{P_0}=\operatorname{tr}(P_0\operatorname{Cov})+\|c-\bar\gamma\|^2_{P_0}$. $\mathcal C$ is a closed ellipsoid in the $P_0$-metric, centered at the posterior mean. For $p=1$ it is the interval $[\bar\gamma-r,\bar\gamma+r]$ with $r=\bigl(1/(P_0(1-\rho_{\rm hi}))-\operatorname{Var}(\pi_\gamma)\bigr)^{1/2}$. By Lemma 3.10, $\mathcal C$ is nonempty and contains the whole segment $[\gamma^\star,\bar\gamma]$, so in particular both the mode and the mean. Membership in the rest of $\mathcal C$ depends on $\operatorname{Cov}(\pi_\gamma)$, which is usually unknown. **The mode–mean segment is the part of $\mathcal C$ that can be verified in practice.**

---

## 3. Common lemmas

**Lemma 3.1 (tilt representation).** With $V_i(\beta)=-q_i(\beta)+\tfrac12\beta^\top\Psi_i^{-1}\beta$,
$$f_i(\beta\mid\gamma)\propto\exp\bigl\{\theta_i^\top\beta-V_i(\beta)\bigr\},\qquad\theta_i=A_i\gamma,$$
and $K_i^{\rm lo}\preceq\nabla^2V_i\preceq K_i^{\rm hi}$ almost everywhere.

*Proof.* Expand $-\tfrac12(\beta-H_i\gamma)^\top\Psi_i^{-1}(\beta-H_i\gamma)$. $\square$

**Lemma 3.2 (covariance sandwich).** For every $\gamma$,
$$(K_i^{\rm hi})^{-1}\preceq\operatorname{Cov}(\beta_i\mid\gamma)\preceq(K_i^{\rm lo})^{-1}.$$

*Proof.* *Upper bound:* Brascamp–Lieb. *Lower bound:* the Cramér–Rao inequality for the location family, $\operatorname{Cov}\succeq(E\nabla^2V_i)^{-1}$, which needs $V_i\in C^1$ with absolutely continuous gradient; this is where $C^1$ junctions are used. $\square$

**Lemma 3.3 (drift Jacobian; strong log-concavity of $\pi_\gamma$).** With $d(\gamma)=E[\gamma'\mid\gamma]$,
$$\nabla d(\gamma)=P_0^{-1}R(\gamma),\qquad R(\gamma)=\sum_iA_i^\top\operatorname{Cov}(\beta_i\mid\gamma)A_i,\qquad R_{\rm lo}\preceq R(\gamma)\preceq R_{\rm hi},$$
$$-\nabla^2\log\pi_\gamma(\gamma)=P_0-R(\gamma)\succeq K_\gamma .$$

*Proof.* $\nabla_\theta E[\beta\mid\theta]=\operatorname{Cov}$ for exponential tilts, and $\theta_i=A_i\gamma$. The second identity follows by differentiating $\log\int e^{q_i}\phi$ twice (Tweedie). Then apply Lemma 3.2. $\square$

**Lemma 3.4 (quantile Lipschitz bound, one dimension).** If $g_s\propto e^{sw-U(w)}$ on $\mathbb R$ with $U''\ge\kappa>0$, the $u$-quantile $w_u(s)$ satisfies $0\le\partial_sw_u\le1/\kappa$.

*Proof.* With $G(w)=\int_{-\infty}^w(\mu_s-t)g_s(t)\,dt$, one has $\partial_sw_u=G(w_u)/g_s(w_u)\ge0$. The function $H=g_s/\kappa-G$ satisfies $H'=g_s\bigl(w-U'_s(w)/\kappa-\mu_s\bigr)$, whose bracket is non-increasing. So $H'$ changes sign at most once, from $+$ to $-$, and since $H(\pm\infty)=0$, $H\ge0$. $\square$

**Lemma 3.5 (one-dimensional projections inherit strong log-concavity).** Let $X$ have density $\propto e^{\theta^\top x-V(x)}$ with $\nabla^2V\succeq K\succ0$, and let $a\in\mathbb R^q$. Then $W=a^\top X$ has a log-concave density with curvature at least $(a^\top K^{-1}a)^{-1}$. If the tilt is along $a$, $\theta=sa$, the law of $W$ is a one-parameter tilt in $s$ of a fixed base.

*Proof.* Marginals of strongly log-concave measures are strongly log-concave with the projected constant (Saumard and Wellner 2014, §3; a consequence of Prékopa–Leindler). With $\theta=sa$, $\theta^\top x=s\,w$. $\square$

**Lemma 3.5′ (projections with cell-wise curvature).** Let $X$ have density $\propto e^{\theta^\top x-V(x)}$, with $V$ convex, $C^1$ and piecewise quadratic on polyhedral cells $\{\mathcal C\}$, with Hessian $\mathbf H_{\mathcal C}\succ0$ on cell $\mathcal C$. For $a\in\mathbb R^q$, the density $h$ of $W=a^\top X$ satisfies, almost everywhere,
$$-(\log h)''(w)\ \ge\ \min_{\mathcal C\ \text{nonempty}}\bigl(a^\top\mathbf H_{\mathcal C}^{-1}a\bigr)^{-1}.$$
The minimum runs only over cells that actually exist.

*Proof.* Choose linear coordinates $x\mapsto(w,v)$ with $w=a^\top x$; the Jacobian is constant. Then $h(w)\propto\int e^{-V(w,v)}dv$, and differentiating under the integral (valid since $\nabla V$ is Lipschitz) gives
$$-(\log h)''(w)=E[V_{ww}\mid w]-\operatorname{Var}(V_w\mid w).$$
The conditional law of $v$ given $w$ has potential $V(w,\cdot)$, whose Hessian is $V_{vv}\succ0$ almost everywhere. Brascamp–Lieb in its pointwise form (extended from $C^2$ to $C^1$ piecewise-quadratic $V$ by smoothing) gives
$$\operatorname{Var}(V_w\mid w)\le E\bigl[V_{wv}V_{vv}^{-1}V_{vw}\mid w\bigr].$$
Hence $-(\log h)''\ge E\bigl[V_{ww}-V_{wv}V_{vv}^{-1}V_{vw}\mid w\bigr]$. The integrand is the Schur complement of the Hessian, which equals $\bigl((\nabla^2V)^{-1}\bigr)_{ww}^{-1}=\bigl(a^\top(\nabla^2V)^{-1}a\bigr)^{-1}$ in the original coordinates. It is constant on each cell, so the expectation is at least the minimum over nonempty cells. $\square$

Applied to $V_i$ of Lemma 3.1, with $\mathbf H_{\mathcal C}=-\nabla^2q_i|_{\mathcal C}+\Psi_i^{-1}$, this **replaces the global bound $(a^\top(K_i^{\rm lo})^{-1}a)^{-1}$ by the cell-wise one**. Lemma 3.4 then applies with that constant.

**Lemma 3.6 (tilts in $W_2$).** In the setting of Lemma 3.5,
$$\mathrm{KL}(\mu_\theta\|\mu_{\theta'})\le\tfrac12\|\theta-\theta'\|^2_{K^{-1}},\qquad W_{2,K}(\mu_\theta,\mu_{\theta'})\le\|\theta-\theta'\|_{K^{-1}},$$
and $E_{\mu_\theta}X-E_{\mu_{\theta'}}X=\bar C(\theta-\theta')$ with $0\preceq\bar C\preceq K^{-1}$.

*Proof.* The KL divergence is the Bregman divergence of the log-partition function, whose Hessian is $\operatorname{Cov}\preceq K^{-1}$. The $W_2$ bound is Talagrand's inequality in the coordinates $K^{1/2}x$. For the mean, $\bar C=\int_0^1\operatorname{Cov}_{\theta'+s(\theta-\theta')}ds$. $\square$

**Lemma 3.7 (second moment about the mode).** If $-\nabla^2\log\pi\succeq K\succ0$ on $\mathbb R^p$, with mode $x^\star$ and mean $\bar x$, then
$$E\|X-x^\star\|_K^2\le p,\qquad\operatorname{Cov}(X)\preceq K^{-1}.$$
Consequently every $c$ on the segment $[x^\star,\bar x]$ satisfies $E\|X-c\|_K^2\le p$.

*Proof.* Integration by parts gives $E[\nabla U(X)^\top(X-x^\star)]=p$, and $\nabla U(x)^\top(x-x^\star)\ge\|x-x^\star\|_K^2$. The covariance bound is Brascamp–Lieb. For $c$ on the segment, $E\|X-c\|_K^2=\operatorname{tr}(K\operatorname{Cov})+\|c-\bar x\|_K^2\le E\|X-x^\star\|_K^2$. $\square$

**Lemma 3.8 (Gaussian smoothing).**
$$\bigl\|\mathcal N(\mu,P_0^{-1})-\mathcal N(\tilde\mu,P_0^{-1})\bigr\|_{TV}=2\Phi\bigl(\tfrac12\|\mu-\tilde\mu\|_{P_0}\bigr)-1\le\frac{\|\mu-\tilde\mu\|_{P_0}}{\sqrt{2\pi}}.$$
$\square$

**Lemma 3.9 (spectral bound from contraction).** If $W_{1}(K(x,\cdot),K(y,\cdot))\le\rho\,\|x-y\|_{P_0}$ for all $x,y$, then $\rho_\star\le\rho$ for the reversible kernel $K$ (Ollivier 2009, Prop. 30). Here $\rho_\star$ is the $L^2(\pi_\gamma)$ spectral radius of $K$ on mean-zero functions. $\square$

**Lemma 3.10 (the posterior mean and mode are admissible).** Under Assumption PQ, or Assumption LC of §10, let $\bar\gamma$ and $\gamma^\star$ be the mean and the mode of $\pi_\gamma$. Then every point of the segment $[\gamma^\star,\bar\gamma]$, in particular $\bar\gamma$ and $\gamma^\star$, belongs to the admissible set $\mathcal C$ of Definition 2.2:
$$E_{\pi_\gamma}\|\tilde\gamma-c\|^2_{P_0}\le\frac p{1-\rho_{\rm hi}}\qquad\text{for all }c\in[\gamma^\star,\bar\gamma].$$

*Proof.* By Lemma 3.3, $-\nabla^2\log\pi_\gamma\succeq K_\gamma=P_0-R_{\rm hi}\succeq(1-\rho_{\rm hi})P_0=:K$. Apply Lemma 3.7 with this $K$:

- *Mean:* $E\|\tilde\gamma-\bar\gamma\|^2_K=\operatorname{tr}(K\operatorname{Cov})\le\operatorname{tr}(KK^{-1})=p$, by Brascamp–Lieb.
- *Mode:* $E\|\tilde\gamma-\gamma^\star\|^2_K\le p$, by integration by parts.
- *Points between:* $E\|\tilde\gamma-c\|^2_K=\operatorname{tr}(K\operatorname{Cov})+\|c-\bar\gamma\|^2_K\le\operatorname{tr}(K\operatorname{Cov})+\|\gamma^\star-\bar\gamma\|^2_K=E\|\tilde\gamma-\gamma^\star\|^2_K\le p$.

Finally, $\|\cdot\|^2_K=(1-\rho_{\rm hi})\|\cdot\|^2_{P_0}$. $\square$

---

## 4. Case 1: scalar $\gamma$, scalar $\beta_i$ ($p=q=1$)

Here
$$a_{\rm hi}=\frac1{P_0}\sum_i\frac{H_i^2}{\Psi_i(1+c_{i,\min}\Psi_i)},\qquad a_{\rm lo}=\frac1{P_0}\sum_i\frac{H_i^2}{\Psi_i(1+c_{i,\max}\Psi_i)},\qquad\kappa_\gamma=P_0(1-a_{\rm hi}).$$

**Theorem 4.1 (pathwise contraction).** Under the synchronous coupling (common uniforms in block 1, drawn by inversion; common normals in block 2),
$$\gamma'-\tilde\gamma'=M(\gamma-\tilde\gamma),\qquad0\le M\le a_{\rm hi}.$$

*Proof.* Lemma 3.4 applied to $\beta_i\mid\gamma$, with $s=A_i\gamma$ and $\kappa=K_i^{\rm lo}$, gives $\beta_i-\tilde\beta_i=D_iA_i(\gamma-\tilde\gamma)$ with $0\le D_i\le(K_i^{\rm lo})^{-1}$. The common noise cancels in block 2. $\square$

**Theorem 4.2 (total variation).** For $k\ge1$ and any reference point $c$ with $E_{\pi_\gamma}(\tilde\gamma-c)^2\le1/\kappa_\gamma$ (in particular any point between the mode and the mean of $\pi_\gamma$),
$$\bigl\|\mathcal L(\gamma^{(k)})-\pi_\gamma\bigr\|_{TV}\le\sqrt{\frac{P_0}{2\pi}}\;a_{\rm hi}^{\,k}\Bigl(|\gamma^{(0)}-c|+\kappa_\gamma^{-1/2}\Bigr).$$

*Proof.* Couple with a stationary copy, then apply Lemma 3.8, Theorem 4.1, Lemma 3.7 and Cauchy–Schwarz. $\square$

**Theorem 4.3 (rate sandwich).** $a_{\rm lo}\le\rho_\star\le a_{\rm hi}$.

*Proof.* *Upper bound:* Theorem 4.1 and Lemma 3.9. *Lower bound:* with $f(\gamma)=\gamma$, $\rho_\star\ge\operatorname{Cov}(X,d(X))/\operatorname{Var}(X)$. Lemma 3.3 gives $d'\ge a_{\rm lo}$, so $\operatorname{Cov}(X,d(X))=\tfrac12E[(X-X')(d(X)-d(X'))]\ge a_{\rm lo}\operatorname{Var}X$. $\square$

---

## 5. Case 2: vector $\gamma$, scalar $\beta_i$ ($p\ge1$, $q=1$)

Here $H_i$ is a $1\times p$ row, and $\tilde A_{\rm hi}=P_0^{-1/2}\bigl[\sum_iH_i^\top H_i/(\Psi_i(1+c_{i,\min}\Psi_i))\bigr]P_0^{-1/2}$.

**Theorem 5.1 (pathwise contraction).** Under the synchronous coupling,
$$\gamma'-\tilde\gamma'=M(\gamma-\tilde\gamma),\qquad0\preceq P_0^{1/2}MP_0^{-1/2}=P_0^{-1/2}\Bigl[\sum_iA_i^\top D_iA_i\Bigr]P_0^{-1/2}\preceq\tilde A_{\rm hi},$$
so $\|\gamma'-\tilde\gamma'\|_{P_0}\le\rho_{\rm hi}\|\gamma-\tilde\gamma\|_{P_0}$.

*Proof.* As for Theorem 4.1. The tilt of $\beta_i$ is one-dimensional, with $\theta_i=A_i\gamma\in\mathbb R$. $\square$

**Theorem 5.2 (total variation).** For $k\ge1$ and any $c$ with $E\|\tilde\gamma-c\|_{P_0}^2\le p/(1-\rho_{\rm hi})$, an ellipsoid containing the segment between the mode and the mean,
$$\bigl\|\mathcal L(\gamma^{(k)})-\pi_\gamma\bigr\|_{TV}\le\frac{\rho_{\rm hi}^{\,k}}{\sqrt{2\pi}}\Bigl(\|\gamma^{(0)}-c\|_{P_0}+\sqrt{\frac p{1-\rho_{\rm hi}}}\Bigr).$$

*Proof.* As for Theorem 4.2, in the $P_0$-norm, with Lemma 3.7 for $K=(1-\rho_{\rm hi})P_0$. $\square$

**Theorem 5.3 (rate sandwich).** $\rho_{\rm lo}\le\rho_\star\le\rho_{\rm hi}$.

*Proof.* *Upper bound:* as before. *Lower bound:* apply the Case 1 argument to each coordinate of $P_0^{1/2}\gamma$ and sum. By Lemma 3.3, $\tfrac12E[(X-X')^\top(\tilde d(X)-\tilde d(X'))]\ge\rho_{\rm lo}\operatorname{tr}\operatorname{Cov}X$, so some coordinate has autocorrelation at least $\rho_{\rm lo}$. $\square$

*Remark.* The rate $\rho_{\rm hi}$ governs the slowest direction. Rates in individual directions would require the matrices $P_0^{1/2}MP_0^{-1/2}$ for different pairs of states to share eigenvectors. They do in the Gaussian case, where they are constant, but not in general.

---

## 6. Case 3: scalar $\gamma$, vector $\beta_i$ ($p=1$, $q\ge1$)

Here $H_i\in\mathbb R^{q\times1}$, so each $A_i=\Psi_i^{-1}H_i$ is a column vector, and
$$a_{\rm hi}=\frac1{P_0}\sum_iA_i^\top(K_i^{\rm lo})^{-1}A_i,\qquad a_{\rm lo}=\frac1{P_0}\sum_iA_i^\top(K_i^{\rm hi})^{-1}A_i,\qquad\kappa_\gamma=P_0(1-a_{\rm hi}).$$

**Key observation.** Block 2 depends on $\beta_i$ only through the scalar $w_i=A_i^\top\beta_i$, and the tilt of $\beta_i$ is along $A_i$: $\theta_i=\gamma A_i$. By Lemma 3.5, $w_i\mid\gamma$ is a one-parameter tilt in $\gamma$ of a base with curvature at least $(A_i^\top(K_i^{\rm lo})^{-1}A_i)^{-1}$.

**Theorem 6.1 (pathwise contraction).** Couple block 1 by drawing each $w_i$ by inversion with a common uniform, and then $\beta_i\mid w_i$ with common randomness. Couple block 2 by common normals. Then
$$\gamma'-\tilde\gamma'=M(\gamma-\tilde\gamma),\qquad0\le M\le a_{\rm hi}.$$

*Proof.* By Lemma 3.4 with $\kappa=(A_i^\top(K_i^{\rm lo})^{-1}A_i)^{-1}$, $w_i-\tilde w_i=D_i(\gamma-\tilde\gamma)$ with $0\le D_i\le A_i^\top(K_i^{\rm lo})^{-1}A_i$. Block 2 gives $\gamma'-\tilde\gamma'=P_0^{-1}\sum_i(w_i-\tilde w_i)$, since the noise is common. The coupling of $\beta_i$ given $w_i$ is irrelevant to $\gamma'$. $\square$

**Theorem 6.2 (total variation).** Theorem 4.2 holds verbatim with the $a_{\rm hi}$ and $\kappa_\gamma$ of this section.

**Theorem 6.3 (rate sandwich).** $a_{\rm lo}\le\rho_\star\le a_{\rm hi}$, by the proof of Theorem 4.3 with Lemma 3.3.

*Remark.* Although $\beta_i$ is multivariate, the rate is governed by a **single direction per group**: the curvature of the pieces along $(K_i^{\rm lo})^{-1/2}A_i$. Curvature in directions to which $\gamma$ does not load affects the rate only through $(K_i^{\rm lo})^{-1}$.

---

## 7. Case 4: vector $\gamma$, vector $\beta_i$ ($p\ge1$, $q\ge1$)

**Two couplings.** Group $i$ influences $\gamma'$ through $A_i^\top\beta_i\in\mathbb R^p$, which ranges over a subspace of dimension $r_i=\operatorname{rank}(A_i)\le\min(p,q)$.

- If $r_i=1$, write $A_i=a_iv_i^\top$. The tilt is one-dimensional in $s_i=v_i^\top\gamma$, and the quantile coupling of Case 3 applies, giving a **pathwise** contraction for that group.
- If $r_i\ge2$, no monotone one-dimensional coupling is available, and the $W_2$-optimal coupling is used instead.

The rate is the same either way.

**Theorem 7.1 ($L^2$ contraction at $\rho_{\rm hi}$).** Couple block 1, for each group, by any coupling with $E\|\beta_i-\tilde\beta_i\|^2_{K_i^{\rm lo}}\le\|A_i(\gamma-\tilde\gamma)\|^2_{(K_i^{\rm lo})^{-1}}$: for example, the quantile coupling when $r_i=1$, or a $W_{2,K_i^{\rm lo}}$-optimal coupling in general. Couple block 2 by common normals. Then
$$E\|\gamma'-\tilde\gamma'\|_{P_0}^2\le\rho_{\rm hi}^2\,\|\gamma-\tilde\gamma\|_{P_0}^2 .$$
If every $r_i=1$, the contraction is pathwise: $\|\gamma'-\tilde\gamma'\|_{P_0}\le\rho_{\rm hi}\|\gamma-\tilde\gamma\|_{P_0}$.

*Proof.* With common noise, $P_0^{1/2}(\gamma'-\tilde\gamma')=P_0^{-1/2}\sum_iA_i^\top\Delta\beta_i$. Stack $v=\bigl((K_i^{\rm lo})^{1/2}\Delta\beta_i\bigr)_{i=1}^J$ and let $B=\bigl[P_0^{-1/2}A_1^\top(K_1^{\rm lo})^{-1/2},\ \dots,\ P_0^{-1/2}A_J^\top(K_J^{\rm lo})^{-1/2}\bigr]$. Then
$$P_0^{1/2}(\gamma'-\tilde\gamma')=Bv,\qquad\|Bv\|^2\le\lambda_{\max}(BB^\top)\,\|v\|^2,\qquad BB^\top=\tilde A_{\rm hi}.$$
So, pathwise, $\|\gamma'-\tilde\gamma'\|^2_{P_0}\le\rho_{\rm hi}\sum_i\|\Delta\beta_i\|^2_{K_i^{\rm lo}}$. Taking expectations and using the coupling property (Lemma 3.6 for the $W_2$-optimal coupling),
$$E\|\gamma'-\tilde\gamma'\|^2_{P_0}\le\rho_{\rm hi}\sum_i\|A_i\Delta\gamma\|^2_{(K_i^{\rm lo})^{-1}}=\rho_{\rm hi}\,\Delta\gamma^\top R_{\rm hi}\,\Delta\gamma\le\rho_{\rm hi}^2\,\|\Delta\gamma\|^2_{P_0}.$$
The rank-one statement is Theorem 5.1 or 6.1 applied group by group. $\square$

*Remark.* An earlier version split $\Delta\beta_i$ into a mean part and a fluctuation and bounded each separately, which gave $\sqrt{\rho_{\rm hi}^2+\omega\,\rho_{\rm hi}}$ with a leverage $\omega$. That bound need not be below $1$: it can exceed $1$ when $\rho_{\rm hi}>0.62$ and one group dominates. The weighted Cauchy–Schwarz step above bounds both parts at once, removes the leverage term, and gives $\rho_{\rm hi}<1$ always.

**Theorem 7.2 (total variation).** For $k\ge1$ and any $c$ with $E\|\tilde\gamma-c\|_{P_0}^2\le p/(1-\rho_{\rm hi})$,
$$\bigl\|\mathcal L(\gamma^{(k)})-\pi_\gamma\bigr\|_{TV}\le\frac{\rho_{\rm hi}^{\,k}}{\sqrt{2\pi}}\Bigl(\|\gamma^{(0)}-c\|_{P_0}+\sqrt{\frac p{1-\rho_{\rm hi}}}\Bigr).$$

*Proof.* Iterate Theorem 7.1, so $E\|\Delta_k\|_{P_0}\le(E\|\Delta_k\|_{P_0}^2)^{1/2}\le\rho_{\rm hi}^k(E\|\Delta_0\|_{P_0}^2)^{1/2}$. By Minkowski's inequality and Lemma 3.7, $(E\|\Delta_0\|^2_{P_0})^{1/2}\le\|\gamma^{(0)}-c\|_{P_0}+\sqrt{p/(1-\rho_{\rm hi})}$. Apply Lemma 3.8. $\square$

**Theorem 7.3 (rate sandwich).** $\rho_{\rm lo}\le\rho_\star\le\rho_{\rm hi}$.

*Proof.* *Upper bound:* $W_1\le W_2$, Theorem 7.1 and Lemma 3.9. *Lower bound:* as in Theorem 5.3. $\square$

*Remarks.*

1. **One rate for all four cases.** Cases 1–3 and Case 4 share the rate $\rho_{\rm hi}=\lambda_{\max}(\tilde A_{\rm hi})$, the rate of the "all groups on their flattest pieces" Gaussian model. The cases differ only in whether the contraction is pathwise (every $r_i=1$) or holds in $L^2$.
2. **Why the pathwise form matters little here.** The TV bound uses $E\|\Delta_k\|_{P_0}$, which the $L^2$ contraction controls just as well.
3. **Always below 1.** $\rho_{\rm hi}<1$ in every case, since $R_{\rm hi}\preceq\sum_iH_i^\top\Psi_i^{-1}H_i=P_0-\Lambda$.

**Remark (the admissible reference set is not a small set).** The set of admissible reference points, $\{c:E\|\tilde\gamma-c\|^2_{P_0}\le p/(1-\rho_{\rm hi})\}$, enters the proofs only through one step: bounding $E\|\gamma^{(0)}-\tilde\gamma^{(0)}\|_{P_0}$, the expected distance from the start to a stationary draw. For any point $c$ whatsoever,
$$E\|\gamma^{(0)}-\tilde\gamma^{(0)}\|_{P_0}\le\|\gamma^{(0)}-c\|_{P_0}+\bigl(E\|c-\tilde\gamma\|^2_{P_0}\bigr)^{1/2}.$$
The admissible set is simply where the second term is at most the worst-case constant $\sqrt{p/(1-\rho_{\rm hi})}$.

- **No minorization is involved.** A small set in the Meyn–Tweedie sense is a set $C$ with $K^m(x,\cdot)\ge\epsilon\,\nu(\cdot)$ for all $x\in C$. Nothing in Theorems 4.2–7.2 uses such a condition, and the bounds hold from every start.
- **What replaces drift and minorization.**
  - *Global contraction* (Theorems 4.1–7.1): two copies started anywhere move together at rate $\rho_{\rm hi}$. This holds everywhere because the curvature floor $C_{i,\min}+\Psi_i^{-1}$ is uniform, so no drift toward a small set is needed.
  - *One-shot Gaussian smoothing* (Lemma 3.8): the final $\gamma$-update turns the coupling distance into a TV bound. For any two means the Gaussian laws overlap, so this step plays the role a minorization plays in a drift–minorization argument, but quantitatively and from every start (compare the one-shot coupling of Roberts and Rosenthal 2002, and of Madras and Sezer 2010).
- **When a small set would be needed.** If the curvature floor vanishes, as for the exact logit, where $\inf_\beta nb''(\beta)=0$, contraction is no longer global. One would then fall back on a drift condition toward a region where contraction holds (Hairer and Mattingly's weak Harris theorem). In that argument a small set, or a contracting set, *is* key.

---

## 8. How the shape of the likelihood enters

All rates depend on the likelihood only through the curvature bounds $C_{i,\min}$ and $C_{i,\max}$, that is, through the flattest and the steepest pieces.

**Scalar $\beta_i$.** For pieces positioned by the rule of the companion document, in group $i$'s standardized units $s_i$ (formal expansion),
$$c_{i,\min}s_i^2\approx(1-\lambda_{0,i})-1.05\,|\varepsilon_i|-0.55\,\delta_i,\qquad c_{i,\max}s_i^2\approx(1-\lambda_{0,i})+1.05\,|\varepsilon_i|-0.55\,\delta_i .$$
Here $\varepsilon_i$ and $\delta_i$ are the standardized skewness and kurtosis of group $i$'s conditional, and $\lambda_{0,i}$ is the prior's share of the precision at the mode. Hence:

- **skewness** widens the sandwich $[\rho_{\rm lo},\rho_{\rm hi}]$;
- **heavier tails** ($\delta_i>0$) shift both ends toward $1$.

To first order, a group's contribution to $a_{\rm hi}$ changes by
$$\frac{H_i^2}{P_0\,s_i^2\,(1+\bar c_i\Psi_i)^2}\bigl(1.05\,|\varepsilon_i|+0.55\,\delta_i\bigr),\qquad\bar c_i=\frac{1-\lambda_{0,i}}{s_i^2}.$$
So the shape matters most when $\bar c_i\Psi_i$ is near $1$, that is, when likelihood and group prior carry comparable precision.

**Vector $\beta_i$.** With the projection construction, $C_{i,\min}=\sum_jc_{ij,\min}x_{ij}x_{ij}^\top$, and each covariate pattern $j$ carries its own shape numbers $\varepsilon_{ij}$ and $\delta_{ij}$ along $x_{ij}$. Skewness and kurtosis lower the curvature floor **direction by direction**. They affect the rate to the extent that those directions align with the loadings: with $(K_i^{\rm lo})^{-1/2}A_i$ in Case 3, and with the leading eigenvectors of $\tilde A_{\rm hi}$ in Cases 2 and 4.

**Gaussian limit.** If every likelihood is Gaussian, $C_{i,\min}=C_{i,\max}$, so $\tilde A_{\rm lo}=\tilde A_{\rm hi}$ is the Gaussian rate matrix and $\rho_\star=\rho_{\rm hi}$, its largest eigenvalue, in agreement with Nygren (2020).

---

## 9. What the piecewise-quadratic likelihood looks like when $\beta_i$ is a vector

Details are in the companion document; this is a summary.

- **A direct grid does not work.** Two quadratics can be joined with matching value and gradient along a hyperplane $\{n^\top\beta=c\}$ only if their Hessians differ by a multiple of $nn^\top$ (companion, Lemma 6.1). So tangent quadratics of a general log-likelihood placed on a $q$-dimensional grid cannot be joined smoothly.
- **Construction by projections.** For a generalized linear model, $\ell_i(\beta)=\sum_j\ell_{ij}(x_{ij}^\top\beta)$, summed over the $m_i$ distinct covariate patterns. Each $\ell_{ij}$ is replaced by a one-dimensional three-piece C¹ function $q_{ij}$, and $q_i(\beta)=\sum_jq_{ij}(x_{ij}^\top\beta)$. The result:
  - is concave and $C^1$;
  - is quadratic on the cells of the arrangement of the $2m_i$ hyperplanes $\{x_{ij}^\top\beta=\tau_{ij}^\pm\}$;
  - has Hessian $-\sum_jc_{ij,k}x_{ij}x_{ij}^\top$ in each cell. Crossing one hyperplane changes a single rank-one term, consistent with the join condition above.
- **Example.** Random intercept and slope with a binary covariate: patterns $x_{i1}=(1,0)$ and $x_{i2}=(1,1)$. The four lines $\beta_0=\tau_{i1}^\pm$ and $\beta_0+\beta_1=\tau_{i2}^\pm$ cut the plane into nine cells. In each cell the Hessian is $-(c_{i1,k}e_1e_1^\top+c_{i2,k'}\mathbf 1\mathbf 1^\top)$, and $C_{i,\min}=c_{i1,\min}e_1e_1^\top+c_{i2,\min}\mathbf 1\mathbf 1^\top\succ0$ because the two patterns span $\mathbb R^2$. If the patterns do not span $\mathbb R^q$, the group prior supplies the remaining curvature: $K_i^{\rm lo}=C_{i,\min}+\Psi_i^{-1}\succ0$ regardless.
- **Sampling.** Along any line the log-density is one-dimensional, $C^1$, concave and piecewise quadratic, with at most $2m_i+1$ pieces. Coordinate or direction updates are therefore exact draws by inversion. Joint exact draws are possible by enumerating cells when $m_i$ is small. (Theorems 6.1 and 7.1 assume exact draws from $f_i(\cdot\mid\gamma)$; inner sweeps add a further contraction factor, which is open.)

### 9.1 Rates for projection-built likelihoods

For $q_i(\beta)=\sum_jq_{ij}(x_{ij}^\top\beta)$, with term $j$ having piece curvatures $c_{ij,k}=n_{ij}b''(t_{ij,k})$, Assumption PQ holds with
$$C_{i,\min}=\sum_jc_{ij,\min}\,x_{ij}x_{ij}^\top,\qquad C_{i,\max}=\sum_jc_{ij,\max}\,x_{ij}x_{ij}^\top,$$
because on every cell $-\nabla^2q_i=\sum_jc_{ij,k(j)}x_{ij}x_{ij}^\top$ and each coefficient lies in $[c_{ij,\min},c_{ij,\max}]$. Hence
$$K_i^{\rm lo}=\sum_jc_{ij,\min}\,x_{ij}x_{ij}^\top+\Psi_i^{-1},$$
and all of §2–§7 applies.

- **$p=1$** (Case 3): pathwise contraction at
$$a_{\rm hi}=\frac1{P_0}\sum_iA_i^\top\Bigl(\sum_jc_{ij,\min}x_{ij}x_{ij}^\top+\Psi_i^{-1}\Bigr)^{-1}A_i .$$
- **$p\ge2$** (Case 4): $L^2$ contraction at $\rho_{\rm hi}$, pathwise for groups whose $A_i$ has rank one.

**Sharpening to the cells that exist.** $C_{i,\min}$ assumes every term is on its flattest piece at once. That combination need not occur in any cell. What the coupling proofs can use instead depends on the coupling.

1. **Rank-one groups** (all groups when $p=1$ or $q=1$; groups in $S_1$ in Case 4). The quantile coupling needs only the curvature of the one-dimensional statistic $a_i^\top\beta_i$, with $a_i$ spanning the range of $A_i$. By Lemma 3.5′ this is at least $\min_{\mathcal C}(a_i^\top\mathbf H_{i,\mathcal C}^{-1}a_i)^{-1}$ over **nonempty** cells. So group $i$'s contribution to the rate can be taken as
$$\max_{\mathcal C\ \text{nonempty}}\ A_i^\top\bigl(\textstyle\sum_jc_{ij,k_{\mathcal C}(j)}x_{ij}x_{ij}^\top+\Psi_i^{-1}\bigr)^{-1}A_i\qquad(p=1),$$
   which is never larger than $A_i^\top(K_i^{\rm lo})^{-1}A_i$. The proofs of Theorems 6.1 and 7.1 go through unchanged with this constant. **Rigorous.**
2. **Groups of rank $\ge2$** (the $L^2$ coupling). Lemma 3.6 uses Talagrand's inequality, which needs one uniform matrix bound $K$ with $\nabla^2V_i\succeq K$ everywhere; a cell-wise maximum cannot be plugged in directly. However, **any** $K_i$ satisfying $K_i\preceq\mathbf H_{i,\mathcal C}$ for every nonempty cell is valid in Lemmas 3.2, 3.6 and 3.7, and in Theorem 7.1. $K_i^{\rm lo}=C_{i,\min}+\Psi_i^{-1}$ is one such choice, but not necessarily the best. Choosing the $K_i$ to minimize $\lambda_{\max}\bigl(P_0^{-1/2}\sum_iA_i^\top K_i^{-1}A_iP_0^{-1/2}\bigr)$ subject to $K_i\preceq\mathbf H_{i,\mathcal C}$ for all nonempty cells is a convex program: $K\mapsto A^\top K^{-1}A$ is matrix-convex, and the constraints are linear matrix inequalities. Its solution is at least as good as $C_{i,\min}$. **Rigorous.**
3. **Which cells are nonempty.** For each assignment of pieces to terms, nonemptiness is a linear feasibility problem in $\beta$. With few patterns, enumerate the assignments. With many, $C_{i,\min}$ is always valid, because it covers all assignments, feasible or not.

The lower end of the sandwich does **not** sharpen cell by cell. It rests on Cramér–Rao, $\operatorname{Cov}\succeq(E\nabla^2V)^{-1}$, and an average of cell Hessians can have a smaller $a^\top(\cdot)^{-1}a$ than every individual cell. What remains valid is any **uniform upper bound** $K_i^{\rm hi}\succeq\mathbf H_{i,\mathcal C}$ over the nonempty cells. $C_{i,\max}+\Psi_i^{-1}$ is one such bound, and the tightest can again be found by a convex program.

**Dependence on the data and on the design.**

- More trials per pattern (larger $n_{ij}$) make $C_{i,\min}$ larger. That group contributes less to the rate, because its $\beta_i$ is pinned by its own data.
- Directions not spanned by the group's patterns get curvature only from $\Psi_i^{-1}$, as in a model without data in those directions.
- Skewness and kurtosis of term $j$ lower $c_{ij,\min}$ along $x_{ij}$ (§8). They raise the rate only insofar as $x_{ij}$ aligns with the loading $(K_i^{\rm lo})^{-1/2}A_i$.

### 9.2 Number of iterations required

**Corollary 9.1 (run length).** Let $\tau\in(0,1)$, and let $c$ be admissible. Set
$$D=\|\gamma^{(0)}-c\|_{P_0}+\sqrt{\frac p{1-\rho_{\rm hi}}},\qquad z_\tau=2\,\Phi^{-1}\Bigl(\frac{1+\tau}2\Bigr).$$
Then $\|\mathcal L(\gamma^{(k)})-\pi_\gamma\|_{TV}\le\tau$ for every
$$k\ \ge\ \max\Bigl\{1,\ \Bigl\lceil\frac{\log(D/z_\tau)}{\log(1/\rho_{\rm hi})}\Bigr\rceil\Bigr\}.$$
The same holds for the joint law of $(\gamma^{(k)},\beta^{(k+1)})$.

*Proof.* Solve the bounded form $2\Phi\bigl(\tfrac12\rho_{\rm hi}^kD\bigr)-1\le\tau$ for $k$. Since $z_\tau\approx\tau\sqrt{2\pi}$ for small $\tau$, the linear form of the bound gives essentially the same $k$. $\square$

**Starting at the posterior mode makes the distance term exactly zero.** By Lemma 3.7 the mode $\gamma^\star$ of $\pi_\gamma$ is admissible, so with $\gamma^{(0)}=c=\gamma^\star$,
$$D=\sqrt{\frac p{1-\rho_{\rm hi}}},\qquad k\ \ge\ \frac{\log\bigl(\sqrt{p/(1-\rho_{\rm hi})}\,/\,z_\tau\bigr)}{\log(1/\rho_{\rm hi})}.$$

- **Everything is known before sampling.** $\rho_{\rm hi}$ comes from the pieces, and $\gamma^\star$ is the mode of the marginal posterior of $\gamma$ under the three-piece model, computed deterministically by Newton's method with closed-form gradient and Hessian (Construction document, Procedure 5.4). The required $k$ is fully determined by the model and the data.
- **Which mode.** $\gamma^\star$ must be the mode of the **marginal** posterior of $\gamma$, with the $\beta_i$ integrated out, not the joint mode of $(\gamma,\beta)$. The posterior mean, or any point between mode and mean, gives the same $D$. A start elsewhere adds $\|\gamma^{(0)}-\gamma^\star\|_{P_0}$, in posterior-precision units.

**Example** (the 20-group model of `Piecewise_Gaussian_GLMM_Specification.md` §9B: $p=1$, $n_i=10$, $\Psi=0.5$, $\Lambda=0.01$, pieces from the three-case rule; computed $\rho_{\rm hi}=a_{\rm hi}=0.5443$). Started at $\gamma^\star$, $D=(1-0.5443)^{-1/2}=1.481$:

| tolerance $\tau$ | iterations $k$ |
|---|---|
| 0.05 | 5 |
| 0.01 | 7 |
| 0.001 | 11 |

Each further factor of $10$ in tolerance costs $\log10/\log(1/\rho_{\rm hi})\approx3.8$ iterations. These counts are guarantees under the worst-case rate. The actual chain may mix faster, since its rate lies in $[\rho_{\rm lo},\rho_{\rm hi}]$.

### 9.3 Justification for using $\rho_{\rm hi}$ with the true (logit) posterior

The theorems above are proved for the sampler whose group likelihoods **are** the piecewise-quadratic functions. They are intended to guide the starting point and run length of the sampler for the **true** posterior. That use rests on two arguments, which are supported numerically but not proved:

1. **A better likelihood approximation.** The three-piece fit tracks the skewed, flatter side of each group likelihood that a single Gaussian at the mode misses. A single-Gaussian rate therefore tends to be optimistic, while the piecewise rate reflects the slower mixing on the flat side.
2. **A joint worst case.** $\rho_{\rm hi}$ places every group on its flattest piece simultaneously. An individual group's true conditional can have less curvature than its flattest piece, but only for part of its mass and not for all groups at once. The joint worst case offsets this, except far from the posterior mass.

The check below illustrates argument 2.

The bound puts **every** group on its flattest piece **simultaneously**. The true logit conditional of an individual $\beta_i$ can have less curvature than its flattest piece (beyond the flat-side tangency point). But that happens one group at a time, and only for part of the conditional's mass. The joint worst case should more than offset it, except far from the posterior mass.

**Check** (positioning_rule/true_rate_check.py; the 20-group model of §9.2, data in Specification §9B; the posterior of $\gamma$ has mean $0.45$ and s.d. $0.36$). The true chain's local contraction factor is the drift slope
$$a_{\rm true}(\gamma)=\frac1{P_0}\sum_i\frac{H_i^2}{\Psi_i^2}\operatorname{Var}_{\rm true}(\beta_i\mid H_i\gamma),$$
computed by quadrature under the exact logit likelihood, compared with $a_{\rm hi}=0.5443$:

| $\gamma$ | distance from posterior mean | $a_{\rm true}(\gamma)$ | groups whose true conditional variance exceeds their piece bound |
|---|---|---|---|
| $-1$ | $4.0$ s.d. | 0.483 | 1 |
| $0$ | $1.2$ s.d. | 0.478 | 0 |
| $0.45$ | $0$ | 0.482 | 0 |
| $1$ | $1.5$ s.d. | 0.491 | 0 |
| $2$ | $4.3$ s.d. | 0.521 | 3 |
| $-3$ | $9.6$ s.d. | 0.539 | 4 |
| $3$ | $7.1$ s.d. | 0.569 | 9 |
| $5$ | $12.6$ s.d. | 0.703 | 12 |

Within about four posterior s.d. the true contraction factor stays below $a_{\rm hi}$, by about $4$–$12\%$, even where individual groups exceed their own bound. The worst case only fails far out, where most groups sit in their flat tails at once. A start at the posterior mode does not visit that region, which supports using $\rho_{\rm hi}$ for the true chain started there. **How much posterior mass lies where the bound is not conservative** (positioning_rule/tail_prob.py):

- $a_{\rm true}(\gamma)$ crosses $a_{\rm hi}$ at $\gamma=2.53$ and $\gamma=-3.14$, that is, $5.7$ and $9.9$ posterior s.d. from the mean.
- The exact (logit) posterior mass beyond these points is $7.4\times10^{-9}$.
- For the **approximate** posterior, strong log-concavity gives the rigorous bound $P(|\gamma-\bar\gamma|\ge r)\le2e^{-\kappa_\gamma r^2/2}=9.9\times10^{-7}$ at $r=5.7$ s.d. This does **not** transfer to the true posterior. The true marginal is log-concave but only $\Lambda$-strongly so, since the logit likelihood has no curvature floor, and $e^{-\Lambda r^2/2}$ is useless here.
- The $7.4\times10^{-9}$ above is a deterministic numerical computation for the true posterior, not a theorem.

**Status.** There is no proof that $\rho_{\rm hi}$ bounds the true chain's rate, and in general it does not: far enough out it fails. What can be said is quantitative but computational. The region where the conservative rate is not conservative carries negligible true posterior mass in this example ($7.4\times10^{-9}$ by quadrature).

## 10. General log-concave likelihoods: comparison with a Gaussian

The theorems of §4–§7 depend on the group likelihoods only through the curvature bounds $C_{i,\min}$ and $C_{i,\max}$, and they hold for **any** finite number of pieces. This gives statements for general likelihoods whose curvature is bounded below by that of a Gaussian likelihood with precision $C_{i,\min}$. We prove them in two ways: directly, and as the limit of piecewise-quadratic approximations, which is the route that motivated them.

**What is and is not compared.** Two kinds of comparison with the Gaussian model must be kept apart:

- *Comparisons of rates* (contraction coefficients and $L^2$ spectral radii). These are proved below: the true model's rate lies between the exact rates of two Gaussian models.
- *Comparisons of distances.* A statement such as $\|\mathcal L(\gamma^{(k)})-\pi_\gamma\|_{TV}\le\|\mathcal L^{\rm G}(\gamma^{(k)})-\pi^{\rm G}_\gamma\|_{TV}$, the actual TV of the true sampler against the actual TV of the Gaussian sampler, is **not** proved. It was attempted earlier in this project without success. The TV bound below is a coupling bound whose rate is the Gaussian rate. Its form, decaying like $\rho^k$ with a stationary-coupling constant, is weaker than the exact Gaussian analysis of Nygren (2020), which gives $a^{2k}$ decay from the center. The direct proof is available precisely because the statement is about a bound, not about the Gaussian chain's actual distance.

### 10.1 Statement

**Assumption LC (curvature bounds).** Each $\ell_i:\mathbb R^q\to\mathbb R$ is concave, with locally Lipschitz gradient, and
$$C_{i,\min}\preceq-\nabla^2\ell_i(\beta)\preceq C_{i,\max}\qquad\text{for almost every }\beta,$$
for symmetric $0\preceq C_{i,\min}\preceq C_{i,\max}$, where $C_{i,\max}$ may be $+\infty$ in some directions (no upper bound). That is, each observation carries **at least** the precision $C_{i,\min}$, and at most $C_{i,\max}$, of a Gaussian observation.

Define $K_i^{\rm lo}$, $\tilde A_{\rm hi}$, $\rho_{\rm hi}$ and $\rho_{\rm lo}$ exactly as in Definition 2.1. **These are the rates of the Gaussian hierarchical models whose group likelihoods are Gaussian with precision $C_{i,\min}$ and $C_{i,\max}$ respectively.** Write $\rho^{\rm G}(C)$ for the convergence rate of the two-block sampler in the Gaussian model with group-likelihood precisions $C=(C_i)$. Then $\rho_{\rm hi}=\rho^{\rm G}(C_{\min})$ and $\rho_{\rm lo}=\rho^{\rm G}(C_{\max})$.

**Theorem 10.1 (comparison of rates with Gaussian models).** Under Assumption LC, for the two-block Gibbs sampler with the true likelihoods $e^{\ell_i}$:

1. **Contraction.** There is a coupling with $E\|\gamma'-\tilde\gamma'\|^2_{P_0}\le\rho^{\rm G}(C_{\min})^2\,\|\gamma-\tilde\gamma\|^2_{P_0}$, pathwise when every $A_i$ has rank one.
2. **Total variation.** Let $\mathcal C$ be the admissible reference set of Definition 2.2, with $\rho_{\rm hi}=\rho^{\rm G}(C_{\min})$ and $\pi_\gamma$ the true marginal posterior. Under Assumption LC, $-\nabla^2\log\pi_\gamma\succeq(1-\rho_{\rm hi})P_0$ (Lemma 3.3), so $\mathcal C$ contains the segment between the mode and the mean of $\pi_\gamma$. Then for every $k\ge1$ and **every $c\in\mathcal C$** (the statement is made only for such $c$),
$$\bigl\|\mathcal L(\gamma^{(k)})-\pi_\gamma\bigr\|_{TV}\le\frac{\rho^{\rm G}(C_{\min})^{\,k}}{\sqrt{2\pi}}\Bigl(\|\gamma^{(0)}-c\|_{P_0}+\sqrt{\frac p{1-\rho^{\rm G}(C_{\min})}}\Bigr).$$
3. **Rate sandwich.** $\rho^{\rm G}(C_{\max})\le\rho_\star\le\rho^{\rm G}(C_{\min})$, where the lower bound requires $C_{i,\max}<\infty$.

In words:

- **Rates (items 1 and 3).** The true sampler's coupling contraction coefficient is at most the exact rate of the Gaussian model whose likelihood precision is the smallest curvature $C_{\min}$. Its asymptotic geometric ($L^2$ spectral) rate lies between the exact rates of the Gaussian models with precisions $C_{\max}$ and $C_{\min}$. These are genuine comparisons, since for a Gaussian model $\rho^{\rm G}$ **is** the exact rate.
- **Distances (item 2).** The TV bound uses the Gaussian rate $\rho^{\rm G}(C_{\min})$, but it is a bound, not a comparison with the Gaussian sampler's actual TV. It is weaker than the Gaussian TV results of Nygren (2020). It decays like $\rho^k$ rather than $a^{2k}$ from the center, and carries the constant $\sqrt{p/(1-\rho)}$ from coupling with an independent stationary copy.

For Gaussian likelihoods, $C_{\min}=C_{\max}$, and the rate statements reduce to the exact Gaussian rate. The TV bound remains weaker than the exact Gaussian TV.

**Corollary 10.2 (start inside the admissible set).**

1. *Best reference point.* Since the bound of Theorem 10.1(2) holds for every $c\in\mathcal C$, it holds with the distance term replaced by $\operatorname{dist}_{P_0}(\gamma^{(0)},\mathcal C)=\min_{c\in\mathcal C}\|\gamma^{(0)}-c\|_{P_0}$:
$$\bigl\|\mathcal L(\gamma^{(k)})-\pi_\gamma\bigr\|_{TV}\le\frac{\rho^{\rm G}(C_{\min})^{\,k}}{\sqrt{2\pi}}\Bigl(\operatorname{dist}_{P_0}(\gamma^{(0)},\mathcal C)+\sqrt{\frac p{1-\rho^{\rm G}(C_{\min})}}\Bigr).$$
2. *Start inside $\mathcal C$.* If $\gamma^{(0)}\in\mathcal C$, in particular if $\gamma^{(0)}$ is the posterior mode, the posterior mean, or any point between them, then taking $c=\gamma^{(0)}$ makes the distance term vanish:
$$\boxed{\ \bigl\|\mathcal L(\gamma^{(k)})-\pi_\gamma\bigr\|_{TV}\le\rho^{\rm G}(C_{\min})^{\,k}\sqrt{\frac{p}{2\pi\,\bigl(1-\rho^{\rm G}(C_{\min})\bigr)}}\qquad(k\ge1).\ }$$
   For $p=1$ this is $\rho^k/\sqrt{2\pi(1-\rho)}$.
3. *Run length.* With $\gamma^{(0)}\in\mathcal C$, the TV distance is at most $\tau$ once
$$k\ \ge\ \max\Bigl\{1,\ \Bigl\lceil\frac{\log\bigl(\sqrt{p/(1-\rho^{\rm G}(C_{\min}))}\,/\,z_\tau\bigr)}{\log\bigl(1/\rho^{\rm G}(C_{\min})\bigr)}\Bigr\rceil\Bigr\},\qquad z_\tau=2\Phi^{-1}\Bigl(\frac{1+\tau}2\Bigr),$$
   using the bounded form $2\Phi(\tfrac12\rho^kD)-1$ of Lemma 3.8 with $D=\sqrt{p/(1-\rho)}$.

*Proof.* Item 1: minimize the bound of Theorem 10.1(2) over $c\in\mathcal C$, a closed set, so the minimum is attained. Item 2: $\operatorname{dist}_{P_0}(\gamma^{(0)},\mathcal C)=0$. Item 3: solve $2\Phi(\tfrac12\rho^kD)-1\le\tau$ for $k$. $\square$

*Remark.* The value of the bound is the same for every start in $\mathcal C$. This is a property of the bound, which caps the spread of $\pi_\gamma$ by the worst case $p/(1-\rho)$, not of the chain: the actual distance from stationarity does depend on where in $\mathcal C$ the chain starts. For Gaussian likelihoods with $p=1$, $\operatorname{Var}(\pi_\gamma)=1/(P_0(1-\rho))$ exactly, so $\mathcal C=\{\bar\gamma\}$ and the bound distinguishes every start. For $p>1$, $\operatorname{tr}(P_0\operatorname{Cov})=\operatorname{tr}\bigl((I-\tilde A)^{-1}\bigr)\le p/(1-\rho)$, with equality only when every eigenvalue of $\tilde A$ equals $\rho$. So even for Gaussian likelihoods $\mathcal C$ is generally a nondegenerate ellipsoid.

**Remark 10.3 (when the rate is below 1, and how close to 1 it can be).**

1. **Always below 1 for a proper prior on $\gamma$.** Since $C_{i,\min}\succeq0$ implies $(K_i^{\rm lo})^{-1}=(C_{i,\min}+\Psi_i^{-1})^{-1}\preceq\Psi_i$, Loewner monotonicity gives
$$\rho^{\rm G}(C_{\min})\ \le\ \rho^{\rm G}(0)=\lambda_{\max}\Bigl(P_0^{-1/2}\sum_iH_i^\top\Psi_i^{-1}H_iP_0^{-1/2}\Bigr)=1-\lambda_{\min}\bigl(P_0^{-1/2}\Lambda P_0^{-1/2}\bigr).$$
   This is strictly below $1$ whenever $\Lambda\succ0$, that is, whenever the population prior on $\gamma$ is proper, for **every** hierarchical model with log-concave group likelihoods. For an improper flat prior ($\Lambda=0$), $\rho^{\rm G}(0)=1$ and the bound says nothing, unless $C_{i,\min}\succ0$ supplies curvature.
2. **Close to 1 for log-concave GLMs with weak priors.** The logit, probit and Poisson likelihoods have curvature tending to $0$ in the tails, so their global floor is $C_{i,\min}=0$, and the best rate Theorem 10.1 gives is $\rho^{\rm G}(0)=1-\lambda_{\min}(P_0^{-1/2}\Lambda P_0^{-1/2})$. With a weak prior, $\Lambda$ is small compared with the information $\sum_iH_i^\top\Psi_i^{-1}H_i$, and $\rho^{\rm G}(0)$ is close to $1$. For scalar $\gamma$, $\rho^{\rm G}(0)=1-\Lambda/P_0$.
3. **Illustration** (the 20-group model of §9.2: $\Lambda=0.01$, $\Psi_i=0.5$, $H_i$ evenly spaced on $[-1,1]$, so $P_0=14.75$; logit likelihoods with $n_i=10$; start in $\mathcal C$; tolerance $\tau=0.01$; Corollary 10.2(3)):

| curvature floor used | rate | iterations for $\tau=0.01$ |
|---|---|---|
| true logit likelihood, $C_{\min}=0$ | $\rho^{\rm G}(0)=1-0.01/14.75=0.99932$ | $10{,}813$ |
| three-piece approximation, $c_{i,\min}=n_iw(t_{i,\rm far})$ | $\rho_{\rm hi}=0.5443$ | $7$ |

   Both numbers are valid bounds for the model they refer to. The first is valid for the true logit model but uninformative. The second is valid for the piecewise model and, by the argument of §9.3, a conservative guide for the true model.
4. **Why the piecewise construction matters.** For log-concave GLMs the global curvature floor is $0$, and with it the bound degenerates toward the prior-only rate. The piecewise approximation replaces the global floor by the curvature of the flattest piece, which is placed where the posterior mass ends (about one to two conditional s.d. from the mode). That turns an uninformative bound into a usable one. The price is that the guarantee becomes a guarantee for the approximating model, with the transfer to the true model argued in §9.3 rather than proved.

**Proposition 10.4 (ridge-modified likelihoods: a floor that does not depend on the prior).** For $\omega\in(0,1]$ and a reference point $m_i$, define the ridge-modified log-likelihood
$$\ell_i^{\omega}(\beta)=(1-\omega)\,\ell_i(\beta)+\omega\,T_{m_i}(\beta),$$
where $T_{m_i}$ is the second-order Taylor polynomial of $\ell_i$ at $m_i$, with Hessian $-C_i^Q$, $C_i^Q=-\nabla^2\ell_i(m_i)$. A natural choice of $m_i$ is the mode of $\beta_i$ given $\eta_i=H_i\gamma^\star$, which is computed deterministically. If $\ell_i$ satisfies Assumption LC with bounds $[C_{i,\min},C_{i,\max}]$, then $\ell_i^{\omega}$ satisfies it with
$$C^{\omega}_{i,\min}=(1-\omega)C_{i,\min}+\omega C_i^Q\ \succeq\ \omega\,C_i^Q,\qquad C^{\omega}_{i,\max}=(1-\omega)C_{i,\max}+\omega C_i^Q .$$
Theorem 10.1 therefore applies to the ridge-modified model with rate
$$\rho^{\rm G}(C^{\omega}_{\min})\ \le\ \rho^{\rm G}(\omega C^Q)=\lambda_{\max}\Bigl(P_0^{-1/2}\sum_iA_i^\top\bigl(\omega C_i^Q+\Psi_i^{-1}\bigr)^{-1}A_iP_0^{-1/2}\Bigr).$$

*Proof.* $-\nabla^2\ell_i^\omega=(1-\omega)(-\nabla^2\ell_i)+\omega C_i^Q$, which is concave with the stated bounds. $\square$

**The rate stays below 1 as the prior becomes flat.** For scalar $\gamma$ and $\beta_i$, letting $\Lambda\to0$ (with $\sum_iH_i^2/\Psi_i>0$),
$$\rho^{\rm G}(\omega C^Q)\ \to\ \frac{\sum_i\dfrac{H_i^2}{\Psi_i}\cdot\dfrac1{1+\omega c_i^Q\Psi_i}}{\sum_i\dfrac{H_i^2}{\Psi_i}}\ \le\ \max_i\frac1{1+\omega c_i^Q\Psi_i}<1:$$
a weighted average of the per-group Gaussian shrinkage factors at curvature $\omega c_i^Q$. Unlike $\rho^{\rm G}(0)=1-\Lambda/P_0$, this does **not** tend to $1$ as the prior weakens. The floor $\omega C_i^Q$ comes from the data, through the curvature at $m_i$, not from the prior.

*Illustration* (`positioning_rule/ridge_rate.py`; the 20-group model of §9.2, same simulated $y_i$; logit, $n_i=10$, $\Psi_i=0.5$; $m_i$ is the mode of $\beta_i$ given $\eta_i=0.445H_i$; curvatures $c_i^Q=10\,w(m_i)\in[1.96,2.50]$; start in $\mathcal C$, $\tau=0.01$):

| $\omega$ | rate, $\Lambda=0.01$ | iterations | rate, $\Lambda=0$ (flat prior) | iterations |
|---|---|---|---|---|
| 0 (true logit) | 0.9993 | 10,813 | 1 | none |
| 0.01 | 0.98810 | 494 | 0.98877 | 526 |
| 0.02 | 0.97713 | 242 | 0.97780 | 249 |
| 0.05 | 0.94565 | 93 | 0.94629 | 94 |
| 0.10 | 0.8975 | 45 | 0.8981 | 45 |
| 0.25 | 0.7786 | 18 | 0.7791 | 18 |
| 0.50 | 0.6379 | 10 | 0.6383 | 10 |
| 1 (Laplace, single Gaussian) | 0.4687 | 6 | 0.4690 | 6 |

**Small $\omega$.** Expanding $(1+\omega c_i^Q\Psi_i)^{-1}\approx1-\omega c_i^Q\Psi_i$ (scalar case),
$$1-\rho^{\rm G}(\omega C^Q)\ \approx\ \frac{\Lambda+\omega\sum_iH_i^2c_i^Q}{P_0}.$$
So a small ridge weight acts like **adding prior precision $\omega\sum_iH_i^2c_i^Q$ to $\gamma$**: a fraction $\omega$ of the information the data carry about $\gamma$, when that information is evaluated at the modes. $\omega=0.01$ is therefore analogous to a 1% prior weight (as with the glmbayes default $\mathrm{pwt}=0.01$), but supplied through the likelihood, so it does not vanish when $\Lambda\to0$. With $\omega=0.01$ in the example: $1-\rho\approx(0.01+0.01\cdot16.73)/14.75=0.0120$, against the exact $0.0119$. Run lengths scale like $\log(\cdot)/(1-\rho)$, so $\omega=0.01$ gives about $500$ iterations for $\tau=0.01$ and about $700$ for $\tau=0.001$. That compares with $10{,}813$ (or no bound at all for a flat prior) without the ridge.

**Interpretation and trade-offs.**

1. **A rigorous rate for a model close to the true one.** The ridge-modified model is a genuine log-concave model with a positive floor. Theorem 10.1 applies to it exactly, for any prior strength, including a flat prior on $\gamma$, with no piecewise construction.
2. **The choice of $\omega$** trades fidelity against the rate. $\omega=0$ is the true model (rate near $1$ for weak priors). $\omega=1$ is the single-Gaussian (Laplace) model, whose rate is optimistic for the true posterior. Intermediate $\omega$ gives a model whose likelihood is pulled partway toward its quadratic approximation. Its posterior differs from the true one by an amount controlled, as in Theorem 5.1 of the construction document, by the differences in the groups' conditional means. To first order in $\omega$, those differences are proportional to $\omega$ times each conditional's skewness (a formal statement, not proved here).
3. **Relation to the piecewise floor.** The three-piece floor $c_{i,\min}=n_iw(t_{i,\rm far})$ corresponds to a data-adaptive ridge weight of about $c_{i,\min}/c_i^Q$, which is about $0.5$–$0.7$ in the single-group examples of the companion (for $n=10$: $1.38/2.23$ at $y=8$, $0.85/1.80$ at $y=10$, with the earlier positioning). The ridge construction makes the floor an explicit, user-chosen fraction of the curvature at the mode.
4. **Sampling the ridge-modified model.** $e^{\ell_i^\omega}\phi(\cdot;\eta,\Psi_i)$ is log-concave, but for the logit it is not a finite mixture of truncated normals. Exact block-1 draws can use likelihood-subgradient accept–reject (Nygren and Nygren 2006). The three-piece approximation applied to $\ell_i^\omega$ is another option: it inherits the floor $\omega C_i^Q$, so its pieces have curvature at least $\omega c_i^Q$.
5. **Relation to earlier work.** This is the ridge-modified model of the earlier ridge-bound documents (`Ridge_TV_Bounds_Symmetric.md`), now combined with the Gaussian-comparison rate of Theorem 10.1.

### 10.2 Proof A (direct)

Every lemma of §3 uses only the curvature bounds, never the piecewise structure.

- *Lemma 3.1* holds with $V_i=-\ell_i+\tfrac12\beta^\top\Psi_i^{-1}\beta$, and $K_i^{\rm lo}\preceq\nabla^2V_i\preceq K_i^{\rm hi}$ almost everywhere.
- *Lemmas 3.2–3.3* hold: Brascamp–Lieb for the upper bound, and Cramér–Rao for the lower bound. Cramér–Rao requires $\nabla V_i$ absolutely continuous, which Assumption LC provides.
- *Lemma 3.4* (quantile) requires only $U''\ge\kappa$ in the sense of distributions.
- *Lemma 3.5* (projection) and *Lemma 3.6* (tilts) require only $\nabla^2V\succeq K$.

The proofs of Theorems 4.1–7.3 then apply verbatim. $\square$

### 10.3 Proof B (limit of piecewise-quadratic approximations)

This is the route that motivated the result. It shows why the finite-$K$ theorems imply the general one.

1. **Approximations that respect the curvature floor.** Take scalar $\beta_i$ for clarity; vector $\beta_i$ follows with the projection construction when $\ell_i$ is a sum over covariate patterns. Let $q_i^{(K)}$ be the $K$-piece tangent construction of the companion document, with tangency points becoming dense: maximal spacing $\to0$ on compacts and the outermost points $\to\pm\infty$. Each piece's curvature is $-\ell_i''(t_{ik})\in[c_{i,\min},c_{i,\max}]$, so
$$c_{i,\min}\le c_{i,\min}^{(K)}\le c_{i,\max}^{(K)}\le c_{i,\max}\qquad\text{for every }K .$$
*Feasibility.* The C¹ switch point between adjacent tangency points lies between them whenever $\ell_i'$ is convex or concave on that interval (Proposition 3.8(1) of the companion, in general form). It is therefore enough to include the finitely many inflection points of $\ell_i'$ among the tangency points; for the logit, the single point $p=\tfrac12$. Then no adjacent pair straddles an inflection, and every switch point is feasible.

Hence $\rho_{\rm hi}^{(K)}\le\rho^{\rm G}(C_{\min})$, $\kappa_\gamma^{(K)}\ge P_0(1-\rho^{\rm G}(C_{\min}))$ and $\rho_{\rm lo}^{(K)}\ge\rho^{\rm G}(C_{\max})$, **uniformly in $K$**, by monotonicity of $c\mapsto(c+1/\Psi)^{-1}$ and of $\lambda_{\max}$ and $\lambda_{\min}$ in the Loewner order.
2. **Convergence of the model.** $q_i^{(K)}\to\ell_i$ uniformly on compacts, since tangent quadratics agree with $\ell_i$ to second order at the tangency points. Moreover $e^{q_i^{(K)}}\phi(\cdot;\eta,\Psi_i)\le e^{T_{\rm env}}\phi$ for an integrable envelope, because all pieces are concave with curvature at least $c_{i,\min}+1/\Psi_i$ after adding the prior. By dominated convergence, for every $\gamma$:
   - $f_i^{(K)}(\cdot\mid\gamma)\to f_i(\cdot\mid\gamma)$ in total variation (Scheffé);
   - the one-step kernels converge, $\|K^{(K)}(\gamma,\cdot)-K(\gamma,\cdot)\|_{TV}\to0$;
   - $\pi_\gamma^{(K)}\to\pi_\gamma$ in TV.
3. **Convergence of the $k$-step laws.** For fixed $k$ and $\gamma^{(0)}$, $\|\mathcal L^{(K)}(\gamma^{(k)})-\mathcal L(\gamma^{(k)})\|_{TV}\le\sum_{j<k}E\|K^{(K)}(\gamma^{(j)},\cdot)-K(\gamma^{(j)},\cdot)\|_{TV}\to0$, by the triangle inequality and dominated convergence.
4. **Passing the bound to the limit.** By the triangle inequality,
$$\|\mathcal L(\gamma^{(k)})-\pi_\gamma\|_{TV}\le\|\mathcal L-\mathcal L^{(K)}\|_{TV}+\|\mathcal L^{(K)}-\pi^{(K)}_\gamma\|_{TV}+\|\pi^{(K)}_\gamma-\pi_\gamma\|_{TV}.$$
The middle term is at most the finite-$K$ bound, which by step 1 is at most the bound of Theorem 10.1(2), for every $K$. The admissible reference points converge as well: the mode and mean of $\pi^{(K)}_\gamma$ converge to those of $\pi_\gamma$. Letting $K\to\infty$ removes the outer terms. The rate statements follow because the uniform constants are inherited by the limit. $\square$

*Why this works.* The finite-$K$ bound never sees $K$: it sees only the extreme curvatures. Tangent pieces never have curvature outside the range of the true likelihood's curvature, so the bound for every approximation is dominated by the same expression evaluated at the Gaussian rate $\rho^{\rm G}(C_{\min})$, and that domination survives the limit. What passes to the limit is the **bound**, with the Gaussian rate in it, not a comparison with the Gaussian sampler's actual distance from stationarity.

### 10.4 Consequences

1. **Any log-concave likelihood ($C_{i,\min}=0$): geometric ergodicity with an explicit rate.** With $C_{i,\min}=0$, the logit, probit and Poisson all qualify, because their curvature tends to $0$ in the tails. Then $K_i^{\rm lo}=\Psi_i^{-1}$ and
$$\rho_\star\le\rho^{\rm G}(0)=\lambda_{\max}\Bigl(P_0^{-1/2}\sum_iH_i^\top\Psi_i^{-1}H_iP_0^{-1/2}\Bigr)=1-\lambda_{\min}\bigl(P_0^{-1/2}\Lambda P_0^{-1/2}\bigr)<1,$$
the rate of the model with **no** likelihood information. This is explicit and always below $1$, but close to $1$ when the prior on $\gamma$ is weak. That is why the piecewise approximation, with its positive floor, is needed for useful run lengths.
2. **Logit lower bound (Böhning).** $-\ell_i''\le n_i/4$, so $\rho_\star\ge\rho^{\rm G}(n/4)$. No logit hierarchical model mixes faster than the Gaussian model in which every group has precision $n_i/4$, its maximal Fisher information.
3. **Local floors.** If the curvature bound holds only on a region $\mathcal R$, Theorem 10.1 holds for the chain while it stays in $\mathcal R$. Turning that into a global statement needs control of the probability of leaving $\mathcal R$ (§9.3; open problem 6). The piecewise construction can be read as a device that turns a local curvature floor, where the posterior mass is, into a global one.
4. **Sharpness of the rate, not of the TV bound.** The **rate** bound is attained: for Gaussian likelihoods with precision $C_{i,\min}$, $\rho_\star=\rho^{\rm G}(C_{\min})$ exactly, so no better rate follows from the curvature floor alone. The **TV bound** is not attained, even for Gaussians (convergence document, comparison in TV_Bound_Three_Piece_Gibbs.md §5).
5. **The stronger statement remains open:** $\|\mathcal L(\gamma^{(k)})-\pi_\gamma\|_{TV}\le\|\mathcal L^{\rm G}(\gamma^{(k)})-\pi^{\rm G}_\gamma\|_{TV}$ for the Gaussian model with precision $C_{\min}$ and a matched start. It would make the exact Gaussian TV formulas of Nygren (2020) directly usable for non-Gaussian likelihoods.

---

## 11. Summary

| case | contraction | rate $\rho$ | TV bound ($k\ge1$; $c$ admissible) | sandwich for $\rho_\star$ |
|---|---|---|---|---|
| 1: $p=q=1$ | pathwise | $a_{\rm hi}=\frac1{P_0}\sum\frac{H_i^2}{\Psi_i(1+c_{i,\min}\Psi_i)}$ | $\sqrt{P_0/2\pi}\,a_{\rm hi}^k(\lvert\gamma^{(0)}-c\rvert+\kappa_\gamma^{-1/2})$ | $[a_{\rm lo},a_{\rm hi}]$ |
| 2: $p\ge1$, $q=1$ | pathwise | $\rho_{\rm hi}=\lambda_{\max}(\tilde A_{\rm hi})$ | $(2\pi)^{-1/2}\rho_{\rm hi}^k(\lVert\gamma^{(0)}-c\rVert_{P_0}+\sqrt{p/(1-\rho_{\rm hi})})$ | $[\rho_{\rm lo},\rho_{\rm hi}]$ |
| 3: $p=1$, $q\ge1$ | pathwise, via the projection $A_i^\top\beta_i$ | $a_{\rm hi}=\frac1{P_0}\sum A_i^\top(K_i^{\rm lo})^{-1}A_i$ | as Case 1 | $[a_{\rm lo},a_{\rm hi}]$ |
| 4: $p,q\ge1$ | pathwise for rank-one groups; $L^2$ otherwise | $\rho_{\rm hi}=\lambda_{\max}(\tilde A_{\rm hi})$ | $(2\pi)^{-1/2}\rho_{\rm hi}^k(\lVert\gamma^{(0)}-c\rVert_{P_0}+\sqrt{p/(1-\rho_{\rm hi})})$ | $[\rho_{\rm lo},\rho_{\rm hi}]$ |

---

## 12. Open problems

1. **$a^{2k}$-type decay from the center.** The Gaussian analysis gives TV decaying like $a^{2k}$ from the mean (Nygren 2020); these bounds give $a^k$. Closing the gap requires controlling the variance deficit rather than the coupling distance.
2. **Pathwise couplings for rank $\ge2$.** A monotone multivariate coupling would make the Case 4 contraction pathwise for every group. It is not needed for the bounds stated here.
3. **Inner sweeps** in block 1 when exact draws of vector $\beta_i$ are replaced by coordinate updates.
4. **Directional rates** without commuting rate matrices.
5. **Random variance components** $\Psi_i$ and $\Lambda$.
6. **Transfer to the true posterior.** (For likelihoods with a positive global curvature floor, §10 settles this: the Gaussian comparison bound applies to the true model directly.) A proof that $\rho_{\rm hi}$, or the run length of Corollary 9.1, is conservative for the true chain, or a rigorous bound on the true posterior mass where it is not. Possible tools: perturbation bounds for contracting chains, and the log-concavity of the true marginal of $\gamma$. Not pursued here.

## References

- Liu, J. S., Wong, W. H. and Kong, A. (1994). Covariance structure of the Gibbs sampler with applications to the comparisons of estimators and augmentation schemes. *Biometrika* 81, 27–40.
- Nygren, K. N. (2020). Multivariate normal distances (manuscript dated 20 July 2020). *(Final title to be filled in.)*
- Ollivier, Y. (2009). Ricci curvature of Markov chains on metric spaces. *Journal of Functional Analysis* 256, 810–864.
- Madras, N. and Sezer, D. (2010). Quantitative bounds for Markov chain convergence: Wasserstein and total variation distances. *Bernoulli* 16, 882–908.
- Roberts, G. O. and Rosenthal, J. S. (2002). One-shot coupling for certain stochastic recursive sequences. *Stochastic Processes and their Applications* 99, 195–208.
- Saumard, A. and Wellner, J. A. (2014). Log-concavity and strong log-concavity: a review. *Statistics Surveys* 8, 45–114.
