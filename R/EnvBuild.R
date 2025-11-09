#' @name EnvelopeBuild 
#' @title 
#' GPU-Accelerated Envelope Construction for Posterior Simulation
#'
#' @details
#' Constructs an enveloping function for posterior simulation using a grid of
#' tangency points. The envelope is used in accept-reject sampling to guarantee
#' iid draws from the posterior distribution. The implementation follows
#' \insertCite{Nygren2006}{glmbayes}, with extensions for GPU acceleration
#' (via OpenCL), dynamic grid optimization, and parallelized evaluation.
#'
#' The envelope is typically built around the posterior mode \eqn{\theta^\star} for a model in standard 
#' form (which in this context means a model with a diagonal posterior precision matrix
#' and prior identity precision matrix - \code{glmb_Standardize_Model}). It uses dimension-specific width parameters 
#' \eqn{\omega_i} derived from the precision
#' matrix. Tangency points are selected per dimension, and the full grid is
#' formed via Cartesian expansion. Negative log-likelihood and gradient values
#' are computed at each grid point, either on CPU or GPU depending on the
#' \code{use_opencl} flag. These values are used to construct a piecewise
#' envelope function that dominates the posterior density.
#' 
#' @section Models in standard form:
#'
#' Following Nygren & Nygren (2006, Section 3.3), the envelope construction
#' assumes the model has been reparameterized into *standard form*. In this form:
#' \itemize{
#'   \item The prior precision matrix is the identity.
#'   \item The posterior precision matrix, evaluated at the posterior mode
#'         \eqn{\theta^\star}, is diagonal.
#' }
#'
#' Example 2 in the paper illustrates the special case of a zero-mean normal
#' prior with identity covariance. In this setting:
#' \itemize{
#'   \item The generalized likelihood-subgradient density at a tangency point
#'         \eqn{\bar\theta} has mean vector \eqn{-c(\bar\theta)} and covariance \eqn{I}.
#'   \item The normalizing integrals over restricted sets factorize across
#'         dimensions. Specifically, for a rectangular set
#'         \eqn{A = \{\theta : l_L \leq \theta \leq l_U\}}, we have
#'         \deqn{
#'           \int_{\theta \in A} q^{\bar{\theta}}(\theta)\, d\theta
#'           \;=\;
#'           \prod_{r=1}^{p}
#'           \Big[
#'             \Phi\!\big(l_{U,r} + c_{r}(\bar{\theta})\big)
#'             - \Phi\!\big(l_{L,r} + c_{r}(\bar{\theta})\big)
#'           \Big],
#'         }
#'         and the truncated expectation in coordinate \eqn{r} is
#'         \deqn{
#'           \mathbb{E}_{\tilde{q}^{\bar{\theta}}}[\theta_r \mid \theta \in A]
#'           =
#'           -\,c_{r}(\bar{\theta})
#'           +
#'           \frac{
#'             \phi\!\big(l_{L,r} + c_{r}(\bar{\theta})\big)
#'             - \phi\!\big(l_{U,r} + c_{r}(\bar{\theta})\big)
#'           }{
#'             \Phi\!\big(l_{U,r} + c_{r}(\bar{\theta})\big)
#'             - \Phi\!\big(l_{L,r} + c_{r}(\bar{\theta})\big)
#'           }.
#'         }
#' }
#'
#' These closed-form expressions explain why \code{EnvelopeBuild} evaluates
#' \code{logU}, \code{loglt}, and \code{logrt} using univariate normal CDFs
#' and densities, rather than numerical integration. The gradients
#' (\code{cbars}) directly determine the shifted means of the restricted
#' densities, and the separability across dimensions makes the grid-based
#' construction computationally tractable.
#'
#' Models with Zellner's \eqn{g}-priors are essentially in standard form,
#' since in the whitened design space both the prior and likelihood precisions
#' are diagonal. Each dimension still needs to be scaled so that the prior
#' precision is exactly the identity matrix. For other models, standard form
#' can be achieved by reparameterization (e.g. via Cholesky of the posterior
#' precision) or by shifting part of the prior quadratic form into the
#' likelihood.

#' @section Construction of restricted subgradient densities:
#'
#' Following Remark 5 in Nygren & Nygren (2006), each unrestricted
#' likelihood-subgradient density \eqn{q_{\bar{\theta}}(\cdot)} can be
#' restricted to a subset \eqn{A \subset \Theta}. The restricted density is
#' defined as
#' \deqn{
#'   \tilde{q}_{\bar{\theta}}(\theta)
#'   =
#'   \frac{ q_{\bar{\theta}}(\theta)\,\mathbf{1}_{\{\theta \in A\}} }
#'        { \int_{\theta' \in A} q_{\bar{\theta}}(\theta')\, d\theta' },
#' }
#' and the corresponding constant is
#' \deqn{
#'   \tilde{a}(\theta_bar)
#'   =
#'   a(\theta_bar)
#'   \int_{\theta \in A} q_{\theta_bar}(\theta)\, d\theta,
#' }
#' where \eqn{a(\theta_bar)} is the global normalizing constant from
#' Theorem 1. For every \eqn{\theta \in A}, the identity
#' \deqn{
#'   \tilde{a}(\theta_bar) \cdot h_{\theta_bar}(\theta) \cdot
#'   \tilde{q}_{\theta_bar}(\theta)
#'   = \pi(\theta \mid y)
#' }
#' holds, ensuring that the restricted densities reproduce the posterior
#' when combined.
#'
#' The constant \eqn{a(\bar{\theta})} is defined in Theorem 1 as
#' \deqn{
#'   a(\bar{\theta}) =
#'   \frac{g(\bar{\theta})\,\mathrm{MGF}\!\big(-c(\bar{\theta})\big)}
#'        {f(y)\,\exp\!\big(-c(\bar{\theta})^{T}\bar{\theta}\big)},
#' }
#' where \eqn{g(\bar{\theta})} is the reference density at the tangency point,
#' \eqn{c(\bar{\theta})} is the subgradient of the log-likelihood, and
#' \eqn{\mathrm{MGF}(-c(\bar{\theta}))} is the moment-generating function of
#' the prior evaluated at \eqn{-c(\bar{\theta})}.
#'
#' The envelope function \eqn{h_{\bar{\theta}}(\theta)} is defined as
#' \deqn{
#'   h_{\bar{\theta}}(\theta) =
#'   \frac{\exp\!\big(-c(\bar{\theta})^{T}\bar{\theta}\big)\,f(y\mid\theta)}
#'        {\exp\!\big(-c(\bar{\theta})^{T}\theta\big)\,g(\bar{\theta})},
#' }
#' and satisfies
#' \deqn{
#'   0 \le h_{\bar{\theta}}(\theta) \le 1
#'   \quad \forall\,\theta \in \Theta,
#'   \qquad
#'   h_{\bar{\theta}}(\bar{\theta}) = 1
#'   \quad \text{if } f(y\mid\bar{\theta}) = g(\bar{\theta}).
#' }
#'
#' In the standardized model (zero-mean normal prior with identity covariance
#' and diagonal posterior precision at the mode), the restricted integral
#' \eqn{\int_{A} q_{\bar{\theta}}(\theta)\, d\theta} factorizes across
#' dimensions and can be evaluated in closed form using normal CDFs:
#' \deqn{
#'   \int_{\theta \in A} q_{\bar{\theta}}(\theta)\, d\theta
#'   =
#'   \prod_{r=1}^{p}
#'   \Big[
#'     \Phi\!\big(l_{U,r} + c_{r}(\bar{\theta})\big)
#'     - \Phi\!\big(l_{L,r} + c_{r}(\bar{\theta})\big)
#'   \Big].
#' }
#'
#' In standardized models, where the prior is \eqn{\mathcal{N}(0, I)}, the
#' moment-generating function simplifies to
#' \deqn{
#'   \mathrm{MGF}\big(-c(\bar{\theta})\big)
#'   = \exp\left(\tfrac{1}{2} c(\bar{\theta})^{T} c(\bar{\theta})\right),
#' }
#' so the constant \eqn{a(\bar{\theta})} from Theorem 1 becomes
#' \deqn{
#'   a(\bar{\theta}) =
#'   \frac{g(\bar{\theta})\,\exp\left(\tfrac{1}{2} c(\bar{\theta})^{T} c(\bar{\theta})\right)}
#'        {f(y)\,\exp\!\big(-c(\bar{\theta})^{T} \bar{\theta}\big)}.
#' }
#'
#' This is why \code{EnvelopeBuild} computes and stores \code{logU},
#' \code{loglt}, and \code{logrt} using univariate normal CDF evaluations.
#' The constants \eqn{\tilde{a}(\bar{\theta})} are then obtained by scaling
#' the global constant with these integrals, and the mixture weights
#' (\code{PLSD}) are normalized accordingly. In practice:
#' \itemize{
#'   \item \code{Set_Grid_C2_pointwise} evaluates restricted densities
#'         at each grid point.
#'   \item \code{LLconst} stores the log of the restricted integrals.
#'   \item \code{setlogP_C2} computes \eqn{\tilde{a}(\bar{\theta})} and
#'         normalizes mixture weights.
#' }
#' 
#' @section Mixture construction and tractable probabilities:
#'
#' Claim 2 in Nygren & Nygren (2006) shows that the posterior density
#' \eqn{\pi(\theta \mid y)} can be expressed as a mixture of restricted
#' likelihood-subgradient densities. Let \eqn{A_1, \dots, A_m} be a partition
#' of the parameter space \eqn{\Theta}, and define
#' \deqn{
#'   \tilde{q}_{\bar{\theta}}(\theta)
#'   =
#'   \sum_{i=1}^{m} \tilde{p}_i \, q^{\bar{\theta}}_{A_i}(\theta),
#'   \qquad
#'   \tilde{p}_i=  \frac{ \tilde{a}_i }{ \sum_{j=1}^{k} \tilde{a}_j } 
#'   }
#' 
#' 
#'
#' In Remark 6 of Nygren & Nygren (2006), the mixture weights
#' \eqn{\tilde{p}_i} for each restricted likelihood-subgradient density
#' \eqn{q^{\bar{\theta}_i}_{A_i}} are defined as
#' \deqn{
#'   \tilde{p}_i
#'   =
#'   \frac{
#'     g(\bar{\theta}_i)\,\mathrm{MGF}\big(-c(\bar{\theta}_i)\big)
#'     \int_{\theta \in A_i} q^{\bar{\theta}_i}(\theta)\, d\theta
#'     \,/\, \exp\big(-c(\bar{\theta}_i)^{T} \bar{\theta}_i\big)
#'   }{
#'     \sum_{j=1}^{k}
#'     g(\bar{\theta}_j)\,\mathrm{MGF}\big(-c(\bar{\theta}_j)\big)
#'     \int_{\theta \in A_j} q^{\bar{\theta}_j}(\theta)\, d\theta
#'     \,/\, \exp\big(-c(\bar{\theta}_j)^{T} \bar{\theta}_j\big)
#'   }.
#' }
#' This expression reflects the full normalization of the mixture, where each
#' \eqn{\tilde{p}_i} is proportional to the restricted constant
#' \eqn{\tilde{a}_i(\bar{\theta}_i)} from Remark 5, and the denominator sums
#' over all such constants across the partition. The resulting mixture
#' \eqn{\tilde{q}^{\bar{\theta}}(\theta)} is a valid approximation to the
#' posterior density \eqn{\pi(\theta \mid y)}.
#'
#' Remark 6 emphasizes that these mixture weights are tractable to compute
#' in standardized models. When the prior is \eqn{\mathcal{N}(0, I)} and the
#' posterior precision is diagonal at the mode, each integral
#' \eqn{\int_{A_i} q^{\bar{\theta}}(\theta)\, d\theta} factorizes across
#' dimensions and can be evaluated using normal CDFs:
#' \deqn{
#'   \int_{\theta \in A_i} q^{\bar{\theta}}(\theta)\, d\theta
#'   =
#'   \prod_{r=1}^{p}
#'   \Big[
#'     \Phi\!\big(l_{U,r}^{(i)} + c_{r}(\bar{\theta})\big)
#'     - \Phi\!\big(l_{L,r}^{(i)} + c_{r}(\bar{\theta})\big)
#'   \Big],
#' }
#' where \eqn{l_{L}^{(i)}} and \eqn{l_{U}^{(i)}} are the bounds defining
#' region \eqn{A_i}.
#'
#' This tractability is central to the envelope construction. It allows
#' \code{EnvelopeBuild} to compute:
#' \itemize{
#'   \item \code{LLconst}, which stores the log of each restricted integral
#'         \eqn{\log \int_{A_i} q^{\bar{\theta}}(\theta)\, d\theta}.
#'   \item \code{PLSD}, which stores the normalized mixture weights
#'         \eqn{\tilde{p}_i}.
#' }
#' These quantities are used to construct the combined density
#' \eqn{\tilde{q}_{\bar{\theta}}(\theta)} and to evaluate the envelope
#' approximation to the posterior. Because all components are normalized and
#' tractable, the mixture is both valid and computationally efficient.
#' @section Log-scale properties of the envelope function:
#'
#' The envelope function \eqn{h_{\bar{\theta}}(\theta)} is defined in
#' Theorem 1 as
#' \deqn{
#'   h_{\bar{\theta}}(\theta) =
#'   \frac{\exp\!\big(-c(\bar{\theta})^{T} \bar{\theta}\big)\,f(y \mid \theta)}
#'        {\exp\!\big(-c(\bar{\theta})^{T} \theta\big)\,g(\bar{\theta})}.
#' }
#' When \eqn{g(\bar{\theta}) = f(y \mid \bar{\theta})}, this simplifies to
#' \deqn{
#'   h_{\bar{\theta}}(\theta) =
#'   \exp\!\big( c(\bar{\theta})^{T}(\theta - \bar{\theta}) \big)
#'   \cdot \frac{f(y \mid \theta)}{f(y \mid \bar{\theta})}.
#' }
#' Taking logarithms yields
#' \deqn{
#'   \log h_{\bar{\theta}}(\theta)
#'   =
#'   c(\bar{\theta})^{T}(\theta - \bar{\theta})
#'   + \log f(y \mid \theta) - \log f(y \mid \bar{\theta}),
#' }
#' which is tractable as long as the log-likelihood \eqn{\log f(y \mid \theta)}
#' is. In particular, if the log-likelihood is concave or piecewise affine,
#' then \eqn{\log h_{\bar{\theta}}(\theta)} inherits that structure and can be
#' efficiently evaluated across grid regions.
#'
#' This tractability is central to the envelope construction: it allows
#' \code{EnvelopeBuild} to evaluate the envelope function pointwise using
#' \code{Set_Grid_C2_pointwise}, and ensures that the resulting approximation
#' remains bounded between 0 and 1. At the tangency point, we recover
#' \deqn{
#'   h_{\bar{\theta}}(\bar{\theta}) = 1,
#' }
#' confirming that the envelope touches the posterior density exactly.
#' The key inequality that ensures envelope dominance follows from the
#' subgradient inequality for concave functions. If \eqn{\log f(y \mid \theta)}
#' is concave and \eqn{c(\bar{\theta})} is a subgradient at \eqn{\bar{\theta}},
#' then
#' \deqn{
#'   \log f(y \mid \theta)
#'   \le \log f(y \mid \bar{\theta}) + c(\bar{\theta})^{T}(\theta - \bar{\theta}),
#' }
#' which implies
#' \deqn{
#'   \log h_{\bar{\theta}}(\theta)
#'   \le 0,
#'   \qquad
#'   h_{\bar{\theta}}(\theta) \le 1.
#' }
#' This inequality guarantees that the envelope function dominates the posterior
#' density pointwise, as required by Theorem 1. Equality holds at the tangency
#' point \eqn{\theta = \bar{\theta}}, where
#' \deqn{
#'   h_{\bar{\theta}}(\bar{\theta}) = 1.
#' }
#' @section Use of the envelope during sampling:
#'
#' The functions \code{rnnorm_reg_std_cpp()} and \code{rnnorm_reg_std_cpp_parallel()}
#' use the envelope to generate posterior samples via rejection sampling. Although not exported,
#' these functions are called internally by \code{rnnorm_reg_cpp()}, which in turn is invoked by
#' the user-facing function \code{rNormal_reg()}. Together, these routines implement
#' envelope-based sampling for generalized linear models with log-concave likelihood functions
#' and multivariate normal priors.
#'
#' The envelope provides a mixture of restricted likelihood-subgradient densities,
#' each defined over a region \eqn{A_i}, with associated mixture weights
#' \eqn{\tilde{p}_i} stored in \code{PLSD}. The sampling proceeds as follows:
#'
#' \enumerate{
#'   \item A region index \eqn{J(i)} is drawn from the discrete distribution
#'         defined by \code{PLSD}.
#'   \item A candidate \eqn{\theta_i} is drawn from the restricted density
#'         \eqn{q^{\bar{\theta}_{J(i)}}_{A_{J(i)}}}, using the normal CDF bounds
#'         \code{loglt} and \code{logrt}, and subgradient vector \code{cbars}.
#'         Simulation for each dimension uses the internal C++ function
#'         \code{ctrnorm_cpp()}, which explicitly uses these inputs.
#'   \item The log-likelihood \eqn{\log f(y \mid \theta_i)} is computed and
#'         stored in \code{testll[0]} using the appropriate likelihood function
#'         \code{f2}.
#' }
#'
#' The acceptance test is performed using the inequality
#' \deqn{
#'   \log(U_2) \le \mathrm{LLconst}[J(i)] + \mathrm{cbars}[J(i), ]^{T} \theta_i
#'   + \log f(y \mid \theta_i),
#' }
#' which is equivalent to
#' \deqn{
#'   \log(U_2) \le \log f(y \mid \theta_i) - \left( \log f(y \mid \bar{\theta}_{J(i)}) - c(\bar{\theta}_{J(i)})^{T}(\theta_i - \bar{\theta}_{J(i)}) \right),
#' }
#' 
#' where:
#' \itemize{
#'   \item \code{LLconst[J(i)]} stores the precomputed quantity
#'         \eqn{-\log f(y \mid \bar{\theta}_{J(i)}) - c(\bar{\theta}_{J(i)})^{T} \bar{\theta}_{J(i)}},
#'         computed during envelope construction via \code{setlogP_C2()}.
#'   \item \code{cbars[J(i), ]} is the precomputed subgradient vector
#'         \eqn{c(\bar{\theta}_{J(i)})}, extracted via \code{cbars(J(i), _)}.
#'         It defines the exponential tilt direction used to evaluate the envelope.
#'   \item \code{testll[0]} is the log-likelihood at the candidate draw
#'         \eqn{\theta_i}, evaluated using the model specified by
#'         \code{family} and \code{link}.
#'   \item \eqn{-\log(U_2)} is the threshold from a uniform draw
#'         \eqn{U_2 \sim \mathrm{Unif}(0,1)}.
#' }
#'
#' The right-hand side of this inequality is always non-positive, and equals zero
#' when \eqn{\theta_i = \bar{\theta}_{J(i)}}. This reflects the fact that the envelope
#' is tangent to the log-likelihood at each \eqn{\bar{\theta}_j}, and lies above it elsewhere.
#'
#' This procedure guarantees that accepted samples are drawn from the posterior
#' \eqn{\pi(\theta \mid y)}. The envelope ensures bounded rejection probability,
#' and the mixture structure allows efficient sampling across regions. The output
#' \code{out} contains accepted draws, and \code{draws} records the number of
#' attempts per sample.
#'
#' The components returned by \code{EnvelopeBuild()} are used in specific steps of the
#' sampling procedure as follows:
#' \itemize{
#'   \item \code{PLSD} is used to randomly select a region index \eqn{J(i)} from the envelope mixture.
#'   \item \code{loglt} and \code{logrt} define the truncated normal bounds for each dimension,
#'         used together with \code{cbars} to generate candidate values \eqn{\theta_i}.
#'   \item \code{cbars} provides the subgradient vectors \eqn{c(\bar{\theta}_j)} used both for
#'         candidate generation and for computing the acceptance test.
#'   \item \code{LLconst} stores precomputed constants used in the acceptance inequality,
#'         avoiding recomputation of posterior terms at tangency points.
#'   \item \code{logU} stores the per-dimension log-density contributions for each region,
#'         computed during envelope setup. These values are summed to produce \code{logP},
#'         which determines the mixture weights \code{PLSD}.
#'   \item \code{logP} contains the total log-probabilities for each grid component,
#'         which are normalized to form the mixture weights \code{PLSD}.
#'   \item \code{thetabars} stores the tangency points \eqn{\bar{\theta}_j} used to define
#'         subgradients and region-specific densities.
#'   \item \code{GridIndex} encodes the sampling type (tail, center, line) used for each
#'         dimension and region, guiding how each coordinate is simulated.
#' }
#' @section Algorithmic steps (linked to theory):
#'
#' The implementation of \code{EnvelopeBuild} follows the envelope construction
#' in Nygren & Nygren (2006) for models in standard form (See Section 3-3.3). 
#' Each computational step corresponds to a theoretical guarantee:
#'
#' 1. **Compute width parameters \eqn{\omega_i} from the diagonal precision matrix.**  In particular,
#' let \eqn{\theta^{\ast}} denote the unique posterior mode. For each dimension \eqn{i},
#'  define
#' \deqn{
#'   \omega_{i} :=
#'   \frac{\sqrt{2} - \exp\!\big(-1.20491 - 0.7321\,\sqrt{0.5 - \partial^{2}\log f(\theta^{\ast}\mid y)/\partial\theta_{i}^{2}}\big)}
#'        {\sqrt{1 - \partial^{2}\log f(\theta^{\ast}\mid y)/\partial\theta_{i}^{2}}}.
#' }
#' 
#'    As seen from the above, the widths \eqn{\omega_i} are derived from the local curvature of the 
#'    log-likelihood at the posterior mode. This ensures that the three-interval construction per
#'    dimension below yields an envelope whose efficiency does not deteriorate with sample size.
#'
#' 2. **Use the width parameters to construct intervals around the posterior mode \eqn{\theta^\star}.**  Specifically,
#' we set
#' \deqn{
#'   \ell_{i,1} = \theta^{\ast}_{i} - 0.5\,\omega_{i}, \quad
#'   \ell_{i,2} = \theta^{\ast}_{i} + 0.5\,\omega_{i},
#' }
#' and construct three intervals per dimension:
#' \deqn{
#'   A_{i,1} = (-\infty,\ell_{i,1}), \quad
#'   A_{i,2} = [\ell_{i,1},\ell_{i,2}], \quad
#'   A_{i,3} = (\ell_{i,2},\infty).
#' }
#'
#' For each dimension \eqn{i}, let \eqn{J_{i} = \{1,2,3\}} and define
#' \eqn{J = \prod_{i=1}^{p} J_{i}}, which has \eqn{3^{p}} elements. Each
#' \eqn{j \in J} is a vector \eqn{(j_{1},\ldots,j_{p})}, and we define
#' \deqn{
#'   A^{\ast}_{j} = \prod_{i=1}^{p} A_{i,j_{i}}.
#' }
#' The collection \eqn{A^{\ast} = \{A^{\ast}_{j} : j \in J\}} forms a partition of \eqn{\Theta}.
#'
#'    
#' 3. **For each member of the partition, select tangency points \eqn{\theta^\star \pm \omega_i}.**
#' 
#'  For each \eqn{j \in J}, define index sets
#' \deqn{
#'   C_{j1} = \{i : j_{i} = 1\}, \quad
#'   C_{j2} = \{i : j_{i} = 2\}, \quad
#'   C_{j3} = \{i : j_{i} = 3\}.
#' }
#' 
#' The tangency points \eqn{\bar{\theta}_{j}} are then defined componentwise by
#' \deqn{
#'   \bar{\theta}_{j,i} =
#'   \begin{cases}
#'     \theta^{\ast}_{i} - \omega_{i}, & i \in C_{j1}, \\
#'     \theta^{\ast}_{i},              & i \in C_{j2}, \\
#'     \theta^{\ast}_{i} + \omega_{i}, & i \in C_{j3}.
#'   \end{cases}
#' }
#'   
#'    The tangency points are hence chosen so that the envelope touches the log-likelihood at
#'    representative points in each interval, guaranteeing dominance and tightness.
#'
#' 4. **Build the full grid of tangency points (Cartesian product across dimensions).**  
#'    
#'    The Cartesian product of per-dimension partitions yields the \eqn{3^p} restricted
#'    densities described in the paper, ensuring coverage of the full parameter space.
#'
#' 5. **Evaluate negative log-likelihood and gradients at each grid point to construct the 
#'    likelihood subgradient densities and to facilitate accept rejection sampling**  
#'      
#'    The subgradients \eqn{c(\bar{\theta})} are part of the definitions of Likelihood subgradient densities 
#'    (see definition 2 below) while both the subgradients and the negative log-likelihoods (through the function 
#'    \eqn{h_{\bar{\theta}}(.)}) are used during the accept-reject procedure. 
#'    CPU and GPU routines compute these values efficiently.
#'    - On CPU: via \code{f2_*} and \code{f3_*} routines.
#'    - On GPU: via \code{f2_f3_opencl}, which computes both together
#'    
#' 6. **Call \code{Set_Grid_C2_pointwise} to evaluate restricted multivariate normal log-densities.**  
#'    (Claim 2; Remark 5).  
#'    Each restricted density corresponds to a subset of the partition, normalized
#'    as in Remark 5.
#'
#' 7. **Call \code{setlogP_C2} to compute component log-probabilities and constants.**  
#'    (Remark 6).  
#'    The constants \eqn{\tilde{a}} and mixture weights \eqn{\tilde{p}_i} are computed
#'    explicitly as in Remark 6, ensuring that the mixture envelope is properly normalized.
#'
#' 8. **Normalize probabilities (\code{PLSD}) and optionally sort grid components.**  
#'    (Claim 2).  
#'    Normalization ensures that the mixture of restricted densities forms a valid
#'    dominating density for the posterior. Sorting is an implementation detail to
#'    improve sampling efficiency.
#' 
#' 
#' 
#' 
#' @section Formal definition and key claim:
#' \strong{Definition 2.} A probability density function \eqn{q(\cdot)} is a
#' generalized likelihood-subgradient probability density for a posterior density
#' \eqn{\pi(\cdot\mid y)} with prior density \eqn{\pi(\cdot)} and likelihood
#' function \eqn{f(y\mid\cdot)} at a point \eqn{\bar{\theta} \in \Theta} if there
#' exists a subgradient \eqn{c(\bar{\theta})} for the negative of the log of a
#' function \eqn{g} at \eqn{\bar{\theta}} such that:
#' \itemize{
#'   \item (a) \eqn{g(\cdot)} bounds \eqn{f(y\mid\cdot)} from above,
#'   \item (b) \eqn{\mathrm{MGF}(-c(\bar{\theta})) = \int_{\Theta} \exp\!\big(-c(\bar{\theta})^{T}\theta\big)\,\pi(\theta)\,d\theta} is finite,
#'   \item (c) \eqn{\forall\,\theta \in \Theta:~ q(\theta) = \exp\!\big(-c(\bar{\theta})^{T}\theta\big)\,\pi(\theta)\,/\,\mathrm{MGF}\!\big(-c(\bar{\theta})\big)}.
#' }
#'
#' \emph{Special cases:}
#' \itemize{
#'   \item If \eqn{g(\theta)} is the likelihood function, then we call \eqn{q(\theta)} a
#'         likelihood-subgradient density. Log-concave likelihood functions are ubiquitous
#'         in statistical modeling. Models with log-concave likelihood functions include
#'         the Poisson and logit regression models, as well as some survival models.
#'   \item If \eqn{g(\theta) = \bar{f} \ge f(y\mid\theta)} for all \eqn{\theta}, then we are in the
#'         Bayesian context discussed earlier and the prior is a generalized likelihood-subgradient density.
#' }
#'
#' The Appendix provides a more detailed discussion on the existence of likelihood-subgradient densities.
#' These are guaranteed to exist at every point if the prior is a finite mixture of multivariate normals
#' and the likelihood function is log-concave. The generalized likelihood-subgradient density is particularly
#' simple in the case of finite mixtures of multivariate normal priors.
#'
#' \strong{Claim 1.} Suppose that the prior \eqn{\pi(\cdot)} in Definition 2 is a finite mixture
#' of multivariate normals \eqn{\sum_{i=1}^{k} p_{i}\,\pi_{i}(\cdot\mid \mu_{i}, \Sigma_{i})}. Let
#' \eqn{c(\bar{\theta})} be a subgradient for \eqn{-\ln(g(\cdot))} at \eqn{\bar{\theta}}. Then
#' \eqn{\mathrm{MGF}(-c(\bar{\theta}))} is finite and is in the form:
#' \deqn{
#'   \mathrm{MGF}\!\big(-c(\bar{\theta})\big)
#'   = \sum_{i=1}^{k} p_{i}\,\exp\!\Big(-c(\bar{\theta})^{T}\mu_{i} + \tfrac{1}{2}\,c(\bar{\theta})^{T}\Sigma_{i}\,c(\bar{\theta})\Big).
#' }
#'
#' The corresponding generalized likelihood-subgradient density is also a mixture of multivariate normals in the form:
#' \deqn{
#'   q(\theta) = \sum_{i=1}^{k} \tilde{p}_{i}\,\pi_{i}\!\big(\theta \mid \tilde{\mu}_{i}, \Sigma_{i}\big),
#' }
#' where
#' \deqn{
#'   \tilde{p}_{i}
#'   = \frac{p_{i}\,\exp\!\Big(-c(\bar{\theta})^{T}\mu_{i} + \tfrac{1}{2}\,c(\bar{\theta})^{T}\Sigma_{i}\,c(\bar{\theta})\Big)}
#'          {\mathrm{MGF}\!\big(-c(\bar{\theta})\big)}
#'   \quad\text{and}\quad
#'   \tilde{\mu}_{i} = \mu_{i} - \Sigma_{i}\,c(\bar{\theta}).
#' }
#'
#' @section Theorem 1 (envelope dominance and equality at tangency):
#' Let \eqn{q_{\bar{\theta}}(\cdot)} be a generalized likelihood-subgradient
#' density as in Definition 2. Define:
#' \deqn{
#'   a(\bar{\theta}) =
#'   \frac{g(\bar{\theta})\,\mathrm{MGF}\!\big(-c(\bar{\theta})\big)}
#'        {f(y)\,\exp\!\big(-c(\bar{\theta})^{T}\bar{\theta}\big)}
#' }
#' and
#' \deqn{
#'   h_{\bar{\theta}}(\theta) =
#'   \frac{\exp\!\big(-c(\bar{\theta})^{T}\bar{\theta}\big)\,f(y\mid\theta)}
#'        {\exp\!\big(-c(\bar{\theta})^{T}\theta\big)\,g(\bar{\theta})}.
#' }
#'
#' Then
#' \deqn{
#'   a(\bar{\theta})\,q_{\bar{\theta}}(\theta)
#'   \;\ge\;
#'   a(\bar{\theta})\,h_{\bar{\theta}}(\theta)\,q_{\bar{\theta}}(\theta)
#'   \;=\;
#'   \pi(\theta\mid y),
#' }
#' and
#' \deqn{
#'   0 \le h_{\bar{\theta}}(\theta) \le 1
#'   \quad \forall\,\theta \in \Theta.
#' }
#'
#' Finally, if \eqn{f(y\mid\bar{\theta}) = g(\bar{\theta})}, then
#' \deqn{
#'   h_{\bar{\theta}}(\bar{\theta}) = 1.
#' }
#'
#' @section Claim 2 (mixtures over a partition) and Remark 5 (restricted densities):
#' \strong{Claim 2.} Let \eqn{A_{1},A_{2},\ldots,A_{k}} be a finite partition of \eqn{\Theta} and
#' let \eqn{\tilde{q}_{1}(\cdot),\tilde{q}_{2}(\cdot),\ldots,\tilde{q}_{k}(\cdot)} be associated restricted densities
#' such that there exist associated constants \eqn{\tilde{a}_{1},\tilde{a}_{2},\ldots,\tilde{a}_{k}} and functions
#' \eqn{\tilde{h}_{1}(\cdot),\tilde{h}_{2}(\cdot),\ldots,\tilde{h}_{k}(\cdot)} satisfying the following:
#' \itemize{
#'   \item (a) \eqn{\theta \in A_{i} \Rightarrow \tilde{a}_{i}\,\tilde{h}_{i}(\theta)\,\tilde{q}_{i}(\theta) = \pi(\theta \mid y)}, \eqn{i=1,2,\ldots,k}.
#'   \item (b) \eqn{\theta \in A_{i} \Rightarrow 0 \le \tilde{h}_{i}(\theta) \le 1}, \eqn{i=1,2,\ldots,k}.
#' }
#' Define a new density by
#' \deqn{
#'   \tilde{q}(\theta) = \sum_{i=1}^{k} \tilde{p}_{i}\,\tilde{q}_{i}(\theta),
#' }
#' where \eqn{\tilde{p}_{i} = \tilde{a}_{i} \big/ \big(\sum_{j=1}^{k} \tilde{a}_{j}\big)}.
#' Let \eqn{\tilde{a} = \sum_{j=1}^{k} \tilde{a}_{j}}, and let \eqn{\tilde{h}(\theta)} be a function satisfying
#' \eqn{\theta \in A_{i} \Rightarrow \tilde{h}(\theta) = \tilde{h}_{i}(\theta)}, \eqn{i=1,2,\ldots,k}.
#' Then
#' \deqn{
#'   \tilde{a}\,\tilde{q}(\theta) \;\ge\; \tilde{a}\,\tilde{h}(\theta)\,\tilde{q}(\theta) \;=\; \pi(\theta \mid y)
#' }
#' and
#' \deqn{
#'   0 \le \tilde{h}(\theta) \le 1 \quad \forall\,\theta \in \Theta.
#' }
#' In other words, a bounding function for the full space can be constructed by combining bounding
#' functions for the individual elements of the partition. Our main interest is in mixtures of restricted
#' generalized likelihood-subgradient densities. \emph{Remark 5} shows how to use generalized
#' likelihood-subgradient densities to construct restricted densities of the form required for Claim 2.
#'
#' \strong{Remark 5.} Let \eqn{q_{\bar{\theta}}(\cdot)} be a generalized likelihood-subgradient density as in Theorem 1.
#' Define a restricted density \eqn{\tilde{q}_{\bar{\theta}}(\cdot)} and a corresponding constant \eqn{\tilde{a}(\bar{\theta})}
#' on a set \eqn{A \subset \Theta} by
#' \deqn{
#'   \tilde{q}_{\bar{\theta}}(\theta) =
#'   \frac{q_{\bar{\theta}}(\theta)}{\int_{\theta \in A} q_{\bar{\theta}}(\theta)\,d\theta}
#' }
#' and
#' \deqn{
#'   \tilde{a}(\bar{\theta}) =
#'   a(\bar{\theta}) \int_{\theta \in A} q_{\bar{\theta}}(\theta)\,d\theta.
#' }
#' Then for every \eqn{\theta \in A}, we have
#' \deqn{
#'   \tilde{a}(\bar{\theta})\,h_{\bar{\theta}}(\theta)\,\tilde{q}_{\bar{\theta}}(\theta) = \pi(\theta \mid y).
#' }
#'
#'
#' @section Remark 6 (explicit mixture weights and constants):
#' Suppose that the restricted densities in Claim 2 are restricted generalized
#' likelihood-subgradient densities as in Remark 5. Then both the overall
#' constant \eqn{\tilde{a}} and the mixture probabilities \eqn{\tilde{p}_{i}}
#' admit explicit formulas:
#'
#' \deqn{
#'   \tilde{a}
#'   = \frac{1}{f(y)}
#'     \sum_{j=1}^{k}
#'       g(\bar{\theta}_{j})\,
#'       \mathrm{MGF}\!\big(-c(\bar{\theta}_{j})\big)\,
#'       \int_{\theta \in A_{j}} q_{\bar{\theta}_{j}}(\theta)\,d\theta\,
#'       \exp\!\big(-c(\bar{\theta}_{j})^{T}\bar{\theta}_{j}\big)
#' }
#'
#' and
#'
#' \deqn{
#'   \tilde{p}_{i}
#'   =
#'   \frac{
#'     g(\bar{\theta}_{i})\,
#'     \mathrm{MGF}\!\big(-c(\bar{\theta}_{i})\big)\,
#'     \int_{\theta \in A_{i}} q_{\bar{\theta}_{i}}(\theta)\,d\theta\,
#'     \exp\!\big(-c(\bar{\theta}_{i})^{T}\bar{\theta}_{i}\big)
#'   }{
#'     \sum_{j=1}^{k}
#'       g(\bar{\theta}_{j})\,
#'       \mathrm{MGF}\!\big(-c(\bar{\theta}_{j})\big)\,
#'       \int_{\theta \in A_{j}} q_{\bar{\theta}_{j}}(\theta)\,d\theta\,
#'       \exp\!\big(-c(\bar{\theta}_{j})^{T}\bar{\theta}_{j}\big)
#'   }.
#' }
#'
#' The expression for \eqn{\tilde{a}} provides insight into how the partition
#' \eqn{\{A_{j}\}} and the positioning of the tangency points
#' \eqn{\bar{\theta}_{j}, j=1,2,\ldots,k}, affect its value.
#' In fact, the optimal placement of the tangencies for each element of the
#' partition should satisfy the property that the tangency points are the
#' expectations of the resulting restricted likelihood-subgradient densities.

#' @section Example 2 (standard normal prior, restricted set):
#' Suppose that the prior \eqn{\pi(\cdot)} in Definition 2 is a
#' \eqn{p}-dimensional multivariate normal density with mean vector 0
#' and variance-covariance matrix \eqn{I}, the identity matrix.
#' Let \eqn{q_{\bar{\theta}}(\cdot)} be a generalized likelihood-subgradient
#' density at \eqn{\bar{\theta}}.
#'
#' It is straightforward to verify that \eqn{q_{\bar{\theta}}(\cdot)} has
#' mean vector \eqn{-c(\bar{\theta})} and variance-covariance matrix \eqn{I}.
#'
#' Define a restricted set
#' \deqn{
#'   A = \{\theta \in \Theta : \ell^{L} \le \theta \le \ell^{U}\},
#' }
#' for some vectors \eqn{\ell^{L}} and \eqn{\ell^{U}}.
#'
#' Then
#' \deqn{
#'   \int_{\theta \in A} q_{\bar{\theta}}(\theta)\,d\theta
#'   = \prod_{r=1}^{p} \Big[ \Phi(\ell^{U}_{r} + c_{r}(\bar{\theta}))
#'                          - \Phi(\ell^{L}_{r} + c_{r}(\bar{\theta})) \Big],
#' }
#' where \eqn{\Phi(\cdot)} denotes the standard normal cumulative distribution function.
#'
#' The expectation under the restricted density is
#' \deqn{
#'   \mathbb{E}_{\tilde{q}_{\bar{\theta}}}[\theta]
#'   = -c(\bar{\theta})
#'     - \lambda\!\big(\ell^{L} + c(\bar{\theta}),\,\ell^{U} + c(\bar{\theta})\big),
#' }
#' where \eqn{\lambda(\cdot)} is a vector-valued function with
#' \eqn{r}th component given by Mills' ratio:
#' \deqn{
#'   \lambda_{r}
#'   = \frac{\varphi(\ell^{L}_{r} + c_{r}(\bar{\theta}))
#'          - \varphi(\ell^{U}_{r} + c_{r}(\bar{\theta}))}
#'          {\Phi(\ell^{U}_{r} + c_{r}(\bar{\theta}))
#'          - \Phi(\ell^{L}_{r} + c_{r}(\bar{\theta}))},
#' }
#' with \eqn{\varphi(\cdot)} the standard normal density.
#'
#' This example illustrates that in the standard normal prior case,
#' the generalized likelihood-subgradient density remains normal with
#' shifted mean \eqn{-c(\bar{\theta})}, and that restriction to a box
#' set \eqn{A} yields closed-form expressions for both the normalizing
#' constant and the expectation.
#' @section Remarks on sampling from restricted normals:
#' \strong{Remark 7.} Sampling from the restricted normal densities in
#' Example 2 can be implemented using the inverse-transform method
#' (see, e.g., Fishman 1999).
#'
#' \strong{Remark 8.} In many applications, the inverse-transform method
#' of sampling from the restricted density in Example 2 will require
#' evaluating the cumulative normal distribution function (or its logarithm)
#' in the extreme tail of a normal distribution. Accurate computation in
#' this regime requires numerical procedures with uniformly small relative
#' errors. Authors presenting such procedures include Hart (1957, 1966)
#' and Bryc (2002).
#' 
#' 
#' @section Theorem 2 (log-concave univariate models with normal priors):
#' In the univariate case with a normal prior (variance = 1) and normal data,
#' the partitioning approach requires careful positioning of the intervals and
#' corresponding restricted likelihood-subgradient densities.
#'
#' Empirical investigation shows:
#' \itemize{
#'   \item A single optimally positioned likelihood-subgradient density deteriorates
#'         in performance as the number of data points increases.
#'   \item The same deterioration occurs for optimally positioned two-interval partitions.
#'   \item Remarkably, the optimal three-interval partition does not suffer this deterioration:
#'         the enveloping function remains a close approximation even as the sample size grows.
#' }
#'
#' Define the posterior mode \eqn{\theta^{\ast}} and set
#' \deqn{
#'   \omega :=
#'   \frac{\sqrt{2} - \exp\!\big(-1.20491 - 0.7321\,(0.5 - \partial^{2}\log f(\theta^{\ast}\mid y)/\partial\theta^{2})\big)}
#'        {1 - \partial^{2}\log f(\theta^{\ast}\mid y)/\partial\theta^{2}}.
#' }
#'
#' Then define the partition points
#' \deqn{
#'   \ell_{1} := \theta^{\ast} - 0.5\,\omega, \quad
#'   \ell_{2} := \theta^{\ast} + 0.5\,\omega,
#' }
#' and the three regions
#' \deqn{
#'   A_{1} := (-\infty,\ell_{1}), \quad
#'   A_{2} := [\ell_{1},\ell_{2}], \quad
#'   A_{3} := (\ell_{2},\infty).
#' }
#' The tangency points are chosen as
#' \deqn{
#'   \bar{\theta}_{1} = \theta^{\ast} - \omega, \quad
#'   \bar{\theta}_{2} = \theta^{\ast}, \quad
#'   \bar{\theta}_{3} = \theta^{\ast} + \omega.
#' }
#'
#' Note: if the data represent \eqn{N} observations from a normal density with unit variance,
#' then \eqn{-\partial^{2}\log f(\theta^{\ast}\mid y)/\partial\theta^{2} = N}.
#'
#' \strong{Theorem 2.} Consider this normal data model and let \eqn{\tilde{a}^{\ast}(N)}
#' denote the value of \eqn{\tilde{a}} at sample size \eqn{N}. Then
#' \deqn{
#'   \lim_{N \to \infty} \tilde{a}^{\ast}(N) = \frac{2}{\sqrt{\pi}}.
#' }
#'
#' \emph{Proof sketch.} By symmetry, \eqn{\tilde{a}_{1}(\bar{\theta}_{1}) = \tilde{a}_{3}(\bar{\theta}_{3})}.
#' It follows that
#' \deqn{
#'   \lim_{N \to \infty} \tilde{a}^{\ast}(N)
#'   = \lim_{N \to \infty} \big(\tilde{a}_{2}(\bar{\theta}_{2}) + 2\,\tilde{a}_{3}(\bar{\theta}_{3})\big)
#'   = \frac{1}{\sqrt{\pi}} + \frac{2}{2\sqrt{\pi}}
#'   = \frac{2}{\sqrt{\pi}},
#' }
#' where the second-to-last equality follows from Claims A.1 and A.4 in the Appendix.
#' @section Log-concave models with multivariate normal priors:
#' To ensure that the quality of the enveloping function remains good in the
#' multivariate case, models are first reparameterized into a standard form.
#'
#' \strong{Definition 3.} A probability model with a multivariate normal prior
#' and log-concave likelihood function is in \emph{standard form} if:
#' \itemize{
#'   \item (a) the prior variance-covariance matrix is the identity matrix, and
#'   \item (b) the Hessian of the log-posterior density evaluated at the posterior mode
#'         is a diagonal matrix.
#' }
#'
#' \emph{Remarks on reparameterization:}
#' \itemize{
#'   \item \strong{Remark 11.} If the log-likelihood function is concave and twice
#'         continuously differentiable, then a Cholesky decomposition of the posterior
#'         precision at the unique posterior mode can be used to reparameterize the model
#'         so that the posterior precision at the mode is diagonal.
#'   \item \strong{Remark 12.} For any positive definite matrix \eqn{P}, there exists a
#'         positive definite diagonal matrix \eqn{D} such that \eqn{P - D} is also positive definite.
#'   \item \strong{Remark 13.} Let \eqn{P} and \eqn{D} be as above. Then the following two
#'         models have the same posterior density:
#'         \enumerate{
#'           \item A model with prior mean vector \eqn{\mu}, prior precision matrix \eqn{P},
#'                 and log-likelihood function \eqn{LL(\theta)}.
#'           \item A model with prior mean vector \eqn{\mu}, prior precision matrix \eqn{D},
#'                 and log-likelihood function
#'                 \deqn{
#'                   LL^{\ast}(\theta) = -\tfrac{1}{2}(\theta - \mu)^{T}(P - D)(\theta - \mu) + LL(\theta).
#'                 }
#'         }
#'   \item \strong{Remark 14.} If a model has a multivariate normal prior with a diagonal
#'         variance-covariance matrix and the posterior precision at the posterior mode is
#'         also diagonal, then the model can be reparameterized into standard form.
#'   \item \strong{Remark 15.} If a probability model with a multivariate normal prior and
#'         log-concave likelihood function has a twice continuously differentiable log-likelihood,
#'         then it can be reparameterized into standard form.
#' }
#'
#' \emph{Multivariate partition construction:}
#'
#' Let \eqn{\theta^{\ast}} denote the unique posterior mode. For each dimension \eqn{i},
#' define
#' \deqn{
#'   \omega_{i} :=
#'   \frac{\sqrt{2} - \exp\!\big(-1.20491 - 0.7321\,(0.5 - \partial^{2}\log f(\theta^{\ast}\mid y)/\partial\theta_{i}^{2})\big)}
#'        {1 - \partial^{2}\log f(\theta^{\ast}\mid y)/\partial\theta_{i}^{2}}.
#' }
#'
#' Then set
#' \deqn{
#'   \ell_{i,1} = \theta^{\ast}_{i} - 0.5\,\omega_{i}, \quad
#'   \ell_{i,2} = \theta^{\ast}_{i} + 0.5\,\omega_{i},
#' }
#' and construct three intervals per dimension:
#' \deqn{
#'   A_{i,1} = (-\infty,\ell_{i,1}), \quad
#'   A_{i,2} = [\ell_{i,1},\ell_{i,2}], \quad
#'   A_{i,3} = (\ell_{i,2},\infty).
#' }
#'
#' For each dimension \eqn{i}, let \eqn{J_{i} = \{1,2,3\}} and define
#' \eqn{J = \prod_{i=1}^{p} J_{i}}, which has \eqn{3^{p}} elements. Each
#' \eqn{j \in J} is a vector \eqn{(j_{1},\ldots,j_{p})}, and we define
#' \deqn{
#'   A^{\ast}_{j} = \prod_{i=1}^{p} A_{i,j_{i}}.
#' }
#' The collection \eqn{A^{\ast} = \{A^{\ast}_{j} : j \in J\}} forms a partition of \eqn{\Theta}.
#'
#' For each \eqn{j \in J}, define index sets
#' \deqn{
#'   C_{j1} = \{i : j_{i} = 1\}, \quad
#'   C_{j2} = \{i : j_{i} = 2\}, \quad
#'   C_{j3} = \{i : j_{i} = 3\}.
#' }
#'
#' The tangency points \eqn{\bar{\theta}_{j}} are then defined componentwise by
#' \deqn{
#'   \bar{\theta}_{j,i} =
#'   \begin{cases}
#'     \theta^{\ast}_{i} - \omega_{i}, & i \in C_{j1}, \\
#'     \theta^{\ast}_{i},              & i \in C_{j2}, \\
#'     \theta^{\ast}_{i} + \omega_{i}, & i \in C_{j3}.
#'   \end{cases}
#' }
#'
#' \strong{Remark 16.} The mixture-generalized likelihood-subgradient density
#' that results from this construction is a mixture of restricted multivariate
#' normal densities, for which a straightforward sampling procedure exists.
#' 
#' 
#' @section Subgradient density formulation:
#' Each grid component corresponds to a tilted multivariate normal density,
#' normalized using the moment-generating function (MGF). In the single-point
#' case, centered at the posterior mode \eqn{\theta^\star}, the density is:
#' \deqn{
#' f(\theta) = \frac{1}{(2\pi)^{p/2} |A|^{-1/2} \cdot \text{MGF}_A(c)} \exp\left( -\frac{1}{2} (\theta - \mu)^T A (\theta - \mu) + c^T (\theta - \theta^\star) \right)
#' }
#' where:
#' - \eqn{A} is the precision matrix,
#' - \eqn{\mu} is the prior mean vector,
#' - \eqn{c} is the gradient of the log-likelihood at \eqn{\theta^\star},
#' - \eqn{\text{MGF}_A(c)} is the moment-generating function:
#' \deqn{
#' \text{MGF}_A(c) = \exp\left( \frac{1}{2} c^T A^{-1} c \right)
#' }
#'
#' This closed-form density dominates the posterior locally and is used when
#' \code{Gridtype = 1}. For richer envelopes, multiple such components are
#' constructed at tangency points \eqn{\theta_j}, each with its own gradient
#' \eqn{c_j}, and combined into a mixture:
#' \deqn{
#' f_{\text{env}}(\theta) = \sum_{j=1}^{K} p_j f_j(\theta)
#' }
#' where the weights \eqn{p_j} are computed using log-CDF differences and constants:
#' \deqn{
#' \log p_j = \log \Phi(U_j) - \log \Phi(L_j) - \text{NegLL}_j + \text{LLconst}_j
#' }
#' @section Algorithmic steps:
#' 1. Compute width parameters \eqn{\omega_i} from the diagonal of the precision matrix.
#' 2. Construct intervals around the posterior mode \eqn{\theta^\star}.
#' 3. Select tangency points at the mode and at \eqn{\theta^\star \pm \omega_i}.
#' 4. Build the full grid of tangency points (Cartesian product across dimensions).
#' 5. Evaluate negative log-likelihood and gradients at each grid point:
#'    - On CPU: via \code{f2_*} and \code{f3_*} routines.
#'    - On GPU: via \code{f2_f3_opencl}, which computes both in parallel.
#' 6. Call \code{Set_Grid_C2_pointwise} to evaluate restricted multivariate normal
#'    log-densities in parallel.
#' 7. Call \code{setlogP_C2} to compute component log-probabilities and constants.
#' 8. Normalize probabilities (\code{PLSD}) and optionally sort grid components
#'    by probability if \code{sortgrid = TRUE}.
#' @section Gridtype logic:
#' The \code{Gridtype} argument controls how many tangency points are used per dimension:
#' - 1: Threshold rule. If \eqn{1 + a_i \le 2/\sqrt{\pi}}, use a single-point envelope at the mode;
#'      otherwise use three points.
#' - 2: Dynamic optimization via \code{EnvelopeOpt}, which balances grid build cost and
#'      expected acceptance rate. Grid size is scaled by \code{n} and the number of
#'      OpenCL cores when GPU is enabled.
#' - 3: Always use three points per dimension.
#' - 4: Always use a single point (mode only).
#' @section Supported families and links:
#' The following families and link functions are supported:
#' - Binomial: logit, probit, cloglog
#' - Quasibinomial: logit, probit
#' - Poisson: log
#' - Quasipoisson: log
#' - Gamma: log
#' - Gaussian: identity
#'
#' GPU acceleration (\code{use_opencl = TRUE}) is available for all of the above
#' except Gaussian, which is always evaluated on CPU.
#' @section GPU acceleration:
#' When \code{use_opencl = TRUE}, likelihood and gradient evaluations are
#' offloaded to the GPU using OpenCL. This can substantially reduce runtime for
#' high-dimensional models or large grids. Results are mathematically equivalent
#' to the CPU version, but small numerical differences may occur due to
#' floating-point arithmetic. If reproducibility across hardware is critical,
#' prefer the CPU path.
#'
#' If OpenCL support was not detected at compile time, the flag is ignored and
#' the CPU implementation is used. Diagnostic messages are printed when
#' \code{verbose = TRUE}.
#' @section Verbose output:
#' When \code{verbose = TRUE}, the function prints:
#' - Grid type, number of draws, OpenCL usage, and detected core count.
#' - Grid size after expansion.
#' - Time-stamped messages when entering the grid loop, starting likelihood
#'   evaluations, starting gradient evaluations, and invoking GPU kernels.
#' - Messages when setting grid values, computing log-probabilities, and sorting.
#' 
#' Any constants needed by the sampling are added to a list and returned.
#'
#' @param bStar     Point at which envelope should be centered (typically posterior mode).
#' @param A         Diagonal precision matrix for the log-likelihood in standard form.
#' @param y         A vector of observations of length \code{m}.
#' @param x         A design matrix of dimension \code{m × p}.
#' @param mu        A vector giving the prior means of the variables.
#' @param P         Prior precision matrix of the variables (positive-definite).
#' @param alpha     Offset vector.
#' @param wt        A vector of weights.
#' @param family    Family for the envelope: \code{binomial}, \code{quasibinomial}, \code{poisson}, \code{quasipoisson}, or \code{Gamma}.
#' @param link      Link function ("logit", "probit", "cloglog" for binomial; "log" for Poisson/Gamma).
#' @param Gridtype  Method to determine the number of subgradient densities in the grid.
#' @param n         Number of draws from the posterior (used for grid sizing).
#' @param n_envopt Effective sample size passed to EnvelopeOpt for grid construction.
#'   Defaults to match `n`. Larger values encourage tighter envelopes.
#' @param sortgrid  Logical; if \code{TRUE}, sort the envelope descending by component probability.
#' @param use_opencl Logical; if \code{TRUE}, use OpenCL for gradient evaluations.
#' @param verbose   Logical; if \code{TRUE}, print progress messages.
#'
#' @param GridIndex A matrix indicating, for each grid component, whether the component
#'   lies in the left tail, center, or right tail of the density. Rows correspond to
#'   grid components; columns correspond to standardized variables.
#' @param cbars     A matrix containing the subgradient of the (adjusted) negative log-likelihood
#'   at each grid component.
#' @param Lint      A matrix storing the lower and upper bounds for each grid component,
#'   depending on whether sampling is from the left, center, or right.
#'
#' @param logP      A matrix (typically two columns) with information for each grid component.
#'   The first column usually holds the output from \code{Set_Grid()}, corresponding to
#'   the restricted normal density.
#' @param NegLL     A vector of negative log-likelihood evaluations at each grid component.
#' @param G3        A matrix of tangency points used in the grid.
#'    
#' @return
#' \describe{
#'
#'   \item{\code{EnvelopeBuild()}}{A list of envelope components used for accept-reject sampling:
#'     \describe{
#'       \item{\code{GridIndex}}{Integer matrix encoding sampling type (tail, center, line) per dimension and region.}
#'       \item{\code{thetabars}}{Matrix of tangency points \eqn{\bar{\theta}_j} for each grid region.}
#'       \item{\code{cbars}}{Matrix of subgradients \eqn{c(\bar{\theta}_j)} of the negative log-likelihood at tangency.}
#'       \item{\code{loglt}}{Matrix of log left-tail probabilities per dimension and region.}
#'       \item{\code{logrt}}{Matrix of log right-tail probabilities per dimension and region.}
#'       \item{\code{logU}}{Matrix of selected per-dimension log-density contributions (tail/center) for each region.}
#'       \item{\code{logP}}{Matrix of total log-probabilities per region (first column); used to derive mixture weights.}
#'       \item{\code{PLSD}}{Vector of normalized mixture weights over grid regions used to draw region indices.}
#'       \item{\code{LLconst}}{Vector of acceptance-test constants per region used in the inequality for rejection sampling.}
#'     }
#'   }
#'
#'   \item{\code{Set_Grid()}}{A list of matrices computed for grid-based log-density evaluation:
#'     \describe{
#'       \item{\code{Down}}{Lower bounds for truncated-normal evaluation per dimension and region.}
#'       \item{\code{Up}}{Upper bounds for truncated-normal evaluation per dimension and region.}
#'       \item{\code{lglt}}{Log left-tail probabilities (from \eqn{(-\infty, \mathrm{Up}]}) per dimension and region.}
#'       \item{\code{lgrt}}{Log right-tail probabilities (from \eqn{[\mathrm{Down}, \infty)}) per dimension and region.}
#'       \item{\code{lgct}}{Log central-interval probabilities (from \eqn{[\mathrm{Down}, \mathrm{Up}]}) per dimension and region.}
#'       \item{\code{logU}}{Selected log-probability per grid cell based on \code{GridIndex} (tail or center).}
#'       \item{\code{logP}}{Matrix with row-wise sums of \code{logU} (first column) used to form mixture weights.}
#'     }
#'   }
#'
#'   \item{\code{setlogP()}}{A list with updated mixture-weight and acceptance constants:
#'     \describe{
#'       \item{\code{logP}}{Input \code{logP} with its second column populated by the log of unnormalized visit probabilities per region (mixture denominators).}
#'       \item{\code{LLconst}}{Vector of acceptance constants \eqn{-\log f(y \mid \bar{\theta}_j) - c(\bar{\theta}_j)^{T}\bar{\theta}_j} used in the accept-reject test.}
#'     }
#'   }
#' }
#' 
#'  @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt


#' @usage EnvelopeBuild(bStar,A,y,x,mu,P,alpha,wt,family = "binomial",link = "logit", 
#' Gridtype = 2L,n = 1L,n_envopt=NULL,sortgrid = FALSE,use_opencl = FALSE,verbose = FALSE)
#' @rdname EnvelopeBuild
#' @export
EnvelopeBuild <- function(
    bStar, A, y, x, mu, P, alpha, wt,
    family     = "binomial",
    link       = "logit",
    Gridtype   = 2L,
    n          = 1L,
    n_envopt   = NULL,       # effective sample size for EnvelopeOpt
    sortgrid   = FALSE,
    use_opencl = FALSE,
    verbose    = FALSE
) {
  # normalize n_envopt: if not supplied, fall back to n
  if (is.null(n_envopt)) {
    n_envopt <- n
  }
  
  # validate: must be a single non‑NA integer >= 1
  if (length(n_envopt) != 1L || is.na(n_envopt) || n_envopt < 1) {
    stop("`n_envopt` must be a positive integer scalar.")
  }
  # coerce safely to integer
  n_envopt <- as.integer(n_envopt)
  
  if (family == "gaussian") {
    return(.EnvelopeBuild_Ind_Normal_Gamma(
      bStar, A, y, x, mu, P, alpha, wt,
      family = family, link = link,
      Gridtype = Gridtype, n = n,    
      n_envopt  = n_envopt, 
      sortgrid = sortgrid,
      use_opencl = use_opencl,
      verbose    = verbose
    ))
  }
  
  .EnvelopeBuild_cpp(
    bStar, A, y, x, mu, P, alpha, wt,
    family    = family,
    link      = link,
    Gridtype  = Gridtype,
    n         = n,
    n_envopt  = n_envopt,
    sortgrid  = sortgrid,
    use_opencl = use_opencl,
    verbose    = verbose
  )
}



#' @usage Set_Grid(GridIndex, cbars, Lint)
#' @rdname EnvelopeBuild
#' @export
Set_Grid <- function(GridIndex, cbars, Lint) {
  .Set_Grid_cpp(GridIndex, cbars, Lint)
}



#' @usage setlogP(logP, NegLL, cbars, G3)
#' @rdname EnvelopeBuild
#' @export
#' @keywords internal
setlogP <- function(logP, NegLL, cbars, G3) {
  .setlogP_cpp(logP, NegLL, cbars, G3)
}