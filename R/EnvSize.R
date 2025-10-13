#' Envelope Sizing and Optimization
#'
#' \code{EnvelopeSize()} is the high‑level entry point that constructs
#' per‑dimension grids and expected draw counts, while \code{EnvelopeOpt()}
#' performs the adaptive optimization used when \code{Gridtype = 2}.
#'
#' These functions implement the grid sizing logic used in envelope construction
#' for rejection sampling. They make use of the theory described in
#' \insertCite{Nygren2006}{glmbayes} and the general implementation outlined in
#' \insertCite{glmbayesSimmethods}{glmbayes}. 
#' 
#' @param a Numeric vector of diagonal precisions for the log‑likelihood
#'   (posterior precision is \eqn{1 + a_i}).
#' @param G1 Numeric matrix of candidate grid points (3 × l1).
#' @param Gridtype Integer code controlling grid sizing logic:
#'   \itemize{
#'     \item 1 = static threshold test
#'     \item 2 = adaptive optimization via \code{EnvelopeOpt()}
#'     \item 3 = always three‑point grid
#'     \item 4 = always single‑point grid
#'   }
#' @param n Integer; number of posterior draws to generate (used for grid sizing).
#' @param n_envopt Integer; effective sample size passed to \code{EnvelopeOpt}.
#'   Defaults to -1, which means "use n".
#' @param use_opencl Logical; if \code{TRUE}, attempt GPU acceleration.
#' @param verbose Logical; if \code{TRUE}, print progress messages.
#'
#' @param a1 Numeric vector of diagonal elements of the data precision matrix
#'   (used by \code{EnvelopeOpt}).
#' @param core_cnt Integer; number of OpenCL cores or parallel workers available
#'   (default 1). When >1, envelope build cost is scaled down to reflect
#'   parallel construction.
#'
#' @section Gridtype Logic and Candidates per Draw:
#'
#' The envelope sizing logic follows the analysis of \insertCite{Nygren2006}{glmbayes}.
#'
#' \describe{
#'
#'   \item{Gridtype 1: Static Threshold}{
#'     For each dimension \eqn{i}, if
#'     \eqn{\sqrt{1 + a_i} \leq 2/\sqrt{\pi} \approx 1.128379},
#'     then a single tangent at the posterior mode suffices.
#'     Expected candidates per draw in that dimension:
#'     \eqn{\sqrt{1 + a_i}}.  
#'
#'     Otherwise, a symmetric three‑point envelope is used at
#'     \eqn{\{\theta^\star_i - \omega_i, \theta^\star_i, \theta^\star_i + \omega_i\}},
#'     with expected candidates per draw bounded above by
#'     \eqn{2/\sqrt{\pi}}.
#'   }
#'
#'   \item{Gridtype 2: Adaptive Optimization}{
#'     Each dimension is assigned either a single‑point or three‑point envelope
#'     by minimizing
#'     \deqn{T_\mathrm{total}(g_i) = T_\mathrm{build}(g_i) + T_\mathrm{sample}(n, acc_i(g_i)).}
#'
#'     The optimizer balances build cost (grows with number of tangents) against
#'     sampling cost (decreases as acceptance improves).  
#'
#'     Expected candidates per draw:
#'     \eqn{\prod_j \text{scaleest}[i,j]}, where each factor is either
#'     \eqn{\sqrt{1+a_j}} (single‑point) or \eqn{2/\sqrt{\pi}} (three‑point),
#'     depending on the optimization outcome.
#'   }
#'
#'   \item{Gridtype 3: Always Three‑Point}{
#'     Every dimension uses a symmetric three‑point envelope.  
#'     Expected candidates per draw:
#'     \deqn{\left(\tfrac{2}{\sqrt{\pi}}\right)^k}
#'     for \eqn{k} dimensions, as shown in Theorem 3 of
#'     \insertCite{Nygren2006}{glmbayes}.
#'   }
#'
#'   \item{Gridtype 4: Always Single‑Point}{
#'     Every dimension uses a single tangent at the posterior mode.  
#'     Expected candidates per draw:
#'     \deqn{\prod_{i=1}^k \sqrt{1 + a_i}}
#'     (Example 1 in \insertCite{Nygren2006}{glmbayes}).
#'   }
#'
#' }
#'
#' @details
#' \code{EnvelopeSize()} returns the constructed grid (\code{G2}),
#' index vectors (\code{GIndex1}), expected draw count (\code{E_draws}),
#' and the per‑dimension grid index.  
#'
#' \code{EnvelopeOpt()} implements the adaptive optimization used in
#' \code{Gridtype = 2}, ranking dimensions by posterior variance and
#' promoting them to three‑point tangents when the tradeoff is favorable.
#'
#' @return
#' \describe{
#'   \item{\code{EnvelopeSize()}}{A list with components \code{G2}, \code{GIndex1},
#'   \code{E_draws}, and \code{gridindex}.}
#'   \item{\code{EnvelopeOpt()}}{An integer vector of length \eqn{l1} with entries
#'   1 (single‑point) or 3 (three‑point).}
#' }
#'
#' @seealso \code{\link{EnvelopeBuild}}, \code{\link{EnvelopeSort}}
#'
#' @references
#' \insertAllCited{}
#' @example inst/examples/Ex_EnvelopeOpt.R
#' @export
#' @keywords internal

#' @usage EnvelopeSize(a, G1, Gridtype = 2L, n = 1000L, n_envopt = -1,
#'                     use_opencl = FALSE, verbose = FALSE)
#' @rdname EnvelopeSize
#' @export
EnvelopeSize <- function(a,
                         G1,
                         Gridtype   = 2L,
                         n          = 1000L,   # <-- updated default
                         n_envopt   = -1,
                         use_opencl = FALSE,
                         verbose    = FALSE) {
  .EnvelopeSize(a, G1, Gridtype, n, n_envopt, use_opencl, verbose)
}


#' @usage EnvelopeOpt(a1,n,core_cnt=1L)
#' @rdname EnvelopeSize
#' @export


EnvelopeOpt<-function(a1,n,core_cnt=1L){
  
  core_cnt <- as.integer(core_cnt)
  if (is.na(core_cnt) || core_cnt < 1L) core_cnt <- 1L
  
  
  a1rank<-rank(1/(1+a1))
  l1<-length(a1)
  
  dimcount<-matrix(0,(l1+1),l1)
  scaleest<-matrix(0,(l1+1),l1)
  intest<-c(1:(l1+1))
  slopeest<-c(1:(l1+1))
  
  dimcount[1,]<-diag(diag(l1))
  scaleest[1,]<-sqrt(1+a1)
  slopeest[1]<-prod(scaleest[1,])
  
  for(i in 2:(l1+1)){
    dimcount[i,]<-dimcount[i-1,]
    scaleest[i,]<-scaleest[i-1,]
    for(j in 1:l1){
      if(a1rank[j]==i-1){ 
        dimcount[i,j]<-3
        scaleest[i,j]<-2/sqrt(pi) 
      }
    }
##    intest[i]<-3^(i-1)
    intest[i]<-(3^(i-1))
    slopeest[i]<-prod(scaleest[i,])
  }
  evalest<-(intest/core_cnt)+n*slopeest
  minindex<-0
  for(j in 1:(l1+1)){if(evalest[j]==min(evalest)){minindex<-j}}
  
##  message("Estimated draws per Acceptance: ", slopeest[minindex])
  
  dimcount[minindex,]
  
}

