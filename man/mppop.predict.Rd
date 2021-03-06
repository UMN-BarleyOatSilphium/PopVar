% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mppop_predict.R
\name{mppop.predict}
\alias{mppop.predict}
\alias{mppop_predict2}
\title{Predict genetic variance and genetic correlations in multi-parent populations
using a deterministic equation.}
\usage{
mppop.predict(
  G.in,
  y.in,
  map.in,
  crossing.table,
  parents,
  n.parents = 4,
  tail.p = 0.1,
  self.gen = 10,
  DH = FALSE,
  models = c("rrBLUP", "BayesA", "BayesB", "BayesC", "BL", "BRR"),
  n.core = 1,
  ...
)

mppop_predict2(
  M,
  y.in,
  marker.effects,
  map.in,
  crossing.table,
  parents,
  n.parents = 4,
  tail.p = 0.1,
  self.gen = 10,
  DH = FALSE,
  models = c("rrBLUP", "BayesA", "BayesB", "BayesC", "BL", "BRR"),
  n.core = 1,
  ...
)
}
\arguments{
\item{G.in}{See \code{G.in} in \code{\link[PopVar]{pop.predict}}.}

\item{y.in}{See \code{y.in} in \code{\link[PopVar]{pop.predict}}.}

\item{map.in}{See \code{map.in} in \code{\link[PopVar]{pop.predict}}.}

\item{crossing.table}{See \code{crossing.table} in \code{\link[PopVar]{pop.predict}}.}

\item{parents}{See \code{parents} in \code{\link[PopVar]{pop.predict}}.}

\item{n.parents}{Integer number of parents per cross. May be 2 or 4. If \code{crossing.table} is passed,
this argument is ignored.}

\item{tail.p}{See \code{tail.p} in \code{\link[PopVar]{pop.predict}}.}

\item{self.gen}{The number of selfing generations in the potential cross. Can be an integer or \code{Inf} for
recombinant inbreds. Note: \code{self.gen = 1} corresponds to an F2 population.}

\item{DH}{Indicator if doubled-haploids are to be induced after the number of selfing generations indicated by
\code{self.gen}. For example, if \code{self.gen = 0} and \code{DH = TRUE}, then doubled-haploids are assumed
to be induced using gametes from F1 plants.}

\item{models}{See \code{models} in \code{\link[PopVar]{pop.predict}}.}

\item{n.core}{Number of cores for parallelization. Parallelization is supported
only on a Linux or Mac OS operating system; if working on a Windows system, the function
is executed on a single core.}

\item{...}{Additional arguments to pass depending on the choice of \code{model}.}

\item{M}{A Matrix of marker genotypes of dimensions \code{nLine} x \code{nMarker}, coded as
-1, 0, and 1.}

\item{marker.effects}{A data frame of marker effects. The first column should include the marker name and
subsequent columns should include the marker effects. Supercedes \code{y.in} if passed.}
}
\description{
Predicts the genotypic mean, genetic variance, and
usefulness criterion (superior progeny mean) in a set of multi-parent populations
using marker effects and a genetic map. If more than two traits are specified,
the function will also return predictions of the genetic correlation in the population
and the correlated response to selection.
}
\details{
Predictions are based on the deterministic equations specified by Allier et al. (2019).

The \code{mppop.predict} function takes similarly formatted arguments as the \code{\link[PopVar]{pop.predict}} function
in the \code{PopVar} package. For the sake of simplicity, we also include the \code{mppop_predict2} function, which
takes arguments in a format more consistent with other genomewide prediction packages/functions.

If you select a \code{model} other than "rrBLUP", you must specify the following additional arguments:
\itemize{
\item{\code{nIter}: See \code{\link[PopVar]{pop.predict}}}.
\item{\code{burnIn}: See \code{\link[PopVar]{pop.predict}}}.
}
}
\references{
Allier, A., L. Moreau, A. Charcosset, S. Teyssèdre, and C. Lehermeier, 2019 Usefulness Criterion and Post-selection Parental
Contributions in Multi-parental Crosses: Application to Polygenic Trait Introgression. G3 (Bethesda) 9: 1469–1479.
https://doi.org/https://doi.org/10.1534/g3.119.400129
}
