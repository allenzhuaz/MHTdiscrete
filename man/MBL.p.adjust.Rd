% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FDRSD.R
\name{MBL.p.adjust}
\alias{MBL.p.adjust}
\title{The adjusted p-values for Modified Benjamini-Liu (BL) step-down FDR controlling procedure.}
\usage{
MBL.p.adjust(p, p.set, alpha, make.decision)
}
\arguments{
\item{p}{numeric vector of p-values (possibly with \code{\link[base]{NA}}s). Any other R is coerced by \code{\link[base]{as.numeric}}. Same as in \code{\link[stats]{p.adjust}}.}

\item{p.set}{a list of numeric vectors, where each vector is the vector of all attainable p-values containing the available p-value for the corresponding hypothesis.}

\item{alpha}{significant level used to compare with adjusted p-values to make decisions, the default value is 0.05.}

\item{make.decision}{logical; if  \code{TRUE}, then the output include the decision rules compared adjusted p-values with significant level \eqn{\alpha}}
}
\value{
A numeric vector of the adjusted p-values (of the same length as \code{p}).
}
\description{
The function for calculating the adjusted p-values based on original available p-values and all attaianble p-values.
}
\note{
The MBL procedure for discrete data controls FDR under the specific dependence assumption where the joint distribution of statistics from true nulls are independent of the joint distribution of statistics from false nulls.
}
\examples{
p <- c(pbinom(1,8,0.5),pbinom(1,5,0.75),pbinom(1,6,0.6))
p.set <-list(pbinom(0:8,8,0.5),pbinom(0:5,5,0.75),pbinom(0:6,6,0.6))
MBL.p.adjust(p,p.set)
}
\references{
Benjamini, Y., and Liu, W. (1999).
A step-down multiple hypotheses testing procedure that controls the false discovery rate under independence.
 \emph{Journal of Statistical Planning and Inference}, \strong{82}: 163-170.
}
\seealso{
\code{\link{MBH.p.adjust}},  \code{\link{MBY.p.adjust}}
}
\author{
Yalin Zhu
}
