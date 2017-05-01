#' The adjusted p-values for Sidak single-step FWER controlling procedure.
#'
#' The function for calculating the adjusted p-values based on original available p-values.
#'
#' @usage
#' Sidak.p.adjust(p, alpha, make.decision)
#' @param p numeric vector of p-values (possibly with \code{\link[base]{NA}}s). Any other R is coerced by \code{\link[base]{as.numeric}}. Same as in \code{\link[stats]{p.adjust}}.
#' @param alpha significant level used to compare with adjusted p-values to make decisions, the default value is 0.05.
#' @param make.decision logical; if  \code{TRUE}, then the output include the decision rules compared adjusted p-values with significant level \eqn{\alpha}
#' @return
#' A numeric vector of the adjusted p-values (of the same length as \code{p}) if  \code{make.decision = FALSE}, or a list including original p-values, adjusted p-values and decision rules if \code{make.decision = TRUE}.
#' @seealso  \code{\link[stats]{p.adjust}}.
#' @author Yalin Zhu
#' @examples
#' p <- c(pbinom(1,8,0.5),pbinom(1,5,0.75),pbinom(1,6,0.6))
#' Sidak.p.adjust(p)
#' @export

Sidak.p.adjust <- function(p, alpha = 0.05, make.decision = FALSE){
  adjP <- 1-(1-p)^(length(p))
  if (make.decision==FALSE){
    return(adjP)
  } else{
    return(list(method= "Sidak", significant.level = alpha, Result = data.frame(raw.p = p, adjust.p=adjP, decision=ifelse(adjP<=alpha, "reject","accept"))))
  }
}

#' The adjusted p-values for Modified Bonferroni single-step FWER controlling procedure.
#'
#' The function for calculating the adjusted p-values based on original available p-values and all attaianble p-values
#'
#' @usage
#' MBonf.p.adjust(p, p.set, alpha, make.decision)
#' @param p numeric vector of p-values (possibly with \code{\link[base]{NA}}s). Any other R is coerced by \code{\link[base]{as.numeric}}. Same as in \code{\link[stats]{p.adjust}}.
#' @param p.set a list of numeric vectors, where each vector is the vector of all attainable p-values containing the available p-value for the corresponding hypothesis.
#' @param alpha significant level used to compare with adjusted p-values to make decisions, the default value is 0.05.
#' @param make.decision logical; if  \code{TRUE}, then the output include the decision rules compared adjusted p-values with significant level \eqn{\alpha}
#' @note The attainable p-value refers to the element of domain set of p-value for the corresponding hypothesis. For continuous test statistics, the p-value under true null are uniform distributed in (0,1), thus the p-values are attainable everywhere between 0 and 1. But for discrete test statistics, the p-value can only take finite values bewtween 0 and 1, that is the attainable p-values for discrete case are finite and countable, so we can assign them in a finite list \code{p.set}.
#' @return
#' A numeric vector of the adjusted p-values (of the same length as \code{p}) if  \code{make.decision = FALSE}, or a list including original p-values, adjusted p-values and decision rules if \code{make.decision = TRUE}.
#' @seealso \code{\link{Tarone.p.adjust}},  \code{\link{MixBonf.p.adjust}},  \code{\link[stats]{p.adjust}}.
#' @author  Yalin Zhu
#' @examples
#' p <- c(pbinom(1,8,0.5),pbinom(1,5,0.75),pbinom(1,6,0.6))
#' p.set <-list(pbinom(0:8,8,0.5),pbinom(0:5,5,0.75),pbinom(0:6,6,0.6))
#' MBonf.p.adjust(p,p.set)
#' @export

MBonf.p.adjust <- function(p, p.set, alpha = 0.05, make.decision = FALSE){
  m <- length(p)
  adjP <- numeric(m); pCDF <- matrix(NA,m,m)
  for(i in 1:m){
    for(j in 1:m){
      pCDF[i,j] <- max(p.set[[j]][p.set[[j]] <= p[i]],0)
    }
    adjP[i] <- min(sum(pCDF[i,]),1)
  }
  if (make.decision==FALSE){
  return(adjP)
  } else{
    return(list(method= "Modified Bonferroni", significant.level = alpha, Result = data.frame(raw.p = p, adjust.p=adjP, decision=ifelse(adjP<=alpha, "reject","accept"))))
  }
}

#' The adjusted p-values for Mixed Bonferroni single-step FWER controlling procedure.
#'
#' The function for calculating the adjusted p-values based on original available p-values and the attaianble p-values for the discrete test statistics.
#'
#' @usage
#' MixBonf.p.adjust(pc, pd, pd.set, alpha, make.decision)
#' @param pc numeric vector of the available p-values (possibly with \code{\link[base]{NA}}s) for the continuous test statistics. Any other R is coerced by \code{\link[base]{as.numeric}}. Same as in \code{\link[stats]{p.adjust}}.
#' @param pd numeric vector of the available p-values (possibly with \code{\link[base]{NA}}s) for the discrete test statistics. Any other R is coerced by \code{\link[base]{as.numeric}}. Same as in \code{\link[stats]{p.adjust}}.
#' @param pd.set a list of numeric vectors, where each vector is the vector of all attainable p-values containing the available p-value for the corresponding hypothesis for discrete data.
#' @param alpha significant level used to compare with adjusted p-values to make decisions, the default value is 0.05.
#' @param make.decision logical; if  \code{TRUE}, then the output include the decision rules compared adjusted p-values with significant level \eqn{\alpha}
#' @return
#' A numeric vector of the adjusted p-values (of the same length as \code{p}) if  \code{make.decision = FALSE}, or a list including original p-values, adjusted p-values and decision rules if \code{make.decision = TRUE}.
#' @seealso \code{\link{Tarone.p.adjust}},  \code{\link{MBonf.p.adjust}},  \code{\link[stats]{p.adjust}}.
#' @note The arguments include three parts, the available p-values need to be reorganized in advance. Gather all available p-values for continuous data as \code{pc}, and all available p-values for discrete data as \code{pd}. The attainable p-value refers to the element of domain set of p-value for the corresponding hypothesis for discrete test statistics, the p-value can only take finite values bewtween 0 and 1, that is, the attainable p-values for discrete case are finite and countable, so we can assign them in a finite list \code{pd.set}. The function returns the  adjusted p-values with the first part for continuous data of the same length as \code{pc}, and second part for discrete data of the same length as \code{pd}
#' @author  Yalin Zhu
#' @export
MixBonf.p.adjust <- function(pc, pd, pd.set, alpha = 0.05, make.decision = FALSE){
  mc <- length(pc); md <- length(pd); m <- mc+md
  p <- c(pc,pd)
  adjP <- numeric(m); pCDF <- matrix(NA,m,md)
  for(i in 1:m){
    for(j in 1:md){
      pCDF[i,j] <- max(pd.set[[j]][pd.set[[j]] <= p[i]],0)
    }
    adjP[i] <- min(mc*p[i]+sum(pCDF[i,]),1)
  }
  if (make.decision==FALSE){
    return(adjP)
  } else{
    return(list(method= "Mixed Bonferroni", significant.level = alpha, Result = data.frame(raw.p = p, adjust.p=adjP, decision=ifelse(adjP<=alpha, "reject","accept"))))
  }
}

#' The adjusted p-values for Tarone's single-step FWER controlling procedure.
#'
#' The function for calculating the adjusted p-values based on original available p-values and all attaianble p-values.
#'
#' @usage
#' Tarone.p.adjust(p, p.set, alpha, make.decision)
#' @param p numeric vector of p-values (possibly with \code{\link[base]{NA}}s). Any other R is coerced by \code{\link[base]{as.numeric}}. Same as in \code{\link[stats]{p.adjust}}.
#' @param p.set a list of numeric vectors, where each vector is the vector of all attainable p-values containing the available p-value for the corresponding hypothesis.
#' @param alpha significant level used to compare with adjusted p-values to make decisions, the default value is 0.05.
#' @param make.decision logical; if  \code{TRUE}, then the output include the decision rules compared adjusted p-values with significant level \eqn{\alpha}
#' @return
#' A numeric vector of the adjusted p-values (of the same length as \code{p}) if  \code{make.decision = FALSE}, or a list including original p-values, adjusted p-values and decision rules if \code{make.decision = TRUE}.
#' @seealso \code{\link{MBonf.p.adjust}},  \code{\link{MixBonf.p.adjust}},  \code{\link[stats]{p.adjust}}.
#' @author Yalin Zhu
#' @references
#' Tarone, R. E. (1990).
#' A modified Bonferroni method for discrete data.
#' \emph{Biometrics}, \strong{46}: 515-522.
#'
#' @examples
#' p <- c(pbinom(1,8,0.5),pbinom(1,5,0.75),pbinom(1,6,0.6))
#' p.set <-list(pbinom(0:8,8,0.5),pbinom(0:5,5,0.75),pbinom(0:6,6,0.6))
#' Tarone.p.adjust(p,p.set)
#' @export

Tarone.p.adjust <- function(p, p.set, alpha = 0.05, make.decision = FALSE){
  m <- length(p); adjP <- numeric(m)
  minP <- sort(sapply(p.set,min))
  for (j in 1:m){
    for (i in 1:m ){
      if(p[j]>=max(minP)){q <- m}
      else if (p[j]>=minP[i] & p[j]<minP[i+1]){q <- i}
    }
    adjP[j] <- min(1,q*p[j])
  }
  if (make.decision==FALSE){
    return(adjP)
  } else{
    return(list(method= "Tarone", significant.level = alpha, Result = data.frame(raw.p = p, adjust.p=adjP, decision=ifelse(adjP<=alpha, "reject","accept"))))
  }
}

