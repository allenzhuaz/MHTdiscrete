#' The adjusted p-values for Modified Benjamini-Hochberg (BH) step-up FDR controlling procedure.
#'
#' The function for calculating the adjusted p-values based on original available p-values and all attaianble p-values.
#'
#' @usage
#' MBH.p.adjust(p, p.set, alpha, make.decision)
#' @param p numeric vector of p-values (possibly with \code{\link[base]{NA}}s). Any other R is coerced by \code{\link[base]{as.numeric}}. Same as in \code{\link[stats]{p.adjust}}.
#' @param p.set a list of numeric vectors, where each vector is the vector of all attainable p-values containing the available p-value for the corresponding hypothesis.
#' @param alpha significant level used to compare with adjusted p-values to make decisions, the default value is 0.05.
#' @param make.decision logical; if  \code{TRUE}, then the output include the decision rules compared adjusted p-values with significant level \eqn{\alpha}
#' @return
#' A numeric vector of the adjusted p-values (of the same length as \code{p}).
#' @seealso \code{\link{MBY.p.adjust}},  \code{\link{MBL.p.adjust}}
#' @author Yalin Zhu
#' @references
#' Benjamini, Y., and Hochberg, Y. (1995).
#'  Controlling the false discovery rate: a practical and powerful approach to multiple testing.
#'  \emph{Journal of the Royal Statistical Society Series B}, \strong{57}: 289-300.
#' @examples
#' p <- c(pbinom(1,8,0.5),pbinom(1,5,0.75),pbinom(1,6,0.6))
#' p.set <-list(pbinom(0:8,8,0.5),pbinom(0:5,5,0.75),pbinom(0:6,6,0.6))
#' MBH.p.adjust(p,p.set)
#' @export

MBH.p.adjust <- function(p,p.set, alpha = 0.05, make.decision = FALSE){
  o <- order(p); ro <- order(o); m <- length(p)
  sort.p <- p[o]; sort.p.set <- p.set[o]
  adjP <- numeric();pCDF <- matrix(NA,m,m)
  for(i in m:1){
    for(j in 1:m){
      pCDF[i,j] <- max(sort.p.set[[j]][sort.p.set[[j]] <= sort.p[i]],0)
    }
    c <- sum(pCDF[i,1:m])/i
    adjP[i] <- ifelse(i==m,c,min(adjP[i+1],c))
  }
  if (make.decision==FALSE){
    return(adjP[ro])
  } else{
    return(list(method= "Modified Benjamini-Hochberg (BH)", significant.level = alpha, Result = data.frame(raw.p = p, adjust.p=adjP[ro], decision=ifelse(adjP[ro]<=alpha, "reject","accept"))))
  }
}

#' The adjusted p-values for Gilbert-Tarone-BH step-up FDR controlling procedure.
#'
#' The function for calculating the adjusted p-values based on original available p-values and all attaianble p-values.
#'
#' @usage
#' GTBH.p.adjust(p, p.set, alpha, make.decision)
#' @param p numeric vector of p-values (possibly with \code{\link[base]{NA}}s). Any other R is coerced by \code{\link[base]{as.numeric}}. Same as in \code{\link[stats]{p.adjust}}.
#' @param p.set a list of numeric vectors, where each vector is the vector of all attainable p-values containing the available p-value for the corresponding hypothesis.
#' @param alpha significant level used to compare with adjusted p-values to make decisions, the default value is 0.05.
#' @param make.decision logical; if  \code{TRUE}, then the output include the decision rules compared adjusted p-values with significant level \eqn{\alpha}
#' @return
#' A numeric vector of the adjusted p-values (of the same length as \code{p}).
#' @seealso \code{\link{GTBY.p.adjust}},  \code{\link{MBH.p.adjust}},  \code{\link{MBY.p.adjust}}
#' @author Yalin Zhu
#' @references
#' Gilbert, P. B. (2005).
#' A modified false discovery rate multiple-comparisons procedure for discrete data, applied to human immunodeficiency virus genetics.
#' \emph{Journal of the Royal Statistical Society: Series C (Applied Statistics)}, \strong{54}: 143-158.
#'
#' Benjamini, Y., and Hochberg, Y. (1995).
#'  Controlling the false discovery rate: a practical and powerful approach to multiple testing.
#'  \emph{Journal of the Royal Statistical Society Series B}, \strong{57}: 289-300.
#'
#' @examples
#' p <- c(pbinom(1,8,0.5),pbinom(1,5,0.75),pbinom(1,6,0.6))
#' p.set <-list(pbinom(0:8,8,0.5),pbinom(0:5,5,0.75),pbinom(0:6,6,0.6))
#' GTBH.p.adjust(p,p.set)
#' @export

GTBH.p.adjust <- function(p,p.set, alpha = 0.05, make.decision = FALSE){
  o <- order(p); ro <- order(o); m <- length(p); adjP <- numeric(m)
  sort.p <- p[o]; sort.p.set <- p.set[o]
  minP <- sort(sapply(sort.p.set,min))
  for (j in m:1){
    for (i in 1:m ){
      if(sort.p[j]>=max(minP)){q <- m}
      else if (sort.p[j]>=minP[i] & sort.p[j]<minP[i+1]){q <- i}
    }
    adjP[j] <-  ifelse(j==q,sort.p[j],min(adjP[j+1],q*sort.p[j]/j))
  }
  if (make.decision==FALSE){
    return(adjP[ro])
  } else{
    return(list(method= "Gilbert-Tarone-BH", significant.level = alpha, Result = data.frame(raw.p = p, adjust.p=adjP[ro], decision=ifelse(adjP[ro]<=alpha, "reject","accept"))))
  }
}



#' The adjusted p-values for Modified Benjamini-Yekutieli (BY) step-up FDR controlling procedure.
#'
#' The function for calculating the adjusted p-values based on original available p-values and all attaianble p-values.
#'
#' @usage
#' MBY.p.adjust(p, p.set, alpha, make.decision)
#' @param p numeric vector of p-values (possibly with \code{\link[base]{NA}}s). Any other R is coerced by \code{\link[base]{as.numeric}}. Same as in \code{\link[stats]{p.adjust}}.
#' @param p.set a list of numeric vectors, where each vector is the vector of all attainable p-values containing the available p-value for the corresponding hypothesis.
#' @param alpha significant level used to compare with adjusted p-values to make decisions, the default value is 0.05.
#' @param make.decision logical; if  \code{TRUE}, then the output include the decision rules compared adjusted p-values with significant level \eqn{\alpha}
#' @return
#' A numeric vector of the adjusted p-values (of the same length as \code{p}).
#' @seealso  \code{\link{MBH.p.adjust}},  \code{\link{MBL.p.adjust}}
#' @author Yalin Zhu
#' @references
#' Benjamini, Y., and Yekutieli, D. (2001).
#' The control of the false discovery rate in multiple testing under dependency.
#' \emph{Annals of Statistics}, \strong{29}: 1165-1188.
#' @examples
#' p <- c(pbinom(1,8,0.5),pbinom(1,5,0.75),pbinom(1,6,0.6))
#' p.set <-list(pbinom(0:8,8,0.5),pbinom(0:5,5,0.75),pbinom(0:6,6,0.6))
#' MBY.p.adjust(p,p.set)
#' @export


MBY.p.adjust <- function(p,p.set, alpha = 0.05, make.decision = FALSE){
  o <- order(p); ro <- order(o); m <- length(p); C <- sum(1/c(1:m))
  sort.p <- p[o]; sort.p.set <- p.set[o]
  adjP <- numeric();pCDF <- matrix(NA,m,m)
  for(i in m:1){
    for(j in 1:m){
      pCDF[i,j] <- max(sort.p.set[[j]][sort.p.set[[j]] <= sort.p[i]],0)
    }
    c <- min(1,sum(pCDF[i,1:m])*C/i)
    adjP[i] <- ifelse(i==m,c,min(adjP[i+1],c))
  }
  if (make.decision==FALSE){
    return(adjP[ro])
  } else{
    return(list(method= "Modified Benjamini-Yekutieli (BY)", significant.level = alpha, Result = data.frame(raw.p = p, adjust.p=adjP[ro], decision=ifelse(adjP[ro]<=alpha, "reject","accept"))))
  }
}


#' The adjusted p-values for Gilbert-Tarone-BY step-up FDR controlling procedure.
#'
#' The function for calculating the adjusted p-values based on original available p-values and all attaianble p-values.
#'
#' @usage
#' GTBY.p.adjust(p, p.set, alpha, make.decision)
#' @param p numeric vector of p-values (possibly with \code{\link[base]{NA}}s). Any other R is coerced by \code{\link[base]{as.numeric}}. Same as in \code{\link[stats]{p.adjust}}.
#' @param p.set a list of numeric vectors, where each vector is the vector of all attainable p-values containing the available p-value for the corresponding hypothesis.
#' @param alpha significant level used to compare with adjusted p-values to make decisions, the default value is 0.05.
#' @param make.decision logical; if  \code{TRUE}, then the output include the decision rules compared adjusted p-values with significant level \eqn{\alpha}
#' @return
#' A numeric vector of the adjusted p-values (of the same length as \code{p}).
#' @seealso \code{\link{GTBH.p.adjust}},  \code{\link{MBH.p.adjust}},  \code{\link{MBY.p.adjust}}
#' @author Yalin Zhu
#' @references
#' Gilbert, P. B. (2005).
#' A modified false discovery rate multiple-comparisons procedure for discrete data, applied to human immunodeficiency virus genetics.
#' \emph{Journal of the Royal Statistical Society: Series C (Applied Statistics)}, \strong{54}: 143-158.
#'
#' Benjamini, Y., and Yekutieli, D. (2001).
#' The control of the false discovery rate in multiple testing under dependency.
#' \emph{Annals of Statistics}, \strong{29}: 1165-1188.
#'
#' @examples
#' p <- c(pbinom(1,8,0.5),pbinom(1,5,0.75),pbinom(1,6,0.6))
#' p.set <-list(pbinom(0:8,8,0.5),pbinom(0:5,5,0.75),pbinom(0:6,6,0.6))
#' GTBY.p.adjust(p,p.set)
#' @export

GTBY.p.adjust <-  function(p,p.set, alpha = 0.05, make.decision = FALSE){
  o <- order(p); ro <- order(o); m <- length(p); adjP <- numeric(m)
  sort.p <- p[o]; sort.p.set <- p.set[o]
  minP <- sort(sapply(sort.p.set,min))
  for (j in m:1){
    for (i in 1:m ){
      if(sort.p[j]>=max(minP)){q <- m}
      else if (sort.p[j]>=minP[i] & sort.p[j]<minP[i+1]){q <- i}
    }
    C <- sum(1/c(1:q))
    adjP[j] <-  ifelse(j==q,sort.p[j],min(adjP[j+1],q*C*sort.p[j]/j))
  }
  if (make.decision==FALSE){
    return(adjP[ro])
  } else{
    return(list(method= "Gilbert-Tarone-BY", significant.level = alpha, Result = data.frame(raw.p = p, adjust.p=adjP[ro], decision=ifelse(adjP[ro]<=alpha, "reject","accept"))))
  }
}
