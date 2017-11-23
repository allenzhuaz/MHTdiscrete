#' The adjusted p-values for Modified Holm step-down FWER controlling procedure.
#'
#' The function for calculating the adjusted p-values based on original available p-values and all attaianble p-values.
#'
#' @usage
#' MHolm.p.adjust(p, p.set, alpha, make.decision)
#' @param p numeric vector of p-values (possibly with \code{\link[base]{NA}}s). Any other R is coerced by \code{\link[base]{as.numeric}}. Same as in \code{\link[stats]{p.adjust}}.
#' @param p.set a list of numeric vectors, where each vector is the vector of all attainable p-values containing the available p-value for the corresponding hypothesis.
#' @param alpha significant level used to compare with adjusted p-values to make decisions, the default value is 0.05.
#' @param make.decision logical; if  \code{TRUE}, then the output include the decision rules compared adjusted p-values with significant level \eqn{\alpha}
#' @return
#' A numeric vector of the adjusted p-values (of the same length as \code{p}).
#' @seealso \code{\link{TH.p.adjust}}, \code{\link[stats]{p.adjust}}.
#' @author Yalin Zhu
#' @references
#' Zhu, Y., & Guo, W. (2017).
#' Familywise error rate controlling procedures for discrete data
#' \emph{arXiv preprint} arXiv:1711.08147.
#'
#' Holm, S. (1979).
#' A simple sequentially rejective multiple test procedure.
#' \emph{Scandinavian Journal of Statistics}, \strong{6}: 65-70.
#' @examples
#' p <- c(pbinom(1,8,0.5),pbinom(1,5,0.75),pbinom(1,6,0.6))
#' p.set <-list(pbinom(0:8,8,0.5),pbinom(0:5,5,0.75),pbinom(0:6,6,0.6))
#' MHolm.p.adjust(p,p.set)
#' ## Compare with the traditional Holm adjustment
#' p.adjust(p,method = "holm")
#' ## Compare with the Tarone-Holm adjustment
#' TH.p.adjust(p,p.set)
#' @export

MHolm.p.adjust <- function(p,p.set, alpha = 0.05, make.decision = FALSE){
  o <- order(p); ro <- order(o); m <- length(p)
  sort.p <- p[o]; sort.p.set <- p.set[o]
  adjP <- numeric(m); pCDF <- matrix(NA,m,m)
  for(i in 1:m){
    for (j in i:m){
      pCDF[i,j] <- max(sort.p.set[[j]][sort.p.set[[j]] <= sort.p[i]],0)
    }
    c <- min(1,sum(pCDF[i,i:m]))
    adjP[i] <- ifelse(i==1,c,max(adjP[i-1],c))
  }
  if (make.decision==FALSE){
    return(adjP[ro])
  } else{
    return(list(method= "Modified Holm", significant.level = alpha, Result = data.frame(raw.p = p, adjust.p=adjP[ro], decision=ifelse(adjP[ro]<=alpha, "reject","accept"))))
  }
}

#' The adjusted p-values for Tarone-Holm step-down FWER controlling procedure.
#'
#' The function for calculating the adjusted p-values based on original available p-values and all attaianble p-values.
#'
#' @usage
#' TH.p.adjust(p, p.set, alpha, make.decision)
#' @param p numeric vector of p-values (possibly with \code{\link[base]{NA}}s). Any other R is coerced by \code{\link[base]{as.numeric}}. Same as in \code{\link[stats]{p.adjust}}.
#' @param p.set a list of numeric vectors, where each vector is the vector of all attainable p-values containing the available p-value for the corresponding hypothesis.
#' @param alpha significant level used to compare with adjusted p-values to make decisions, the default value is 0.05.
#' @param make.decision logical; if  \code{TRUE}, then the output include the decision rules compared adjusted p-values with significant level \eqn{\alpha}
#' @return
#' A numeric vector of the adjusted p-values (of the same length as \code{p}).
#' @seealso \code{\link{MHolm.p.adjust}}, \code{\link[stats]{p.adjust}}.
#' @author Yalin Zhu
#' @references
#' Hommel, G., & Krummenauer, F. (1998).
#' Improvements and modifications of Tarone's multiple test procedure for discrete data.
#' \emph{Biometrics}, \strong{54}: 673-681.
#'
#' Holm, S. (1979).
#' A simple sequentially rejective multiple test procedure.
#' \emph{Scandinavian Journal of Statistics}, \strong{6}: 65-70.
#'
#' @examples
#' p <- c(pbinom(1,8,0.5),pbinom(1,5,0.75),pbinom(1,6,0.6))
#' p.set <-list(pbinom(0:8,8,0.5),pbinom(0:5,5,0.75),pbinom(0:6,6,0.6))
#' TH.p.adjust(p,p.set)
#' @export


TH.p.adjust <- function(p,p.set, alpha = 0.05, make.decision = FALSE){
  o <- order(p); ro <- order(o); m <- length(p)
  sort.p <- p[o]; sort.p.set <- p.set[o]; adjP <- numeric(m)
  j <- 1
  while (j <= m){
    minP <- sort(sapply(sort.p.set[j:m],min))
    for (i in 1:(m-j+1)){
      if (sort.p[j]>=max(minP)){q <- m-j+1}
      else if (sort.p[j]>=minP[i] & sort.p[j]<minP[i+1]){q <- i}
    }
    c <- min(1,q*sort.p[j])
    adjP[j] <- ifelse(j==1,c,max(adjP[j-1],c))
    j <- j+1
  }
  if (make.decision==FALSE){
    return(adjP[ro])
  } else{
    return(list(method= "Tarone-Holm", significant.level = alpha, Result = data.frame(raw.p = p, adjust.p=adjP[ro], decision=ifelse(adjP[ro]<=alpha, "reject","accept"))))
  }
}



