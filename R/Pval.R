#' Calculating p-values for discrete data
#'
#' The function for calculating the original available p-values and all attaianble p-values for the corresponding hypothesis.
#'
#' @usage
#' getPval(raw.data, test.type, alternative)
#' @param raw.data original data set with count number for treatment group and study group. The data set type could be \code{\link[base]{matrix}} or \code{\link[base]{data.frame}}.
#' @param test.type there are two discrete test available now, must be one of \code{"FET"} for Fisher's Exact Test and \code{"BET"} for Binomial Exact Test.
#' @param alternative indicates the alternative hypothesis and must be one of \code{"two.sided"}, \code{"greater"} or \code{"less"}.
#' @return
#' A numeric vector of the adjusted p-values (of the same length as \eqn{p}).
#' @author Yalin Zhu
#' @references
#' Zhu, Y., & Guo, W. (2017).
#' Familywise error rate controlling procedures for discrete data
#' \emph{arXiv preprint} arXiv:1711.08147.
#'
#' Clopper, C. J. & Pearson, E. S. (1934).
#' The use of confidence or fiducial limits illustrated in the case of the binomial.
#' \emph{Biometrika}, \strong{26}: 404-413.
#'
#' Fisher, R. A. (1922).
#' On the Interpretation of \eqn{\chi^2} from Contingency Tables, and the Calculation of P.
#' \emph{Journal of the Royal Statistical Society}, \strong{85}: 87-94.
#'
#' @examples
#'  ## Using Fisher's Exact Test to get the avaiable and attainablep-values
#'  # import raw data set as data.frame type
#'  df <-  data.frame(X1=c(4, 2, 2, 13, 6, 8, 4, 0, 1), N1 = rep(148, 9),
#'  	X2 = c(0, 0, 1, 3, 2, 1, 2, 2, 2), N2 = rep(132, 9))
#'  # obtain the avaiable p-values and attainable p-values using two-sided Fisher's Exact Test
#'  getPval(raw.data=df, test.type = "FET",alternative = "two.sided")
#'  # store the avaiable p-values
#' p <- getPval(raw.data=df, test.type = "FET",alternative = "two.sided")[[1]]; p
#'  # store the attainable p-values
#' p.set <- getPval(raw.data=df, test.type = "FET",alternative = "two.sided")[[2]]; p.set
#' @export

getPval <- function(raw.data, test.type = c("FET", "BET"), alternative = c("two.sided","greater","less")){
  if (class(raw.data)!="data.frame" & class(raw.data)!="matrix"){
    stop("The raw data type has to be 'data.frame' or 'matrix'")
  }
  else{
   test.type <- match.arg(test.type)
   alternative <- match.arg(alternative)
    p <- c(); p.set <- list()
  if (test.type=="FET"){
    for (i in 1:nrow(raw.data)){
      p[i] <- stats::fisher.test(matrix(c(raw.data[i,1],raw.data[i,2]-raw.data[i,1],raw.data[i,3],raw.data[i,4]-raw.data[i,3]),2,2))[[1]]
      s<-c()
      for (j in 0:(raw.data[i,1]+ raw.data[i,3])){
      s[j+1] <- stats::fisher.test(matrix(c(j,raw.data[i,2]-j,raw.data[i,1]+raw.data[i,3]-j,raw.data[i,4]-(raw.data[i,1]+raw.data[i,3]-j)),2,2), alternative = alternative)[[1]]
    }
    p.set[[i]] <- sort(s)
  }
 }
 if (test.type=="BET"){
   for (i in 1:nrow(raw.data)){
     p[i] <- stats::binom.test(c(raw.data[i,1],raw.data[i,2]), p = 0.5, alternative = alternative)[[3]]
     s<-c()
     for (j in 0:(raw.data[i,1]+ raw.data[i,2])){
       s[j+1] <- stats::binom.test(c(j, raw.data[i,1]+raw.data[i,2]-j), p = 0.5, alternative = alternative)[[3]]
     }
     p.set[[i]] <- sort(s)
   }
  }
  return(list(available.pvalue=p, attainable.pvalue=p.set))
 }
}
