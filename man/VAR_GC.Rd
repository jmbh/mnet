\name{VAR_GC}
\alias{VAR_GC}
\title{
Function to test whether parameters of two VAR(1) models are different
}
\description{
Function to test whether parameters of two VAR(1) models are different
}
\usage{
VAR_GC(data, vars, dayvar, beepvar,
       groups, test = "parametric",
       nP = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A n x p data matrix.
}
  \item{vars}{
An integeger vector indicating the column numbers of the variables that should be modeled in the pair of VAR models.
}
  \item{dayvar}{
Integener indicating the column number of the variable specifying the day of the measurement point. This is designed for EMA studies. If not applicable, this argument can be left unspecified.
}
  \item{beepvar}{
Integener indicating the column number of the variable specifying the number of the measurement occasion during a specific day. The full set of measurement occasions has to be consecutive and increasing sequence of integers (e.g., 1, 2,..., 5). This is designed for EMA studies. If not applicable, this argument can be left unspecified.
}
  \item{groups}{
Integer indiciating the column number of the group variable. The values of the group variable have to be 1 and 2.
}
  \item{test}{
If \code{test = "parmetric"} p-values are computed using the t-distribution, using the Welch-Satterthwaite equation to approximate degrees of freedom (df). This is the default. Alternatively, \code{test = "permutation"} uses a permutation to sample from a distribution in which the null hypothesis that no group differences exist is true. The permutation uses \code{nP} permuted datasets.
}
  \item{nP}{
The number of permuted datasets used if  \code{test = "permutation"}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}

\value{

\item{Call}{
Returns the call of the function
}

\item{phi_diff}{
A p x p matrix of differences in phi coefficients (Group 1 - Group 2) in the empirical data. The test-statistics.
}

\item{phi_pval}{
A p x p matrix with pvalues corresponding to \code{phi_diff}.
}

\item{int_diff}{
A p numeric vector of differences in intercepts (Group 1 - Group 2).
}

\item{int_pvar}{
A p numeric vector of pvalues corresponding to \code{int_diff}
}

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Jonas Haslbeck <jonashaslbeck@protonmail.com>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

\donttest{

library(mlVAR) # for simulateVAR() function

# Specify Model
p <- 4
A1 <- diag(p) * 0.8
A2 <- diag(p) * 0.8
A2[2,1] <- 0.7

# Simulate datasets
Nt <- 500
set.seed(13)  # for reproducibility
data1_x <- simulateVAR(A1, means=rep(0, p), Nt = Nt, residuals=.1)
data2_x <- simulateVAR(A2, means=rep(0, p), Nt = Nt, residuals=.1)

# Add beep and day vars
dayvar1 <- dayvar2 <- rep(1:(Nt/5), each=5)
beepvar1 <- beepvar2 <- rep(1:5, Nt/5)

# Add grouping var
groups1 <- rep(1, Nt)
groups2 <- rep(2, Nt)

# Combine
data1 <- data.frame(cbind(dayvar1, beepvar1, groups1, data1_x))
data2 <- data.frame(cbind(dayvar2, beepvar2, groups2, data2_x))
colnames(data1) <- colnames(data2) <- c("dayvar", "beepvar", "groups", paste0("V", 1:4))
data <- rbind(data1, data2)


# Call
out <- VAR_GC(data = data,
              vars = 4:7,
              dayvar = 1,
              beepvar = 2,
              groups = 3)

round(out$phi_pval, 2)
round(out$phi_pval[2,1], 2) # worked!

}


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
