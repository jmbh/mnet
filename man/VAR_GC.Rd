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
       groups, boot = FALSE,
       nBoot = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{vars}{
%%     ~~Describe \code{vars} here~~
}
  \item{dayvar}{
%%     ~~Describe \code{dayvar} here~~
}
  \item{beepvar}{
%%     ~~Describe \code{beepvar} here~~
}
  \item{groups}{
%%     ~~Describe \code{groups} here~~
}
  \item{boot}{
NOT IMPLEMENTED YET: If \code{boot = TRUE} p-values are computed using bootstrapping instead of normal theory.
}
  \item{nBoot}{
%%     ~~Describe \code{nBoot} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
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
