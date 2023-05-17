\name{mlVAR_GC}
\alias{mlVAR_GC}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Testing Group Differences in mlVAR
}
\description{
Does inference on group differences in mlVAR models using a permutation test. The function is a wrapper around the function \code{mlVAR()} from the \code{mlVAR} package.
}
\usage{
mlVAR_GC(data, vars, idvar, dayvar = NULL,
         beepvar = NULL, groups, test = "permutation",
         estimator, contemporaneous, temporal, nCores = 1,
         nP = 500, saveModels = FALSE, verbose = FALSE,
         pbar = TRUE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{The data.frame containing the data (rows=cases, columns=variables) including the variables to be modeled, a unique id per subject, and potentially additional variables about the measurement procedure (see \code{dayvar} and \code{dayvar} below).
}
  \item{vars}{
A character vector with the column names of the variables that should be included in the mlVAR models.
}
  \item{idvar}{
The name of the column containing the unique subject id.
}
  \item{dayvar}{
Optional: The name of the column indicating assessment day. Adding this argument makes sure that the first measurement of a day is not regressed on the last measurement of the previous day. IMPORTANT: only add this if the data has multiple observations per day. See also \code{mlVAR()}.
}
  \item{beepvar}{
Optional: The name of the column indicating measurement occasion per day. Adding this argument will cause non-consecutive beeps to be treated as missing! See also \code{mlVAR()}.
}
  \item{groups}{
The name of the column indicating group membership.
}
  \item{test}{
The type of test used. Defaults to \code{test = "permutation"}, which performs a permutation test. Alternatively, \code{test = "parametric"} performs a parametric test based on confidence intervals. For the latter, a (vectorcritical p-value has to be specified.
}
  \item{estimator}{
The estimator to be used. \code{"lmer"} for sequential univariate multi-level estimation, \code{"Mplus"} for multivariate Bayesian estimation (requires Mplus), and \code{"lm"} for fixed effects estimation. See also \code{mlVAR()}.
}
  \item{contemporaneous}{
How should the contemporaneous networks be estimated? These networks are always estimated post-hoc by investigating the residuals of the temporal models. \code{"correlated"} and \code{"orthogonal"} run second multi-level models in which the networks are estimated using node-wise estimation. \code{"fixed"} and \code{"unique"} simply correlate the residuals, either by computing one network for all subjects (fixed) or a single network per per subject. See also \code{mlVAR()}.
}
  \item{temporal}{
How should the temporal effects be estimated? \code{"correlated"} estimates correlated random effects, \code{"orthogonal"} estimates non-correlated random effects and \code{"fixed"} estimates a model in which only the intercept is random. Defaults to \code{"correlated"} when the number of variables is less than 6 and \code{"orthogonal"} otherwise. \code{"unique"} uses lm to estimate an unique model for each subject. See also \code{mlVAR()}.
}
  \item{nCores}{
Number of cores to use in computation. See also \code{mlVAR()}.
}
  \item{nP}{
Number of permutations used to create the sampling distribution under the null hypothesis for the permutation test. Defaults to \code{nP = 500}.
}
  \item{saveModels}{
If true, the pair of mlVAR objects for each permutations are saved in the output. Defaults to \code{saveModels = FALSE} since the resulting output object can become very large.
}
  \item{verbose}{
If true, the the verbose option in foreach, which is used in side \code{mlVAR_GC()}, is switched on. Defaults to \code{verbose = FALSE}.
}
  \item{pbar}{
If true, a progress bar is shown. Defaults to \code{pbar = TRUE}.
}
}
\details{
The output consists of a list with the following entries:
}


\value{
Returns a list with the following entries:

  \item{EmpDiffs}{A list containing five matrices with the empirically observed differences (Group 1 - Group 2) in the between-person network, the fixed and random effects variances of the temporal network, and the fixed and random effects variances of the contemporaneous network.}
  \item{Pval}{Same structure as \code{TrueDiffs} but provides p-values for the corresponding observed group differences. Note that for the parametric test, p-values are only available for fixed lagged effects, fixed contemporaneous (residual) effects, and the between network effects. This is because \code{mlVAR} only provides standard errors for those types of parameters.}
  \item{SampDist}{A list of arrays, containing the sampling distributions for between-person network, the fixed and random effects variances of the temporal network, and the fixed and random effects variances of the contemporaneous network.}
  \item{Models}{If \code{saveModels = TRUE}, this contains a list with all pairs of mlVAR models for each permutation.}
  \item{Runtime}{The runtime of the function in minutes.}
}


\references{
Epskamp, S., Waldorp, L. J., Mõttus, R., & Borsboom, D. (2018). The Gaussian graphical model in cross-sectional and time-series data. Multivariate behavioral research, 53(4), 453-480.
}
\author{
Jonas Haslbeck <jonashaslbeck@protonmail.com>
}

\examples{

\donttest{

# Use simulated example data loaded with package
head(ExampleData)

# Call Permutation test
out <- mlVAR_GC(data = ExampleData,
                vars = c("V1", "V2", "V3"),
                idvar = "id",
                groups = "group",
                nCores = 2, # choose max cores possible on your machine
                nP = 2) # Should be more in practice, see paper!

# P-values for the five parameter types:
out$Pval
# e.g., we see that there is a small p-value for
out$Pval$Phi_mean[1,2]
# which is a parameter with a true group difference
# (which we know since this is simulated data)

# The observed group differences (i.e., the test statistics)
# Can be found in:
out$EmpDiffs

# Specifically, the difference is: group 1 - group 2
# For example:
out$EmpDiffs$Phi_mean[1,2,]
# The true group difference was -0.4

}

}







