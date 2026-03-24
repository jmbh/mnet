# jonashaslbeck@protonmail; March 24, 2026

# ------------------------------------------------------------
# -------- Function to Process mlVAR Outputs -----------------
# ------------------------------------------------------------

Process_mlVAR <- function(object1,
                          object2,
                          contemporaneous = "orthogonal",
                          temporal = "orthogonal",
                          empirical = TRUE) {

  # Number of vars:
  p <- ncol(object1$results$Gamma_Omega_mu$mean)

  # a) Between network (using function from mlVAR package)

  # For very low number of subjects close to the boundary of identifiability
  # it is possible that lme4 in mlVAR estimates zero variances for random intercepts
  # which then does not allow one to estimate the between person network
  # for these cases we set the differences to zero here; later, when calculating
  # p-values we will just exclude those cases

  # # Here we use partial correlations [OLD]:
  # check1 <- any(is.na(object1$results$Omega_mu$pcor$mean))
  # if(check1) {
  #   btw_1 <- matrix(NA, p, p)
  # } else {
  #   btw_1 <- mlVAR::getNet(object1, "between", nonsig="show")
  # }
  # check2 <- any(is.na(object2$results$Omega_mu$pcor$mean))
  # if(check2) {
  #   btw_2 <- matrix(NA, p, p)
  # } else {
  #   btw_2 <- mlVAR::getNet(object2, "between", nonsig="show")
  # }
  # btw_diff <- btw_1 - btw_2
  # if(empirical & any(c(check1, check2))) warning("Random intercept variance was estimated to be zero for some variablesin mlVAR(). Therefore, no between-network can be obtained.")

  # We use averaged regression coefficients, to use the same parameter as in the parametric test
  # Note: also computing the partial correlations requires averaging regression coefficients
  check1 <- any(is.na(object1$results$Gamma_Omega_mu$mean))
  if(check1) {
    btw_1 <- matrix(NA, p, p)
  } else {
    btw_1 <- object1$results$Gamma_Omega_mu$mean
    btw_1 <- (btw_1 + t(btw_1)) / 2
  }
  check2 <- any(is.na(object2$results$Gamma_Omega_mu$mean))
  if(check2) {
    btw_2 <- matrix(NA, p, p)
  } else {
    btw_2 <- object2$results$Gamma_Omega_mu$mean
    btw_2 <- (btw_2 + t(btw_2)) / 2
  }
  btw_diff <- btw_1 - btw_2
  if(empirical & any(c(check1, check2))) warning("Random intercept variance was estimated to be zero for some variablesin mlVAR(). Therefore, no between-network can be obtained.")

  # b.1) VAR: fixed effects
  phi_fix_1 <- object1$results$Beta$mean
  phi_fix_2 <- object2$results$Beta$mean
  phi_fix_diff <- phi_fix_1 - phi_fix_2

  # b.2) VAR: RE sds
  phi_RE_sd_1 <- object1$results$Beta$SD
  phi_RE_sd_2 <- object2$results$Beta$SD
  phi_RE_sd_diff <- phi_RE_sd_1 - phi_RE_sd_2

  # c.1) Contemp: fixed effects
  Gam_fix_1 <- object1$results$Gamma_Theta$mean
  Gam_fix_1 <- (Gam_fix_1 + t(Gam_fix_1)) / 2 # Apply AND-rule
  Gam_fix_2 <- object2$results$Gamma_Theta$mean
  Gam_fix_2 <- (Gam_fix_2 + t(Gam_fix_2)) / 2 # Apply AND-rule
  Gam_fix_diff <- Gam_fix_1 - Gam_fix_2

  # c.2) Contemp: RE sds
  Gam_RE_sd_1 <- object1$results$Gamma_Theta$SD
  Gam_RE_sd_1 <- (Gam_RE_sd_1 + t(Gam_RE_sd_1)) / 2 # Apply AND-rule
  Gam_RE_sd_2 <- object2$results$Gamma_Theta$SD
  Gam_RE_sd_2 <- (Gam_RE_sd_2 + t(Gam_RE_sd_2)) / 2 # Apply AND-rule
  Gam_RE_sd_diff <- Gam_RE_sd_1 - Gam_RE_sd_2


  ## Prepare output list
  # Now both the differences and the raw parameters
  outlist <- list("between" = list(btw_1, btw_2),
                  "phi_fix" = list(phi_fix_1, phi_fix_2),
                  "phi_RE_sd" = list(phi_RE_sd_1, phi_RE_sd_2),
                  "gam_fix" = list(Gam_fix_1, Gam_fix_2),
                  "gam_RE_sd" = list(Gam_RE_sd_1, Gam_RE_sd_2))

  return(outlist)

} # eoF




