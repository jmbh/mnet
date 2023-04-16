# jonashaslbeck@protonmail; April 5, 2023

# ------------------------------------------------------------
# -------- Function to Process mlVAR Outputs -----------------
# ------------------------------------------------------------


Process_mlVAR <- function(object1,
                          object2,
                          contemporaneous = "orthogonal",
                          temporal = "orthogonal") {

  # Number of vars:
  p <- ncol(object1$results$Gamma_Omega_mu$mean)

  # a) Between network (using function from mlVAR package)
  btw_1 <- mlVAR::getNet(object1, "between", nonsig="show")
  btw_2 <- mlVAR::getNet(object2, "between", nonsig="show")
  btw_diff <- btw_1 - btw_2

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


  outlist <- list("diff_between" = btw_diff,
                  "diff_phi_fix" = phi_fix_diff,
                  "diff_phi_RE_sd" = phi_RE_sd_diff,
                  "diff_gam_fix" = Gam_fix_diff,
                  "diff_gam_RE_sd" = Gam_RE_sd_diff)

  return(outlist)

} # eoF




