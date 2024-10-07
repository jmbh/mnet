
VAR_GC <- function(data,
                   vars,
                   dayvar,
                   beepvar,
                   groups,
                   boot=FALSE,  # TO ADD
                   nBoot=1000) {

  # --- Input Checks ---
  # TODO: Add mode, match at least mlVAR_GC()

  # (3) Is the grouping variable specified properly?
  v_groups <- as.numeric(unlist(data[, groups]))
  if(any(!(v_groups %in% 1:2))) stop("Groups need to be specified by a vector of 1s and 2s referring to the two groups.")

  # --- Get Info ---
  data1 <- data[data[, groups]==1, ]
  data2 <- data[data[, groups]==2, ]

  n1 <- nrow(data1)
  n2 <- nrow(data2)

  data1x <- as.matrix(data1[, vars])
  data2x <- as.matrix(data2[, vars])

  p <- ncol(data1x)

  # --- Compute which rows are usable for each dataset ---
  v_pdb1 <- f_pdb(data1[, dayvar], data1[, beepvar])
  v_pdb2 <- f_pdb(data2[, dayvar], data2[, beepvar])

  # --- Estimate VAR model ---
  phi1 <- phi2 <- array(NA, dim=c(p, p, 2)) # estimates + SE
  int1 <- int2 <- matrix(NA, p, 2) # estimates + SE

  for(j in 1:p) {
    # Fit Model 1
    mod_j1 <- lm(data1x[which(v_pdb1), j] ~ data1x[which(v_pdb1)-1, ])
    coefs1 <- coef(mod_j1)
    SEs1 <- coef(summary(mod_j1))[, "Std. Error"]
    phi1[j, , 1] <- coefs1[-1]
    phi1[j, , 2] <- SEs1[-1]
    int1[j, 1] <- coefs1[1]
    int1[j, 2] <- SEs1[1]
    # Fit Model 2
    mod_j2 <- lm(data2x[which(v_pdb1), j] ~ data2x[which(v_pdb1)-1, ])
    coefs2 <- coef(mod_j2)
    SEs2 <- coef(summary(mod_j2))[, "Std. Error"]
    phi2[j, , 1] <- coefs2[-1]
    phi2[j, , 2] <- SEs2[-1]
    int2[j, 1] <- coefs2[1]
    int2[j, 2] <- SEs2[1]
  } # end for: j

  # --- Get p-value for all comparisons ---
  # Phi matrix
  phi_diff <- (phi1[, , 1] - phi2[, , 1]) / sqrt(phi1[, , 2]^2 + phi2[, , 2]^2)
  df <- (phi1[, , 2]^2 + phi2[, , 1]^2)^2 / ((phi1[, , 2]^4 / (n1 - 1)) + (phi2[, , 1]^4 / (n2 - 1))) # Degrees of freedom (approximation using Welch-Satterthwaite equation)
  phi_pval <- 2* pt(-abs(phi_diff), df=df)

  # Intercepts
  int_diff <- (int1[, 1] - int2[, 1]) / sqrt(int1[, 2]^2 + int2[, 2]^2)
  df <- (int1[, 2]^2 + int2[, 1]^2)^2 / ((int1[, 2]^4 / (n1 - 1)) + (int1[,1]^4 / (n2 - 1))) # Degrees of freedom (approximation using Welch-Satterthwaite equation)
  int_pval <- 2* pt(-abs(int_diff), df=df)

  # --- Bootstrap Sampling Distribution ---
  if(boot) {

  } # end if

  # --- Compile Output ---
  outlist <- list("phi_diff" = phi_diff,
                  "phi_pval" = phi_pval,
                  "int_diff" = int_diff,
                  "int_pval" = int_pval)
  return(outlist)

} # eoF
