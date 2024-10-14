
VAR_GC <- function(data,
                   vars,
                   dayvar,
                   beepvar,
                   groups,
                   test="parametric",  # TO ADD
                   nP=100) {

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

  # TODO: Allow for possibility that either beepvar/dayvar or both are not specified

  # --- Estimate VAR model ---
  # This is needed as the test-statistic for both the parametric and the bootstrap test

  phi1 <- phi2 <- array(NA, dim=c(p, p, 2)) # estimates + SE
  int1 <- int2 <- matrix(NA, p, 2) # estimates + SE

  # Loop variables
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

  # Get differences
  phi_diff_emp <- phi1[, , 1] - phi2[, , 1]
  int_diff_emp <- int1[, 1] - int2[, 1]

  # ---------- Parametric Test ----------
  if(test == "parametric") {
    # --- Get p-value for all comparisons ---
    # Phi matrix
    phi_tstat <- (phi1[, , 1] - phi2[, , 1]) / sqrt(phi1[, , 2]^2 + phi2[, , 2]^2)
    df <- (phi1[, , 2]^2 + phi2[, , 1]^2)^2 / ((phi1[, , 2]^4 / (n1 - 1)) + (phi2[, , 1]^4 / (n2 - 1))) # Degrees of freedom (approximation using Welch-Satterthwaite equation)
    phi_pval <- 2* pt(-abs(phi_tstat), df=df)

    # Intercepts
    int_tstat <- (int1[, 1] - int2[, 1]) / sqrt(int1[, 2]^2 + int2[, 2]^2)
    df <- (int1[, 2]^2 + int2[, 1]^2)^2 / ((int1[, 2]^4 / (n1 - 1)) + (int1[,1]^4 / (n2 - 1))) # Degrees of freedom (approximation using Welch-Satterthwaite equation)
    int_pval <- 2* pt(-abs(int_tstat), df=df)

  } # end if: parametric

  # --- Permutation Test ---
  if(test == "permutation") {

    # Get number of predictable timepoints per dataset
    n1_pd <- sum(v_pdb1)
    n2_pd <- sum(v_pdb2)
    n_pd <- n1_pd + n2_pd

    # Combine datasets
    v_pdb_cmb <- c(v_pdb1, v_pdb2)
    data_cmb_x <- rbind(data1x, data2x)

    # Storage
    phi_diff <- array(NA, dim=c(p, p, nP))
    int_diff <- matrix(NA, nP, p)

    for(b in 1:nP) {

      # Loop variables
      for(j in 1:p) {

        # We sample from a null distribution with E[no differences]

        ## Fit Model 1
        perm_ind_1 <- which(v_pdb_cmb)[sample(1:n_pd, size=n1_pd, replace=TRUE)]
        mod_j1 <- lm(data_cmb_x[perm_ind_1, j] ~ data_cmb_x[perm_ind_1-1, ])
        coefs1 <- coef(mod_j1)
        # Fit Model 2
        perm_ind_2 <- which(v_pdb_cmb)[sample(1:n_pd, size=n2_pd, replace=TRUE)]
        mod_j2 <- lm(data_cmb_x[perm_ind_2, j] ~ data_cmb_x[perm_ind_2-1, ])
        coefs2 <- coef(mod_j2)
        # Save Differences
        phi_diff[j, , b] <- coefs1[-1] - coefs2[-1]
        int_diff[b, j] <- coefs1[1] - coefs2[1]

      } # end: variables

    } # end: boot

    # ----- Compute p-values ------

    ## Phi
    # Create new array to compute p-values efficiently
    phi_diff_calc <- array(NA, dim=c(p, p, nP, 2))
    phi_diff_calc[, , , 1] <- phi_diff
    phi_diff_calc[, , 1:nP, 2] <- phi_diff_emp
    # compute p-vals
    phi_pval <- apply(phi_diff_calc, 1:2, function(x) {
      sum(abs(x[, 1]) > abs(x[, 2]))
    } )/nP

    ## Intercepts
    # Create array to comp p-vals
    int_diff_calc <- array(NA, dim=c(nP, p, 2))
    int_diff_calc[, , 1] <- int_diff
    int_diff_calc[, , 2] <- matrix(int_diff_emp, nrow=nP, ncol=4, byrow = TRUE)
    # Compute pvals
    int_pval <- apply(int_diff_calc, 2, function(x) {
      sum(abs(x[, 1]) > abs(x[, 2]))
    } )/nP

  } # end if

  # --- Compile Output ---
  outlist <- list("Call"=list("data" = data,
                              "vars" = vars,
                              "dayvar" = dayvar,
                              "beepvar" = beepvar,
                              "groups" = groups,
                              "test" = test,
                              "nP" = nP),
                  "phi_diff" = phi_diff_emp,
                  "phi_pval" = phi_pval,
                  "int_diff" = int_diff_emp,
                  "int_pval" = int_pval)
  return(outlist)

} # eoF











