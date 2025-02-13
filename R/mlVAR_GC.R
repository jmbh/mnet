# jonashaslbeck@protonmail; May 17th, 2024

# ------------------------------------------------------------
# -------- Function for Permutation Test ---------------------
# ------------------------------------------------------------

# Inputs:
# - two datasets with same structure + same inputs as for mlVAR()
# - number of samples in permutation test

# Output:
# - sampling distribution for each parameter
# - test statistic based on the mlVAR models estimated on the two datasets

mlVAR_GC <- function(data, # data including both groups
                     vars, # variables to be included in mlVAR (same in both data sets)
                     idvar, # variable indicating the nesting/subject id (same in both data sets)
                     dayvar = NULL,
                     beepvar = NULL,
                     groups, # indicates which case belongs to which group
                     test = "permutation", # can also be "parametric"
                     partest = "Welch",
                     paired = FALSE,
                     estimator, # same as in ml
                     contemporaneous, # same as in ml
                     temporal, # same as in ml
                     scale, # is data scaled?
                     # trueMeans, # allow specifying true means for developmental purposes
                     nCores = 1,
                     nP = 500, # number of samples in permutation test
                     saveModels = FALSE, # if TRUE, all models are saved; defaults to FALSE to save memory
                     saveEmpModels = FALSE,
                     verbose = FALSE, # if TRUE, verbose mode is activated in foreach
                     pbar = TRUE # if TRUE, a progress bar is being shown
) {

  timer_total <- proc.time()[3]

  # ------ Input Checks -----

  # (1) Are the data numerical?
  v_class <- apply(data[, vars], 2, function(x) class(x))
  if(any( !(v_class %in% c("numeric", "integer")) )) stop("Modeled variables need to be provided as integer or numeric variables.")

  # (2) Can we find the all variables?
  check_coln1 <- c(vars, idvar) %in% colnames(data)
  if(any(!check_coln1)) stop("Specified variable names could not be found in the data.")

  # (3) Is the grouping variable specified properly?
  v_groups <- as.numeric(unlist(data[, groups]))
  if(any(!(v_groups %in% 1:2))) stop("Groups need to be specified by a vector of 1s and 2s referring to the two groups.")

  # (4) Are valid tests selected?
  if(!(test %in% c("permutation", "parametric"))) stop('The available tests are "permutation" and "parametric".')

  # (5) No paired parametric test
  if(test == "parametric" & paired == TRUE) stop("No parametric test available for paired samples. Use the paired permutation test instead.")

  # (6) Specification of "beepvar" and "dayvar"
  if(!is.null(beepvar) & is.null(dayvar)) stop("If you use only one time vector, use 'dayvar' and not 'beepvar')")

  data1 <- data[v_groups==1, ]
  data2 <- data[v_groups==2, ]

  # Get IDs
  ids1 <- sapply(data1[, idvar], as.character)
  ids2 <- sapply(data2[, idvar], as.character)
  v_ids <- c(ids1, ids2)
  v_u_ids <- unique(v_ids)
  u_ids1 <- unique(ids1)
  u_ids2 <- unique(ids2)
  n_subj_G1 <- length(u_ids1)
  n_subj_G2 <- length(u_ids2)
  n_subj <- n_subj_G1 + n_subj_G2

  # (6) Independent samples test: Are IDs unique across datasets? [required for indepdendent samples]
  if(paired == FALSE) {
    v_intersec <- intersect(u_ids1, u_ids2)
    if(length(v_intersec) > 0) stop("IDs need to be unique across two datasets.")
  }

  # (7) Paired samples test: do we have measurements in both "groups" form each person?
  if(paired == TRUE) {
    ind_same_IDs <- all(sort(u_ids1) == sort(u_ids2))
    if(!ind_same_IDs) stop("Data for each ID is required in both 'groups'.")
  }

  # ------ Collect passed down arguments -----
  if(missing(estimator)) estimator <- "default"
  if(missing(contemporaneous)) contemporaneous <- "orthogonal"
  if(missing(temporal)) temporal <- "orthogonal"
  if(missing(nCores)) nCores <- 1
  if(missing(scale)) scale <- TRUE
  # if(missing(trueMeans)) trueMeans <- NULL

  # Copy call
  Call <- list("vars" = vars,
               "idvar" = idvar,
               "dayvar" = dayvar,
               "beepvar" = beepvar,
               "groups" = groups,
               "test" = test,
               "estimator" = estimator,
               "contemporaneous" = contemporaneous,
               "temporal" = temporal,
               "scale" = scale,
               "nCores" = nCores,
               "nP" = nP,
               "saveModels" = saveModels,
               "verbose" = verbose,
               "pbar" = pbar)


  # ------ Get Basic Info -----

  # number of variables
  p <- length(vars)

  # Combine in list/single matrix for shorter code below
  l_data <- list(data1, data2)
  m_data_cmb <- rbind(data1, data2)

  # Get Number of subjects
  v_Ns <- c(length(u_ids1), length(u_ids2))
  totalN <- sum(v_Ns)

  if(pbar) pb <- txtProgressBar(0, nP+1, style = 3) else pb <- NULL


  # ------ Loop Over Permutations -----

  if(test == "permutation") {

    # Storage
    if(saveModels) l_out_mods <- list(vector("list", length = nP),
                                      vector("list", length = nP))

    # Setup cores, if multi-core
    if(nCores > 1) {
      cl <- makeCluster(nCores, outfile = "")
      registerDoParallel(cl)
    }

    # Define a function that suppresses package startup messages
    quiet_library <- function(pkg) {
      suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    }

    # Call foreach:
    out_P <- foreach(b = 1:nP,
                     .packages = c("mlVAR"),
                     .export = c("m_data_cmb", "vars", "idvar", "estimator",
                                 "contemporaneous", "temporal", "totalN", "v_Ns",
                                 "v_ids", "pb", "pbar", "dayvar", "beepvar", "paired",
                                 "quiet_library"),
                     .verbose = verbose) %dopar% {

                       # Suppress starting message of mnet
                       quiet_library("mnet")

                       # --- Make permutation ---

                       # For independent samples
                       if(paired == FALSE) {
                         # Here we randomly distribute persons across the group, keeping with the original group sizes

                         # This is done in a way that keeps the size in each group exactly the same as in the real groups
                         v_ids_rnd <- v_u_ids[sample(1:totalN, size=totalN, replace=FALSE)]
                         v_ids_1 <- v_ids_rnd[1:v_Ns[1]]
                         v_ids_2 <- v_ids_rnd[(v_Ns[1]+1):totalN]

                         # Split data based on permutations
                         data_h0_1 <- m_data_cmb[v_ids %in% v_ids_1, ]
                         data_h0_2 <- m_data_cmb[v_ids %in% v_ids_2, ]
                         l_data_h0 <- list(data_h0_1, data_h0_2)
                       } # end if

                       # For dependent samples
                       if(paired == TRUE) {
                         # here we go through each participant, and switch the group label (or not), determined by a coin flip

                         unique_IDs <- unique(data[, idvar])
                         n_unique_IDs <- length(unique_IDs)

                         # Create copies of dataset
                         data_perm <- data

                         for(i in 1:n_unique_IDs) {
                           flip <- sample(c(FALSE, TRUE), size=1, prob = c(0.5, 0.5))
                           if(flip) data_perm[data_perm[, idvar] == unique_IDs[i], groups] <- abs(data_perm[data_perm[, idvar] == unique_IDs[i], groups] -3)
                         } # end: permuting

                         # Collect two permuted datasets in list
                         l_data_h0 <- list(data_perm[data_perm[, groups]==1, ],
                                           data_perm[data_perm[, groups]==2, ])

                       } # end if: paired?


                       # --- Fit mlVAR models ---

                       # Output list
                       l_pair_b <- list()
                       l_models <- list()

                       for(j in 1:2) {

                         # TODO: make this variable specification of dayvar/beepvar less hacky

                         # Neither daybvar/beepvar are specified
                         if(is.null(dayvar) & is.null(beepvar)) {
                           l_pair_b[[j]] <- mlVAR(data = l_data_h0[[j]],
                                                  vars = vars,
                                                  idvar = idvar,
                                                  estimator = estimator,
                                                  contemporaneous = contemporaneous,
                                                  temporal = temporal,
                                                  scale = scale,
                                                  nCores = 1, # we use parallelization across resamples (not nodes in nodewise estimation here)
                                                  verbose = FALSE,
                                                  lags = 1)
                           # Only dayvar is specified
                         } else if (!is.null(dayvar) & is.null(beepvar)) {
                           l_pair_b[[j]] <- mlVAR(data = l_data_h0[[j]],
                                                  vars = vars,
                                                  idvar = idvar,
                                                  estimator = estimator,
                                                  contemporaneous = contemporaneous,
                                                  temporal = temporal,
                                                  scale = scale,
                                                  nCores = 1,
                                                  dayvar = dayvar, # now also provided
                                                  # beepvar = beepvar,
                                                  verbose = FALSE,
                                                  lags = 1)
                           # Both dayvar and beepvar are specified
                         } else {

                           l_pair_b[[j]] <- mlVAR(data = l_data_h0[[j]],
                                                  vars = vars,
                                                  idvar = idvar,
                                                  estimator = estimator,
                                                  contemporaneous = contemporaneous,
                                                  temporal = temporal,
                                                  scale = scale,
                                                  nCores = 1,
                                                  dayvar = dayvar, # now also provided
                                                  beepvar = beepvar,
                                                  verbose = FALSE,
                                                  lags = 1)


                         } # end if: dayvar specified

                         if(saveModels) l_models[[j]] <- l_pair_b[[j]]

                       } # end loop: J=2 groups

                       # All differences are: Group 1 - Group 2
                       diffs_b <- Process_mlVAR(object1 = l_pair_b[[1]],
                                                object2 = l_pair_b[[2]],
                                                empirical = FALSE)

                       outlist_b <- list("diff_between" = diffs_b$diff_between,
                                         "diff_phi_fix" = diffs_b$diff_phi_fix,
                                         "diff_phi_RE_sd" = diffs_b$diff_phi_RE_sd,
                                         "diff_gam_fix" = diffs_b$diff_gam_fix,
                                         "diff_gam_RE_sd" = diffs_b$diff_between,
                                         "Models" = l_models)

                       # Set progress bar
                       if(pbar) setTxtProgressBar(pb, b)

                       return(outlist_b)

                     } # end foreach: over permutations

    # Close down cores, if multi-core
    if(nCores>1) stopCluster(cl)

    # ------ Loop results into objects for Sampling Distribution -----

    # Collect sampling distributions (for now) for:
    # a) between-person partial correlations
    # b.1) VAR/phi fixed effects
    # b.2) VAR/phi random effects sds
    # c.1) Contemp./Gamma fixed effects
    # c.2) Contemp./Gamma random effects sds
    # Create Storage
    a_between <- array(NA, dim=c(p, p, nP))
    a_phi_fixed <- array(NA, dim=c(p, p, nP))
    a_phi_RE_sd <- array(NA, dim=c(p, p, nP))
    a_gam_fixed <- array(NA, dim=c(p, p, nP))
    a_gam_RE_sd <- array(NA, dim=c(p, p, nP))

    for(b in 1:nP) {

      # Fill into arrays
      a_between[, , b] <- out_P[[b]]$diff_between
      a_phi_fixed[, , b] <- out_P[[b]]$diff_phi_fix[, , 1] # TODO: adapt also to higher order lags
      a_phi_RE_sd[, , b] <- out_P[[b]]$diff_phi_RE_sd[, , 1] # TODO: adapt also to higher order lags
      a_gam_fixed[, , b] <- out_P[[b]]$diff_gam_fix
      a_gam_RE_sd[, , b] <- out_P[[b]]$diff_gam_RE_sd

      if(saveModels) l_out_mods[[b]] <- out_P[[b]]$Models

    } # end for: loop in permutations

    if(saveModels) l_out_ret <- l_out_mods else l_out_ret <- NULL

  } # end if: test=permutation



  # ------ Create Test Statistics -----

  l_out_emp <- list()

  for(j in 1:2) {

    # browser()

    # TODO: make this variable specification of dayvar/beepvar less hacky
    if(is.null(dayvar)) {
      # if(is.null(trueMeans)) {

      # neither dayvar/beepvar are specified
      l_out_emp[[j]] <-  mlVAR(data = l_data[[j]],
                               vars = vars,
                               idvar = idvar,
                               estimator = estimator,
                               contemporaneous = contemporaneous,
                               temporal = temporal,
                               scale = scale,
                               nCores = 1,
                               verbose = FALSE,
                               lags = 1)

      # Only dayvar is specified
    } else if (!is.null(dayvar) & is.null(beepvar)) {

      # Fit mlVAR
      l_out_emp[[j]] <-  mlVAR(data = l_data[[j]],
                               vars = vars,
                               idvar = idvar,
                               estimator = estimator,
                               contemporaneous = contemporaneous,
                               temporal = temporal,
                               scale = scale,
                               nCores = 1,
                               dayvar = dayvar,
                               # beepvar = beepvar,
                               verbose = FALSE,
                               lags = 1)

      # Both are specified
    } else {

      # Fit mlVAR
      l_out_emp[[j]] <-  mlVAR(data = l_data[[j]],
                               vars = vars,
                               idvar = idvar,
                               estimator = estimator,
                               contemporaneous = contemporaneous,
                               temporal = temporal,
                               scale = scale,
                               nCores = 1,
                               dayvar = dayvar,
                               beepvar = beepvar,
                               verbose = FALSE,
                               lags = 1)


    } # end if: dayvar specified


    # } else {
    #
    #   trueMeans_j <- trueMeans[trueMeans$id %in% l_data[[j]]$id, c(idvar, vars)]
    #
    #   # Fit mlVAR; note: true means spec currently only possible with mlVAR dev
    #   l_out_emp[[j]] <-  mlVAR(data = l_data[[j]],
    #                            vars = vars,
    #                            idvar = idvar,
    #                            estimator = estimator,
    #                            contemporaneous = contemporaneous,
    #                            temporal = temporal,
    #                            scale = scale,
    #                            trueMeans = trueMeans_j,
    #                            nCores = 1,
    #                            verbose = FALSE,
    #                            lags = 1)
    #
    # } # end if: true means specified?

    # } else {

    # if(is.null(trueMeans)) {



    # } else {
    #
    #   # Get true means for group J
    #   trueMeans_j <- trueMeans[trueMeans$id %in% l_data[[j]]$id, c(idvar, vars)]
    #
    #   # Fit mlVAR; note: true means spec currently only possible with mlVAR dev
    #   l_out_emp[[j]] <-  mlVAR(data = l_data[[j]],
    #                            vars = vars,
    #                            idvar = idvar,
    #                            estimator = estimator,
    #                            contemporaneous = contemporaneous,
    #                            temporal = temporal,
    #                            scale = scale,
    #                            trueMeans = trueMeans_j,
    #                            nCores = 1,
    #                            dayvar = dayvar,
    #                            beepvar = beepvar,
    #                            verbose = FALSE,
    #                            lags = 1)
    #
    # } # end if: true means specified?



  } # Loop: 2 groups

  # Final step
  if(pbar)  setTxtProgressBar(pb, nP + 1)
  runtime <- proc.time()[3] - timer_total


  # --- Matrices with True differences ---

  diffs_true <- Process_mlVAR(object1 = l_out_emp[[1]],
                              object2 = l_out_emp[[2]])


  # ------ Compute p-values based on permutation test -----

  if(test == "permutation") {

    # b.1) VAR: fixed effects
    m_pval_phi_fix <- matrix(NA, p, p)
    for(i in 1:p) for(j in 1:p) m_pval_phi_fix[i,j] <- mean(abs(a_phi_fixed[i,j,]) > abs(diffs_true$diff_phi_fix[i,j,]))

    # b.2) VAR: RE sds
    m_pval_phi_RE_sd <- matrix(NA, p, p)
    for(i in 1:p) for(j in 1:p) m_pval_phi_RE_sd[i,j] <- mean(abs(a_phi_RE_sd[i,j,])>abs(diffs_true$diff_phi_RE_sd[i,j,]))

    # c.1) Contemp: fixed effects
    m_pval_gam_fixed <- matrix(NA, p, p)
    for(i in 2:p) for(j in 1:(i-1)) m_pval_gam_fixed[i,j] <- mean(abs(a_gam_fixed[i,j,])>abs(diffs_true$diff_gam_fix[i,j]))

    # c.2) Contemp: RE sds
    m_pval_gam_RE_sd <- matrix(NA, p, p)
    for(i in 2:p) for(j in 1:(i-1)) m_pval_gam_RE_sd[i,j] <- mean(abs(a_gam_RE_sd[i,j,])>abs(diffs_true$diff_gam_RE_sd[i,j]), na.rm=TRUE)
    # Especially in conditions with low SNR, Random effects cannot be estimated and mlVAR/lme4 returns errors

    # a) Between
    m_pval_btw <- matrix(NA, p, p)
    for(i in 2:p) for(j in 1:(i-1)) m_pval_btw[i,j] <- mean(abs(a_between[i,j,])>abs(diffs_true$diff_between[i,j]), na.rm=TRUE)
    # Note: the na.rm=TRUE in the means is there so we can exclude the rare cases in which between-networks
    #       have not been calculated, because random effects variances of intercepts were set to zero

    # ------ Create Output List -----

    # Save the two mlVAR models estimated on the empirical data?
    if(saveEmpModels==FALSE) l_out_emp <- NULL

    outlist <- list("Call" = Call,
                    "EmpDiffs" = list("Lagged_fixed" = diffs_true$diff_phi_fix,
                                      "Lagged_random" = diffs_true$diff_phi_RE_sd,
                                      "Contemp_fixed" = diffs_true$diff_gam_fix,
                                      "Contemp_random" = diffs_true$diff_gam_RE_sd,
                                      "Between" = diffs_true$diff_between),
                    "Pval" = list("Lagged_fixed" = m_pval_phi_fix,
                                  "Lagged_random" = m_pval_phi_RE_sd,
                                  "Contemp_fixed" = m_pval_gam_fixed,
                                  "Contemp_random" = m_pval_gam_RE_sd,
                                  "Between" = m_pval_btw),
                    "SampDist" = list("Lagged_fixed" = a_phi_fixed,
                                      "Lagged_random" = a_phi_RE_sd,
                                      "Contemp_fixed" = a_gam_fixed,
                                      "Contemp_random" = a_gam_RE_sd,
                                      "Between" = a_between),
                    "Models" = l_out_ret,
                    "EmpModels" = l_out_emp,
                    "Runtime_min" = runtime / 60)

  } # end if: permutation

  # ------ Compute p-values based on standard errors [parametric test] -----
  if(test == "parametric") {

    # --- Lagged effects ---
    bet_1 <- l_out_emp[[1]]$results$Beta$mean[, , 1]
    bet_1se <- l_out_emp[[1]]$results$Beta$SE[, , 1]
    bet_2 <- l_out_emp[[2]]$results$Beta$mean[, , 1]
    bet_2se <- l_out_emp[[2]]$results$Beta$SE[, , 1]
    # Compute degrees of freedom
    if(partest == "Student") df = n_subj - 2
    if(partest == "Welch") df = ((bet_1se^2 / n_subj_G1 + bet_2se^2 / n_subj_G2)^2) /
      ((bet_1se^4 / (n_subj_G1^2 * (n_subj_G1 - 1))) + (bet_2se^4 / (n_subj_G2^2 * (n_subj_G2 - 1))))
    # Compute t-statistic
    t_stat <- abs(bet_1-bet_2)/ sqrt(bet_1se^2 + bet_2se^2)
    # get pvalues
    m_pval_phi_fix <- pt(t_stat, df = df, lower.tail = FALSE) * 2 # times two to make 2-sided

    # --- Contemporaneous effects ---
    # Average across nodewise reg to get estimates and SEs
    bet_1 <- (l_out_emp[[1]]$results$Gamma_Theta$mean + t(l_out_emp[[1]]$results$Gamma_Theta$mean)) / 2
    bet_1se <- (l_out_emp[[1]]$results$Gamma_Theta$SE + t(l_out_emp[[1]]$results$Gamma_Theta$SE)) / 2
    bet_2 <-  (l_out_emp[[2]]$results$Gamma_Theta$mean + t(l_out_emp[[2]]$results$Gamma_Theta$mean)) / 2
    bet_2se <- (l_out_emp[[2]]$results$Gamma_Theta$SE + t(l_out_emp[[2]]$results$Gamma_Theta$SE)) / 2
    # Compute degrees of freedom
    if(partest == "Student") df = n_subj-2
    if(partest == "Welch") df = ((bet_1se^2 / n_subj_G1 + bet_2se^2 / n_subj_G2)^2) /
      ((bet_1se^4 / (n_subj_G1^2 * (n_subj_G1 - 1))) + (bet_2se^4 / (n_subj_G2^2 * (n_subj_G2 - 1))))
    # Compute t-statistic
    t_stat <- abs(bet_1-bet_2)/ sqrt(bet_1se^2 + bet_2se^2)
    # get pvalues
    m_pval_gam_fixed <- pt(t_stat, df = df, lower.tail = FALSE) * 2 # times two to make 2-sided
    m_pval_gam_fixed[upper.tri(m_pval_gam_fixed)] <- NA
    diag(m_pval_gam_fixed) <- NA

    # --- Between ---
    # Average across nodewise reg to get estimates and SEs
    bet_1 <- (l_out_emp[[1]]$results$Gamma_Omega_mu$mean + t(l_out_emp[[1]]$results$Gamma_Omega_mu$mean)) / 2
    bet_1se <- (l_out_emp[[1]]$results$Gamma_Omega_mu$SE + t(l_out_emp[[1]]$results$Gamma_Omega_mu$SE)) / 2
    bet_2 <- (l_out_emp[[2]]$results$Gamma_Omega_mu$mean + t(l_out_emp[[2]]$results$Gamma_Omega_mu$mean)) / 2
    bet_2se <- (l_out_emp[[2]]$results$Gamma_Omega_mu$SE + t(l_out_emp[[2]]$results$Gamma_Omega_mu$SE)) / 2
    # Compute degrees of freedom
    if(partest == "Student") df = n_subj-2
    if(partest == "Welch") df = ((bet_1se^2 / n_subj_G1 + bet_2se^2 / n_subj_G2)^2) /
      ((bet_1se^4 / (n_subj_G1^2 * (n_subj_G1 - 1))) + (bet_2se^4 / (n_subj_G2^2 * (n_subj_G2 - 1))))
    # Compute t-statistic
    t_stat <- abs(bet_1-bet_2)/ sqrt(bet_1se^2 + bet_2se^2)
    # get p-values
    m_betw_sign <- pt(t_stat, df = df, lower.tail = FALSE) * 2 # times two to make 2-sided
    m_betw_sign[upper.tri(m_pval_gam_fixed)] <- NA
    diag(m_betw_sign) <- NA


    # ------ Create Output List -----

    # Save the two mlVAR models estimated on the empirical data?
    if(saveEmpModels==FALSE) l_out_emp <- NULL


    outlist <- list("Call" = Call,
                    "EmpDiffs" = list("Lagged_fixed" = diffs_true$diff_phi_fix,
                                      "Lagged_random" = diffs_true$diff_phi_RE_sd,
                                      "Contemp_fixed" = diffs_true$diff_gam_fix,
                                      "Contemp_random" = diffs_true$diff_gam_RE_sd,
                                      "Between" = diffs_true$diff_between),
                    "Pval" = list("Lagged_fixed" = m_pval_phi_fix,
                                  "Contemp_fixed" = m_pval_gam_fixed,
                                  "Between" = m_betw_sign),
                    "EmpModels" = l_out_emp,
                    "Runtime_min" = runtime / 60)

  } # end if: test=parametric


  # ------ Return Output -----

  class(outlist) <- c("list", "mlVAR_GC")

  return(outlist)

} # eoF





