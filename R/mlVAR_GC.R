# jonashaslbeck@protonmail; April 17, 2023

# ------------------------------------------------------------
# -------- Function for Permutation Test ---------------------
# ------------------------------------------------------------

# Inputs:
# - two datasets with same structure + same inputs as for mlVAR()
# - number of samples in permutation test

# Output:
# - sampling distribution for each parameter
# - test statistic based on the mlVAR models estimated on the two datasets

mlVAR_GC <- function(data1, # dataset of group 1
                     data2, # dataset of group 2
                     vars, # variables to be included in mlVAR (same in both data sets)
                     idvar, # variable indicating the nesting/subject id (same in both data sets)
                     dayvar = NULL,
                     beepvar = NULL,
                     estimator, # same as in ml
                     contemporaneous, # same as in ml
                     temporal, # same as in ml
                     nCores = 1,
                     nP = 500, # number of samples in permutation test
                     saveModels = FALSE, # if TRUE, all models are saved; defaults to FALSE to save memory
                     verbose = FALSE, # if TRUE, verbose mode is activated in foreach
                     pbar = TRUE # if TRUE, a progress bar is being shown
) {


  # ------ Input Checks -----

  # (1) Are the data numerical?
  v_class <- apply(cbind(data1[, vars], data2[, vars]), 2, function(x) class(x))
  if(any( !(v_class %in% c("numeric", "integer")) )) stop("Modeled variables need to be provided as integer or numeric variables.")

  # (2) Can we find the all variables?
  check_coln1 <- c(vars, idvar) %in% colnames(data1)
  if(any(!check_coln1)) stop("Specified variable names could not be found in dataset 1.")
  check_coln2 <- c(vars, idvar) %in% colnames(data2)
  if(any(!check_coln2)) stop("Specified variable names could not be found in dataset 2.")

  # (3) Are IDs unique across datasets?
  # Get IDs
  ids1 <- sapply(data1[, idvar], as.character)
  ids2 <- sapply(data2[, idvar], as.character)
  v_ids <- c(ids1, ids2)
  v_u_ids <- unique(v_ids)
  u_ids1 <- unique(ids1)
  u_ids2 <- unique(ids2)

  v_intersec <- intersect(u_ids1, u_ids2)
  if(length(v_intersec) > 0) stop("IDs need to be unique across two datasets.")


  # ------ Collect passed down arguments -----
  if(missing(estimator)) estimator <- "default"
  if(missing(contemporaneous)) contemporaneous <- "orthogonal"
  if(missing(temporal)) temporal <- "orthogonal"
  if(missing(nCores)) nCores <- 1


  # ------ Get Basic Info -----

  # number of variables
  p <- length(vars)

  # Combine in list/single matrix for shorter code below
  l_data <- list(data1, data2)
  m_data_cmb <- rbind(data1, data2)

  # Get Number of subjects
  v_Ns <- c(length(u_ids1), length(u_ids2))
  totalN <- sum(v_Ns)


  # ------ Loop Over Permutations -----

  # Storage
  if(saveModels) l_out_mods <- list(vector("list", length = nP),
                                    vector("list", length = nP))

  # Setup cores, if multi-core
  if(nCores > 1) {
    cl <- makeCluster(nCores, outfile = "")
    registerDoParallel(cl)
  }

  if(pbar) pb <- txtProgressBar(0, nP+1, style = 3)

  timer_total <- proc.time()[3]

  # Call foreach:
  out_P <- foreach(b = 1:nP,
                   .packages = c("mlVAR", "mnet"),
                   .export = c("m_data_cmb", "vars", "idvar", "estimator",
                               "contemporaneous", "temporal", "totalN", "v_Ns",
                               "v_ids"),
                   .verbose = verbose) %dopar% {

                     # --- Make permutation ---
                     # This is done in a way that keeps the size in each group exactly the same as in the real groups
                     v_ids_rnd <- v_u_ids[sample(1:totalN, size=totalN, replace=FALSE)]
                     v_ids_1 <- v_ids_rnd[1:v_Ns[1]]
                     v_ids_2 <- v_ids_rnd[(v_Ns[1]+1):totalN]

                     # Split data based on permutations
                     data_h0_1 <- m_data_cmb[v_ids %in% v_ids_1, ]
                     data_h0_2 <- m_data_cmb[v_ids %in% v_ids_2, ]
                     l_data_h0 <- list(data_h0_1, data_h0_2)

                     # --- Fit mlVAR models ---

                     # Output list
                     l_pair_b <- list()
                     l_models <- list()

                     for(j in 1:2) {

                       # TODO: make this variable specification of dayvar/beepvar less hacky
                       if(is.null(dayvar)) {
                         l_pair_b[[j]] <- mlVAR(data = l_data_h0[[j]],
                                                vars = vars,
                                                idvar = idvar,
                                                estimator = estimator,
                                                contemporaneous = contemporaneous,
                                                temporal = temporal,
                                                nCores = 1, # we use parallelization across resamples (not nodes in nodewise estimation here)
                                                verbose = FALSE,
                                                lags = 1) # TODO: later allow also higher order lags (see also below)
                       } else {
                         l_pair_b[[j]] <- mlVAR(data = l_data_h0[[j]],
                                                vars = vars,
                                                idvar = idvar,
                                                estimator = estimator,
                                                contemporaneous = contemporaneous,
                                                temporal = temporal,
                                                nCores = 1,
                                                dayvar = dayvar, # now also provided
                                                beepvar = beepvar,
                                                verbose = FALSE,
                                                lags = 1) # TODO: later allow also higher order lags (see also below)
                       } # end if: dayvar specified

                       if(saveModels) l_models[[j]] <- l_pair_b[[j]]

                     } # end loop: J=2 groups

                     # browser()

                     # All differences are: Group 1 - Group 2
                     diffs_b <- Process_mlVAR(object1 = l_pair_b[[1]],
                                              object2 = l_pair_b[[2]])

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


  # ------ Create Test Statistics -----

  l_out_emp <- list()

  for(j in 1:2) {

    # TODO: make this variable specification of dayvar/beepvar less hacky
    if(is.null(dayvar)) {
      l_out_emp[[j]] <-  mlVAR(data = l_data[[j]],
                               vars = vars,
                               idvar = idvar,
                               estimator = estimator,
                               contemporaneous = contemporaneous,
                               temporal = temporal,
                               nCores = 1,
                               verbose = FALSE,
                               lags = 1)
    } else {
      l_out_emp[[j]] <-  mlVAR(data = l_data[[j]],
                               vars = vars,
                               idvar = idvar,
                               estimator = estimator,
                               contemporaneous = contemporaneous,
                               temporal = temporal,
                               nCores = 1,
                               dayvar = dayvar,
                               beepvar = beepvar,
                               verbose = FALSE,
                               lags = 1)
    } # end if: dayvar specified

  } # Loop: 2 groups

  # Final step
  if(pbar)  setTxtProgressBar(pb, nP + 1)
  runtime <- proc.time()[3] - timer_total


  # --- Matrices with True differences ---

  diffs_true <- Process_mlVAR(object1 = l_out_emp[[1]],
                              object2 = l_out_emp[[2]])


  # ------ Compute p-values -----

  # a) between
  m_pval_btw <- matrix(NA, p, p)
  for(i in 2:p) for(j in 1:(i-1)) m_pval_btw[i,j] <- mean(abs(a_between[i,j,])>abs(diffs_true$diff_between[i,j]), na.rm=TRUE)
  # Note: the na.rm=TRUE in the means is there so we can exclude the rare cases in which between-networks
  #       have not been calculated, because random effects variances of intercepts were set to zero

  # b.1) VAR: fixed effects
  m_pval_phi_fix <- matrix(NA, p, p)
  for(i in 1:p) for(j in 1:p) m_pval_phi_fix[i,j] <- mean(abs(a_phi_fixed[i,j,])>abs(diffs_true$diff_phi_fix[i,j,]))

  # b.2) VAR: RE sds
  m_pval_phi_RE_sd <- matrix(NA, p, p)
  for(i in 1:p) for(j in 1:p) m_pval_phi_RE_sd[i,j] <- mean(abs(a_phi_RE_sd[i,j,])>abs(diffs_true$diff_phi_RE_sd[i,j,]))

  # c.1) Contemp: fixed effects
  m_pval_gam_fixed <- matrix(NA, p, p)
  for(i in 2:p) for(j in 1:(i-1)) m_pval_gam_fixed[i,j] <- mean(abs(a_gam_fixed[i,j,])>abs(diffs_true$diff_gam_fix[i,j]))

  # c.2) Contemp: RE sds
  m_pval_gam_RE_sd <- matrix(NA, p, p)
  for(i in 2:p) for(j in 1:(i-1)) m_pval_gam_RE_sd[i,j] <- mean(abs(a_gam_RE_sd[i,j,])>abs(diffs_true$diff_gam_RE_sd[i,j]))


  # ------ Create Output List -----

  if(saveModels) l_out_ret <- l_out_mods else l_out_ret <- NULL

  outlist <- list("TrueDiffs" = list("Between" = diffs_true$diff_between,
                                     "Phi_mean" = diffs_true$diff_phi_fix,
                                     "Phi_sd" = diffs_true$diff_phi_RE_sd,
                                     "Gam_mean" = diffs_true$diff_gam_fix,
                                     "Gam_sd" = diffs_true$diff_gam_RE_sd),
                  "Pval" = list("Between" = m_pval_btw,
                                "Phi_mean" = m_pval_phi_fix,
                                "Phi_sd" = m_pval_phi_RE_sd,
                                "Gam_mean" = m_pval_gam_fixed,
                                "Gam_sd" = m_pval_gam_RE_sd),
                  "SampDist" = list("Between" = a_between,
                                    "Phi_mean" = a_phi_fixed,
                                    "Phi_sd" = a_phi_RE_sd,
                                    "Gam_mean" = a_gam_fixed,
                                    "Gam_sd" = a_gam_RE_sd),
                  "Models" = l_out_ret,
                  "Runtime_min" = runtime / 60)

  # ------ Return Output -----

  return(outlist)


} # eoF





