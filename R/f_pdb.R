
# Takes dayvar and beepvar, and returns for a lag-1 VAR model whether any given time point can be predicted

f_pdb <- function(dayvar, beepvar) {
  n <- length(dayvar)
  v_pdb <- rep(NA, n)
  v_pdb[1] <- FALSE
  for(i in 2:n) {
    day_eq <- dayvar[i] == dayvar[i-1]
    beep_eq <- beepvar[i] == (beepvar[i-1]+1)
    v_pdb[i] <- ifelse(day_eq & beep_eq, TRUE, FALSE)
  }
  return(v_pdb)
} #eoF

