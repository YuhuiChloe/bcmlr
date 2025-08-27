bcmlr_model_select <- function(data, init, thinning, prior, alpha_f = 0.1, threshold = 0.5, min_size = 30, 
                               max_num_cp = 5, num_iter = 5000, num_warmup = 2500, print_progress = FALSE){
  # Select the number of changepoints
  old_out = bcmlr(data, init = init, num_CP = max_num_cp, prior = prior, alpha_f = alpha_f, thinning = thinning,
              min_size = min_size, num_iter = num_iter, num_warmup = num_warmup,
              print_outputs = FALSE, print_progress = print_progress, model_selection = TRUE)
  boolean = (old_out$AUC < threshold) # TRUE if less than threshold, FALSE otherwise
  rej = rowSums(!boolean) # Count number of FALSE in each row
  tb = table(factor(rej, levels = seq_len(max_num_cp)))/length(rej) 
  num_cp = as.integer(names(tb)[which.max(tb)]) # the posterior mode is taken as in the CP estimate
  probs = as.numeric(tb)
  names(probs) = seq_len(max_num_cp) 
  # Refit with all the data: 
  out = bcmlr(data, init = init, num_CP = num_cp, prior = prior, thinning = 1, min_size = min_size, num_iter = num_iter, num_warmup = num_warmup,
              print_outputs = FALSE, print_progress = print_progress, model_selection = FALSE)
  out$num_cp = num_cp # estimate of the number of CPs
  out$max_cps = old_out$Kappa # the max number of CPs from the "old output" before refitting on all the data
  out$num_cp_dist = rej # posterior distribution of the number of CPs
  out$probs_num_cp = probs 
  return(out)
}
