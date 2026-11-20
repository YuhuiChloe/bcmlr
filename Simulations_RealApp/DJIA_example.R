# Run after installing the updated bcmlr package.
# The 15000-iteration runs are analysis settings, not a short software test.

library(bcmlr)
data("DJIA_networks", package = "bcmlr")
X <- DJIA_networks$net_stats
date_list <- DJIA_networks$date_list

set.seed(2026)
fit <- bcmlr_model_select(data = X, max_num_cp = 10,
                          min_size = 10, alpha_f = 0.1,
                          num_iter = 15000, num_warmup = 10000,
                          print_progress = TRUE)
print(fit$num_cp)
print(fit$probs_num_cp)
print(bcmlr::plot_histogram(fit$num_cp_dist))

if (fit$num_cp > 0) {
  print(date_list[fit$Kappa_mode])
  print(bcmlr::plot_changepoints(fit$Kappa, X = X, dates = date_list))
  print(bcmlr::plot_beta_diff_CIs(fit$Beta, coef_labels = colnames(X)))
  print(bcmlr::plot_log_prob_diff(fit$P, fit$Kappa_mode))
} else {
  message("Zero changepoints selected; conditional difference plots are skipped.")
}

# Optional alternative: condition on three changepoints.
# set.seed(2026)
# fixed_fit <- bcmlr(data = X, num_CP = 3, min_size = 10,
#                    num_iter = 15000, num_warmup = 10000)
# print(bcmlr::plot_changepoints(fixed_fit$Kappa, X = X, dates = date_list))
# print(bcmlr::plot_beta_diff_CIs(fixed_fit$Beta, coef_labels = colnames(X)))
# print(bcmlr::plot_log_prob_diff(fixed_fit$P, fixed_fit$Kappa_mode))
