# PART 1. 
# How CPDstergm authors implemented their real data applications on the DIJA data from the ecp package
# https://github.com/allenkei/CPDstergm_demo/blob/main/experiment_real.R

#%%%%%%%%%%%%%%%%%%%%#
# DJIA Data from ecp #
#%%%%%%%%%%%%%%%%%%%%#
library(ecp)
data(DJIA)
market <- DJIA$market

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Convert to adjacency matrices #
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
date_vec <- DJIA$dates[1:1138]
rownames(market) <- date_vec
start <- which(date_vec == '2007-01-01')
end <- which(date_vec == '2010-01-04')
date_range <- start:end
df <- list()
for(i in 1:length(date_range)){
  temp <- market[date_range[i]:(date_range[i]+3),]
  temp <- ifelse(cor(temp)< 0, 1, 0)
  diag(temp) <- 0
  df[[i]] <- temp
}

#%%%%%%%%%%%%%%%%%%%%%%%%#
# Get network statistics #
#%%%%%%%%%%%%%%%%%%%%%%%%#
library(ergm)
net_stats = matrix(0, nrow = length(df), ncol = 3)
node_degree <- rep(0, dim(df[[1]])[1])
for (t in seq_len(length(df))) {
  mat <- df[[t]]
  node_degree <- node_degree + rowSums(mat)
}; rm(mat)
node_status <- ifelse(node_degree > median(node_degree), "H", "M")
net_stats <- matrix(0, nrow = length(df), ncol = 3)
for (i in seq_along(df)) {
  yt <- network(df[[i]], directed = FALSE)
  # attach nodal attribute
  set.vertex.attribute(yt, "node_status", node_status)
  net_stats[i, ] <- summary(
    yt ~ edges + triangles + nodematch("node_status")
  )
}
head(net_stats)
# sum(is.na(net_stats))


# PART 2. 
# Fit the bcmlr method, visualize the results, and fit e.divisive on the same data. 

#%%%%%%%%%%%#
# Fit bcmlr #
#%%%%%%%%%%%#
X = as.matrix(net_stats)
out = bcmlr_model_select(data = X, init = "even", prior_beta = "Gaussian", prior_kappa = "default", 
                         thinning = 5, min_size = 10,
                         alpha_f = 0.1, threshold = 0.5, max_num_cp = 10, 
                         num_iter = 15000, num_warmup = 10000, print_progress = TRUE, print_outputs = TRUE)
str(out)

#%%%%%%%%%%%%%%%%%%%%#
# Plot the numCPs #
#%%%%%%%%%%%%%%%%%%%%#
hist(out$num_cp_dist,
     breaks = seq(-0.5, 10.5, by = 1),
     prob = TRUE,
     labels = TRUE,
     ylim = c(0, 0.4),
     col = "lightgrey",
     main = "Posterior distribution of the number of fitted changepoint that correspond to true changepoints",
     xlab = "Number of changepoints",
     ylab = "Posterior probability",
     xaxt = "n")
axis(1, at = 0:10)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Get a list of days: 2004-09-15 to 2005-05-04 #
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
date_list = DJIA$dates[start:end]
# Check the dates when BCMLR predicted a change.
date_list[out$Kappa_mode]

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# CPs AGAINST THE SERIES: ggplot #
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

net_stats <- as.data.frame(net_stats)
covariate_names <- c("Edges", "Triangles", "Homophily")
colnames(net_stats) <- covariate_names

# 1) CP summary (95% CI)
lbs   <- apply(out$Kappa, 2, quantile, probs = 0.025)
ubs   <- apply(out$Kappa, 2, quantile, probs = 0.975)
k_ci  <- tibble(xmin = lbs, xmax = ubs)
k_mode <- tibble(x = out$Kappa_mode)

# 2) Weekly dates
dates <- seq.Date(from = as.Date("2007-01-01"),
                  to   = as.Date("2010-01-04"),
                  by   = "week")
# Make sure length matches
stopifnot(length(dates) == nrow(net_stats))

# 3) Long data with dates
DJIA_net <- net_stats |>
  mutate(date = dates) |>
  select(date, Edges, Triangles, Homophily) |>
  pivot_longer(-date, names_to = "series", values_to = "value")
# Fix factor order
DJIA_net$series <- factor(DJIA_net$series,
                          levels = covariate_names)
# Convert CP summaries (t indices) into dates
k_ci   <- k_ci  |> mutate(xmin = dates[xmin], xmax = dates[xmax])
k_mode <- k_mode |> mutate(x = dates[x])
k_true = 
k_true <- k_true |> mutate(x = dates[x])

# 4) Plot
ggplot(DJIA_net, aes(x = date, y = value)) +
  geom_rect(data = k_ci, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax),
            ymin = -Inf, ymax = Inf,
            fill = "blue", alpha = 0.18, color = NA) +
  geom_vline(data = k_mode,
             aes(xintercept = x),
             linetype = "dashed", linewidth = 0.7, color = "blue") +
  geom_line(linewidth = 0.5, color = "black") +
  facet_wrap(~ series, ncol = 1, scales = "free_y") +
  labs(
    x = "Date",
    y = NULL,
    title = "DJIA network statistics with posterior changepoint modes and CIs",
    subtitle = "Blue dashed lines = posterior modes; Blue bands = 95% credible intervals"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(size = 11, face = "bold"),
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 10),
    axis.title.x = element_text(size = 11),
    panel.spacing = unit(6, "pt")
  ) +
  scale_x_date(
    expand = expansion(mult = c(0.01, 0.02)),
    date_breaks = "6 months",     # tick every ~6 months
    date_labels = "%b %Y"         # format: "Jan 2007"
  )

#%%%%%%%%%%%%%%%%%%%%#
# FIT another method #
#%%%%%%%%%%%%%%%%%%%%#
X = as.matrix(net_stats)
est = e.divisive(X, min.size = 10)$estimates 
length(est)
est
