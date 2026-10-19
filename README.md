# bcmlr

Generalized Bayesian multiple changepoint analysis with multinomial log loss.
`bcmlr` fits a specified number of changepoints; `bcmlr_model_select` selects a
count using held-out AUC confidence limits and refits the full series. Posterior
summaries describe changepoint locations and the covariates distinguishing
neighboring segments.

**Development version.** This package corrects substantive defects in the
original scripts. Existing paper results have not been revalidated with these
corrections. See [AUDIT.md](AUDIT.md) for evidence and numerical implications,
and [RELEASE_CHECKLIST.md](RELEASE_CHECKLIST.md) for release decisions. The
original scripts are retained in `legacy/`, outside the installed package.

## Installation

Install this reviewed development snapshot from its local directory:

```r
install.packages(c("remotes", "roxygen2", "testthat"))
remotes::install_local("path/to/bcmlr", dependencies = NA)
library(bcmlr)
```

After the package files have been committed and pushed to the repository:

```r
remotes::install_github("YuhuiChloe/bcmlr")
library(bcmlr)
```

The existing upstream revision is a collection of scripts, so the second command
requires the package update to be published first. Core installation does not
require Python, `ergm`, `LearnBayes`, `MASS`, or `PGdensity`. `ecp` is optional for
reconstructing the example networks or using its initializer; `reticulate` and
Python `ruptures` are optional for `kcp` and `multirank` initialization.

## Load the DJIA example

```r
djia <- load_djia()
X <- djia$net_stats
dates <- djia$dates
dim(X) # 158 observations by 3 network statistics
colnames(X) # Edges, Triangles, Homophily
```

The small RDS from the original repository is included with the package. It
contains network statistics, not adjacency matrices or posterior draws. The
158 networks represent four-week rolling windows of weekly returns for
29 stocks, with weekly date labels from January 1, 2007 through January 4,
2010. Negative correlations indicate edges. Homophily counts edges whose stocks
share the risk-orientation label defined by their cumulative degrees.

The original `ecp::DJIA` dates are reverse chronological. The construction uses
rows `i:(i+3)` in that native order and labels the network with `dates[i]`, just
as the original analysis does. Thus these windows extend toward earlier dates
in the source; the output network series itself is ordered chronologically.

The construction follows [Kei et al. (2025)](https://github.com/allenkei/CPDstergm_demo)
using `DJIA` from the [ecp package](https://cran.r-project.org/package=ecp).
To reconstruct the networks and statistics:

```r
install.packages("ecp")
rebuilt <- prepare_djia()
stopifnot(isTRUE(all.equal(unname(rebuilt$net_stats), unname(X))))
length(rebuilt$networks) # 158
dim(rebuilt$networks[[1]]) # 29 by 29
```

`prepare_djia()` checks that every window is available and rejects undefined
correlations. Its direct edge, triangle, and homophily counts implement the same
undirected statistics as the original `ergm` formulas.

## Select the number of changepoints

```r
set.seed(20261009)
selection <- bcmlr_model_select(
  X, max_num_cp = 10, min_size = 10,
  num_iter = 15000, num_warmup = 10000,
  print_progress = TRUE
)
selection$num_cp
selection$probs_num_cp
plot_histogram(selection)
```

These use the DJIA settings in the paper: Gaussian coefficient prior with
covariance `3 * I`, AUC confidence level 90%, threshold 0.5, and every fifth
observation held out. The default `alpha_f = 0.1` preserves the historical code
and DJIA settings; the methods text states a general default of 0.05. Specify
`alpha_f = 0.05` for 95% AUC confidence intervals. `thinning` denotes observation
holdout spacing, **not** thinning of the MCMC draws.

The histogram corresponds to the quantity displayed in paper Figure 2: the
number of fitted changepoints supported by the AUC criterion in each posterior
draw. It is not a reversible-jump posterior distribution over different models.
The returned fit already refits all observations with the modal count, including
a valid one-segment fit when that count is zero.

## Fit a specified number of changepoints

For an analysis using the selected count, the refit is already in `selection`.
The paper conditions its DJIA interpretation on three changepoints. To make the
same conditioning choice explicitly:

```r
set.seed(20261010)
fit <- bcmlr(
  X, num_CP = 3, min_size = 10,
  num_iter = 15000, num_warmup = 10000,
  print_progress = TRUE
)
dates[fit$Kappa_mode]
```

The seeds above are new workflow seeds, not the seeds of the published analysis.
The corrected sampler may select another count or estimate different locations.
The saved historical posterior draws and original RNG state were not supplied,
so exact reproduction of the printed figures is not claimed.

## Plot the DJIA posterior summaries

All plotting functions return `ggplot` objects and their numerical summaries
are accessible through `$data`.

```r
# Figure 1: observed statistics, changepoint modes and 95% credible intervals
p1 <- plot_changepoints(fit, dates = dates)

# Figure 2: supported changepoint counts from the selection stage
p2 <- plot_histogram(selection)

# Figure 3: neighboring coefficient differences
p3 <- plot_beta_diff_CIs(fit)

# Figure 4: posterior log-probability differences and pointwise 95% intervals
p4 <- plot_log_prob_diff(fit)

print(p1)
print(p2)
print(p3)
print(p4)
ggplot2::ggsave("DJIA-log-probability-differences.pdf", p4,
                width = 10, height = 3.5)
```

The fourth plot computes `x_i^T (beta_next - beta_current)` for every retained
draw, using the stored standardized covariates. This is exactly the difference
between the two log classifier probabilities. Its posterior mean is generally
different from the log ratio of posterior mean probabilities. Each panel uses
the two neighboring segments determined by the supplied changepoint point
estimates. The bands are pointwise credible intervals.

The coefficient and log-probability panels can also be combined:

```r
install.packages("patchwork")
patchwork::wrap_plots(p3, p4, ncol = 1) +
  patchwork::plot_annotation(tag_levels = "a")
```

Every interval conditions on the fitted count. If the count distribution is
diffuse, refit the plausible counts separately and compare them, as discussed
in the paper. The interpretation functions report that no contrasts exist for
a zero-changepoint fit; `plot_changepoints()` still displays the observed series.

## Main arguments

| Argument | Meaning and default |
|---|---|
| `data` | Numeric observations in rows; vector, matrix or numeric data frame |
| `num_CP` | Fixed nonnegative count for `bcmlr` |
| `max_num_cp` | Upper bound for selection; 5 |
| `init` | `"even"`; also validated numeric boundaries or optional external methods |
| `prior_beta` | `"Gaussian"` or `"horseshoe"`; Gaussian by default |
| `prior_kappa` | Segment-length prior (`"default"`) or a joint prior mass function |
| `sd_beta` | Gaussian **variance**; 3, despite its historical name |
| `min_size` | Minimum segment length in original observations; 30 |
| `thinning` | Holdout spacing; 1 for fixed-count fits, 5 for selection |
| `alpha_f` | AUC confidence-interval error probability; 0.1 |
| `threshold` | Lower AUC confidence-limit cutoff; 0.5 |
| `num_iter`, `num_warmup` | Total iterations and discarded warmup; 5000 and 2500 |
| `store_prob` | Store probability draws; TRUE |

For DeLong intervals, selection requires `min_size >= 2 * thinning`. Both
classes must contain at least two held-out observations. Original-series
boundaries are restricted to training indices during selection, then unrestricted
in the full-data refit. Some tightly packed boundary configurations may therefore
be infeasible on the selection grid even if the simple total-length check passes.

## Returned objects and compatibility

Let `S = num_iter - num_warmup`, `p` be the number of covariates, and
`J = num_CP + 1`.

| Field | Shape or meaning |
|---|---|
| `Kappa` | `S x num_CP` draws, on the original time scale |
| `Kappa_mode` | Marginal posterior modes; ties choose the smaller index |
| `Beta` | `S x p x J`; the final class is zero |
| `P` | `S x N x J`, including held-out rows; optional |
| `Beta_mean`, `P_mean` | Posterior means; `P_mean` is optional |
| `data`, `X`, `center`, `scale` | Raw and standardized data and preprocessing |
| `num_cp_dist`, `probs_num_cp` | Selection counts and their empirical probabilities |
| `selection_AUC`, `max_cps` | AUC lower limits and upper-bound fit modes |

Single-changepoint output now follows the same convention as multiple-change
output. The original binary script instead stored a coefficient for the later
class against a zero earlier class. To reuse that old coefficient matrix, make
the sign conversion explicit:

```r
# Only for a legacy ONE-changepoint result:
B_old <- as.matrix(old_fit$Beta)
B <- array(0, c(nrow(B_old), ncol(B_old), 2))
B[, , 1] <- -B_old
plot_beta_diff_CIs(B)
```

For a historical full-data multiple-changepoint fit with array-valued `Beta`:

```r
plot_changepoints(old_fit$Kappa, S = X, dates = dates,
                  kappa_mode = old_fit$Kappa_mode)
plot_beta_diff_CIs(old_fit$Beta, coef_names = colnames(X))
# Valid only if the historical fit standardized ALL rows of X with scale():
plot_log_prob_diff(old_fit, X = scale(X))
```

For a historical holdout fit, reconstruct the training-only centering and
scaling instead. Do not pass raw statistics to `plot_log_prob_diff(X = ...)`.
Marginal modes need not form a valid joint partition. The log-difference plot
requires ordered boundaries and lets you supply a representative joint draw
explicitly via `kappa_mode` if necessary.


## Citation

Wang, Y., Thomas, A. M., and Jauch, M. (2026). *A generalized Bayesian approach
to multiple changepoint analysis*. Manuscript dated March 12, 2026.

See `citation("bcmlr")`. Bibliographic publication details should be updated
when available. The upstream repository had no license at the audited revision;
a release license and data redistribution permissions remain to be confirmed.


### Installation

If `remotes::install_github()` fails with **HTTP error 401: Bad credentials**,
R may be using an invalid or expired GitHub personal access token saved on
your computer. Try either solution below.

#### Option 1: Replace the saved token

[Generate a new token in your own GitHub account](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/managing-your-personal-access-tokens),
then run:

```r
install.packages("gitcreds")  # Only if not already installed
gitcreds::gitcreds_set()
```

Choose **Replace these credentials** and enter the new token when prompted.
Then retry installation:

```r
remotes::install_github(
  "YuhuiChloe/bcmlr",
  dependencies = NA,
  upgrade = "never"
)
```

#### Option 2: Install from the public ZIP archive

This method avoids using the saved GitHub credential:

```r
remotes::install_url(
  "https://github.com/YuhuiChloe/bcmlr/archive/refs/heads/main.zip",
  dependencies = NA,
  upgrade = "never"
)
```

GitHub generates this ZIP automatically from the `main` branch. R downloads
and extracts it into temporary storage, then installs the package.
No manual ZIP download is required.

After either installation succeeds, load the package:

```r
library(bcmlr)
```

To check the installed package location and the current R session's
temporary directory:

```r
# Installed package location
find.package("bcmlr")

# Temporary download directory
tempdir()
```

Temporary installation files may already have been removed. The installed
package remains available through `library(bcmlr)`.
