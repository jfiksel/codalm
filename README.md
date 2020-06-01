
# The codalm R package

The codalm R package implements the methods as described in [A Transformation-free Linear Regression for Compositional Outcomes and Predictors](https://arxiv.org/abs/2004.07881) (Fiksel, Zeger, and Datta, 2020)

## Installation

You can install the current version of codalm for github using the `remotes` package.

``` r
remotes::install_github("jfiksel/codalm")
```

## Example

We will start by loading the `codalm` R package, as well as the `ggtern` package,
which we will use to access the data.

```r
library(codalm)
library(ggtern)
```

We will now load in the data from the ggtern package. We will be analyzing how
two different methods (image analysis or microscopic inspection) estimate the composition
of 30 white blood cells. The format that we need is for both compositions to be in matrices,
with one row per observation. We will also normalize the rows of these matrices to ensure 
that they sum to 1, although the `codalm` function would also take care of this for us.

```r
data("WhiteCells", package = 'ggtern')
image <- subset(WhiteCells, Experiment == "ImageAnalysis")
image_mat <- as.matrix(image[,c("G", "L", "M")])
microscopic <- subset(WhiteCells, Experiment == "MicroscopicInspection")
microscopic_mat <- as.matrix(microscopic[,c("G", "L", "M")])

image_mat  <- image_mat  / rowSums(image_mat)
microscopic_mat <- microscopic_mat / rowSums(microscopic_mat)
```

To estimate the coefficient matrix B, we can use the `codalm` function.

```r
B_est <- codalm(yout = microscopic_mat, ypred = image_mat)
B_est
```

We can also use the bootstrap to estimate 95% confidence intervals. We will specify
the `accelerate` argument as TRUE to take advantage of the Squarem algorithm for bootstrapping. We will
only use 50 bootstrap iterations as an example (R = 50), but is recommended to do more.

```r
B_ci <- coda_lm_ci(yout = microscopic_mat, ypred = image_mat,
           accelerate = TRUE, R = 50, ci_type = "perc", conf = .95)
B_ci$ci_L
B_ci$ci_U
```

Finally, we will do a permutation test for linear independence. Again, we will only do 50
permutations as an example, but in practice this number should be higher.

```r
indep_test_pval <- codalm_indep_test(yout = microscopic_mat, ypred = image_mat, accelerate = TRUE,
                                     nperms = 50, init.seed = 123)
indep_test_pval
```

