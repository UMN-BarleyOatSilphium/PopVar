
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PopVar

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/PopVar)](https://cran.r-project.org/package=PopVar)

<!-- badges: end -->

The main attribute of ‘PopVar’ is the prediction of genetic variance in
bi-parental populations, from which the package derives its name.
‘PopVar’ contains a set of functions that use phenotypic and genotypic
data from a set of candidate parents to 1) predict the mean, genetic
variance, and superior progeny value of all, or a defined set of
pairwise bi-parental crosses, and 2) perform cross-validation to
estimate genome-wide prediction accuracy of multiple statistical models.
More details are available in Mohammadi, Tiede, and Smith (2015). Crop
Sci. <doi:10.2135/cropsci2015.01.0030>. A dataset ‘think\_barley.rda’ is
included for reference and examples.

## Installation

You can install the released version of PopVar from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("PopVar")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("UMN-BarleyOatSilphium/PopVar")
```

## Examples

Below are some example uses of the functions in `PopVar`:

``` r
# Load the package
library(PopVar)

# Load the example data
data("think_barley", package = "PopVar")
```

### Predictions using simulated populations

The code below simulates a single population of 1000 individuals for
each of 150 crosses. For the sake of speed, the marker effects are
predicted using RR-BLUP and no cross-validation is performed.

``` r
out <- pop.predict(G.in = G.in_ex, y.in = y.in_ex, map.in = map.in_ex,
                   crossing.table = cross.tab_ex,
                   nInd = 1000, nSim = 1, 
                   nCV.iter = 1, models = "rrBLUP")
```

The function returns a list, one element of which is called
`predictions.` This element is itself a list of matrices containing the
predictions for each trait. They can be combined as such:

``` r
predictions1 <- lapply(X = out$predictions, FUN = function(x) {
  x1 <- as.data.frame(apply(X = x, MARGIN = 2, FUN = unlist), stringsAsFactors = FALSE)
  cbind(x1[,c("Par1", "Par2")], sapply(X = x1[,-1:-2], as.numeric)) 
})

# Display the first few lines of the predictions for grain yield
knitr::kable(head(predictions1$Yield_param.df))
```

| Par1     | Par2      | midPar.Pheno | midPar.GEBV |   pred.mu | pred.mu\_sd | pred.varG | pred.varG\_sd | mu.sp\_low | mu.sp\_high | low.resp\_FHB | low.resp\_DON | low.resp\_Height | high.resp\_FHB | high.resp\_DON | high.resp\_Height | cor\_w/\_FHB | cor\_w/\_DON | cor\_w/\_Height |
| :------- | :-------- | -----------: | ----------: | --------: | ----------: | --------: | ------------: | ---------: | ----------: | ------------: | ------------: | ---------------: | -------------: | -------------: | ----------------: | -----------: | -----------: | --------------: |
| FEG66-08 | MN97-125  |       99.200 |   100.65607 | 100.70196 |          NA |  7.370914 |            NA |   96.02586 |    105.4725 |      23.06119 |      23.58167 |         79.00802 |       25.21990 |       25.14274 |          75.73304 |    0.3023413 |    0.2889576 |     \-0.3522443 |
| MN99-71  | FEG90-31  |          NaN |    99.17635 |  99.19179 |          NA |  5.060506 |            NA |   95.29299 |    103.0734 |      21.57387 |      23.69024 |         77.78220 |       24.81191 |       24.82185 |          75.63849 |    0.4807210 |    0.1349519 |     \-0.1822760 |
| MN96-141 | FEG183-52 |          NaN |   101.28761 | 101.19122 |          NA |  5.561630 |            NA |   97.17308 |    105.1518 |      23.51747 |      23.75884 |         80.00274 |       27.40918 |       28.47112 |          75.76285 |    0.7308630 |    0.7395142 |     \-0.5590697 |
| MN99-02  | FEG183-52 |          NaN |   100.78451 | 100.69620 |          NA | 10.623268 |            NA |   95.11679 |    106.3322 |      25.50111 |      24.06731 |         75.22959 |       26.25233 |       26.71408 |          73.74650 |    0.1343738 |    0.4162794 |     \-0.1575685 |
| FEG99-10 | FEG148-56 |          NaN |    98.68505 |  98.67753 |          NA |  4.172511 |            NA |   95.12261 |    102.1896 |      19.30170 |      19.94683 |         85.95189 |       22.74299 |       24.00809 |          80.30345 |    0.5882959 |    0.6392855 |     \-0.5836813 |
| MN99-62  | MN01-46   |      105.775 |   102.29724 | 102.41770 |          NA |  5.030651 |            NA |   98.59222 |    106.0020 |      26.07156 |      28.33927 |         73.21921 |       26.31851 |       28.37572 |          74.14646 |    0.1360669 |  \-0.1075501 |       0.3622568 |

### Predictions using deterministic equations

Generating predictions via simulated populations can become
computationally burdensome when many thousands or hundreds of thousands
of crosses are possible. Fortunately, deterministic equations are
available to generate equivalent predictions in a fraction of the time.
These equations are provided in the `pop.predict2` and `pop_predict2`
functions.

The `pop.predict2` function takes arguments in the same format as
`pop.predict`. We have eliminated the arguments for marker filtering and
imputation and cross-validation, as the `pop.predict2` function does not
support these steps. (You may continue to conduct cross-validation using
the `x.val` function.) Therefore, the genotype data input for
`pop.predict2` **must not contain any missing data**. Further, these
predictions assume fully inbred parents, so marker genotypes must only
be coded as -1 or 1. The data `G.in_ex_imputed` contains genotype data
that is formatted properly.

Below is an example of using the `pop.predict2`
function:

``` r
out2 <- pop.predict2(G.in = G.in_ex_imputed, y.in = y.in_ex, map.in = map.in_ex,
                     crossing.table = cross.tab_ex, model = "rrBLUP")

knitr::kable(head(subset(out2, trait == "Yield")))
```

|    | parent1  | parent2   | trait |  pred\_mu | pred\_varG | pred\_musp\_low | pred\_musp\_high | cor\_w\_FHB | cor\_w\_DON | cor\_w\_Yield | cor\_w\_Height | pred\_cor\_musp\_low\_FHB | pred\_cor\_musp\_low\_DON | pred\_cor\_musp\_low\_Yield | pred\_cor\_musp\_low\_Height | pred\_cor\_musp\_high\_FHB | pred\_cor\_musp\_high\_DON | pred\_cor\_musp\_high\_Yield | pred\_cor\_musp\_high\_Height |
| -- | :------- | :-------- | :---- | --------: | ---------: | --------------: | ---------------: | ----------: | ----------: | ------------: | -------------: | ------------------------: | ------------------------: | --------------------------: | ---------------------------: | -------------------------: | -------------------------: | ---------------------------: | ----------------------------: |
| 3  | FEG66-08 | MN97-125  | Yield | 100.41166 |   7.768361 |        95.52021 |         105.3031 |   0.3412831 |   0.3303378 |            NA |    \-0.4004154 |                  22.27035 |                  23.76433 |                          NA |                     78.84012 |                   24.96750 |                   25.91641 |                           NA |                      75.11838 |
| 7  | MN99-71  | FEG90-31  | Yield |  98.98546 |   5.298104 |        94.94591 |         103.0250 |   0.4959745 |   0.2776666 |            NA |    \-0.1612311 |                  21.86563 |                  23.41697 |                          NA |                     77.04669 |                   24.85938 |                   25.18159 |                           NA |                      75.75153 |
| 11 | MN96-141 | FEG183-52 | Yield | 101.07137 |   5.571654 |        96.92884 |         105.2139 |   0.6132612 |   0.6961326 |            NA |    \-0.4603399 |                  24.36629 |                  23.67733 |                          NA |                     79.59345 |                   27.88028 |                   28.40806 |                           NA |                      75.16009 |
| 15 | MN99-02  | FEG183-52 | Yield | 100.82967 |   9.499736 |        95.42053 |         106.2388 |   0.0249962 |   0.4211875 |            NA |    \-0.1929275 |                  26.15174 |                  23.73807 |                          NA |                     74.76672 |                   26.29180 |                   26.60344 |                           NA |                      72.86270 |
| 19 | FEG99-10 | FEG148-56 | Yield |  98.01593 |   4.121802 |        94.45292 |         101.5789 |   0.5959335 |   0.6401323 |            NA |    \-0.5080571 |                  19.67552 |                  19.58928 |                          NA |                     85.81854 |                   23.38585 |                   24.23917 |                           NA |                      80.69997 |
| 23 | MN99-62  | MN01-46   | Yield | 102.51094 |   4.673196 |        98.71709 |         106.3048 |   0.0755600 | \-0.0827658 |            NA |      0.4010907 |                  26.07342 |                  28.37124 |                          NA |                     72.67125 |                   26.20580 |                   28.16102 |                           NA |                      74.51340 |

Note that the output of `pop.predict2` is no longer a list, but a data
frame containing the combined predictions for all traits.

The formatting requirements of `G.in` for `pop.predict` and
`pop.predict2` are admittedly confusing. Marker genotype data is
commonly stored as a *n* x *p* matrix, where *n* is the number of
entries and *p* the number of markers. The function `pop_predict2`
accommodates this general marker data storage. Here is an
example:

``` r
out3 <- pop_predict2(M = G.in_ex_mat, y.in = y.in_ex, map.in = map.in_ex,
                     crossing.table = cross.tab_ex, model = "rrBLUP")

knitr::kable(head(subset(out2, trait == "Yield")))
```

|    | parent1  | parent2   | trait |  pred\_mu | pred\_varG | pred\_musp\_low | pred\_musp\_high | cor\_w\_FHB | cor\_w\_DON | cor\_w\_Yield | cor\_w\_Height | pred\_cor\_musp\_low\_FHB | pred\_cor\_musp\_low\_DON | pred\_cor\_musp\_low\_Yield | pred\_cor\_musp\_low\_Height | pred\_cor\_musp\_high\_FHB | pred\_cor\_musp\_high\_DON | pred\_cor\_musp\_high\_Yield | pred\_cor\_musp\_high\_Height |
| -- | :------- | :-------- | :---- | --------: | ---------: | --------------: | ---------------: | ----------: | ----------: | ------------: | -------------: | ------------------------: | ------------------------: | --------------------------: | ---------------------------: | -------------------------: | -------------------------: | ---------------------------: | ----------------------------: |
| 3  | FEG66-08 | MN97-125  | Yield | 100.41166 |   7.768361 |        95.52021 |         105.3031 |   0.3412831 |   0.3303378 |            NA |    \-0.4004154 |                  22.27035 |                  23.76433 |                          NA |                     78.84012 |                   24.96750 |                   25.91641 |                           NA |                      75.11838 |
| 7  | MN99-71  | FEG90-31  | Yield |  98.98546 |   5.298104 |        94.94591 |         103.0250 |   0.4959745 |   0.2776666 |            NA |    \-0.1612311 |                  21.86563 |                  23.41697 |                          NA |                     77.04669 |                   24.85938 |                   25.18159 |                           NA |                      75.75153 |
| 11 | MN96-141 | FEG183-52 | Yield | 101.07137 |   5.571654 |        96.92884 |         105.2139 |   0.6132612 |   0.6961326 |            NA |    \-0.4603399 |                  24.36629 |                  23.67733 |                          NA |                     79.59345 |                   27.88028 |                   28.40806 |                           NA |                      75.16009 |
| 15 | MN99-02  | FEG183-52 | Yield | 100.82967 |   9.499736 |        95.42053 |         106.2388 |   0.0249962 |   0.4211875 |            NA |    \-0.1929275 |                  26.15174 |                  23.73807 |                          NA |                     74.76672 |                   26.29180 |                   26.60344 |                           NA |                      72.86270 |
| 19 | FEG99-10 | FEG148-56 | Yield |  98.01593 |   4.121802 |        94.45292 |         101.5789 |   0.5959335 |   0.6401323 |            NA |    \-0.5080571 |                  19.67552 |                  19.58928 |                          NA |                     85.81854 |                   23.38585 |                   24.23917 |                           NA |                      80.69997 |
| 23 | MN99-62  | MN01-46   | Yield | 102.51094 |   4.673196 |        98.71709 |         106.3048 |   0.0755600 | \-0.0827658 |            NA |      0.4010907 |                  26.07342 |                  28.37124 |                          NA |                     72.67125 |                   26.20580 |                   28.16102 |                           NA |                      74.51340 |

#### Benchmarking and comparisons

The code below compares the functions `pop.predict` and `pop.predict2`
with respect to computation time and results:

``` r
time1 <- system.time({
  capture.output(pop.predict.out <- pop.predict(
    G.in = G.in_ex_imputed, y.in = y.in_ex, map.in = map.in_ex, crossing.table = cross.tab_ex,
    nInd = 1000, nSim = 1, nCV.iter = 1, models = "rrBLUP"))
})

time2 <- system.time({pop.predict2.out <- pop.predict2(
  G.in = G.in_ex_imputed, y.in = y.in_ex, map.in = map.in_ex,
  crossing.table = cross.tab_ex,model = "rrBLUP")})

c(pop.predict = time1[[3]], pop.predict2 = time2[[3]])
##  pop.predict pop.predict2 
##        27.28         0.64
```

Plot
results

``` r
predictions1 <- lapply(X = pop.predict.out$predictions, FUN = function(x) {
  x1 <- as.data.frame(apply(X = x, MARGIN = 2, FUN = unlist), stringsAsFactors = FALSE)
  cbind(x1[,c("Par1", "Par2")], sapply(X = x1[,-1:-2], as.numeric))
})

pop.predict.out1 <- predictions1$Yield_param.df[,c("Par1", "Par2", "pred.varG")]
pop.predict2.out1 <- subset(pop.predict2.out, trait == "Yield", c(parent1, parent2, pred_varG))

toplot <- merge(pop.predict.out1, pop.predict2.out1, by.x = c("Par1", "Par2"),
                by.y = c("parent1", "parent2"))

plot(pred.varG ~ pred_varG, toplot,
     xlab = "pop.predict2", ylab = "pop.predict",
     main = "Comparsion of predicted genetic variance")
```

<img src="man/figures/README-compare2-1.png" width="100%" />

### Multi-parent populations
