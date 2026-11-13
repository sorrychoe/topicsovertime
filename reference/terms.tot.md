# Most probable words per topic

Most probable words per topic

## Usage

``` r
# S3 method for class 'tot'
terms(x, n = 10L, ...)
```

## Arguments

- x:

  A fitted `"tot"` model.

- n:

  Number of words per topic.

- ...:

  Unused.

## Value

A character matrix with `n` rows and one column per topic.

## Examples

``` r
set.seed(1)
sim <- simulate_tot(n_docs = 100, n_topics = 3)
fit <- tot(sim$x, sim$time, k = 3, n_iter = 50)
terms(fit, 5)
#>      Topic1 Topic2 Topic3
#> [1,] "w100" "w6"   "w6"  
#> [2,] "w26"  "w26"  "w75" 
#> [3,] "w12"  "w97"  "w30" 
#> [4,] "w87"  "w33"  "w26" 
#> [5,] "w84"  "w51"  "w86" 
```
