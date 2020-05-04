# distinct: a method for differential analyses via hierarchical permutation tests

`distinct` is a statistical method to perform differential testing between two or more groups of distributions; differential testing is performed via non-parametric permutation tests on the cumulative distribution functions (cdfs) of each sample.
`distinct` is a general and flexible tool: due to its fully non-parametric nature, which makes no assumptions on how the data was generated, it can be applied to a variety of datasets.
It is particularly suitable to perform differential state analyses on single cell data (i.e., differential analyses within sub-populations of cells), such as single cell RNA sequencing (scRNA-seq) and high-dimensional flow or mass cytometry (HDCyto) data.

At present, covariates are not allowed, and only 2-group comparisons are implemented.
In future releases (i.e., soon), we will allow for covariates and for differential testing between more than 2 groups.

A pre-print will follow in the coming months.

***

**`distinct` is still work in progress. Major updates are expected until summer.**

***

## Bioconductor installation 
`distinct` is available on [Bioconductor](https://bioconductor.org/packages/distinct) and can be installed with the command:
``` r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("distinct")
```

## Devel installation from github
To install the package from github, use `devtools` (available [here](https://github.com/hadley/devtools)):
``` r
devtools::install_github("SimoneTiberi/distinct")
```

To install the package jointly with its vignette remove `--no-build-vignettes` from `build_opts`:
``` r
devtools::install_github("SimoneTiberi/distinct", build_opts = c("--no-resave-data", "--no-manual"))
```

## Vignette
The vignette illustrating how to use the package can be accessed on the 
[Bioconductor website](https://www.bioconductor.org/packages/release/bioc/vignettes/distinct/inst/doc/distinct.pdf)
or from R via:
``` r
vignette("distinct")
```
or
``` r
browseVignettes("distinct")
```
