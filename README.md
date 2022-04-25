# distinct: a method for differential analyses via hierarchical permutation tests

<img src="inst/extdata/distinct.png" width="200" align="right"/> 

`distinct` is a statistical method to perform differential testing between two or more groups of distributions; differential testing is performed via non-parametric permutation tests on the cumulative distribution functions (cdfs) of each sample.
`distinct` is a general and flexible tool: due to its fully non-parametric nature, which makes no assumptions on how the data was generated, it can be applied to a variety of datasets.
It is particularly suitable to perform differential state analyses on single cell data (i.e., differential analyses within sub-populations of cells), such as single cell RNA sequencing (scRNA-seq) and high-dimensional flow or mass cytometry (HDCyto) data.
The method also allows for nuisance covariates (such as batch effects).

>  Simone Tiberi, Helena L Crowell, Pantelis Samartsidis, Lukas M Weber, and Mark D Robinson (2020).
>
> distinct: a novel approach to differential distribution analyses.
>
> *bioRxiv* doi: [10.1101/2020.11.24.394213](https://doi.org/10.1101/2020.11.24.394213)

## Bioconductor installation 
`distinct` is available on [Bioconductor](https://bioconductor.org/packages/distinct) and can be installed with the command:
``` r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("distinct")
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
