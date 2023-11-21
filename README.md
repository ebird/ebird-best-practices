Best Practices for Using eBird Data

<!-- badges: start -->
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

This repository contains the source for the [Best Practices for Using eBird Data](https://ebird.github.io/ebird-best-practices/)
guide, which is a supplement to *Analytical guidelines to increase the value of community science data: An example using eBird data to estimate species distributions* ((Johnston et al. 2021)(https://onlinelibrary.wiley.com/doi/10.1111/ddi.13271)). 

## Installation

The R packages used in this guide can be installed with:

``` r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("ebird/ebird-best-practices")
```

Please cite this guide as:

> Strimas-Mackey, M., W.M. Hochachka, V. Ruiz-Gutierrez, O.J. Robinson, E.T. Miller, T. Auer, S. Kelling, D. Fink, A. Johnston. 2023. Best Practices for Using eBird Data. Version 2.0. https://ebird.github.io/ebird-best-practices/. Cornell Lab of Ornithology, Ithaca, New York.
