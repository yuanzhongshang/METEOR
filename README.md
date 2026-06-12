# METEOR

METEOR (Multivariate Estimation Tool for Exposure-Outcomes with Robustness), is an R package for efficient statistical inference of multi-outcomes mendelian randomization analysis. METEOR utilizes a set of correlated SNPs, self-adaptively accounts for the sample structure of both exposure and outcomes, the uncertainty that these correlated SNPs may exhibit multiple pleiotropic effects. The term 'self-adaptive' represents that METEOR is able to automatically infer the sample structure and the probability that a SNP has specific pleiotropic effect from the data at hand. METEOR places the inference of the causal effects into a likelihood-framework and relies on a scalable sampling-based algorithm to obtain calibrated $p$-values.

# Installation

It is easy to install the development version of METEORPLE package using the 'devtools' package.

```{r eval=FALSE}
#install.packages("devtools")
library(devtools)
install_github("Liye222/METEOR")
```

# Usage

The main function in the package is METEOR, you can find the instructions by `'?METEOR'`.

```{r eval=FALSE}
library(METEOR)

?METEOR
```

# Quick Start

See [Tutorial](https://Liye222.github.io/METEOR/) for detailed documentation and examples.

# Development

This R package is developed by Liye Zhang.
