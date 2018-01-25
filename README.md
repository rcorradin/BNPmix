# BNPmix
Bayesian Nonparametric Mixture models - an efficient  based C++package to model density distributions in R.

## Prerequisite and how to install
The BNPmix package is based on RcppArmadillo tools, the GitHub version is a source pack, if not present and on some machine (i.e. Windows) it is necessary to install the 'Rtools34' package to compile.

```
install.packages('Rtools34')
```

To install the BNPmix package, if it is not present, it is necessary to install previously the 'devtools' package:

```
install.packages('devtools')
```

Finally, it is possible to get the BNPmix package installing by the source code on GitHub:

```
library(devtools)
install_github("rcorradin/BNPmix")
```

