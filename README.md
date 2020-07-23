# HDdiffindiff
Doubly Robust Semiparametric Difference-in-Differences Estimators with High Dimensional Data
If you use this package, please cite our article below:

[Ning, Peng and Tao (2020), "Doubly Robust Semiparametric Difference-in-Differences Estimators with High Dimensional Data"]()

To install this package
```R
library(devtools) 
install_version("flare","1.6.0")
install_github("psdsam/HDdiffindiff")
library(HDdiffindiff)
```

The functions in this package are:
highdimdiffindiff_crossfit: Main function to implement the Doubly Robust Semiparametric Difference-in-Differences Estimators with sample splitting
highdimdiffindiff_crossfit_inside: called by highdimdiffindiff_crossfit to compute and debias the estimators. 
Examplehighdimdiffindiff: Implement a simple example to demonstrate the method. 

Example usage:
```R
f = Examplehighdimdiffindiff(n0=1000)
```


