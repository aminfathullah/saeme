# Package `saeme` (Small Area Estimation with Measurement Error)

This is an R package for small area estimation when covariates subject 
to sampling errors. 

This package is currently including function for:
* Fay-Herriot with Measurement Model<sup>[1]</sup>


Installation
----
1. Install `devtools` package

    ```r
    install.packages('devtools')
    ```
2. Install `saeme` package
    ```r
    devtools::install_github('aminfathullah/saeme')
    ```
3. Load `saeme` package
    ````r
    library(saeme)
    ````

Usage Example
----
Estimating mean body mass index

```r
data(bmi)

#calculate SAE with measurement error
sae_me <- FHme(y = bmi$y, x = cbind(1, bmi$x), vardir = bmi$mse_y, C = cbind(0, bmi$mse_x))

sae_me
summary(sae_me)
plot(sae_me)
```

References
----
[1] Ybarra, L. M., & Lohr, S. L. (2008). Small area estimation when auxiliary information is measured with error. Biometrika, 919-931.
