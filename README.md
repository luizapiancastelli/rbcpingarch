
`rbcpingarch` is an R package providing code to simulate from and fit
the Bivariate Conditional Poisson INGARCH(1,1) process proposed in
[Piancastelli, Barreto-Souza and
Ombao](https://arxiv.org/abs/2011.08799). Our package links to `rstan`,
so proper setup of C++ compiler is in order for usage. We recommend
testing this via installation of `Rcpp` and `RcppArmadillo`. Once this
is done, the following command line installs `rbcpingarch` from its
GitHub repository.

``` r
# install.packages("devtools")
devtools::install_github("luizapiancastelli/rbcpingarch")
```

### BCP-INGARCH simulation

The next chunk of code exemplifies using the `rBCPINGARCH` for
simulating a realisation of the BCP-INGARCH(1,1) process. This function
takes as input the model parameters with `A` being a (2x2) matrix of
parameters related to the previous conditional means, `B` a (2x2) matrix
of INGARCH parameters on the previous counts, `omega` the intercept
vector, `phi` the correlation parameter and `n`, and the trajectory
lenght. A brief description of these arguments is included in a call to
`help(rBCPINGARCH)`, and a full model specification is in Definition 2.1
of the paper’s preprint linked above.

The output of `rBCPINGARCH` is a (`n`x 2) integer matrix, with one
realisation illustrated in the code below.

``` r
A = diag(c(0.3, 0.1)) #Diagonal (2x2) matrix
B = diag(c(0.2, 0.3)) #Diagonal (2x2) matrix, can be non diagonal as long as stationarity conditions satisfied
omega = c(2,1) 
phi = 0.3 
n = 500  #Time series length

Y = rBCPINGARCH(A, B, omega, phi, n)
head(Y)
#>      [,1] [,2]
#> [1,]    5    1
#> [2,]    2    1
#> [3,]    3    0
#> [4,]    4    0
#> [5,]    1    0
#> [6,]    4    4
cor(Y)
#>           [,1]      [,2]
#> [1,] 1.0000000 0.6239631
#> [2,] 0.6239631 1.0000000

Y_df= data.frame('y1' =Y[,1], 'y2' = Y[,2], 't' = 1:nrow(Y))
Y_df = pivot_longer(Y_df, cols = starts_with('y'), names_to = 'serie', values_to = 'count')

ggplot(Y_df, aes(x = t, y = count, color = serie))+
  geom_line()+theme_bw()+
  labs(y = 'Count', x = 'Time')
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

### Model fit

The bivariate time series of confirmed viral hepatites counts in the
nearby Brazilian cities of Goiania and Brasilia explored in
[Piancastelli, Barreto-Souza and
Ombao](https://arxiv.org/abs/2011.08799) is originally collected from
the open [DATASUS platform](https://datasus.saude.gov.br), included in
`rBCPINGARCH` under `data("Y_hep")`. As described in Section 6 of the
manuscript, the study period of this version of the data is from 2001 to
2018.

This data is used in the code chunk below to provide guidelines on using
the `fit_BCP_INGARCH` to fit a BCP-INGARCH(1,1) process via conditional
maximum likelihood estimation. Its arguments are `Y` the observed
bivariate count data, and boolean values `A.diag`, `B.diag` specifying
the configuration of the `A` and `B` matrices as diagonal (`TRUE`) or
non-diagonal (`FALSE`). The diagonal and non-diagonal model fits
analysed in the original manuscript are reproduced in the code below.

``` r
data("Y_hep") #Loads the Brazilian hepatites count data set

fit_diag = fit_BCP_INGARCH(Y_hep, A.diag = TRUE, B.diag = TRUE) #Diagonal model
fit_diag$par #maximised parameter values
#>        A[1]        A[2]        B[1]        B[2]  omega[1,1]  omega[2,1] 
#> 0.465619498 0.481809019 0.429948649 0.383853540 2.310480454 6.518914405 
#>         phi 
#> 0.009613826
fit_diag$se  #asymptotic standard errors
#>         A11         A22         B11         B22      omega1      omega2 
#> 0.056709886 0.039787482 0.042086278 0.026293151 0.578115241 0.909545488 
#>         phi 
#> 0.001149342
fit_diag$loglik 
#> [1] 46115.42

fit_nd = fit_BCP_INGARCH(Y_hep, A.diag = TRUE, B.diag = FALSE) #Non-diagonal model
fit_nd$par
#>        A[1]        A[2]      B[1,1]      B[2,1]      B[1,2]      B[2,2] 
#> 0.496876989 0.447653650 0.398157521 0.018631135 0.001979496 0.399368538 
#>  omega[1,1]  omega[2,1]         phi 
#> 2.216294309 7.011972425 0.009599605
fit_nd$se
#>         A11         A22         B11         B21         B12         B22 
#> 0.059010457 0.043599866 0.046492743 0.031799857 0.007077866 0.029602968 
#>      omega1      omega2         phi 
#> 0.612004557 1.133671734 0.001154630
fit_nd$loglik
#> [1] 46117.22
```

The output is a named list containing the estimated model parameters
`fit$par`, their asymptotic standard errors `fit$se`, maximed
log-likelihood value `fit$loglik`, the fitted process `fit$lambda`, and
the initial values used by the algorithm `fit$initial_values`, among
others.

Selecting between the nested diagonal and non-diagonal model fits can be
done via information criterion with the `model_information` function.
This returns the AIC and BIC of the alternative models, both indicating
the diagonal case for the hepatites count data.

``` r
model_information(fit_diag)
#>       AIC       BIC 
#> -92216.83 -92193.21
model_information(fit_nd)
#>       AIC       BIC 
#> -92216.44 -92186.06
```
