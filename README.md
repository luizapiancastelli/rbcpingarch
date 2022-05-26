
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rbcpingarch

<!-- badges: start -->

<!-- badges: end -->

The goal of rbcpingarch is to …

## Installation

You can install the released version of rbcpingarch from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("rbcpingarch")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("luizapiancastelli/rbcpingarch")
```

  - 
### Data simulation

``` r
A = matrix(c(0.3, 0, 0, 0.1), ncol =2, byrow = TRUE)
B = matrix(c(0.2, 0, 0, 0.3), ncol =2, byrow = TRUE)
omega = c(1,1)
phi = 0.3

Y = rBCPINGARCH(A, B, omega, phi, 500)
head(Y)
#>      [,1] [,2]
#> [1,]    1    1
#> [2,]    2    2
#> [3,]    1    2
#> [4,]    0    1
#> [5,]    2    1
#> [6,]    1    2
cor(Y)
#>          [,1]     [,2]
#> [1,] 1.000000 0.556566
#> [2,] 0.556566 1.000000

Y_df = data.frame('y1' =Y[,1], 'y2' = Y[,2], 't' = 1:nrow(Y))
ggplot(Y_df, aes(x = t))+
  geom_line(aes(y = y1), color = 'red') +  geom_line(aes(y = y2), color = 'blue')+theme_bw()+
  labs(y = 'Simulated count time series data', x = 'Time')
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

### Model fit

    #>        A[1]        A[2]        B[1]        B[2]  omega[1,1]  omega[2,1] 
    #> 0.465619498 0.481809019 0.429948649 0.383853540 2.310480454 6.518914405 
    #>         phi 
    #> 0.009613826
    #>         A11         A22         B11         B22      omega1      omega2 
    #> 0.056709886 0.039787482 0.042086278 0.026293151 0.578115241 0.909545488 
    #>         phi 
    #> 0.001149342
    #> [1] 46115.42
    #>        A[1]        A[2]      B[1,1]      B[2,1]      B[1,2]      B[2,2] 
    #> 0.496876989 0.447653650 0.398157521 0.018631135 0.001979496 0.399368538 
    #>  omega[1,1]  omega[2,1]         phi 
    #> 2.216294309 7.011972425 0.009599605
    #>         A11         A22         B11         B21         B12         B22 
    #> 0.059010457 0.043599866 0.046492743 0.031799857 0.007077866 0.029602968 
    #>      omega1      omega2         phi 
    #> 0.612004557 1.133671734 0.001154630
    #> [1] 46117.22
