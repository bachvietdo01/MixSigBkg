mvking <- setRefClass("mvking",
                    fields = list(
                      mu = "numeric",
                      Sigma = "matrix",
                      eta = "numeric",
                      d = "numeric"
                    ),
                    methods = list(
                      initialize = function(mu, Sigma, eta = 1.5){
                        d <<- length(mu)
                        
                        testthat::expect_true(d > 1)
                        testthat::expect_true(eta > 0)
                        testthat::expect_true(d == dim(Sigma)[2])
                        testthat::expect_true(matrixcalc::is.positive.definite(Sigma))
                        
                        mu <<- mu
                        Sigma <<- Sigma
                        eta <<- eta
                      },
                      rking = function(n) {
                        samp <- matrix(rep(0, d * n), ncol = d)
                        
                        for(i in 1:n) {
                          Psi <-  CholWishart::rInvWishart(1, df = 2*eta - 1, 
                                                           Sigma = Sigma)[,,1]
                          samp[i, ]<- mvtnorm::rmvnorm(n = 1, mean =  mu, 
                                                       sigma = Psi) 
                        }
                        
                        return(samp)
                      },
                      dking = function(x) {
                        testthat::expect_true(length(x) == d)
                        
                        pi <- 3.14159265358979323846
                        Gamma <- solve(Sigma)
                        
                        # compute 1 + trace(.)
                        c <- 1 + sum(diag(Gamma %*% (x - mu) %*% t(x - mu)))
                        c <- (c^eta) * det(Gamma)^0.5 * pi^(0.5 * d) 
                        c <- c * CholWishart::mvgamma(eta - 0.5, p = 2)
                        c <- 1. / c * CholWishart::mvgamma(eta, p = 2)
                        
                        return(c)
                      }
                    ))