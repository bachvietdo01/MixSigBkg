Mixture <- setRefClass("Mixture",
                       fields = list(
                         background = "list",
                         signal = "list",
                         Ks = "integer",
                         Kb = "integer",
                         D = "integer",
                         X = "matrix",
                         N_total = "integer",
                         s = "numeric",
                         sz = "numeric",
                         bz = "numeric",
                         rho = "numeric",
                         sig_kappa0 = "numeric",
                         sig_nu0 = "numeric",
                         sig_m0 = "numeric",
                         sig_S0 = "matrix",
                         sig_alpha = "numeric",
                         alpha0 = "numeric",
                         bg_m0 = "numeric",
                         bg_s0.2 = "numeric",
                         bg_nu0 = "numeric",
                         bg_kappa0 = "numeric",
                         bg_alpha = "numeric",
                         lower_bg = "numeric",
                         upper_bg = "numeric"),

                       methods = list(
                         initialize = function(Ks, Kb, D, X, sz, bz, s, 
                                               sig_kappa0 = 1, 
                                               sig_nu0 = D+2, 
                                               sig_m0 = rep(0, D), 
                                               sig_S0 = diag(rep(1, D)), 
                                               sig_alpha = 1.0 / Ks,
                                               alpha0 = 0.5,
                                               bg_m0 = rep(0, D),
                                               bg_s0.2 = 1.0,
                                               bg_nu0 = D + 2,
                                               bg_kappa0 = 1,
                                               bg_alpha = 1.0 / Kb,
                                               type = 'FMM'){
                           
                           Ks <<- as.integer(Ks)
                           Kb <<- as.integer(Kb)
                           D <<- as.integer(D)
                           N_total <<- nrow(X)
                           s <<- s
                           bz <<- bz
                           sz <<- sz
                           
                           X <<- X
                           alpha0 <<- alpha0
                           sig_alpha <<- sig_alpha

                           sig_kappa0 <<- sig_kappa0
                           sig_nu0 <<- sig_nu0
                           sig_m0 <<- sig_m0
                           sig_S0 <<- sig_S0
                           
                           bg_alpha <<- bg_alpha
                           bg_kappa0 <<- bg_kappa0
                           bg_nu0 <<- bg_nu0
                           bg_m0 <<- bg_m0
                           bg_s0.2 <<- bg_s0.2
                           
                           rho <<- rbeta(1, alpha0, alpha0)
                           
                           sig_init_pars <- list(D = D, kappa0 = sig_kappa0, nu0 = sig_nu0, m0 = sig_m0, S0 = sig_S0)
                           
                           if(type == 'DPSignal') {
                             signal[[1]] <<- DPSignal$new(K = Ks, D = D, X = X, s = s, 
                                                      z = sz, kappa0 = sig_kappa0, alpha = sig_alpha,
                                                      nu0 = sig_nu0, m0 = sig_m0, S0 = sig_S0)
                           } else(type == 'MFMSignal') {
                             signal[[1]] <<- MFMSignal$new(K = Ks, D = D, X = X, s = s, 
                                                          z = sz, kappa0 = sig_kappa0, alpha = sig_alpha,
                                                          nu0 = sig_nu0, m0 = sig_m0, S0 = sig_S0)
                           } else {
                             sig_init_pars <- list(D = D, kappa0 = sig_kappa0, nu0 = sig_nu0, 
                                                   m0 = sig_m0, s0.2 = sig_s0.2)
                             signal[[1]] <<- FMMBackground$new(K = Kb, X = X, s = s, z = bz,
                                                                   alpha = bg_alpha,
                                                                   ComponentClass = IMVNComponent,
                                                                   init_pars = bg_init_pars)
                           }
                           
                           lower_bg <<- apply(X, 2, min)
                           upper_bg <<- apply(X, 2, max)
                           
                         },
                         update_rho = function() {
                           N1 <- sum(s == 1)
                           rho <<- rbeta(1, alpha0 + N1, alpha0 + N_total - N1)
                         },
                        collapsed_gibbs = function(){
                           for(i in 1:N_total){
                             s_probs <- rep(0, 2)
                             #_probs[1] <- log(rho) + signal[[1]]$get_loglik(X[i,])
                             s_probs[2] <- log(1 - rho) - log(prod(upper_bg - lower_bg))
                             
                             if(! is.na(s_probs[1])) {
                               s_probs = softmax(s_probs)
                             } else {
                               s_probs= c(0,1)
                             }
                             s[i] <<- rbinom(1, 1, s_probs[1])
                           }
                          
                          #update_rho()
                          rho <<- 0.75

                          signal[[1]]$set_s(s)
                          signal[[1]]$collapsed_gibbs()
                          if(type != 'FMM') {
                            signal[[1]]$merge_split()
                            signal[[1]]$tidy_up_clusters()
                          }

                          # sz <<- signal[[1]]$z
                        }
                       ))


