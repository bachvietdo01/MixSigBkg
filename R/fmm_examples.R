library(ggplot2)
setwd('~/Documents/MixSigBkg/R/')

Rcpp::sourceCpp('./helpers.cpp')
source('_component.R')
source('_mixture.R')
source('_fmmsignal.R')
source('helpers.R')
source('visualisation.R')

data <- data_gmm

# plot data
ggplot() + geom_point(data = data.frame(x = data$X[,1], y = data$X[,2]), aes(x, y)) 

# prepare to run Gibbs
X = as.matrix(data$X)
Ks = 2
Kb = 1
s = rep(1, nrow(X))
s = data_kmm$s

sz = kmeans(x = X, Ks)$cluster
sz = data_kmm$z

bz = rep(0, nrow(X))

# kmm_init_pars = get_KingComponent_example_init_pars(2)
# kmm_init_pars$Gamma0 = kmm_init_pars$Gamma0 * 900
# 
# m = Mixture(Ks = Ks, Kb = Kb, D = 2, X = X, s = s, sz = sz, bz = bz,
#             sig_componentClass = KingComponent, sig_init_pars = kmm_init_pars)

gmm_init_pars = get_GMVNComponentexample_init_pars(2)
m = Mixture(Ks = Ks, Kb = Kb, D = 2, X = X, s = s, sz = sz, bz = bz,
            sig_componentClass = GMVNComponent, sig_init_pars = gmm_init_pars)

# run Gibbs
niters = 20
for(i in 1:niters) {
  if (i == 1 || i %% 5 == 0 ) {
    print(paste(c('iter:', i)))
    print(plot_2D_MM_signal(m))
    for(k in 1:Ks) {
      print(paste(c('Signal cluster k size:', m$signal[[1]]$components[[k]]$N)))
    }
  }
  
  m$gibbs()
}


ggplot() + 