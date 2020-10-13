library(ggplot2)
setwd('~/Documents/MixSigBkg/R/')

Rcpp::sourceCpp('./helpers.cpp')
source('_component.R')
source('_mixture.R')
source('_dpsignal.R')
source('_mfmsignal.R')
source('_dpsignal.R')
source('helpers.R')
source('visualisation.R')

data = data_gmm
# plot data
ggplot(data = data.frame(x = data$X[,1], y = data$X[,2])) + geom_point(aes(x, y)) 

# prepare to run Gibbs
X = as.matrix(data$X)
Ks = 10
Kb = 1
s = rep(1, nrow(X))

sz = kmeans(x = X, Ks)$cluster
bz = rep(0, nrow(X))

init_pars <- get_GMVNComponent_example_init_pars(2)

m = Mixture(Ks = Ks, Kb = Kb, D = 2, X = X, s = s, sz = sz, bz = bz,
            sig_componentClass = GMVNComponent,
            sig_init_pars = init_pars, type = 'MFMSignal')

# run Gibbs
niters = 50
for(i in 1:niters) {
  if (i == 1 || i %% 5 == 0 ) {
    print(paste(c('iter:', i)))
    print(paste(c('Number of clusters', m$signal[[1]]$K)))
    print(plot_2D_MM_signal(m))
    Ks <- m$signal[[1]]$K
    for(k in 1:Ks) {
      print(paste(c('Signal cluster k size:', m$signal[[1]]$components[[k]]$N)))
    }
  }
  
  m$gibbs()
}
