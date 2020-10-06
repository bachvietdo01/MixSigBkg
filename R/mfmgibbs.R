library(ggplot2)
setwd('~/Downloads/MFMproj/gmm_background/')

Rcpp::sourceCpp('./helpers.cpp')
source('_component.R')
source('_mixture.R')
source('_dpsignal.R')
source('_mfmsignal.R')
source('_dpsignal.R')
source('helpers.R')
source('visualisation.R')

spatial = data
data = NULL
data$X = spatial[1:1000,]

# plot data
ggplot() + geom_point(data = data.frame(data$X), aes(X, Y)) 

# prepare to run Gibbs
X = as.matrix(data$X)
Ks = 10
Kb = 2
s = rep(1, nrow(X))

sz = kmeans(x = X, Ks)$cluster
bz = rep(1, nrow(X))

m = Mixture(Ks = Ks, Kb = Kb, D = 2, X = X, s = s, sz = sz, bz = bz)

# run Gibbs
niters = 1000
for(i in 1:niters) {
  if (i == 1 || i %% 5 == 0 ) {
    print(paste(c('iter:', i)))
    print(paste(c('Number of clusters', m$signal[[1]]$K)))
    print(plot_2D_MM_signal(m))
    for(k in 1:Kb) {
      print(paste(c('Bg cluster k size:', m$signal[[1]]$components[[k]]$N)))
    }
  }
  
  m$collapsed_gibbs()
}
