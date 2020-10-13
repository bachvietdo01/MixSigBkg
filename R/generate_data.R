source('dists.R')

gen_mix_gmm_and_unif <- function(N, K, Umin, Umax, mu, Sigma, rho, pi) {
  d = length(Umin)
  
  X = matrix(rep(0, d* N), ncol = d)
  s = rep(0, N)
  z = rep(0, N)
  
  for (n in 1: N) {
    sn = rbinom(1, 1, rho) # sample s
    s[n] = sn
    
    if (sn == 0) { # sample background from a uniform distribution
      for (p in 1 :d) {
        X[n,p] = runif(1, min = Umin[p], max = Umax[p])
      } 
    } else { # sample signal from a mixture of k multivariate normal
      zn = sample(1:K, 1, prob = pi)
      z[n] = zn
      # Get mu and Sigma for kth normal component
      muk = mu[, zn]
      Sigmak = Sigma[, , zn]
      x = rep(0, d)  
      while(TRUE) {
        ofb = FALSE
        x = mvtnorm::rmvnorm(1, mean = muk, sigma = Sigmak)
        ofb = any((x < Umin) | (x > Umax)) #outside of allowed region
        if(! ofb) {
          X[n,] = x
          break
        }
      }
    }
  }
  
  return(list(X = X, s = s, z = z))
}  

gen_mix_king_and_unif <- function(N, K, Umin, Umax, mu, Sigma, eta, rho, pi) {
  d = length(Umin)
  
  X = matrix(rep(0, d* N), ncol = d)
  s = rep(0, N)
  z = rep(0, N)
  
  kings = list()
  
  for(k in 1:K) {
    muk = mu[, k]
    Sigmak = Sigma[, , k]
    kings[[k]] = mvking$new(muk, Sigmak, eta)
  }
  
  for (n in 1: N) {
    sn = rbinom(1, 1, rho) # sample s
    s[n] = sn
    
    if (sn == 0) { # sample background from a uniform distribution
      for (p in 1 :d) {
        X[n,p] = runif(1, min = Umin[p], max = Umax[p])
      } 
    } else { # sample signal from a mixture of k multivariate normal
      zn = sample(1:K, 1, prob = pi)
      z[n] = zn
      # Get mu and Sigma for kth normal component
      x = rep(0, d)  
      while(TRUE) {
        ofb = FALSE
        x = kings[[zn]]$rking(1)
        ofb = any((x < Umin) | (x > Umax)) #outside of allowed region
        if(! ofb) {
          X[n,] = x
          break
        }
      }
    }
  }
  
  return(list(X = X, s = s, z = z))
}  


# Construct initial parameters
N = 1000;
d = 2;
K = 2;

# Construct mu and Sigma
I = diag(d);
mu0 = c(-10, 10);

mu = matrix(rep(0, d * K), nrow = 2)
mu[,1] = c(-20, 0)
mu[,2] = c(20, 0)

sigma.2 = 3; # know varaince
v = c();
for (k in 1: K) {
  v = c(v, as.vector(sigma.2 * I))
}
Sigma = array(v, c(d, d, K))

# Construct Umin and Umax
Umin = c()
Umax = c()
for (p in 1:d) {
  Umax = c(Umax, max(mu[p,]) + 3 * sigma.2) 
  Umin = c(Umin, min(mu[p,]) - 3 * sigma.2);
}

# Set rho and pi
rho = 0.75;
pi = c(0.65, 0.35)

data_gmm = list()
output_gmm = gen_mix_gmm_and_unif(N, K, Umin, Umax, mu, Sigma, rho, pi)
data_gmm$X = output_gmm$X
data_gmm$s = output_gmm$s
data_gmm$z = output_gmm$z
data_gmm$mu = mu
data_gmm$rho = rho
data_gmm$pi = pi
data_gmm$Sigma = Sigma
data_gmm$lower_bg = Umin
data_gmm$upper_bg = Umax
data_gmm$N = N
data_gmm$K = K
data_gmm$d = d

ggplot() + geom_point(data = data.frame(data_gmm$X), aes(X1, X2))

data_kmm = data_gmm
eta = 1.5
data_kmm$eta = eta

output_kmm = gen_mix_king_and_unif(N, K, Umin, Umax, mu, Sigma, eta, rho, pi)

data_kmm$X = output_kmm$X
data_kmm$s = output_kmm$s
data_kmm$z = output_kmm$z

ggplot() + geom_point(data = data.frame(data_kmm$X), aes(X1, X2))


rm(list=ls()[! ls() %in% c("data_gmm","data_kmm")])