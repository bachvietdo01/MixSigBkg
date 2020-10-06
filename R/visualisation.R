plot_2D_MM_signal = function(obj, z = NULL){
  X = obj$X
  if(is.null(z)){
    z = factor(obj$signal[[1]]$z)
  } else{
    z = factor(z)
  }
  df = data.frame(X1 = X[, 1], X2 = X[, 2], z = z)
  ggplot(df, aes(X1, X2, col = z)) +
    geom_point() +
    theme_bw() +
    theme(legend.position = "none")
}

plot_2D_MM_background = function(obj, z = NULL){
  X = obj$X
  if(is.null(z)){
    z = factor(obj$background[[1]]$z)
  } else{
    z = factor(z)
  }
  df = data.frame(X1 = X[, 1], X2 = X[, 2], z = z)
  ggplot(df, aes(X1, X2, col = z)) +
    geom_point() +
    theme_bw() +
    theme(legend.position = "none")
}

