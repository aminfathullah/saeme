library(expm)
library(MASS)

estimator_FHme <- function(y, X, psi, C, iter)
{
  m <- length(y)
  X <- as.matrix(X)
  p <- dim(X)[2]
  i <- 0

  #initialize the wi
  w <- rep(1,m)
  sigma2_cap <- 0

  while(i < iter)
  {
    # estimate beta
    b_cap <- beta_cap(y, X, C, w)

    # estimate sigma
    sigma2_cap <- sigma_cap(y, X, b_cap, psi, C)

    #new wi
    w <- sapply(1:m, function(i)
    {
      w_i <- 1/(sigma2_cap + psi[i] + t(b_cap)%*%diag(C[i,], nrow = p)%*%b_cap)
      return(w_i)
    })
    i <- i + 1
  }

  # estimate gamma
  g_cap <- sapply(1:m, function(i)
  {
    mse_ri <- sigma2_cap + t(b_cap)%*%diag(C[i,], nrow = p)%*%b_cap

    return(mse_ri/(mse_ri+psi[i]))
  })

  y_me <- sapply(1:m, function(i)
  {
    y_me_i <- g_cap[i]*y[i] + (1 - g_cap[i])*t(as.matrix(X[i,]))%*%b_cap
    return(y_me_i)
  })

  hasil <- list('y_me' = y_me, 'gamma' = g_cap, 'sigma' = sigma2_cap, 'beta' = b_cap, 'psi' = psi)

  return(hasil)
}

beta_cap <- function(y, X, C, w)
{
  m <- length(y)
  p <- dim(X)[2]

  ## compute wXy
  wXy <- Reduce('+', lapply(1:m, function(i)
  {
    X_i <- as.matrix(X[i,])
    wXy_i <- w[i]*X_i*y[i]
    return(wXy_i)
  }))
  ## compute wXX'
  wXX <- Reduce('+', lapply(1:m, function(i)
  {
    X_i <- as.matrix(X[i,])
    wXX_i <- w[i]*(X_i%*%t(X_i))
    return(wXX_i)
  }))
  ## compute wc
  wC <- Reduce('+', lapply(1:m, function(i)
  {
    C_i <- diag(C[i,], nrow = p)
    wC_i <- w[i]*C_i
    return(wC_i)
  }))
  b_cap <- 0
  if (try(class(solve(wXX - wC)), silent = T) == "matrix")
  {
    #if solvable
    b_cap <- solve(wXX - wC) %*% wXy
  }else
  {
    G <- wXX
    root_G <- sqrtm(G)

    inv_G <- 0
    if (sum(eigen(G)$value > 0) == p)
    {
      print('Matrix G is not solvable and it\'s positive definite (see reference)')
      inv_G <- solve(root_G)
    }else
    {
      print('Matrix G is not solvable and it\'s not positive definite (see reference)')
      inv_G <- ginv(root_G)
    }
    eigen_GwCG <- eigen(inv_G%*%wC%*%inv_G)
    orth_P <- eigen_GwCG$vector
    diag_lambda <- diag(eigen_GwCG$values, nrow = p)
    diag_D <- diag(sapply(1:p, function(j)
    {
      Djj <- 0
      if((1 - diag_lambda[j,j]) > 1/m )
        Djj <- 1/(1 - diag_lambda[j,j])

      return(Djj)
    }), nrow = p)

    b_cap <- inv_G%*%orth_P%*%diag_D%*%t(orth_P)%*%inv_G%*%wXy
  }
  return(b_cap)
}

sigma_cap <- function(y, X, b_cap, psi, C)
{
  m <- length(y)
  p <- dim(X)[2]

  first_part <- (m - p)^-1
  second_part <- sum(sapply(1:m, function(i)
  {
    X_i <- as.matrix(X[i,])
    tmp <-  (y[i]-t(X_i)%*%b_cap)^2 - psi[i] - t(b_cap)%*%diag(C[i,], nrow = p)%*%b_cap
    return(tmp)
  }))

  sigma2_cap <- first_part * second_part

  if(sigma2_cap < 0)
    sigma2_cap <- 0
  return(sigma2_cap)
}

jack_knife <- function(y, X, psi, C, j, iter)
{
  m <- length(y)
  X <- as.matrix(X)
  p <- dim(X)[2]
  i <- 0

  #initialize the wi
  w <- rep(1,m)
  sigma2_cap <- 0

  while(i < iter)
  {
    # estimate beta
    b_cap <- beta_cap(y[-j], X[-j,], C[-j,], w[-j])

    # estimate sigma
    sigma2_cap <- sigma_cap(y[-j], X[-j,], b_cap, psi[-j], C[-j,])

    w <- sapply(1:m, function(i)
    {
      w_i <- 1/(sigma2_cap + psi[i] + t(b_cap)%*%diag(C[i,], nrow = p)%*%b_cap)
      return(w_i)
    })
    i <- i + 1
  }

  # estimate gamma
  g_cap <- sapply(1:m, function(i)
  {
      mse_ri <- sigma2_cap + t(b_cap)%*%diag(C[i,], nrow = p)%*%b_cap

      return(mse_ri/(mse_ri+psi[i]))
  })

  y_me <- sapply(1:m, function(i)
  {
      y_me_i <- g_cap[i]*y[i] + (1 - g_cap[i])*t(as.matrix(X[i,]))%*%b_cap
      return(y_me_i)
  })

  hasil <- list('y_me' = y_me, 'gamma' = g_cap, 'beta' = b_cap)
  return(hasil)
}

mse_FHme <- function(y, x, mse_y, mse_x, iter)
{
  X <- as.matrix(x)
  psi <- mse_y
  C <- as.matrix(mse_x)
  m <- length(y)
  p <- dim(X)[2]

  #with full data
  y_me <- estimator_FHme(y, X, psi, C, iter)

  #jackknife
  jack <- lapply(1:m, function(j) jack_knife(y, X, psi, C, j, iter))

  m1 <- sapply(1:m, function(i)
  {
    gam_psi_i <- y_me$gamma[i]*psi[i]

    jack_psi_gamma_i <- ((m-1)/m) * sum(sapply(1:m, function(j)
    {
      jack_gam_i <- jack[[j]]$gamma[i]
      return(gam_psi_i - jack_gam_i*psi[i])
    }))
    return(gam_psi_i + jack_psi_gamma_i)
  })

  m2 <- sapply(1:m, function(i)
  {
    jack_y_me <- ((m-1)/m) * sum(sapply(1:m, function(j)
    {
      return((jack[[j]]$y_me[i] - y_me$y_me[i])^2)
    }))
    return(jack_y_me)
  })
  mse <- m1 + m2

  return(list("mse"=mse))
}

FHme <- function(y, x, vardir, C, iter = 2)
{
  X <- as.matrix(x)
  psi <- as.vector(vardir)
  C <- as.matrix(C)
  m <- length(y)

  if (!identical(dim(X), dim(C)))
    return("X  and MSE of X must be in same length")
  if (length(y) != length(psi))
    return("Y and MSE of Y must be in same length")
  if (m != length(X[,1]))
    return("Y and X must be same length")

  est <- estimator_FHme(y, X, psi, C, iter)

  mse_beta <- mse_FHme(y, x, vardir, C, iter)

  sae <- est
  sae$mse <- mse_beta$mse
  sae$call <- match.call()
  sae$b_delete <- mse_beta$beta
  class(sae) <- 'FHme'
  return(sae)
}

print.FHme <- function(x, ...)
{
  sae <- as.data.frame(x[-c(3,4,5,7)])

  cat("Call:\n")
  print(x$call)
  cat("\nBeta:\n")
  print(as.vector(x$beta))
  cat("\nEsatimated values:\n")
  print(sae)
}

summary.FHme <- function(object, ...)
{
  sum.sae <- list('call' = object$call, 'mse' = object$mse)
  class(sum.sae) <- 'summary.FHme'
  return(sum.sae)
}

print.summary.FHme <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nMSE: \n")
  print(summary(x$mse))
}

plot.FHme <- function(x, ...)
{
  maks <- max(c(x$mse, x$psi))

  plot(x$psi, x$mse, xlim = c(0,maks), ylim = c(0,maks), xlab = 'Direct estimation', ylab = 'Small area estimation', main = 'Mean Square Error')
}
