# MCMC modules
# generate discretized prior of sigma
generate.sigma.prior = function( p, n.points = 999, alpha=10, beta=0 ){
  sigma.value = seq(0.001,0.999,length.out = n.points)
  #take care of small number of species
  if( alpha/p >= 1/2 ){
    alpha = 0.1*p
  }
  tmp = c(0,pbeta( sigma.value, alpha/p, 1/2 + beta - alpha/p ))
  sigma.prior = sapply( 1:n.points, function(x) tmp[x+1]-tmp[x] )
  sigma.prior[length(sigma.prior)] = sigma.prior[length(sigma.prior)] + 1-sum(sigma.prior)
  return( list( sigma.prior = sigma.prior, sigma.value = sigma.value ) )
}

# univariate M-H with normal approx
get.mode = function( data, sigma, T.aug, mu, s ){
  (mu/s^2 + sqrt((mu/s^2)^2+8*data*(2*sigma*T.aug+1/s^2)))/2/(2*sigma*T.aug+1/s^2)
}
get.s = function( Q, data, sigma, T.aug, mu, s ){
  1/sqrt(2*data/Q^2+2*sigma*T.aug+1/s^2)
}
rej.u = function( Q.prev, Q.prop, data, sigma, T.aug, mu, s, mu.prop, s.prop ){
  2*data*log(Q.prop*(Q.prop>0)) - 
    sigma*T.aug*Q.prop^2*(Q.prop>0) + 
    dnorm( Q.prop, mu, s, log = T ) + 
    dnorm( Q.prev, mu.prop, s.prop, log = T ) -
    2*data*log(Q.prev) + sigma*T.aug*Q.prev^2 - 
    dnorm( Q.prev, mu, s, log = T ) - 
    dnorm( Q.prop, mu.prop, s.prop, log = T )
}

# gibbs step for sampling sigma
sigma.gibbs = function( sigma.old, T.aug, Q, sum.species, sigma.value, sigma.prior, sv.log, sp.log, hold = F ){
  #   browser()
  n = length( T.aug )
  p = length( sigma.old )
  if( hold ){
    sigma.old
  }
  else{
    Q.pos = Q*(Q>0)
    tmp.weights = Q.pos^2%*%T.aug
    vapply( 1:p, function(x) {
      log.w = sv.log*sum.species[x] - sigma.value*tmp.weights[x] + sp.log
      w = exp( log.w - max( log.w )  )
      sample( sigma.value, 1, replace = T, prob = w )
    }, 0 )
  }
}

# Gibbs step for sampling augmented variable T
T.aug.gibbs = function( T.aug.old, sigma, Q, sum.sample, hold = F){
  n = length( T.aug.old )
  if( hold ){
    T.aug.old
  }
  else{
    rgamma( n, shape = sum.sample, rate = sigma%*%(Q^2*(Q>0)) )
  }
}

# Gibbs step for sampling Q
Q.gibbs.ind = function( x, T.aug, Sigma.cond ){
  
  #decompress data
  reads = x[1]
  mu.cond = x[2]
  sigma = x[3]
  tmp.Q = x[4]
  
  if( reads==0 ){
    mu.new = mu.cond/(1+2*sigma*T.aug*Sigma.cond^2)
    sigma.new = 1/sqrt(1/Sigma.cond^2+2*sigma*T.aug)
    
    #In real data, this p1 can be really small
    p1.log = pnorm( 0, mu.cond, Sigma.cond, log.p = T ) - pnorm( 0, mu.new, sigma.new, lower.tail = F, log.p = T ) +
      dnorm( 0, mu.new, sigma.new, log = T ) - dnorm( 0, mu.cond, Sigma.cond, log = T )
    p1 = exp( p1.log )
    if( p1 == Inf ){
      p.ratio = 1
    }
    else{
      p.ratio = p1/(1+p1)
    }
    
    if( runif(1)<p.ratio ){
      truncnorm::rtruncnorm( 1, b = 0, mean = mu.cond, sd = Sigma.cond )
    }
    else{
      truncnorm::rtruncnorm( 1, a = 0, mean = mu.new, sd = sigma.new )
    }
  }
  
  #metropolis-hastings for non-zero observation
  else{
    Q.mode = get.mode( reads, sigma, T.aug, mu.cond, Sigma.cond )
    Q.sd = get.s( Q.mode, reads, sigma, T.aug, mu.cond, Sigma.cond )
    Q.prop = rnorm( 1, Q.mode, Q.sd )
    
    if( -rexp(1) < rej.u( tmp.Q, Q.prop, reads, sigma, 
                          T.aug, mu.cond, Sigma.cond, Q.mode, Q.sd ) ){
      Q.prop
    }
    else{
      tmp.Q
    }
  }
}

# Gibbs step for sampling Qij with nij = 0
Q.gibbs.vec.0 = function( para.0 ){
  n.0 = nrow( para.0 )
  mu.new.0 = para.0[,1]/(1+2*para.0[,2]*para.0[,3]*para.0[,4]^2)
  sigma.new.0 = 1/sqrt(1/para.0[,4]^2+2*para.0[,2]*para.0[,3])
  p1.log = pnorm( 0, para.0[,1], para.0[,4], log.p = T ) - pnorm( 0, mu.new.0, sigma.new.0, lower.tail = F, log.p = T ) +
    dnorm( 0, mu.new.0, sigma.new.0, log = T ) - dnorm( 0, para.0[,1], para.0[,4], log = T )
  p1 = exp( p1.log )
  p.ratio = p1/(1+p1)
  p.ratio[is.na(p.ratio)] = 1
  
  label = ( runif( n.0 ) < p.ratio )
  #construct upper and lower vectors
  lower = rep( -Inf, n.0 )
  lower[!label] = 0
  upper = rep( Inf, n.0 )
  upper[label] = 0
  #construct mean and sd vectors
  mean = mu.new.0
  mean[label] = para.0[label,1]
  sd = sigma.new.0
  sd[label] = para.0[label,4]
  
  truncnorm::rtruncnorm( n.0, a = lower, b = upper, mean = mean, sd = sd )
}

# Metropolis step for sampling Qij with nij >0 0
Q.gibbs.vec.n0 = function( para.n0, tmp.Q.n0 ){
  Q.mode = get.mode( para.n0[,1], para.n0[,3], para.n0[,4], para.n0[,2], para.n0[,5] )
  Q.sd = get.s( Q.mode, para.n0[,1], para.n0[,3], para.n0[,4], para.n0[,2], para.n0[,5] )
  n.n0 = length( Q.mode )
  Q.prop = rnorm( n.n0, Q.mode, Q.sd )
  
  #metropolis hastings
  rej = rej.u( tmp.Q.n0, Q.prop, 
               para.n0[,1], para.n0[,3], para.n0[,4], para.n0[,2], para.n0[,5], 
               Q.mode, Q.sd )
  labels = ( -rexp(n.n0) < rej )
  
  Q.ret = tmp.Q.n0
  Q.ret[labels] = Q.prop[labels]
  Q.ret
}

#' Gibbs sampler for Depedent Dirichlet factor model.
#' 
#' @param data Required. A count matrix with species in rows and biological samples in column.
#' @param hyper Required. A list of hyper-parameters in the priors.
#' @param start A list of starting values of model parameters. Default is \code{NA} and the starting values
#'  will be generated automatically based on the input data. See Details for the required fields of the list.
#' @param save.path A string contains the path to save the MCMC results. 
#' For example, \code{save.path="~/sim"} will save results of the ith iteration to ~/sim_i.rds.
#' Default is \code{NA} and a temp directory will be assigned.
#' @param save.obj A list of model parameters that will be saved, default is all parameters.
#' @param burnin A number between 0 and 1. Fraction of burn-in samples. Default is 0.2.
#' @param thinning A positive integer. The MCMC results will be saved every \code{thinning} iterations. 
#' Default is 5.
#' @param step A positive integer. The total number of MCMC iterations. Default is 1000.
#' @param step.disp A positive integer. A message will report the number of iterations finished 
#' every \code{step.disp} iterations. Default is 10. 
#' @return A list with two fields: \code{running.time} and \code{save.path}. 
#' \code{running.time} is the total amount of time for finishing the MCMC simulation. 
#' \code{save.path} is the path to the saved results.
#' @details The Dependent Dirichlet process factor model assumes the observed data is distributed
#' according to a multinomial distrition for each biological sample, conditioning on the probabilities
#' of species, which is assumed to follow a Dependent Dirichlet processes a priori. 
#' The model has two major parts. \code{sigma, Q} directly specify the probabilities of species in each
#' biological samples and \code{X, Y, er, delta, phi} specify the hyper-prior on \code{Q}. \code{T.aug} is
#' an auxilary parameter and does not have direct interpretation. More details on model and prior 
#' specification can be found in Ren et. al. (2016).
#' 
#' Users are required to provide a list containing the values of the hyper-parameters. The list must contain
#' fields as following. 
#' \itemize{
#'  \item \code{alpha} and \code{beta}: hyper-parameters in the prior of \code{sigma}. \code{sigma} follows a 
#' Poisson process on \eqn{(0,1)} with intensity \eqn{\alpha\sigma^{-1-\beta}\cdot(1-\sigma)^{-1/2+\beta}} a priori.
#'  \item \code{a.er} and \code{b.er}: hyper-parameters in the prior of \code{er}. \code{er} follows an inverse
#'   gamma distribution \eqn{1/Gamma(a.er,b.er)} a priori.
#'  \item \code{m}: number of factors the model is assumed to have. Usually a value much smaller than the number
#'  of biological samples.
#'  \item \code{nv}: hyper-parameter of the prior of \code{phi}. \code{phi} follows \eqn{Gamma(nv/2,nv/2)} a priori.
#'  \item \code{a1} and \code{a2}: hyper-parameters in the prior of \code{delta}. The prior of \code{delta[1]} is \eqn{Gamma(a1,1)}
#'  and the prior of \code{delta[2:m]} is independent \eqn{Gamma(a2,1)}. \code{a2} must not be smaller than 1 in order to
#'  achieve factor shrinkage Bhattacharya et al. (2011).
#' }
#' 
#' If the users want to specify the starting values for the model parameters, they can pass a list with fields
#' \code{sigma, Q, T.aug, X, Y, er, delta, phi} to the function augment \code{start}. Assume there are
#' \code{n} biological samples and \code{p} species in \code{data}. Each field is specified as following:
#' \itemize{
#'   \item \code{sigma}: a vector with \code{p} components and the starting values have to be in \eqn{(0,1)}.
#'   \item \code{Q}: a matrix with \code{p} rows and \code{n} columns and the starting values of \code{Q} 
#'   should be positive if the corresponding cell in \code{data} is positive.
#'   \item \code{T.aug}: a vector with \code{n} positive components. 
#'   \item \code{X}: a \eqn{m\times p} matrix.
#'   \item \code{Y}: a \eqn{m\times n} matrix. 
#'   \item \code{er}: a positive scalar. 
#'   \item \code{delta}: a positive vector with \code{m} components.
#'   \item \code{phi}: a \eqn{n\times m} positive matrix. 
#' }
#' 
#' @examples
#' my.hyper = list( nv = 3, a.er = 1, b.er = 0.3, a1 = 3, 
#'                  a2 = 4, m = 10, alpha = 10, beta = 0 )
#' my.sim = SimDirFactorBlock( 1e6, n = 22, p = 68, m = 3, my.hyper, K = 2 )
#' DirFactor( my.sim$data[[1]], my.hyper, save.obj = c("Y", "er"), step.disp = 100 )
#' @export
DirFactor <- function( data, hyper, start = NA, save.path = NA, 
                       save.obj = c("sigma", "Q", "T.aug", "X", "Y", "er", "delta", "phi"), 
                       burnin = 0.2, thinning = 5, step = 1000, step.disp = 10 ){
  n = ncol( data )
  p = nrow( data )
  sum.species = rowSums( data )
  sum.sample = colSums( data )
  
  sigma.dist = generate.sigma.prior( p, hyper$alpha, hyper$beta )
  sigma.value = sigma.dist$sigma.value
  sigma.prior = sigma.dist$sigma.prior
  sv.log = log( sigma.value )
  sp.log = log( sigma.prior )
  
  #initialization
  if(is.na(start)){
    sigma.old = sample( sigma.value, size = p, replace = T, prob = sigma.prior ) 
    T.aug.old = colSums( data )
    Q.old = matrix( 0.5, nrow = p, ncol = n )
    X.old = matrix( rnorm( hyper$m*p ), nrow = hyper$m )
    Y.old = matrix( rnorm( hyper$m*n ), nrow = hyper$m )
    er.old = 1/rgamma( 1, shape = hyper$a.er, rate = hyper$b.er )
    delta.old = c( rgamma( 1, shape = hyper$a1, rate = 1 ), 
                   rgamma( hyper$m-1, shape = hyper$a2, rate = 1 ) )
    phi.old = matrix( rgamma( hyper$m*n, shape = hyper$nv/2, rate = hyper$nv/2 ), nrow = n )
  }else{
    sigma.old = as.vector( start$sigma )
    Q.old = start$Q
    T.aug.old = as.vector( start$T.aug )
    X.old = start$X
    Y.old = start$Y
    er.old = start$er
    delta.old = as.vector(start$delta)
    phi.old = start$phi
  }
  
  all.cache = vector("list", length = length(save.obj))
  names(all.cache) = save.obj
  
  if( is.na(save.path) ){
    save.path = paste(normalizePath(tempdir(), winslash = "/"), "sim", sep="/")
  }
  
  t.t = proc.time()
  for( iter in 1:step ){
    if( iter %% step.disp == 0 )
      message( sprintf("%s iterations finished", iter) )
    
    #sample sigma
    #this is a major time consuming step
    if( 'sigma' %in% save.obj ){
      all.cache$sigma = sigma.old
    }
    sigma.old = sigma.gibbs( sigma.old, T.aug.old, Q.old, sum.species, 
                             sigma.value, sigma.prior,sv.log, sp.log )
    
    #sample T
    if( 'T.aug' %in% save.obj ){
      all.cache$T.aug = T.aug.old
    }
    T.aug.old = T.aug.gibbs( T.aug.old, sigma.old, Q.old, sum.sample )
    
    #sample Q
    #matrix of all parameters
    Q.para.all = cbind( c(data), c(t(X.old)%*%Y.old), rep( sigma.old, n ), c(Q.old),
                        rep( T.aug.old, each = p ), rep( sqrt(er.old), n*p ) )
    if( 'Q' %in% save.obj ){
      all.cache$Q = Q.old
    }
    Q.labels = (Q.para.all[,1]==0)
    Q.0 = Q.gibbs.vec.0( Q.para.all[Q.labels, c(2,3,5,6)] )
    Q.n0 = Q.gibbs.vec.n0( Q.para.all[!Q.labels, c(1,2,3,5,6)], Q.para.all[!Q.labels, 4] )
    #save it into proper locations
    Q.tmp = rep( 0, nrow( Q.para.all ) )
    Q.tmp[Q.labels] = Q.0
    Q.tmp[!Q.labels] = Q.n0
    #convert it back to matrix
    Q.old = matrix( Q.tmp, nrow = p )
    
    #sample Y
    #this is a major time consuming step
    #get tau
    tau.tmp = cumprod( delta.old )
    if( "Y" %in% save.obj ){
      all.cache$Y = Y.old
    }
    Y.old = vapply( 1:n, function(x){
      Sigma.Y = solve( X.old%*%t(X.old)/er.old + diag( phi.old[x,]*tau.tmp ) )
      mu.Y = Sigma.Y%*%X.old%*%Q.old[,x]/er.old
      mvtnorm::rmvnorm( 1, mu.Y, sigma = Sigma.Y )
    }, rep( 0, nrow( Y.old ) ) )
    
    #sample er
    #prior is gamma(a.er,b.er)
    if( "er" %in% save.obj ){
      all.cache$er = er.old
    }
    er.old = 1/rgamma( 1, shape = n*p/2+hyper$a.er, 
                       rate = sum((Q.old-t(X.old)%*%Y.old)^2)/2+hyper$b.er )
    
    #sample X
    Sigma.X = solve( Y.old%*%t(Y.old)/er.old + diag( hyper$m ) )
    if( "X" %in% save.obj ){
      all.cache$X = X.old
    }
    X.mean = Sigma.X%*%Y.old%*%t(Q.old)/er.old
    X.old = X.mean + t( mvtnorm::rmvnorm( p, sigma = Sigma.X ) )
    
    #sample phi
    if( "phi" %in% save.obj ){
      all.cache$phi = phi.old
    }
    phi.old = matrix( rgamma( n*hyper$m, shape = (hyper$nv+1)/2, rate = c( Y.old^2*tau.tmp + hyper$nv )/2 ), 
                      nrow = n, byrow = T )
    
    #sample delta
    delta.tmp = delta.old
    if( "delta" %in% save.obj ){
      all.cache$delta = delta.old 
    }
    for( k in 1:hyper$m ){
      if( k == 1 ){
        delta.prime = cumprod( delta.tmp )/delta.tmp[1]
        delta.tmp[k] = rgamma( 1, shape = n*hyper$m/2 + hyper$a1, 
                               rate = 1 + sum(colSums( t(Y.old^2)*phi.old )*delta.prime)/2
        )
      }
      else{
        delta.prime = cumprod( delta.tmp )[-(1:(k-1))]/delta.tmp[k]
        delta.tmp[k] = rgamma( 1, shape = n*(hyper$m-k+1)/2 + hyper$a2, 
                               rate = 1 + sum( colSums( t(Y.old^2)*phi.old )[-(1:(k-1))]*delta.prime)/2
        )
      }
    }
    delta.old = delta.tmp
    
    if( iter > burnin*step & iter%%thinning == 0 )
      saveRDS( all.cache[save.obj], file = paste( save.path, iter, sep = "_") )
  }
  t.t.1 = proc.time()
  return( list(running.time = t.t.1 - t.t, save.path = save.path) )
}
