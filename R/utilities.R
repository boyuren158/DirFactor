#' @import ggplot2

mt.trace = function( A ){
  sum(diag(A))
}

rv.coef = function( A, B ){
  mt.trace(A%*%B)/sqrt(mt.trace(A%*%A)*mt.trace(B%*%B))
}

gg.heatmap <- function( cormt, x.lab="", y.lab="", title="" ){
  plot.dat = reshape2::melt( cormt )
  ggplot( dat = plot.dat, aes( x=Var1, y=Var2, fill=value ) ) + geom_tile(color = "black") + 
    scale_fill_gradient2( low="steelblue", mid = "white", high = "red", limits=c(-1, 1) ) + 
    theme_bw() + coord_fixed() + 
    theme( axis.title = element_text(face="bold"), plot.title = element_text(face="bold"),
           axis.text = element_blank(), axis.ticks = element_blank(),
           panel.border = element_blank(), panel.grid = element_blank() ) + 
    xlab(x.lab) + ylab(y.lab) + ggtitle( title )
}

gg.color.hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

#' Plot ordination results with credible regions based on posterior simulations from DirFactor.
#' 
#' @param all.cov Posterior samples of normalized Gram matrix/distance matrix between biological samples.
#' @param n.dim Total number of axes visualized. Default is 2.
#' @param labels Labels associated with each biological samples. Default is NA. If not NA, labels will be
#' annote onto the figure.
#' @param types Types of each biological samples. Default is NA. If not NA, the contour corresponds to each
#' biological samples will be colored based on the type of this sample.
#' @param dist A logical flag indicates the types of input. Default is FALSE, meaning each element in \code{all.cov}
#' is a normalized Gram matrix. If set to TRUE, each element in \code{all.cov} is a distance matrix.
#' @note The STATIS method is implemented by \code{DistatisR} (Abdi et. al. (2005)).
#' @return A list contains all pairwise ordination plots between the first \code{n.dim} axes.
#' @export
PlotStatis = function( all.cov, n.dim = 2, labels = NA, types = NA, dist = F ){
  all.cov.array = array( unlist( all.cov), dim = c(dim(all.cov[[1]]), length(all.cov)) )
  all.statis.res = DistatisR::distatis( all.cov.array, Distance = dist, nfact2keep = n.dim )
  compromise.ev = eigen(all.statis.res$res4Splus$Splus)$values
  compromise.prop = compromise.ev[1:n.dim]/sum(compromise.ev)
  compromise.coord = apply( compromise.distatis$res4Splus$PartialF, 2, rbind )
  n.rep = nrow( compromise.coord )/nrow(all.cov[[1]])
  
  for( axis.i in 1:(n.dim-1) ){
    for( axis.j in (axis.i):n.dim ){
      plot.data = data.frame( x = all.coord[,axis.i], y = all.coord[,axis.j], 
                              group = rep( 1:nrow(all.cov[[1]]), n.rep) )
      contr = ggplot() + geom_density2d (data = plot.data, aes( x=x, y=y, group = group) ) +
        theme_bw()
      
      if(!is.na(types)){
        types.color = gg.color.hue( length(levels(types) ) )
        plot.data$types = rep(types, n.rep)
        contr = contr + geom_density2d(aes(color = types))
      }
      if(!is.na(labels)){
        x.annote = tapply( plot.data$x, plot.data$group, mean )
        y.annote = tapply( plot.data$y, plot.data$group, mean )
        annote.data = data.frame( x=x.annote, y=y.annote, 
                                  labels = as.character( labels ) )
        contr = contr + geom_point( data = annote.data, aes( x=x, y=y )  ) +
          with(annote.data, annotate(geom="text", x = x+0.01 , y = y, label = labels, size = 8) )
        if(!is.na(types)){
          annote.data$color = types.color[as.numeric(types)]
          contr = contr + geom_point( aes(color = color) ) + with(annote.data,annotate(color=color))
        }
      }
      contr = contr + xlab(sprintf("Compromise axis %d (%.2f%%)", axis.i, compromise.prop[axis.i]*100)) +
        ylab(sprintf("Compromise axis %d (%.2f%%)", axis.j, compromise.prop[axis.j]*100))
      plot.list = list( plot.list, contr )
    }
  }
}

#' Calculate Rhat statistic and generate traceplots for MCMC results
#' 
#' @param lsMCMC A list of MCMC simulation results. \code{lsMCMC>=2} is required for the Rhat statistic
#' calculation.
#' @param start The iteration number of the first observation.
#' @param end The iteration number of the last observation.
#' @param thin The thinning interval between consecutive observations.
#' @param n.eig Number of eigen values the diagnosis will consider. Default is 1.
#' @note
#' The diagnosis is carried out using the MCMC samples of the first \code{n.eig} eigenvalues of the
#' normalized between-sample Gram matrix.
#' @return A matrix contains the Rhat statistics. The first column is the point estimates and the second
#' column is the upper confidence limits. \code{n.eig} traceplots will also be generated for the first 
#' \code{n.eig} eigenvalues.
#' @export
ConvDiagnosis = function( lsMCMC, start, end, thin, n.eig = 1 ){
  mcmc.eig = lapply( lsMCMC, function(indMCMC){
    coda::mcmc( t( matrix( sapply( indMCMC, function(indIter){
      eig.res = eigen( cov2cor( t(indIter$Y)%*%indIter$Y + diag(rep(indIter$er,ncol(indIter$Y)) ) ) )
      eig.res$values[1:n.eig]
    }), nrow = n.eig ) ), start = start, end = end, thin = thin )
  } )
  #traceplots
  par(mfrow = c(n.eig,1) )
  coda::traceplot( mcmc.eig, smooth = F, type = "l" )
  par(mfrow = c(1,1) )
  #calculate Rhat statistics
  return( coda::gelman.diag( mcmc.eig, multivariate = F )$psrf )
}


#' Generate synthetic OTU abundance table from DirFactor model with block diagonal 
#' between-sample Gram matrix.
#' 
#' @param vcounts Required. A vector of total counts per biological sample. Example: 
#' \code{vcounts=c(1e5,1e6)}.
#' @param n Required. Number of biological samples in the simualted table.
#' @param p Required. Number of species in the simulated dataset.
#' @param m Required. True number of factors in the simulated dataset.
#' @param hyper Required. A list of hyper-parameters in the priors. See Details in \link{DirFactor} for the
#' fields in this list.
#' @param K Number of diagonal blocks in the Gram matrix of biological samples. Default is 1 (not block diagonal).
#' @param a Random weights of species are generated by \code{sigma*Q^a*(Q>0)}. Default is 2.
#' @return A list contains the following fields:
#' \itemize{
#'   \item \code{sigma}: a vector with \code{p} components of simulated \code{sigma}.
#'   \item \code{Q}: a \eqn{p*n} matrix of normal latent variables.
#'   \item \code{X.tru}: a \eqn{m*p} matrix of species latent factors.
#'   \item \code{Y.tru}: a \eqn{m*n} matrix of biological sample latent factors.
#'   \item \code{er}: a scalar of the simulated pure error.
#'   \item \code{data}: a list with the same size of \code{vcounts}. \code{data[[i]]} is the simulated
#'   dataset with total counts \code{vcounts[[i]]}.
#'   }
#' 
#' @examples
#' SimDirFactor( vcounts=c(1e5,1e6), n = 22, p = 68, m = 3, hyper = my.hyper )
#' @export
SimDirFactorBlock = function( vcounts, n, p, m, hyper, K = 1, a = 2 ){
  sigma.dist = generate.sigma.prior( p, alpha = hyper$alpha, beta = hyper$beta )
  sigma.value = sigma.dist$sigma.value
  sigma.prior = sigma.dist$sigma.prior
  sigma = sample( sigma.value, p, replace = T, sigma.prior )
  
  X = matrix( rnorm( p*m ), nrow = m )
  
  #split the m factors into K blocks
  factor.indx.group = split( 1:m, cut( 1:m, breaks = K ) )
  pop.indx.group = split( 1:n, cut( 1:n, breaks = K ) )
  Y = matrix( 0, nrow = m, ncol = n )
  
  for( x in 1:K ){
    raw = matrix( rnorm( length(factor.indx.group[[x]])*length(pop.indx.group[[x]]) ), 
                  ncol = length(pop.indx.group[[x]]) )
    Y[factor.indx.group[[x]],pop.indx.group[[x]]] = raw
  }
  
  er = 1/rgamma( 1, hyper$a.er, hyper$b.er )
  Q = t( apply( t(Y)%*%X, 2, function(x) rnorm( length(x), mean = x, sd = sqrt(er) ) ) )
  final.weights = sigma*(Q*(Q>=0))^a
  data = lapply( vcounts, function(counts) apply( final.weights, 2, function(x) rmultinom( 1, counts, prob = x ) ) )
  
  return( list( sigma = sigma, Q = Q, data = data, X.tru = X, Y.tru = Y, er = er ) )
}

#' Generate synthetic OTU abundance table from DirFactor model with block diagonal-compound symmetric 
#' between-sample Gram matrix.
#' 
#' @inheritParams SimDirFactorBlock
#' @param Y.corr A scalar between -1 to 1. Correlation between biological samples belong to the same 
#' block.
#' @return A list contains the following fields:
#' \itemize{
#'   \item \code{sigma}: a vector with \code{p} components of simulated \code{sigma}.
#'   \item \code{Q}: a \eqn{p*n} matrix of normal latent variables.
#'   \item \code{X.tru}: a \eqn{m*p} matrix of species latent factors.
#'   \item \code{Y.tru}: a \eqn{m*n} matrix of biological sample latent factors.
#'   \item \code{er}: a scalar of the simulated pure error.
#'   \item \code{data}: a list with the same size of \code{vcounts}. \code{data[[i]]} is the simulated
#'   dataset with total counts \code{vcounts[[i]]}.
#'   }
#' 
#' @examples
#' SimDirFactorSym( vcounts=c(1e5,1e6), n = 22, p = 68, m = 3, hyper = my.hyper )
#' @export
SimDirFactorSym = function( vcounts, n, p, m, hyper, K = 2, a = 2, Y.corr = 0.5 ){
  sigma.dist = generate.sigma.prior( p, alpha = hyper$alpha, beta = hyper$beta )
  sigma.value = sigma.dist$sigma.value
  sigma.prior = sigma.dist$sigma.prior
  sigma = sample( sigma.value, p, replace = T, sigma.prior )
  
  X = matrix( rnorm( p*m ), nrow = m )
  #split the m factors into K blocks
  factor.indx.group = split( 1:m, cut( 1:m, breaks = K ) )
  pop.indx.group = split( 1:n, cut( 1:n, breaks = K ) )
  Y = matrix( 0, nrow = m, ncol = n )
  
  for( x in 1:K ){
    cor.mt = diag( 1, length(pop.indx.group[[x]]) )
    cor.mt[upper.tri(cor.mt)] = Y.corr
    cor.mt[lower.tri(cor.mt)] = Y.corr
    raw = mvtnorm::rmvnorm( length(factor.indx.group[[x]]), mean=rep(1,length(pop.indx.group[[x]])), sigma=cor.mt )
    Y[factor.indx.group[[x]],pop.indx.group[[x]]] = raw
  }
  
  er = 1/rgamma( 1, hyper$a.er, hyper$b.er)
  Q = t( apply( t(Y)%*%X, 2, function(x) rnorm( length(x), mean = x, sd = sqrt(er) ) ) )
  final.weights = sigma*(Q*(Q>=0))^a
  data = lapply( vcounts, function(counts) apply( final.weights, 2, function(x) rmultinom( 1, counts, prob = x ) ) )
  
  return( list( data = data, sigma = sigma, Q = Q, X.tru=X, Y.tru = Y, er = er ) )
}

#' Generate synthetic OTU abundance table from DirFactor model with two biological sample clusters.
#' 
#' @param strength Required. A positive scalar. Determine the distance between the centroids of two clusters. 
#' Larger value means larger separation.
#' @inheritParams SimDirFactorBlock
#' @return A list contains the following fields:
#' \itemize{
#'   \item \code{sigma}: a vector with \code{p} components of simulated \code{sigma}.
#'   \item \code{Q}: a \eqn{p*n} matrix of normal latent variables.
#'   \item \code{X.tru}: a \eqn{m*p} matrix of species latent factors.
#'   \item \code{Y.tru}: a \eqn{m*n} matrix of biological sample latent factors.
#'   \item \code{er}: a scalar of the simulated pure error.
#'   \item \code{data}: a list with the same size of \code{vcounts}. \code{data[[i]]} is the simulated
#'   dataset with total counts \code{vcounts[[i]]}.
#'   }
#' 
#' @examples
#' SimDirFactorContour( strength = 1, vcounts=c(1e5,1e6), n = 22, p = 68, m = 3, hyper = my.hyper )
#' @export
SimDirFactorContour = function( strength, vcounts, n, p, m, hyper ){
  sigma.dist = generate.sigma.prior( p, alpha = hyper$alpha, beta = hyper$beta )
  sigma.value = sigma.dist$sigma.value
  sigma.prior = sigma.dist$sigma.prior
  sigma = sample( sigma.value, p, replace = T, sigma.prior )
  
  X = matrix( rnorm( p*m ), nrow = m )
  #split the m factors into 2 blocks
  pop.indx.group = split( 1:n, cut( 1:n, breaks = 2 ) )
  Y = matrix( 0, nrow = m, ncol = n )
  
  #split two groups
  #distance between two centroids is 2*strength
  for( x in 1:2 ){
    Y[,pop.indx.group[[x]]] = (2*x - 3)*strength + 
      matrix( rnorm( m*length(pop.indx.group[[x]]) ), 
              ncol = length(pop.indx.group[[x]]) )
  }
  
  er = 1/rgamma( 1, hyper$a.er, hyper$b.er)
  Q = t( apply( t(Y)%*%X, 2, function(x) rnorm( length(x), mean = x, sd = sqrt(er) ) ) )
  final.weights = sigma*(Q*(Q>=0))
  data = lapply( vcounts, function(counts) apply( final.weights, 2, function(x) rmultinom( 1, counts, prob = x ) ) )
  
  return( list( sigma = sigma, Q = Q, data = data,
                X.tru=X, Y.tru = Y, er = er ) )
}