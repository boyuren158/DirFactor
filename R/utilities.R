#' @import ggplot2

ListtoArray = function( res.all ){
  array( unlist(res.all), dim = c(dim(res.all[[1]]), length(res.all)) )
}

mt.trace = function( A ){
  sum(diag(A))
}

#' Calculate RV coefficients between two correlation matrices
#' @param A, B two correlation matrices.
#' @export
rv.coef = function( A, B ){
  mt.trace(A%*%B)/sqrt(mt.trace(A%*%A)*mt.trace(B%*%B))
}

gg.heatmap <- function( cormt, x.lab="", y.lab="", title="" ){
  plot.dat = reshape2::melt( cormt )
  ggplot( data = plot.dat, aes( x=Var1, y=Var2, fill=value ) ) + geom_tile(color = "black") + 
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
#' @param levels.min Only show the minimum level in the contours. Default is FALSE.
#' @param labels Labels associated with each biological samples. Default is NA. If not NA, labels will be
#' annote onto the figure.
#' @param types Types of each biological samples. Default is NA. If not NA, the contour corresponds to each
#' biological samples will be colored based on the type of this sample.
#' @param dist A logical flag indicates the types of input. Default is FALSE, meaning each element in \code{all.cov}
#' is a normalized Gram matrix. If set to TRUE, each element in \code{all.cov} is a distance matrix.
#' @note The STATIS method is implemented by \code{DistatisR} (Abdi et. al. (2005)).
#' @return A list contains all pairwise ordination plots between the first \code{n.dim} axes.
#' @export
PlotStatis = function( all.cov, n.dim = 2, levels.min = FALSE, labels = NA, types = NA, dist = F, ... ){
  all.cov.array = array( unlist( all.cov ), dim = c(dim(all.cov[[1]]), length(all.cov)) )
  all.statis.res = DistatisR::distatis( all.cov.array, Distance = dist, nfact2keep = n.dim )
  compromise.ev = eigen(all.statis.res$res4Splus$Splus)$values
  compromise.prop = compromise.ev[1:n.dim]/sum(compromise.ev)
  compromise.coord = apply( all.statis.res$res4Splus$PartialF, 2, rbind )
  n.rep = nrow( compromise.coord )/nrow(all.cov[[1]])
  
  apply( combn(1:n.dim, 2), 2, function(axis.idxs){
    axis.i = axis.idxs[1]
    axis.j = axis.idxs[2]
    plot.data = data.frame( x = compromise.coord[,axis.i], y = compromise.coord[,axis.j], 
                            group = rep( 1:nrow(all.cov[[1]]), n.rep) )
    
    if(!any(is.na(types))){
      types = as.factor(types)
      types.color = gg.color.hue( length(levels(types) ) )
      plot.data$types = rep(types, n.rep)
      contr = ggplot() + geom_density2d(data = plot.data, aes( x=x, y=y, group = group, color = types)) + 
        scale_color_manual( values = types.color ) + theme_bw()
    }
    else{
      contr = ggplot() + geom_density2d (data = plot.data, aes( x=x, y=y, group = group) ) +
        theme_bw()
    }
    if(levels.min){
      contr.data = ggplot_build(contr)$data[[1]]
      contr.min = contr.data[grep("-001", contr.data$group, fixed = T),]
      contr = ggplot(NULL) + geom_path( data = contr.min, aes(x=x,y=y,group=group,color=colour), ... ) + 
        scale_color_manual( values = sort(unique(contr.min$colour)), labels = unique(types)[order(unique(contr.min$colour))] ) + 
        theme_bw()
    }
    
    if(!any(is.na(labels))){
      x.annote = tapply( plot.data$x, plot.data$group, mean )
      y.annote = tapply( plot.data$y, plot.data$group, mean )
      annote.data = data.frame( x=x.annote, y=y.annote, 
                                labels = as.character( labels ) )
      contr = contr + geom_point( data = annote.data, aes( x=x, y=y )  ) +
        with(annote.data, annotate(geom="text", x = x+0.01 , y = y, label = labels, size = 6) )
    }
    contr + xlab(sprintf("Compromise axis %d (%.2f%%)", axis.i, compromise.prop[axis.i]*100)) +
      ylab(sprintf("Compromise axis %d (%.2f%%)", axis.j, compromise.prop[axis.j]*100))
  } )
}

#' Calculate Rhat statistic and generate traceplots for MCMC results
#' 
#' @param lsMCMC A list of MCMC simulation results. \code{lsMCMC>=2} is required for the Rhat statistic
#' calculation.
#' @param start The iteration number of the first observation.
#' @param end The iteration number of the last observation.
#' @param thin The thinning interval between consecutive observations.
#' @param title The title of each figure.
#' @param If provided, a list contains the simulation truth for the dataset. 
#' Must include two fields \code{Y}, the latent biological sample factors, and
#' \code{er}, the variance of the pure error.
#' @param n.eig Number of eigen values the diagnosis will consider. Default is 1.
#' @param fast.eig Whether to use fast eigen decomposition algorithm. Default is TRUE and \code{eig_sym} from
#' \code{rARPACK} package will be used.
#' @note
#' The diagnosis is carried out using the MCMC samples of the first \code{n.eig} eigenvalues of the
#' normalized between-sample Gram matrix.
#' @return A matrix contains the Rhat statistics. The first column is the point estimates and the second
#' column is the upper confidence limits. \code{n.eig} traceplots will also be generated for the first 
#' \code{n.eig} eigenvalues.
#' @export
ConvDiagnosis = function( lsMCMC, start, end, thin, title, truth = NA, n.eig = 1, fast.eig = TRUE ){
  mcmc.eig = ListtoArray( lapply( lsMCMC, function(indMCMC){
    t( matrix( sapply( indMCMC, function(indIter){
      if( fast.eig ){
        print("yes")
        rARPACK::eigs_sym( cov2cor( t(indIter$Y)%*%indIter$Y + diag(rep(indIter$er,ncol(indIter$Y)) ) ), k = n.eig, which = "LM" )$values
      }
      else{
        eig.res = eigen( cov2cor( t(indIter$Y)%*%indIter$Y + diag(rep(indIter$er,ncol(indIter$Y)) ) ) )
        eig.res$values[1:n.eig]
      }
    }), nrow = n.eig ) )
  } ) )
  
  mcmc.eig.obj = lapply( 1:n.eig, function(x) 
    lapply( 1:ncol(mcmc.eig[,x,]), function(x.i) coda::mcmc(mcmc.eig[,x,x.i],start=start,end=end,thin=thin) ) )
  
  #traceplots
  if(!any(is.na(truth))){
    ev.tru = eigen( cov2cor( t(truth$Y)%*%truth$Y + diag(truth$er, nrow = ncol(truth$Y)) ) )$values[1:n.eig]
  }
  rhat.all = matrix( nrow = n.eig, ncol = 2 )
  for( i in 1:n.eig ){
    rhat = as.vector( coda::gelman.diag( mcmc.eig.obj[[i]], multivariate = F )$psrf )
    rhat.all[i,] = rhat
    coda::traceplot( mcmc.eig.obj[[i]], ylab = sprintf("Eigenvalue %d", i), lty = 1,
                     main = paste( title, " Rhat=", round( rhat[1], digits = 3 ), sep = "" )  )
    if(!any(is.na(truth))){
      abline( h = ev.tru[i], col = "blue", lwd = 2 )
    }
  }
  #calculate Rhat statistics
  return( rhat.all )
}

#' Plot posterior probability of pairwise classification based on posterior 
#' samples of pairwise distances between biological samples.
#' 
#' @param all.dist A list contains posterior samples of pairwise distances
#' between biological samples. Each element is a distance matrix.
#' @param labels A character vector indicates the label for each biological 
#' sample.
#' @return A ggplot figure object. It illustrates the posterior probability
#' of two biological samples being clustered together. The rows and columns 
#' are sorted by labels of biological samples.
#' @export
PlotClustering = function( all.dist, labels = NA ){
  #get posterior clustering probability
  res.cluster = lapply( all.dist, function(x){
    x = as.matrix(x)
    cluster = fpc::pamk( x, krange = seq(2,ncol(x)-1), diss = T )
    cluster.id = cluster$pamobject$clustering
    outer( cluster.id, cluster.id, "==" )
  })
  res.cluster.mt = array( unlist( res.cluster ), dim = c( dim( res.cluster[[1]] ), length( res.cluster ) ) )
  res.cluster.mean = apply( res.cluster.mt, 1:2, mean )
  
  #plot the clustering results
  if( all(!is.na(labels)) ){
    res.mt.order = res.cluster.mean[order( labels ),order( labels )]
  }
  else{
    res.mt.order = res.cluster.mean
  }
  plot.res = reshape2::melt( res.mt.order )
  
  plt.out = ggplot() + geom_tile( data = plot.res, aes( x=Var1, y=Var2, fill = value ), color = "black" ) + 
    scale_fill_gradient( low = "white", high = "red" ) + guides(fill=guide_legend(title="Posterior\nprobability")) + theme_bw()
  
  if( all(!is.na(labels)) ){
    axis.labels = levels( labels )
    axis.pos = cumsum( table( labels ) ) + 0.5
    axis.text.pos = as.vector( axis.pos - table( labels )/2 )
    x.delim = data.frame( x = axis.pos, xend = axis.pos, yend = rep( 0.5, length(axis.pos) ), y = rep( 0, length(axis.pos) ) )
    y.delim = data.frame( y = axis.pos, yend = axis.pos, xend = rep( 0.5, length(axis.pos) ), x = rep( 0, length(axis.pos) ) )
    plt.out = plt.out + annotate( geom="text", x = axis.text.pos, y = 0, label = axis.labels, size = 6, angle = 30, hjust = 1 ) + 
      annotate( geom="text", y = axis.text.pos, x = 0, label = axis.labels, size = 6, angle = 30, hjust = 1 ) + 
      theme_bw() + coord_fixed() + 
      theme( axis.text = element_blank( ),
             axis.title = element_blank(),
             axis.ticks = element_blank(),
             panel.border = element_blank(),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             legend.title = element_text(size = 18,face = "bold"),
             legend.text = element_text(size = 15) ) + 
      geom_segment( data = x.delim, aes(x = x, xend = xend, y = y, yend = yend ) ) + 
      geom_segment( data = y.delim, aes(x = x, xend = xend, y = y, yend = yend ) )
  }
  plt.out
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
#' my.hyper = list( nv = 3, a.er = 1, b.er = 0.3, a1 = 3, 
#'                  a2 = 4, m = 10, alpha = 10, beta = 0 )
#' SimDirFactorBlock( vcounts=c(1e5,1e6), n = 22, p = 68, m = 3, hyper = my.hyper )
#' @export
SimDirFactorBlock = function( vcounts, n, p, m, hyper, K = 1, a = 2 ){
  sigma.dist = generate.sigma.prior( p, alpha = hyper$alpha, beta = hyper$beta )
  sigma.value = sigma.dist$sigma.value
  sigma.prior = sigma.dist$sigma.prior
  sigma = sample( sigma.value, p, replace = T, sigma.prior )
  
  X = matrix( rnorm( p*m ), nrow = m )
  
  #split the m factors into K blocks
  if( K > 1 ){
    factor.indx.group = split( 1:m, cut( 1:m, breaks = K ) )
    pop.indx.group = split( 1:n, cut( 1:n, breaks = K ) ) 
  }else{
    factor.indx.group = list( 1:m )
    pop.indx.group = list( 1:n )
  }
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
#' my.hyper = list( nv = 3, a.er = 1, b.er = 0.3, a1 = 3, 
#'                  a2 = 4, m = 10, alpha = 10, beta = 0 )
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
#' my.hyper = list( nv = 3, a.er = 1, b.er = 0.3, a1 = 3, 
#'                  a2 = 4, m = 10, alpha = 10, beta = 0 )
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

#' Generate synthetic OTU table from DirFactor model using latent factors provided by users.
#' 
#' @inheritParams SimDirFactorBlock
#' @param sigma A length \code{p} vector. \code{p} is the number of species. 
#' Global abundance parameter for species.
#' @param X A \code{m*p} matrix. \code{m} is the number of factors. Latent 
#' factors for species.
#' @param Y A \code{m*n} matrix. \code{n} is the number of biological samples. 
#' Latent factors for biological samples.
#' @param er A positive scalar. Variance of the pure noise.
#' 
#' @return A list with \code{length(vcounts)} elements. The \code{i}th elment is 
#' an OTU table with total counts per biological sample \code{vcounts[i]}.
#' @export
SimDirFactorCustom = function( vcounts, sigma, X, Y, er ){
  Q = t(X)%*%Y + matrix( rnorm( ncol(X)*ncol(Y), sd = sqrt(er) ), nrow = ncol(X) )
  final.weights = sigma*Q^2*(Q>0)
  data = lapply( vcounts, function(counts) apply( final.weights, 2, function(x) rmultinom( 1, counts, prob = x ) ) )
  return( list( data = data ) )
}