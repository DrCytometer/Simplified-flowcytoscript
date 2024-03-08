
# calculates cross-entropy for tsne plots and calls ce.diff.test


ce.diff.test.tsne <- function( 
    orig.data, tsne.data, 
    event.partition, 
    partition.label = NULL, partition.color = NULL, partition.line.type = NULL, 
    base.test = "ks", base.dist = "ks", 
    prob.sample.n = NULL, dendrogram.order.weight = NULL, 
    result = NULL, cdf.figure = NULL, dendrogram.figure = NULL 
)
{
    stopifnot( nrow( orig.data ) == nrow( tsne.data ) && 
        nrow( orig.data ) == length( event.partition ) )
    
    data.n <- nrow( orig.data )
    
    if ( ! is.null( prob.sample.n ) && prob.sample.n < data.n )
        prob.sample.idx <- sample( data.n, prob.sample.n )
    else
        prob.sample.idx <- 1 : data.n
    
    if ( fcs.use.cached.results && 
        file.exists( fcs.ce.diff.tsne.cache.file.path ) )
    {
        cat( "Using cached results for probability\n" )
        
        load( fcs.ce.diff.tsne.cache.file.path )
    }
    else
    {
        cat( "Calculating probability\n" )
        
        # sampling here temporary, until optimizing dist( tsne.dat ) below
        orig.tsne.prob <- calculate.probability.tsne( 
            orig.data[ prob.sample.idx, ], 
            tsne.data[ prob.sample.idx, ] 
        )
        
        save( orig.tsne.prob, file = fcs.ce.diff.tsne.cache.file.path )
    }
    
    cross.entropy.all <- calculate.cross.entropy( orig.tsne.prob$orig, 
        orig.tsne.prob$tsne )
    
    event.partition.all <- event.partition[ prob.sample.idx ]
    
    ce.diff.test( 
        cross.entropy.all, 
        event.partition.all, 
        partition.label, partition.color, partition.line.type, 
        base.test, base.dist, 
        dendrogram.order.weight, 
        result, cdf.figure, dendrogram.figure
    )
}


calculate.probability.tsne <- function( orig.dat, tsne.dat )
{
    orig.dat.n <- nrow( orig.dat )
    tsne.dat.n <- nrow( tsne.dat )
    
    stopifnot( orig.dat.n == tsne.dat.n )
    
    # find nearest neighbors in original space and their distances
    set.seed.here( fcs.seed.base, "find knn for ce test" )
    
    ce.umap.result <- uwot::umap(orig.dat, 
                              n_neighbors = fcs.ce.diff.tsne.perplexity.factor * fcs.tsne.perplexity + 1,
                              n_epochs = 1, 
                              n_threads = fcs.tsne.threads.n,
                              n_sgd_threads = fcs.tsne.threads.n, 
                              batch = TRUE, verbose = FALSE, ret_model = TRUE,
                              ret_nn = TRUE )
    
    orig.dat.knn <- ce.umap.result$nn$euclidean
    
    orig.dat.self.idx <- sapply( 1 : orig.dat.n, function( ri ) {
        ri.idx <- which( orig.dat.knn$idx[ ri, ] == ri )
        ifelse( length( ri.idx ) == 1, ri.idx, NA )
    } )
    
    stopifnot( ! is.na( orig.dat.self.idx ) )
    
    orig.neigh <- t( sapply( 1 : orig.dat.n, function( ri ) 
        orig.dat.knn$idx[ ri, - orig.dat.self.idx[ ri ] ] ) )
    
    orig.dist2 <- t( sapply( 1 : orig.dat.n, function( ri ) 
        orig.dat.knn$dist[ ri, - orig.dat.self.idx[ ri ] ]^2 ) )
    
    # calculate probabilities associated to distances in original space
    
    # Define the error function outside the apply function
    tsne.perplexity.error <- function(ss, dd2) {
      p <- exp(-dd2 / (2 * ss^2))
      if (sum(p) < .Machine$double.eps)
        p <- 1
      p <- p / sum(p)
      p <- p[p > 0]
      2^(-sum(p * log2(p))) - fcs.tsne.perplexity
    }
    
    # Apply the function to each row of orig.dist2
    orig.stdev <- apply(orig.dist2, 1, function(dd2) {
      # Use vectorized functions to find ss.lower and ss.upper
      ss.lower <- min(dd2[dd2 > 0])
      ss.upper <- max(dd2[is.finite(dd2)])
      
      # Adjust ss.lower and ss.upper
      while (tsne.perplexity.error(ss.upper, dd2) < 0) {
        ss.lower <- ss.upper
        ss.upper <- 2 * ss.upper
      }
      while (tsne.perplexity.error(ss.lower, dd2) > 0) {
        ss.upper <- ss.lower
        ss.lower <- ss.lower / 2
      }
      
      # Find the root of the error function
      uniroot(tsne.perplexity.error, dd2, interval = c(ss.lower, ss.upper), 
              tol = (ss.upper - ss.lower) * .Machine$double.eps^0.25)$root
    })
    
    orig.prob <- t( sapply( 1 : orig.dat.n, function( i ) {
        p <- exp( - orig.dist2[ i, ] /  ( 2 * orig.stdev[ i ]^2 ) )
        p / sum( p )
    } ) )
    
    # symmetrize probabilities in original space
    
    neigh_lengths <- vapply(orig.neigh, length, integer(1))
    
    for (i in seq_len(orig.dat.n)) {
      for (j2 in seq_len(neigh_lengths[i])) {
        j <- orig.neigh[i, j2]
        i2 <- match(i, orig.neigh[j, ])
        
        if (!is.na(i2)) {
          if (j > i) {
            sym.prob <- (orig.prob[i, j2] + orig.prob[j, i2]) / 2
            orig.prob[c(i, j), c(j2, i2)] <- sym.prob
          }
        } else {
          orig.prob[i, j2] <- orig.prob[i, j2] / 2
        }
      }
    }
    
    orig.prob <- sweep( orig.prob, 1, rowSums( orig.prob ), "/" )
    
    # get distances in tsne space for closest neighbors in original space
    
    tsne.dist2 <- t( sapply( 1 : tsne.dat.n, function( i )
        sapply( orig.neigh[ i, ], function( j )
            sum( ( tsne.dat[ i, ] - tsne.dat[ j, ] )^2 )
        )
    ) )
    
    # calculate probabilities associated to distances in tsne representation
    
    tsne.prob.factor <- tsne.dat.n / 
        ( 2 * sum( 1 / ( 1 + dist( tsne.dat )^2 ) ) )
    
    tsne.prob <- t( apply( tsne.dist2, 1, function( dd2 ) 
        p <- tsne.prob.factor / ( 1 + dd2 )
    ) )
    
    list( orig = orig.prob, tsne = tsne.prob )
}


calculate.cross.entropy <- function( prim.prob, secd.prob )
{
    prim.prob.n <- nrow( prim.prob )
    secd.prob.n <- nrow( secd.prob )
    
    prim.prob.m <- ncol( prim.prob )
    secd.prob.m <- ncol( secd.prob )
    
    stopifnot( prim.prob.n == secd.prob.n && prim.prob.m == secd.prob.m )
    
    sapply( 1 : prim.prob.n, function( i ) 
        - sum( prim.prob[ i, ] * log( secd.prob[ i, ] ) )
    )
}

