# run clustering

# run clustering via FlowSOM or Phenograph----------
cat(
  "\n
As the first part of the analysis, the script will cluster your cells into groups.

For this, you can use either Phenograph or FlowSOM.

Phenograph is fast, and will automatically determine how many clusters to generate.
Sometimes it overclusters, particularly in cases with lots of homogeneous cells.

If you prefer to use FlowSOM, you'll need to decide how many clusters you want to find.

One strategy that may be helpful is first run Phenograph, assess approximately how
many clusters appear to really be distinct, then repeat the analysis using FlowSOM
with a defined number of clusters. Alternatively, examine the umap plot in figure_umap
to see how many islands there are.
\n"
)

# calculate umap representation---------------

{
  if ( fcs.use.cached.results && file.exists( fcs.umap.cache.file.path ) )
  {
    cat( "Using cached results for umap\n" )
    
    load( fcs.umap.cache.file.path )
  }
  else { 
    set.seed.here( fcs.seed.base, "calculate umap representation" )
    
    cat( "Calculating UMAP\n" )
    
    umap.result <- uwot::umap(dmrd.data, n_neighbors = fcs.tsne.perplexity,
                              n_epochs = fcs.umap.iter.n, 
                              n_threads = fcs.tsne.threads.n,
                              n_sgd_threads = fcs.tsne.threads.n, 
                              batch = TRUE, verbose = TRUE, ret_model = TRUE,
                              ret_nn = TRUE, ret_extra = "sigma")
    
    save(umap.result, file = fcs.umap.cache.file.path )
  }
}

umap.data <- umap.result$embedding
dimnames( umap.data ) <- NULL

dmrd.knn <- umap.result$nn$euclidean

umap.islands.plot <- ggplot(data.frame(umap.data), aes(x = X1, y = X2 )) + 
  geom_scattermore( alpha = fcs.dmrd.color.alpha, 
    pointsize = fcs.umap.figure.point.size*2,
    pixels = c(700,700))+
  theme_bw() + 
  theme( axis.title = element_blank(), 
         axis.text  = element_blank(), 
         axis.ticks = element_blank(), 
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         strip.background = element_blank(),
         plot.background = element_blank() )

ggsave( 
  file.path( fcs.umap.figure.dir, "umap_islands_plot.jpeg" ), 
  umap.islands.plot, 
  width = fcs.dmrd.figure.width, height = fcs.dmrd.figure.height 
)

umap.islands.plot

# select clustering method, run clustering---------------

clustering.method <- menu( c("Phenograph", "FlowSOM"),
                           title = "Please choose your clustering approach.")

cat("\n")

if (clustering.method == 2){
  fcs.cluster.n <- readline( "Please enter the number of clusters to find using FlowSOM: ")
  fcs.cluster.n <- as.numeric(fcs.cluster.n)
  
  cat("\n")
  
  fcs.cluster <- sprintf( "%02d", 1 : fcs.cluster.n )
  
  fcs.cluster.label <- sprintf( "Cluster-%s", fcs.cluster )
  names( fcs.cluster.label ) <- fcs.cluster
  
  fcs.cluster.color <- rep( 
    fcs.color.pool, 
    ceiling( fcs.cluster.n / fcs.color.pool.n ) 
  )[ 1 : fcs.cluster.n ]
  names( fcs.cluster.color ) <- fcs.cluster
  
  fcs.cluster.line.type <- rep( 
    fcs.line.type.pool, 
    ceiling( fcs.cluster.n / fcs.line.type.pool.n ) 
  )[ 1 : fcs.cluster.n ]
  names( fcs.cluster.line.type ) <- fcs.cluster
  
  fcs.flow.som.dim <- 10
  
  # get FlowSOM clusters---------------
  # build som objects--embedsom method
  cat("Clustering data with FlowSOM via EmbedSOM\n")
  set.seed.here( fcs.seed.base, "get flowsom clusters" )
  
  flow.som <- EmbedSOM::SOM(dmrd.data, xdim = fcs.flow.som.dim, 
                            ydim = fcs.flow.som.dim, batch = TRUE,
                            parallel = TRUE, threads = fcs.tsne.threads.n )
  
  # get clusters
  flow.som.mapping <- flow.som$mapping[ , 1 ]
  flow.som.codes <- flow.som$codes
  
  # get clusters from som mapping
  consensus.cluster <- ConsensusClusterPlus( t( flow.som.codes ),
                                             maxK = fcs.cluster.n, reps = 100, pItem = 0.9, pFeature = 1,
                                             clusterAlg = "hc", innerLinkage = "average", finalLinkage = "average",
                                             distance = "euclidean", 
                                             seed = set.seed.here( fcs.seed.base, "get clusters from som mapping" ) )
  
  flow.som.event.cluster <- consensus.cluster[[ fcs.cluster.n ]]$
    consensusClass[ flow.som.mapping ]
  
  # reorder clusters from bigger to smaller
  flow.som.cluster.rank <- 1 + fcs.cluster.n - 
    rank( table( flow.som.event.cluster ), ties.method = "last" )
  flow.som.event.cluster <- flow.som.cluster.rank[ flow.som.event.cluster ]
  names( flow.som.event.cluster ) <- NULL
  
  # set clusters as a factor
  flow.som.event.cluster <- factor( flow.som.event.cluster, 
                                    levels = 1 : fcs.cluster.n )
  levels( flow.som.event.cluster ) <- fcs.cluster
  
  dmrd.event.cluster <- flow.som.event.cluster
  dmrd.event.cluster.n <- as.vector( table( dmrd.event.cluster ) )
  names( dmrd.event.cluster.n ) <- fcs.cluster
  
  length( dmrd.event.cluster )
  table( dmrd.event.cluster )
  
} else {
  # get phenograph clusters---------------
  set.seed.here( fcs.seed.base, "get clusters with Phenograph" )
  
  # fcs.phenograph.engine is set during flowcytoscript_startup.r: "FastPG" if it
  # installed successfully, "Rphenograph" if we fell back to that instead (e.g.
  # because FastPG's compiler needs OpenMP support that isn't set up on this Mac),
  # or NA if neither is available.
  phenograph.engine <- if (exists("fcs.phenograph.engine")) fcs.phenograph.engine else NA_character_
  if (is.na(phenograph.engine)) {
    phenograph.engine <- if (requireNamespace("FastPG", quietly = TRUE)) "FastPG" 
    else if (requireNamespace("Rphenograph", quietly = TRUE)) "Rphenograph" 
    else NA_character_
  }
  
  if (identical(phenograph.engine, "FastPG")) {
    
    cat("Clustering data with Phenograph (via FastPG)\n")
    
    links <- FastPG::rcpp_parallel_jce( dmrd.knn$idx )
    links <- FastPG::dedup_links( links )
    clusters <- FastPG::parallel_louvain( links )
    
    phenograph.event.cluster <- factor(clusters$communities)
    
  } else if (identical(phenograph.engine, "Rphenograph")) {
    
    cat("Clustering data with Phenograph (via Rphenograph fallback,\n
since FastPG isn't available on this system)\n
Note: Rphenograph builds its own nearest-neighbor graph directly from the\n
reduced data rather than reusing the UMAP neighbor graph the way FastPG\n
does, so cluster boundaries may differ slightly from a FastPG run on the\n
same data.\n")
    
    rphenograph.out <- Rphenograph::Rphenograph( dmrd.data, k = fcs.tsne.perplexity )
    phenograph.event.cluster <- factor( igraph::membership( rphenograph.out[[2]] ) )
    
  } else {
    stop("Neither FastPG nor Rphenograph is installed, so Phenograph clustering\n
isn't available. Please re-run flowcytoscript_startup.r, or choose FlowSOM\n
instead when prompted for a clustering method.")
  }
  
  fcs.cluster.n <- length(unique(phenograph.event.cluster))
  phenograph.cluster.rank <- 1 + fcs.cluster.n - 
    rank( table( phenograph.event.cluster ), ties.method = "last" )
  phenograph.event.cluster <- phenograph.cluster.rank[ phenograph.event.cluster ]
  names( phenograph.event.cluster ) <- NULL
  
  fcs.cluster <- sprintf( "%02d", 1 : fcs.cluster.n )
  
  fcs.cluster.label <- sprintf( "Cluster-%s", fcs.cluster )
  names( fcs.cluster.label ) <- fcs.cluster
  
  fcs.cluster.color <- rep( 
    fcs.color.pool, 
    ceiling( fcs.cluster.n / fcs.color.pool.n ) 
  )[ 1 : fcs.cluster.n ]
  names( fcs.cluster.color ) <- fcs.cluster
  
  fcs.cluster.line.type <- rep( 
    fcs.line.type.pool, 
    ceiling( fcs.cluster.n / fcs.line.type.pool.n ) 
  )[ 1 : fcs.cluster.n ]
  names( fcs.cluster.line.type ) <- fcs.cluster
  
  # set clusters as a factor
  phenograph.event.cluster <- factor( phenograph.event.cluster, 
                                      levels = 1 : fcs.cluster.n )
  levels( phenograph.event.cluster ) <- fcs.cluster
  
  # reorder events
  dmrd.event.cluster <- phenograph.event.cluster
  
  dmrd.event.cluster.n <- as.vector( table( dmrd.event.cluster ) )
  names( dmrd.event.cluster.n ) <- fcs.cluster
  
  length( dmrd.event.cluster )
  table( dmrd.event.cluster )
}



## cluster naming------------

cat(
"\n
Next, the script will try to identify and name the cell types present in every cluster.
For this to work best, you'll need to enter three pieces of information:
1) Whether you're using human or mouse cells
2) Which tissues you're using
3) If you've pre-selected only certain cell types, which cells those are.
\n"
)

cat(
"Please consider the automated naming as a guide, and review the names by checking the heatmaps
and density plots for the clusters. If you don't see the correct cell types being identified,  
check the instructions for cluster naming and add your cell type definitions to the database 
spreadsheet.
\n"
)

# select correct databases based on species---------

species.options <- c("Mouse", "Human")

species.used <- menu( c("Mouse", "Human"), 
                      title = "Please select the species your cells come from. ")

cell.type.filename <- paste0(species.options[species.used], "_celltype_database.xlsx")

cell.database <- read_xlsx( file.path(fcs.src.dir, cell.type.filename) )

# restrict cell types based on tissues used---------

cat("\n
To identify tissue-specific cell types, tell the script which tissue sources
you've used. Immune cell types are always included automatically, since immune
cells are found in every tissue, so you only need to select any additional,
tissue-specific sources here.")
tissue.options <- unique(cell.database$Tissue.restricted)
tissue.type <- multi.menu( tissue.options, title = "Select the tissue or tissues your cells come from. ")

tissue.type <- unique(cell.database$Tissue.restricted)[tissue.type]

# combine whatever was selected with "Immune", rather than requiring the user
# to remember to select Immune every time alongside other tissues. If nothing
# was selected, fall back to using every tissue in the database.
if ( length(tissue.type) == 0 ){
  tissue.type <- unique(cell.database$Tissue.restricted)
} else {
  tissue.type <- union("Immune", tissue.type)
}

cell.database <- dplyr::filter(cell.database, Tissue.restricted %in% tissue.type)


# restrict cell ID based on known input-------
selected.cell.type <- multi.menu( unique(cell.database$Parent.cell.type), 
                                  title = "If you have pre-gated on specific cell types, please them here.
                       If you're using all viable cells, select 1.")

selected.cell.type <- unique(cell.database$Parent.cell.type)[selected.cell.type]


# match markers to database----------
source( file.path( fcs.src.dir, "match_markers.R" ) )

marker.database <- read.csv( file.path( fcs.src.dir, "marker_database.csv" ) )

new.marker.names <- match.markers( fcs.channel.label, marker.database )
new.marker.names <- unname( new.marker.names[ fcs.channel.label ] )  # keep original channel order

# fall back to the original channel label wherever no match was found, rather than
# dropping/failing on that channel
unmatched <- new.marker.names == ""
if ( any( unmatched ) ) {
  cat( "No marker match found for the following channel(s) — keeping the original label:\n" )
  print( fcs.channel.label[ unmatched ] )
  new.marker.names[ unmatched ] <- fcs.channel.label[ unmatched ]
}

match.n <- sapply(new.marker.names, length)

if (any(match.n == 0)) {
  cat("These channel labels could not be matched to any known marker name:\n")
  print(fcs.channel.label[match.n == 0])
  stop("Rename these channels to a recognized marker name, or add them as an alias in the marker names csv.")
}
if (any(match.n > 1)) {
  cat("These channel labels matched more than one marker name (ambiguous):\n")
  print(fcs.channel.label[match.n > 1])
}

new.marker.names <- unlist(new.marker.names)


# scale data for cluster matching--------
flow.data.cluster.median <- apply( dmrd.data, 2, tapply, 
                                   dmrd.event.cluster, median )

# z-score each marker across clusters
safe.scale <- function(x){
  marker.mean <- colMeans(x, na.rm = TRUE)
  marker.sd <- apply(x, 2, sd, na.rm = TRUE)
  marker.sd[ is.na(marker.sd) | marker.sd == 0 ] <- 1
  sweep( sweep(x, 2, marker.mean, "-"), 2, marker.sd, "/" )
}

scaled.fcs.cluster.median <- safe.scale(flow.data.cluster.median)

colnames(scaled.fcs.cluster.median) <- new.marker.names

# prepare marker lists and match clusters-----------
source( file.path( fcs.src.dir, "prepare_marker_lists.R" ) )
source( file.path( fcs.src.dir, "flow_cluster_id_score.R" ) )

marker.list <- prepare_marker_lists( file.path( fcs.src.dir, cell.type.filename ), 
                                             tissue.type, selected.cell.type )

auto.fcs.cluster.label <- fcs.cluster.label


# match clusters to cell type database
scaled.id.score <- flow_cluster_id_score(t(scaled.fcs.cluster.median),
                                         marker_pos = marker.list$markers_positive, 
                                         marker_neg = marker.list$markers_negative,
                                         marker_pos_weight = marker.list$markers_positive_weight,
                                         marker_neg_weight = marker.list$markers_negative_weight )

# assign each cluster its best-matching cell type, and separately record how
# confident that call is: the gap between the top and second-best score.
# Higher values mean more confident (the top cell type won by a wider
# margin); a low value means two (or more) cell types matched the cluster
# almost equally well, so the automated label is more likely to be wrong and
# worth checking by eye against the density plots / cluster ID heatmaps.
cluster.id.confidence <- rep(NA_real_, ncol(scaled.id.score))

for (cluster in 1:ncol(scaled.id.score)){
  ranked.scores <- sort(scaled.id.score[,cluster], decreasing = TRUE)
  auto.fcs.cluster.label[cluster] <- names(ranked.scores)[1]
  cluster.id.confidence[cluster] <- if (length(ranked.scores) > 1){
    ranked.scores[1] - ranked.scores[2]
  } else {
    Inf
  }
}

# plot heatmaps of cluster ID scores--------
# flag clusters whose top two candidate cell types were close in score, i.e.
# had a low confidence value. This threshold is a rule of thumb, not a
# statistical cutoff - it exists to draw attention to the clusters most
# worth double-checking manually, not to replace that check.
fcs.cluster.id.low.confidence.threshold <- 1

low.confidence.clusters <- which(cluster.id.confidence < fcs.cluster.id.low.confidence.threshold)

cat("\n")
if (length(low.confidence.clusters) > 0){
  cat("The following clusters had a close call between their top two matching cell types.
The automated label is still the best match found, but it's worth checking these
against the cluster ID heatmaps and density plots before finalizing names:\n")
  for (cluster in low.confidence.clusters){
    ranked.scores <- sort(scaled.id.score[,cluster], decreasing = TRUE)
    cat(sprintf("  Cluster %s: %s (score %.2f) vs %s (score %.2f)\n",
                colnames(scaled.id.score)[cluster],
                names(ranked.scores)[1], ranked.scores[1],
                names(ranked.scores)[2], ranked.scores[2]))
  }
} else {
  cat("All clusters had a clear best-matching cell type (no close calls).\n")
}
cat("\n")

# save the confidence score alongside the automated labels for later
# reference. Higher confidence means more likely to be correct.
cluster.id.confidence.table <- data.frame(
  cluster = colnames(scaled.id.score),
  auto.label = auto.fcs.cluster.label,
  confidence = cluster.id.confidence
)

write.csv(cluster.id.confidence.table,
          file = file.path(fcs.cluster.id.figure.dir, "cluster_id_confidence.csv"),
          row.names = FALSE)

# also save the full matching score matrix (every candidate cell type against
# every cluster), not just the winning label, so the underlying numbers can
# be reviewed alongside the heatmaps below
write.csv(scaled.id.score,
          file = file.path(fcs.cluster.id.figure.dir, "cluster_id_score_matrix.csv"))

# heatmap with the underlying score for each cell type / cluster combination
# printed on top of each cell, in the original (unclustered) row/column order
plot_score_heatmap_with_values <- function(score.matrix, file.name, width = 1000, height = 1000){
  n.row <- nrow(score.matrix)
  n.col <- ncol(score.matrix)
  
  # image() draws row 1 at the bottom of the plot, so flip the row order
  # here to keep row 1 (first cell type) at the top, as in a normal table
  plot.matrix <- score.matrix[n.row:1, , drop = FALSE]
  
  jpeg( file.name, width = width, height = height )
  par( mar = c(8, 14, 2, 4) )
  
  image( x = 1:n.col, y = 1:n.row, z = t(plot.matrix),
         col = fcs.heatmap.palette, axes = FALSE, xlab = "", ylab = "" )
  
  axis( 1, at = 1:n.col, labels = colnames(score.matrix), las = 2, cex.axis = 0.9 )
  axis( 2, at = 1:n.row, labels = rownames(plot.matrix), las = 2, cex.axis = 0.9 )
  
  for (row in 1:n.row){
    for (col in 1:n.col){
      text( col, row, sprintf("%.2f", plot.matrix[row, col]), cex = 0.7 )
    }
  }
  
  box()
  dev.off()
}

plot_score_heatmap_with_values(scaled.id.score,
                               file.path(fcs.cluster.id.figure.dir, "cluster_id_heatmap.jpg"))

jpeg( file.path(fcs.cluster.id.figure.dir, "cluster_id_heatmap_dendro.jpg"), 
      width = 1000, height = 1000 )
heatmap(scaled.id.score, scale = "none",
        margins = c(5,10),
        xlab = "Clusters", ylab = "matching cell types")
dev.off()


# differentiate similar clusters by variable markers----------------
rownames(flow.data.cluster.median) <- auto.fcs.cluster.label
colnames(flow.data.cluster.median) <- fcs.channel.label

# find (positions of) clusters with non-unique names
non.unique.cluster.names <- auto.fcs.cluster.label[duplicated(auto.fcs.cluster.label)]
# sort into groups
non.unique.groups <- unique(non.unique.cluster.names)

# append high variance marker labels to all clusters until unique
for (cluster.group in non.unique.groups) {
  
  # collect clusters needing renaming
  data.temp <- flow.data.cluster.median[ grepl( cluster.group, rownames(flow.data.cluster.median)),]
  
  # find most variable channels
  variances <- apply(data.temp, 2, var)
  sorted.variances <- sort(variances, decreasing = TRUE)
  
  # find min and max (neg and pos) for each channel
  max.expression <- apply(data.temp, 2, max)
  min.expression <- apply(data.temp, 2, min)
  
  for( cluster in 1:nrow(data.temp)){
    # get position of cluster in full list
    position.in.cluster.list <- rownames(flow.data.cluster.median)[
      grepl(cluster.group, rownames(flow.data.cluster.median))][cluster]
    
    # append markers to name for each positive variable up to the number of clusters in the group or max 4
    markers.to.append.n <- ifelse( nrow(data.temp) < 5, nrow(data.temp), 4 )
    
    for( marker in names(sorted.variances)[1:markers.to.append.n]){
      rownames(flow.data.cluster.median)[as.numeric(names(position.in.cluster.list))] <- 
        if ( abs(flow.data.cluster.median[as.numeric(names(position.in.cluster.list)),marker] - max.expression[marker] ) < abs(flow.data.cluster.median[as.numeric(names(position.in.cluster.list)),marker] - min.expression[marker] )){
          paste( rownames(flow.data.cluster.median)[as.numeric(names(position.in.cluster.list))], marker )
        } else{
          rownames(flow.data.cluster.median)[as.numeric(names(position.in.cluster.list))] <- rownames(flow.data.cluster.median)[as.numeric(names(position.in.cluster.list))]
        }
    }
  }
}

# set names as fcs.cluster.label
fcs.cluster.label <- rownames(flow.data.cluster.median)

cat("
    Plotting histograms for each marker for samples...\n")
source( file.path( fcs.src.dir, "plot_data_histograms.R") )
cat("
    Plotting histograms for each marker for clusters...\n")
source( file.path( fcs.src.dir, "plot_cluster_histograms.R") ) 
cat("\n")

# give option to rename clusters--------
happy.with.cluster.names <- 0

while( happy.with.cluster.names !=1){
  clusters.to.rename <- multi.menu( fcs.cluster.label, title = "Do you want to rename any clusters? ")
  clusters.to.rename <- as.numeric( clusters.to.rename )
  
  if (length(clusters.to.rename) > 1){
    for(cluster in clusters.to.rename){
      fcs.cluster.label[cluster] <- readline(paste0("Please enter a new name for the cluster currently called ", fcs.cluster.label[cluster], ": ") )
    }
    
    cat("\n")
    cat(
      "Your clusters will be named as follows: 
      \n")
    cat( paste0(fcs.cluster.label, collapse = "\n"))
    cat( "\n" )
    cat( "\n" )
    
    happy.with.cluster.names <- menu(c("Yes", "No"), title = "Are you happy with the cluster names? ")
    
  } else if (clusters.to.rename == 0){
      happy.with.cluster.names <- 1
      
      cat("\n")
      cat(
        "Your clusters will be named as follows: 
      \n")
      cat( paste0(fcs.cluster.label, collapse = "\n"))
      
    } else{
      for(cluster in clusters.to.rename){
      fcs.cluster.label[cluster] <- readline(paste0("Please enter a new name for the cluster currently called ", fcs.cluster.label[cluster], ": ") )
      }
      
      cat("\n")
      cat(
        "Your clusters will be named as follows: 
        \n")
      cat( paste0(fcs.cluster.label, collapse = "\n"))
      cat( "\n" )
      cat( "\n" )

      happy.with.cluster.names <- menu(c("Yes", "No"), title = "Are you happy with the cluster names? ")
    }

}



# save cluster counts---------------
cat("
    Exporting cluster counts and percentages as spreadsheets...\n")

flow.cluster.count <- sapply( fcs.cluster, function( fc )
  table( dmrd.event.sample[ dmrd.event.cluster == fc ] ) )

stopifnot( rownames( flow.cluster.count ) == flow.sample )
stopifnot( colnames( flow.cluster.count ) == fcs.cluster )

rownames( flow.cluster.count ) <- flow.sample.label
colnames( flow.cluster.count ) <- fcs.cluster.label

write.csv( flow.cluster.count, file = file.path( fcs.cluster.table.dir, 
                                                 sprintf( "%s.csv", fcs.cluster.table.counts ) ) )

flow.cluster.percent <- flow.cluster.count/rowSums(flow.cluster.count)*100

write.csv( flow.cluster.percent, file = file.path( fcs.cluster.table.dir, 
                                                   sprintf( "%s.csv", "cluster_proportions" ) ) )

# plot density distributions of clusters with new names---------------
# plot density distributions by cluster
if( length(clusters.to.rename) > 1){
  cat("
    Plotting histograms for each marker for clusters...\n")
  source( file.path( fcs.src.dir, "plot_cluster_histograms.R") ) 
} else if (clusters.to.rename != 0){
  cat("
    Plotting histograms for each marker for clusters...\n")
  source( file.path( fcs.src.dir, "plot_cluster_histograms.R") ) 
}

setup.end.time <- Sys.time()
setup.time <- round(difftime(setup.end.time, setup.start.time, units='mins'), 2)
setup.time <- as.numeric(setup.time)

cat("\n")

cat(
"Proceed to the next section."
)
