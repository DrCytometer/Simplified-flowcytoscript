# Adapted from scType
# GNU General Public License v3.0 (https://github.com/IanevskiAleksandr/sc-type/blob/master/LICENSE)
# Written by Aleksandr Ianevski <aleksandr.ianevski@helsinki.fi>, June 2021
#
# edited by Oliver Burton <ob240@cam.ac.uk>, August 2023
#
# Functions on this page:
# flow_cluster_id_score: calculate cluster ID scores and assign cell types
#
# @params: cluster.input.data - input transposed thresholded flow expression matrix (rownames - markers, column names - clusters),
#
# @params: marker_pos - list of markers positively expressed in the cell type 
# @params: marker_neg - list of markers that should not be expressed in the cell type (NULL if not applicable)
# @params: marker_pos_weight - optional list (same names/structure as marker_pos) of named numeric
#                              vectors giving a specificity weight per marker (marker name -> weight).
#                              Markers without an entry, or if the whole argument is NULL, default to
#                              weight 1 (i.e. behaves exactly as the original unweighted scoring).
# @params: marker_neg_weight - as above, for marker_neg

flow_cluster_id_score <- function(cluster.input.data, marker_pos, marker_neg = NULL,
                                   marker_pos_weight = NULL, marker_neg_weight = NULL, ...){
  
  # subset to markers found in the data
  names_mkp_cp <- names(marker_pos)
  names_mkn_cp <- names(marker_neg)
  
  marker_pos <- lapply(1:length(marker_pos), function(x){ 
    MarkerToKeep = rownames(cluster.input.data) %in% as.character(marker_pos[[x]])
    rownames(cluster.input.data)[MarkerToKeep]})
  
  marker_neg = lapply(1:length(marker_neg), function(x){ 
    MarkerToKeep = rownames(cluster.input.data) %in% as.character(marker_neg[[x]])
    rownames(cluster.input.data)[MarkerToKeep]})
  
  names(marker_pos) <- names_mkp_cp
  names(marker_neg) <- names_mkn_cp
  
  # align a weight vector to a (post-filtering) marker vector, defaulting any
  # marker missing an explicit weight (or the whole weight list being absent) to 1
  align_weights <- function(markers_by_type, weights_by_type, names_cp){
    out <- lapply(names_cp, function(ct){
      mk <- markers_by_type[[ct]]
      if (length(mk) == 0) return(numeric(0))
      w_supplied <- if (!is.null(weights_by_type)) weights_by_type[[ct]] else NULL
      w <- setNames(rep(1, length(mk)), mk)
      if (!is.null(w_supplied)){
        matched <- intersect(names(w_supplied), mk)
        w[matched] <- w_supplied[matched]
      }
      w
    })
    names(out) <- names_cp
    out
  }
  
  marker_pos_weight <- align_weights(marker_pos, marker_pos_weight, names_mkp_cp)
  marker_neg_weight <- align_weights(marker_neg, marker_neg_weight, names_mkn_cp)
  
  # subselect only with marker genes
  cluster.input.data = cluster.input.data[unique(c(unlist(marker_pos),unlist(marker_neg))), ]
  
  # deal with cases with only one marker
  if (is.null(ncol(cluster.input.data))){
    cluster.input.data <- data.frame(matrix(cluster.input.data,1))
    row.names(cluster.input.data) <- unique(c(unlist(marker_pos),unlist(marker_neg)))
  }
  
  # combine scores
  combined.score = do.call("rbind", lapply(names(marker_pos), function(mkrs){ 
    sapply(1:ncol(cluster.input.data), function(x) {
      pos_mk = marker_pos[[mkrs]]
      neg_mk = marker_neg[[mkrs]]
      w_pos = marker_pos_weight[[mkrs]]
      w_neg = marker_neg_weight[[mkrs]]
      
      mk_pos = cluster.input.data[pos_mk, x]
      mk_neg = cluster.input.data[neg_mk, x] * -1
      
      # weighted sum, normalised by sqrt(sum(weight^2)) so that markers with
      # weight 1 reduce exactly to the original sqrt(n) normalisation; markers
      # down-weighted for being non-specific (shared by many candidate cell
      # types) contribute proportionally less to both the sum and the
      # normalisation
      sum_t1 = sum(w_pos * mk_pos) / sqrt(sum(w_pos^2))
      sum_t2 = sum(w_neg * mk_neg) / sqrt(sum(w_neg^2))
      
      # both sums can legitimately come out as NaN (0/0) when a cell type has
      # no positive, or no negative, markers present in the panel at all -
      # treat that side of the score as contributing nothing, symmetrically
      if(is.na(sum_t1) || is.nan(sum_t1)){
        sum_t1 = 0;
      }
      if(is.na(sum_t2) || is.nan(sum_t2)){
        sum_t2 = 0;
      }
      sum_t1 + sum_t2
    })
  })) 
  
  
  dimnames(combined.score) = list(names(marker_pos), colnames(cluster.input.data))
  combined.score.max <- combined.score[!apply(is.na(combined.score) | combined.score == "", 1, all),] # remove na rows
  
  combined.score.max
}
