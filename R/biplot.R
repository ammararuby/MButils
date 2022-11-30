#' @title Flexibly generate a biplot from the output of a PCA
#'
#' @description Given the output of a PCA, makes biplot of data showing the
#'   magnitude of the factor loadings that contribute the plotted PCs.
#'
#' @usage biplot(pca, top.n, data, color, base_only)
#'
#' @param pca The output of `prcomp()` on the data to be ordinated.
#' @param top.n The number of loadings to plot (ranked by descending magnitude).
#'   Default 5
#' @param data A `data.frame` of sample metadata variables row names matching
#'   those in the PCA output.
#' @param color The sample metadata variable (a column in provided `data`) by
#'   which to color plotted points.  Taken as a string (i.e. 'subj')
#' @param base_only Should only the base PCA be returned, without loadings?
#'   Default FALSE
#'
#' @import ggplot2
#'
#' @return A ggplot object
#'
#' @export
#'

biplot <- function(pca, top.n = 5, data, color, base_only = FALSE){
     # TODO: Maybe return an object that has all the loadings

     # Calculate loadings
     V <- pca$rotation # Eigenvectors
     L <- diag(pca$sdev) # Diag mtx w/sqrts of eigenvalues on diag.

     # % variance explained
     eigs <- pca$sdev^2
     ve.pc1 <- as.character(100*round(eigs[1] / sum(eigs), 3))
     ve.pc2 <- as.character(100*round(eigs[2] / sum(eigs), 3))

     loadings <- V %*% L

     # Get loadings for first 2 PCs and format for plotting
     pythag <- function(a, b){sqrt(a^2 + b^2)}
     loadings.12 <-
          data.frame(loadings[, 1:2]) %>%
          rename(PC1 = X1, PC2 = X2) %>%
          mutate(variable = row.names(loadings)) %>%
          mutate(length = pythag(PC1, PC2),
                 slope = PC2/PC1,
                 ang = atan(slope)*(180/pi))

     # Select the top loadings contributing to each PC to plot
     loadings.plot <- top_n(loadings.12, top.n, wt = length)

     # What quadrant of the plot is the label in?
     q1 <- filter(loadings.plot, PC1 > 0 & PC2 > 0)
     q2 <- filter(loadings.plot, PC1 < 0 & PC2 > 0)
     q3 <- filter(loadings.plot, PC1 < 0 & PC2 < 0)
     q4 <- filter(loadings.plot, PC1 > 0 & PC2 < 0)

     # Add back sample data
     pca.df <-
          data.frame(pca$x) %>%
          rownames_to_column(var = 'row')
     pca.df <- left_join(pca.df, rownames_to_column(data, var = 'row'))

     # Figure out axis limits to customize scale of data
     # Find the largest magnitude datapoint in the first 2 PCs, then add 5% to it for
     # plotting "room"
     limit <- max(abs(pca.df[, c('PC1','PC2')])) +
          0.05*(max(abs(pca.df[, c('PC1','PC2')])))

     pca.plot <-
          ggplot(pca.df, aes(PC1, PC2, color = .data[[color]])) +
          geom_point(alpha = 0.5) +
          coord_equal() +
          xlim(-limit, limit) + ylim(-limit, limit) +
          labs(x = paste0(' PC1 (', ve.pc1, '%)'),
               y = paste0(' PC2 (', ve.pc2, '%)')) +
          theme_bw()

     pca.biplot <-
          pca.plot +
          geom_segment(data = loadings.plot,
                       aes(x = 0, y = 0,
                           xend = PC1, yend = PC2),
                       color = 'black',
                       arrow = arrow(angle = 15,
                                     length = unit(0.1, 'inches')))

     # Then add geom_text quadrant-by-quadrant, aligning text accordingly
     if (dim(q1)[1] != 0) {
          pca.biplot <- pca.biplot +
               geom_text(data = q1, aes(x = PC1, y = PC2,
                                        hjust = 0,
                                        angle = ang,
                                        label=paste0('   ', variable),
                                        fontface = 'bold'),
                         color = 'black', show.legend = FALSE)
     }
     if (dim(q2)[1] != 0) {
          pca.biplot <- pca.biplot +
               geom_text(data = q2, aes(x = PC1, y = PC2,
                                        hjust = 1,
                                        angle = ang,
                                        label=paste0(variable, '   '),
                                        fontface = 'bold'),
                         color = 'black', show.legend = FALSE)
     }
     if (dim(q3)[1] != 0) {
          pca.biplot <- pca.biplot +
               geom_text(data = q3, aes(x = PC1, y = PC2,
                                        hjust = 1,
                                        angle = ang,
                                        label=paste0(variable, '   '),
                                        fontface = 'bold'),
                         color = 'black', show.legend = FALSE)
     }
     if (dim(q4)[1] != 0) {
          pca.biplot <- pca.biplot +
               geom_text(data = q4, aes(x = PC1, y = PC2,
                                        hjust = 0,
                                        angle = ang,
                                        label=paste0('   ', variable),
                                        fontface = 'bold'),
                         color = 'black', show.legend = FALSE)
     }

     # Return plot object
     if (base_only){
          pca.plot
          }
     else {
          pca.biplot
          }
}
