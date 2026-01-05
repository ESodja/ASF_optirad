# Calculates the matrix for seed dispersal based on the 2Dt dispersal kernel
## adapted for culling search function
# TwoDt.matrix <- function(cells, scale, shape){
#     cells <- (cells*5)+1
#     x = 1:cells
#     y = 1:cells
#     center_x = ceiling(cells/2)-0.5
#     center_y = ceiling(cells/2)-0.5
#     p_cell = matrix(rep(0, cells^2), nrow=cells)
#
#     TwoDt <- function(x, y){
#         dist <- sqrt((x-center_x)^2+(y-center_y)^2)
#         ((shape-1)/(pi * (scale^2)))*((1+((dist^2)/(scale^2)))^-shape)
#     }
#
#     for (i in x){
#         for (j in y){
#             p_cell[j,i] <- dblquad(TwoDt, i, i+1, j, j+1)
#         }
#     }
#
# #     as.big.matrix(p_cell,
# #                   backingfile = paste0('twodtmatrix', cells_x, 'x', cells_y, '_', scale, '_', shape),
# #                   backingpath = paste0(root,'/dm/'),
# #                   descriptorfile = paste0('twodtmatrixdesc', cells_x, 'x', cells_y, '_', scale, '_', shape))
#     return(p_cell)
# }


    TwoDt <- function(dist, shape, scale){
        ((shape-1)/(pi * (scale^2)))*((1+((dist^2)/(scale^2)))^-shape)
    }
