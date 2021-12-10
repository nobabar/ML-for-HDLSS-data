library(smacof)

################################################################################
#                                      MDS                                     #
################################################################################

# scale values
placenta_scale = scale(data_placenta[,-1])

# compute distance matrix
placenta_dist = dist(x = placenta_scale)

# perform multi dimensional scaling to X dimension
placenta_mds = mds(delta = placenta_dist , ndim = 2 , type = "ratio")

################################################################################
#                                Quality control                               #
################################################################################

# mds comes with a stress variable
# we will normalize it in order to be able to compare it

dhat_matrix = as.matrix(placenta_mds$dhat)
d_matrix = as.matrix(placenta_mds$confdist)

denominator = sum(dhat_matrix[upper.tri(dhat_matrix)]^2)

p_ij = dhat_matrix[upper.tri(dhat_matrix)]
d_ij = d_matrix[upper.tri(d_matrix)]
nominator = sum((p_ij - d_ij)^2) 

# compute stress induced by the number of dimensions asked, need to be minimal
normalized_stress = nominator/denominator
