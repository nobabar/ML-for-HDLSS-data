library(ade4)
library(adegraphics)

################################################################################
#                                 center reduce                                #
################################################################################

# PCA is sensitive to input scaling, so we will scale the data first
data_placenta[,-1] <- as.data.frame(scale(data_placenta[,-1],
                                          center=TRUE,
                                          scale=TRUE))

################################################################################
#                                      PCA                                     #
################################################################################

PCA = function(data){
  # start a mock pca
  mock_pca <- dudi.pca(data[-1],
                       scannf=FALSE,
                       nf=2)
  # observe the eingen values to determine optimal dimension reduction
  inertia <- inertia.dudi(mock_pca)
  
  # we can determine the optimal number of dimension graphically
  # plot(inertia$tot.inertia$`cum(%)`)
  
  # but here we will use the condition of keeping 90% of the information
  ndim <- min(which(inertia$tot.inertia$`cum(%)`>90))
  
  # perform pca once again, with correct number of final dimensions
  final_pca <- dudi.pca(data[-1],
                        scannf=FALSE,
                        nf=ndim)
  return(final_pca)
}

################################################################################
#                               apply PCA result                               #
################################################################################

comp_pca = function(data){
  # get PCA results
  res_pca <- PCA(data$train)
  
  # modify train data by PCA output
  data$train[,-1] <- res_pca$l1
  
  
  # cross product with pca transformation matrix
  test_trans <- as.matrix(data$test[-1]) %*% as.matrix(pcaAll$c1)
  
  # modify test data by cross product results
  data$test[,-1] <- as.data.frame(scale(test_trans))
  
  return(data)
}
