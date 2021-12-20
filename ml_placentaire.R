# data juggling ----
library(caret, include.only=c('findCorrelation', 'confusionMatrix'))
library(tibble, include.only='rowid_to_column')
library(dplyr)

# neuralnet ----
# library(caret, include.only='confusionMatrix')
library(neuralnet)
# for plotting neuralnet
# library(gamlss.add)

# PCA ----
library(ade4)
library(adegraphics)

# tests ----
library(mltools)
library(pROC)
library(scatterplot3d)

# seed ----
# force seed for all functions
addTaskCallback(function(...) {set.seed(2);TRUE})


# import data ----

data_placenta <- read.csv("data/membrane.placentaire.tsv", sep="\t")
data_placenta$CI2 <- factor(data_placenta$CI2)

# filter data ----

varSupp <- which(apply(data_placenta, 2, var)== 0)
data_placenta <- data_placenta[, -varSupp]

cor_matrix <- cor(data_placenta[,-1])
cor_col <- findCorrelation(cor_matrix, cutoff=0.9)+1
data_placenta <- data_placenta[, -cor_col]

# center reduce ----

# PCA is sensitive to input scaling, so we will scale the data first
data_placenta[,-1] <- as.data.frame(scale(data_placenta[,-1],
                                          center=TRUE,
                                          scale=TRUE))

# export data ----

# write.csv(data_placenta, "data/membrane_filtered.csv", quote=FALSE)

# data_placenta <- read.csv("data/membrane_filtered.csv", header=TRUE)
# rownames(data_placenta) <- data_placenta[,1]
# data_placenta <- data_placenta[,-1]

# separate into train and test ----

train_test = function(data, proportion){
  data <- data %>% rowid_to_column("rowid")
  train_set <- slice_sample(data, prop = proportion)
  test_set <- data %>% anti_join(as.data.frame(train_set), by = "rowid")
  return(list(train=train_set[,-1], test=test_set[,-1]))
}

# create neural network ----

neural_network = function(data, nneurones){
  nn <- neuralnet(CI2 ~ ., data=data$train,
                  hidden=nneurones,
                  stepmax=1e+09)
  predictions <- c(0, 1)[apply(predict(nn, data$test), 1, which.max)]
  return(predictions)
}

# compute confusion matrix ----

confusion_matrix = function(actuals, preds){
  return(table(preds, actuals))
}

complex_cm = function(actuals, preds){
  return(confusionMatrix(as.factor(preds), actuals))
}

# compute accuracy and Matthew coeff ----

perc_match = function(cmatrix){
  return (sum(diag(cmatrix)) / sum(cmatrix) * 100)
}

acc_mc = function(cmatrix){
  acc <- perc_match(cmatrix)
  mc <- mcc(confusionM=matrix(cmatrix, nrow(cmatrix)))
  return(list(acc=acc, mc=mc))
}

# PCA ----

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

# apply PCA result ----

comp_pca = function(data){
  # get PCA results
  res_pca <- PCA(data$train)
  
  # modify train data by PCA output
  data$train <- cbind(data$train$CI2, res_pca$l1)
  names(data$train) <- c("CI2", paste0("C", seq(length(res_pca$l1))))
  
  # cross product with pca transformation matrix
  test_trans <- as.matrix(data$test[-1]) %*% as.matrix(res_pca$c1)
  
  # modify test data by cross product results
  data$test <- cbind(CI2=data$test$CI2, as.data.frame(scale(test_trans)))
  names(data$test) <- c("CI2", paste0("C", seq(length(res_pca$l1))))
  
  return(data)
}

# Results ----

sets <- train_test(data_placenta, .70)

sets <- comp_pca(sets)

preds <- neural_network(sets, 15)

# Optimize number of neurons ----

opt_2hidden = function(data, start=1, stop, step=1){
  acc_list <- c()
  n1_list <- c()
  n2_list <- c()
  for (n1 in seq(start, stop, step)){
    for (n2 in seq(1, (n1+1)/2)){
      preds <- neural_network(data, c(n1, n2))
      cm <- confusion_matrix(data$test$CI2, preds)  
      acc <- acc_mc(cm)$acc
      
      acc_list <- c(acc_list, acc)
      n1_list <- c(n1_list, n1)
      n2_list <- c(n2_list, n2)
      
      # print(paste("n1: ", n1, "    n2: ", n2, "\n acc: ", acc))
    }
  }
  return(list(acc=acc_list, n1=n1_list, n2=n2_list))
}

opt_1hidden = function(data, start=1, stop, step=1){
  acc_list <- c()
  n_list <- c()
  for (n in seq(start, stop, step)){
    preds <- neural_network(data, n)
    cm <- confusion_matrix(data$test$CI2, preds)
    acc <- acc_mc(cm)$acc

    acc_list <- c(acc_list, acc)
    n_list <- c(n_list, n)
    
    print(paste("n: ", n, "    acc: ", acc))
  }
  return(list(acc=acc_list, n=n_list))
}

optim_res <- opt_2hidden(sets, 1, ncol(sets$train)*2)

x <- optim_res$n1
y <- optim_res$n2
z <- optim_res$acc

# plot n neurones ----

# plot(res$n, res$acc)
# scatterplot3d(z, x, y)

fit <- lm(z~x+y)

grid.lines = 40
x.pred <- seq(min(x), max(x), length.out = grid.lines)
y.pred <- seq(min(y), max(y), length.out = grid.lines)
xy <- expand.grid(x = x.pred, y = y.pred)
z.pred <- matrix(predict(fit, newdata = xy), 
                 nrow = grid.lines, ncol = grid.lines)

fitpoints <- predict(fit)
m <- matrix(fitpoints, nrow = length(x), ncol = length(y))

# scatter plot with regression plane

# library(plot3D)
# scatter3D(x, y, z, pch = 19, cex = 1,colvar = NULL, col="red3", 
#           theta = 80, phi = 20, bty="b",
#           xlab = "n1", ylab = "n2", zlab = "ACC",  
#           surf = list(x = x.pred, y = y.pred, z = z.pred, facets = TRUE,
#                       col=ramp.col(col = c("dodgerblue3","seagreen2"),
#                                    n = 300, alpha=0.9),
#                       border="black"),
#           main = "Accuracy of model and number of hidden neurones")


library(plotly)
plot_ly() %>%
  add_trace(x = x, 
            y = y,
            z = z, 
            type = "scatter3d", 
            mode = "markers",
            marker = list(size = 5, color = 'rgb(17, 157, 255)',
                          line = list(color = 'rgb(0, 0, 0)',
                                      width = 2))) %>%
  add_surface(x = x.pred,
              y = y.pred,
              z = z.pred,
              type = "surface",
              colorscale = "Viridis")

# compute with correct number of neurones ----

# preds <- neural_network(sets, c(optim_res$n1[which.max(optim_res$acc)], 
#                                 optim_res$n2[which.max(optim_res$acc)]))

preds <- neural_network(sets, c(68,2))

cm <- confusion_matrix(sets$test$CI2, preds)
results <- acc_mc(cm)

# ROC curve ----

plot(roc(sets$test$CI2, preds, direction="<"), print.auc=TRUE)
