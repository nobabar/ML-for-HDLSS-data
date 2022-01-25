# libraries ----
# data juggling
library(caret, include.only=c('findCorrelation', 
                              'findLinearCombos', 
                              'nearZeroVar'))
library(tibble, include.only='rowid_to_column')
library(dplyr)

# neuralnet
# library(caret, include.only='confusionMatrix')
library(neuralnet)
# for plotting neuralnet
# library(gamlss.add)

# PCA
library(ade4)
library(adegraphics)
library(factoextra)

# tests
library(mltools)
library(pROC)

#plots
library(viridisLite)
library(plot3D)
library(scatterplot3d)
library(plotly)

# seed ----
# force seed for all functions
addTaskCallback(function(...) {set.seed(2);TRUE})


# import data ----

data_placenta <- read.csv("data/membrane.placentaire.tsv", sep="\t")
data_placenta$CI2 <- factor(data_placenta$CI2)

# filter data ----

varSupp <- caret::nearZeroVar(data_placenta)
data_placenta <- data_placenta[, -varSupp]

cor_col <- findCorrelation(cor(data_placenta[,-1]), cutoff=0.8)
data_placenta <- data_placenta[, -cor_col+1]

col_col <- findLinearCombos(data_placenta[,-1])
data_placenta <- data_placenta[, -col_col$remove+1]

# center reduce ----

par(mfrow=c(1, 2))
boxplot(data_placenta[,-1])

# PCA is sensitive to input scaling, so we will scale the data first
data_placenta[,-1] <- as.data.frame(scale(data_placenta[,-1],
                                          center=TRUE,
                                          scale=TRUE))
boxplot(data_placenta[,-1])
dev.off()

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
                  stepmax=1e+06, rep=3)
  predictions <- c(0, 1)[apply(predict(nn, data$test), 1, which.max)]
  return(predictions)
}

# compute confusion matrix ----

confusion_matrix = function(actuals, preds){
  return(table(preds, actuals))
}

# compute accuracy ----

perc_match = function(cmatrix){
  return (sum(diag(cmatrix)) / sum(cmatrix) * 100)
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
  plot(inertia$tot.inertia$`cum(%)`)
  
  # but here we will use the condition of keeping 80% of the information
  ndim <- min(which(inertia$tot.inertia$`cum(%)`>80))
  abline(v=ndim)
  
  # perform pca once again, with correct number of final dimensions
  final_pca <- dudi.pca(data[-1],
                        scannf=FALSE,
                        nf=ndim)
  
  # graph the variables
  fviz_pca_var(final_pca,
               col.var = "contrib", # Color by contributions to the PC
               gradient.cols = viridis(100)
  )
  
  fviz_pca_ind(final_pca,
               col.ind = sets$train$CI2, # color by groups
               palette = c("#00AFBB",  "#FC4E07"),
               addEllipses = TRUE, # Concentration ellipses
               ellipse.type = "norm",
               legend.title = "Groups",
               axes = c(1, 2)
  )
  
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

# Optimize number of neurons ----

opt_1hidden = function(data){
  fitControl <- caret::trainControl(method = "repeatedcv", 
                                    number = 10, 
                                    repeats = 5, 
                                    classProbs = TRUE,
                                    summaryFunction = caret::twoClassSummary)

  nnetGrid <-  expand.grid(size = seq(from = 1, to = ncol(sets$train)*2, by = 1),
                           decay = c(1, 0.5, 0.1, 5e-2, 1e-2, 5e-3, 1e-3))

  levels(data$CI2) <- c("No", "Yes")

  nnetFit <- caret::train(CI2 ~ ., 
                          data = data,
                          method = "nnet",
                          metric = "ROC",
                          trControl = fitControl,
                          tuneGrid = nnetGrid,
                          verbose = FALSE)
  return(nnetFit)
}

# fit1 <- opt_1hidden(sets$train)
# plot(fit1, lwd=2)

opt_2hidden = function(data, start=1, stop, step=1){
  auc_list <- c()
  n1_list <- c()
  n2_list <- c()
  for (n1 in seq(start, stop, step)){
    for (n2 in seq(1, (n1+1)/2)){
      preds <- neural_network(data, c(n1, n2))
      auc <- auc(data$test$CI2, preds)
      
      auc_list <- c(auc_list, auc)
      n1_list <- c(n1_list, n1)
      n2_list <- c(n2_list, n2)
      
      # print(paste("n1: ", n1, "    n2: ", n2, "\n acc: ", acc))
    }
  }
  return(list(auc=auc_list, n1=n1_list, n2=n2_list))
}

optim_res <- opt_2hidden(sets, stop=floor(ncol(sets$train)*1.3))

x <- optim_res$n1
y <- optim_res$n2
z <- optim_res$auc

fit <- lm(z~x+y)
fit1 <- lm(z~x)
fit2 <- lm(z~y)

fitpoints1 <- predict(fit1)
fitpoints2 <- predict(fit2)
plot(z~x)
lines(fitpoints1)
plot(z~y)
lines(fitpoints2)

grid.lines = 40
x.pred <- seq(min(x), max(x), length.out = grid.lines)
y.pred <- seq(min(y), max(y), length.out = grid.lines)
xy <- expand.grid(x = x.pred, y = y.pred)
z.pred <- matrix(predict(fit, newdata = xy),
                 nrow = grid.lines, ncol = grid.lines)


scatter3D(x, y, z, pch = 19, cex = 1, colvar = NULL, col="red3",
          theta = 120, phi = 20, bty="b",
          xlab = "n1", ylab = "n2", zlab = "AUC",
          surf = list(x = x.pred, y = y.pred, z = z.pred, facets = TRUE,
                      col=ramp.col(col = viridis(100),
                                   n = 300, alpha=0.9),
                      border="black"),
          main = "Accuracy of model and number of hidden neurones")

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

max_auc  <- which(z.pred == max(z.pred), arr.ind = TRUE)
opt_n1 <- x.pred[max_auc[2]]
opt_n2 <- y.pred[max_auc[1]]

# compute with correct number of neurones ----

preds <- neural_network(sets, c(opt_n1, opt_n2))


cm <- confusion_matrix(sets$test$CI2, preds)
acc <- perc_match(cm)

# ROC curve ----

plot(roc(sets$test$CI2, preds, direction="<"), print.auc=TRUE)
