library(dplyr)
library(tibble)
library(neuralnet)
library(gamlss.add)
library(mltools)
library(ggplot2)
library(pROC)
library(ade4)
library(adegraphics)

################################################################################
#                                  import data                                 #
################################################################################

data_placenta <- read.csv("membrane.placentaire.tsv", sep = "\t")
data_placenta$CI2 = factor(data_placenta$CI2)

################################################################################
#                                  filter data                                 #
################################################################################

varSupp <- which(apply(data_placenta, 2, var)== 0)
data_placenta <- data_placenta[, -varSupp]

################################################################################
#                                centrer reduire                               #
################################################################################

data_placenta[,-1] <- as.data.frame(scale(data_placenta[,-1], center=TRUE, scale=TRUE))

# varNan <- unique(as.data.frame(which(is.nan(result), arr.ind=T))$col)
# result <- result[, -varNan]

################################################################################
#                         separate into train and test                         #
################################################################################

train_test = function(data, proportion){
  data <- data %>% rowid_to_column("rowid")
  train_set <- slice_sample(data, prop = proportion)
  test_set <- data %>% anti_join(as.data.frame(train_set), by = "rowid")
  return(list(train=train_set[,-1], test=test_set[,-1]))
}

################################################################################
#                             create neural network                            #
################################################################################

neural_network = function(data, nneurones){
  nn <- neuralnet(CI2 ~ ., data=data$train,
                  hidden=nneurones,
                  stepmax=1e+09)
  preds <- c(0, 1)[apply(predict(nn, data$test), 1, which.max)]
  return(preds)
}

################################################################################
#                           compute confusion matrix                           #
################################################################################

confusion_matrix = function(actuals, preds){
  return(table(preds, actuals))
}

complex_confusion_matrix = function(actuals, preds){
  return(confusionMatrix(as.factor(preds), actuals))
}

################################################################################
#                      compute accuracy and Matthew coeff                      #
################################################################################

perc_match = function(cmatrix){
  return (sum(diag(cmatrix)) / sum(cmatrix) * 100)
}

acc_mc = function(cm){
  acc <- perc_match(cm)
  mc <- mcc(confusionM=matrix(cm, nrow(cm)))
  return(list(acc=acc, mc=mc))
}

################################################################################
#                                      ACP                                     #
################################################################################

sets = train_test(data_placenta, .30)

pcaAll <- dudi.pca(sets$train[-1], scannf = FALSE, nf = 2)
inertia <- inertia.dudi(pcaAll)
plot(inertia$tot.inertia$`cum(%)`)

pcaAll <- dudi.pca(data_placenta[-1], scannf = FALSE, nf = 14)
# s.corcircle(pcaAll$co)

sets$train[,-1] <- pcaAll$l1

sets$test[,-1] <- as.data.frame(scale(as.matrix(sets$test[-1]) %*% as.matrix(pcaAll$c1)))

################################################################################
#                                   Results                                    #
################################################################################

preds <- neural_network(sets, c(8, 4))

cm <- confusion_matrix(sets$test$CI2, preds)
results <- acc_mc(cm)

################################################################################
#                                   ROC curve                                  #
################################################################################


plot(roc(sets$test$CI2, preds, direction="<"), print.auc=TRUE)

