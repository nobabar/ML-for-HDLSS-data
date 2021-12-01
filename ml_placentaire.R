library(dplyr)
library(tibble)
library(neuralnet)
library(gamlss.add)
library(mltools)
library(ggplot2)
library(pROC)

################################################################################
#                                  import data                                 #
################################################################################

data_placenta <- read.csv("membrane.placentaire.tsv", sep = "\t")
data_placenta$CI2 = factor(data_placenta$CI2)

################################################################################
#                         separate into train and test                         #
################################################################################

train_test = function(data, proportion){
  data <- data %>% rowid_to_column("rowid")
  train_set <- slice_sample(data, prop = proportion)
  test_set <- data %>% anti_join(as.data.frame(train_set), by = "rowid")
  return(list(train=train_set, test=test_set))
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

sets = train_test(data_placenta, .30)
preds <- neural_network(sets, c(450, 200))

cm <- confusion_matrix(sets$test$CI2, preds)
results <- acc_mc(cm)

################################################################################
#                                   ROC curve                                  #
################################################################################


plot(roc(sets$test$CI2, preds, direction="<"), print.auc=TRUE)

################################################################################
#                                      ACP                                     #
################################################################################

library(ade4)
library(adegraphics)
pcaAll <- dudi.pca(data_placenta[-1], scannf = FALSE, nf = 2)
inertia.dudi(pcaAll)
pcaAll <- dudi.pca(data_placenta[-1], scannf = FALSE, nf = 11)
s.corcircle(pcaAll$co)
