library(caret, include.only = 'confusionMatrix')
library(neuralnet)

# for plotting neuralnet
# library(gamlss.add)

################################################################################
#                             create neural network                            #
################################################################################

neural_network = function(data, nneurones){
  nn <- neuralnet(CI2 ~ ., data=data$train,
                  hidden=nneurones,
                  stepmax=1e+06)
  preds <- c(0, 1)[apply(predict(nn, data$test), 1, which.max)]
  return(preds)
}

################################################################################
#                           compute confusion matrix                           #
################################################################################

confusion_matrix = function(actuals, preds){
  return(table(preds, actuals))
}

complex_cm = function(actuals, preds){
  return(confusionMatrix(as.factor(preds), actuals))
}

################################################################################
#                      compute accuracy and Matthew coeff                      #
################################################################################

perc_match = function(cmatrix){
  return (sum(diag(cmatrix)) / sum(cmatrix) * 100)
}

acc_mc = function(cmatrix){
  acc <- perc_match(cmatrix)
  mc <- mcc(confusionM=matrix(cmatrix, nrow(cmatrix)))
  return(list(acc=acc, mc=mc))
}

