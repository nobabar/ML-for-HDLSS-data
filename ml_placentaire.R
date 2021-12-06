library(mltools)
library(pROC)
library(scatterplot3d)

################################################################################
#                                   Results                                    #
################################################################################

sets = train_test(data_placenta, .30)

preds <- neural_network(sets, c(180, 90))

cm <- confusion_matrix(sets$test$CI2, preds)
results <- acc_mc(cm)

################################################################################
#                                   ROC curve                                  #
################################################################################

plot(roc(sets$test$CI2, preds, direction="<"), print.auc=TRUE)

################################################################################
#                           Optimize number of neurons                         #
################################################################################

opt_2hidden = function(start=1, stop, step=1, nrep=5){
  acc_list = c()
  n1_list = c()
  n2_list = c()
  for (n1 in seq(data, start, step, stop)){
    for (n2 in seq(1, (n1+1)/2)){
      acc=0
      for (iter in seq(nrep)){
        preds <- neural_network(data, c(n1, n2))
        cm <- confusion_matrix(data$test$CI2, preds)
        acc = acc + acc_mc(cm)$acc
      }
      acc = acc/nrep
      acc_list = c(acc_list, acc)
      n1_list = c(n1_list, n1)
      n2_list = c(n2_list, n2)
      
      print(paste("n1: ", n1, "\t n2: ", n2, "\n acc: ", acc))
    }
  }
  return(list(acc=acc_list, n1=n1_list, n2=n2_list))
}


opt_1hidden = function(data, start=1, stop, step=1, nrep=5){
  acc_list = c()
  n_list = c()
  for (n in seq(start, stop, step)){
    acc=0
    for (iter in seq(1, nrep)){
      preds <- neural_network(data, n)
      cm <- confusion_matrix(data$test$CI2, preds)
      acc = acc + acc_mc(cm)$acc
    }
    acc = acc/nrep
    acc_list = c(acc_list, acc)
    n_list = c(n_list, n)
    
    print(paste("n: ", n, "\t acc: ", acc))
  }
  return(list(acc=acc_list, n=n_list))
}

################################################################################
#                                plot n neurones                               #
################################################################################

# scatterplot3d(acc_list_2, n1_list, n2_list)

# plot(acc_list_1, n_list)