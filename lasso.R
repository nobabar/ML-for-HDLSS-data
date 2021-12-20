library(caret, include.only='findCorrelation')
library(tibble, include.only='rowid_to_column')
library(dplyr)
library(GSelection)

# import data ----

data_placenta <- read.csv("data/membrane.placentaire.tsv", sep="\t")
# data_placenta$CI2 <- factor(data_placenta$CI2)

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

train_test = function(data, proportion){
  data <- data %>% rowid_to_column("rowid")
  train_set <- slice_sample(data, prop = proportion)
  test_set <- data %>% anti_join(as.data.frame(train_set), by = "rowid")
  return(list(train=train_set[,-1], test=test_set[,-1]))
}

sets <- train_test(data_placenta, .70)

# ça marche pô
spam_var <- spam.var.rcv(sets$train[,-1], sets$train$CI2,
                         length(sets$train[,-1]))
spam_var <- rcv(sets$train[,-1], sets$train$CI2,
                1, length(sets$train[,-1]), method="spam")

hsic_var <- hsic.var.rcv(sets$train[,-1], sets$train$CI2, length(sets$train[,-1]))

fit <- feature.selection(sets$train[,-1], sets$train$CI2, length(sets$train[,-1]))

preds <- genomic.prediction(sets$test[,-1], spam_var, hsic_var,
                            fit$spam_selected_feature_index,
                            fit$hsic_selected_feature_index,
                            fit$coefficient.spam,
                            fit$coefficient.hsic)