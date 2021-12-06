library(caret, include.only = 'findCorrelation')
library(tibble, include.only = 'rowid_to_column')
library(dplyr)

################################################################################
#                                  import data                                 #
################################################################################

data_placenta <- read.csv("data/membrane.placentaire.tsv", sep="\t")
data_placenta$CI2 = factor(data_placenta$CI2)

################################################################################
#                                  filter data                                 #
################################################################################

varSupp <- which(apply(data_placenta, 2, var)== 0)
data_placenta <- data_placenta[, -varSupp]

cor_matrix <- cor(data_placenta[,-1])
cor_col <- findCorrelation(cor_matrix, cutoff=0.9)+1
data_placenta <- data_placenta[, -cor_col]

################################################################################
#                                  export data                                 #
################################################################################

write.csv(data_placenta, "data/membrane_filtered.csv", quote=FALSE)

################################################################################
#                         separate into train and test                         #
################################################################################

train_test = function(data, proportion){
  data <- data %>% rowid_to_column("rowid")
  train_set <- slice_sample(data, prop = proportion)
  test_set <- data %>% anti_join(as.data.frame(train_set), by = "rowid")
  return(list(train=train_set[,-1], test=test_set[,-1]))
}
