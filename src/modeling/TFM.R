require(caret)
require(Cubist)

setwd("~/Documents/BioInformatics/stage1/data/")

train <- read.delim("./processed_data/csv/ready_to_use_data_train.csv", header = TRUE, sep = ",", quote = "\"", dec = ".",
                    fill = TRUE, comment.char = "#")
test <- read.delim("./processed_data/csv/ready_to_use_data_test.csv", header = TRUE, sep = ",", quote = "\"", dec = ".",
                    fill = TRUE, comment.char = "#")

#hg <- read.delim("./files/CSV/humangenome.csv", header = TRUE, sep = ",", quote = "\"", dec = ".",
#                   fill = TRUE, comment.char = "#")

ctrl <- cubistControl(unbiased = FALSE,
                      rules = 50,
                      extrapolation = 1,
                      sample = 1,
                      seed = 42,
                      label = "outcome")

tfmModelTHSA <- cubist(train[,c('length','polar_count','hydr_count')],train[,c("thsa")],committees=3,control = ctrl)
tfmModelRHSA <- cubist(train[,c('length','polar_count','hydr_count')],train[,c("rhsa")],committees=3,control = ctrl)
tfmModelLHPSA <- cubist(train[,c('length','polar_count','hydr_count')],train[,c("size")],committees=3,control = ctrl)

TFMPredTHSA <- predict(tfmModelTHSA,test[,c('length','polar_count','hydr_count')])
TFMPredRHSA <- predict(tfmModelRHSA,test[,c('length','polar_count','hydr_count')])
TFMPredLHPSA <- predict(tfmModelLHPSA,test[,c('length','polar_count','hydr_count')])

write.csv(TFMPredTHSA,'./predictions/thsa_tfm_prediction.csv')
write.csv(TFMPredRHSA,'./predictions/rhsa_tfm_prediction.csv')
write.csv(TFMPredLHPSA,'./predictions/lhpsa_tfm_prediction.csv')
