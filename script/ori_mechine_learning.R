
ori_ml <- function(data,SplitData,Method,Labelcol,nCrossValidation,Plotimportant=T,Plotcm=T){
  require("tidyverse")
  require("caret")
  require("mlbench")

  draw_confusion_matrix <- function(cm) {
    
    layout(matrix(c(1,1,2)))
    par(mar=c(2,2,2,2))
    plot(c(100, 345), c(300, 450), type = "n", xlab="", ylab="", xaxt='n', yaxt='n')
    title('CONFUSION MATRIX', cex.main=2)
    
    n1 <- colnames(cm$table)[1]
    n2 <- colnames(cm$table)[2]
    
    # create the matrix 
    rect(150, 430, 240, 370, col='#3F97D0')
    text(195, 435, n1, cex=1.5)
    rect(250, 430, 340, 370, col='#F7AD50')
    text(295, 435, n2, cex=1.5)
    text(125, 370, 'Predicted', cex=1.6, srt=90, font=2)
    text(245, 450, 'Actual', cex=1.6, font=2)
    rect(150, 305, 240, 365, col='#F7AD50')
    rect(250, 305, 340, 365, col='#3F97D0')
    text(140, 400, n1, cex=1.5, srt=90)
    text(140, 335, n2, cex=1.5, srt=90)
    
    # add in the cm results 
    res <- as.numeric(cm$table)
    text(195, 400, res[1], cex=1.9, font=2, col='white')
    text(195, 335, res[2], cex=1.9, font=2, col='white')
    text(295, 400, res[3], cex=1.9, font=2, col='white')
    text(295, 335, res[4], cex=1.9, font=2, col='white')
    
    # add in the specifics 
    plot(c(100, 0), c(100, 0), type = "n", xlab="", ylab="", main = "DETAILS", xaxt='n', yaxt='n',cex.main=1.7)
    text(10, 85, names(cm$byClass[1]), cex=1.5, font=2)
    text(10, 70, round(as.numeric(cm$byClass[1]), 3), cex=1.5)
    text(30, 85, names(cm$byClass[2]), cex=1.5, font=2)
    text(30, 70, round(as.numeric(cm$byClass[2]), 3), cex=1.5)
    text(50, 85, names(cm$byClass[5]), cex=1.5, font=2)
    text(50, 70, round(as.numeric(cm$byClass[5]), 3), cex=1.5)
    text(70, 85, names(cm$byClass[6]), cex=1.5, font=2)
    text(70, 70, round(as.numeric(cm$byClass[6]), 3), cex=1.5)
    text(90, 85, names(cm$byClass[7]), cex=1.5, font=2)
    text(90, 70, round(as.numeric(cm$byClass[7]), 3), cex=1.5)
    
    # add in the accuracy information 
    text(30, 35, names(cm$overall[1]), cex=1.8, font=2)
    text(30, 20, round(as.numeric(cm$overall[1]), 3), cex=1.7)
    text(70, 35, names(cm$overall[2]), cex=1.8, font=2)
    text(70, 20, round(as.numeric(cm$overall[2]), 3), cex=1.7)
  }  
  
  data <- na.omit(data)
  data[[Labelcol]] <- as.factor(data[[Labelcol]])
  training.samples <- data[[Labelcol]] %>% 
    createDataPartition(p = SplitData, list = FALSE)
  train.data  <- data[training.samples, ]
  test.data <- data[-training.samples, ]
  
  FM<-as.formula(paste0(Labelcol,"~."))
  set.seed(123)
  model <- train(
    FM, data = train.data, method = Method,
    trControl = trainControl("cv", number = nCrossValidation),
    #importance = TRUE
  )
  
  if(Plotimportant==T){
    gbmImp <- varImp(model, scale = FALSE)
    N <- ncol(train.data)-1
    if(N <10){
      p <- plot(gbmImp, top = N)
    }else{
      p <- plot(gbmImp, top = 10)
    }
    print(p)
  }

 #Plotcm=T
  if(Plotcm==T){
    predicted.classes <- model %>% predict(test.data)
    cm <- confusionMatrix(predicted.classes, as.factor(test.data[[Labelcol]]))
    draw_confusion_matrix(cm)
  }
  
  }


# rm(list = ls())
# SplitData <- 0.8
# data("PimaIndiansDiabetes2", package = "mlbench")
# data <- na.omit(PimaIndiansDiabetes2)
# Labelcol <- "diabetes"
# Method <- "rf" #"svmLinear","glm"
# nCrossValidation <- 10
# Plotimportant=T
# ori_ml(data,SplitData,Method,Labelcol,nCrossValiation=nCrossValidation,Plotimportant=F,Plotcm=T)
# ori_ml(data,SplitData,Method,Labelcol,nCrossValiation=nCrossValidation,Plotimportant=T,Plotcm=F)
















