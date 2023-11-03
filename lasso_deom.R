lasso_model <- function(data,label,train_ratio=0.5,savepath="lasso.RData"){
  ### data is the matrix data set for trainning and testing
  ### label for data
  ## train_ratio used for split train data from data
  
  library(caret)
  library(glmnet)
  library(pROC)
  train_indices <- createDataPartition(label, p = train_ratio, list = FALSE)
  train_data <- data[train_indices, ]
  test_data <- data[-train_indices, ]
  
  #####lasso
  {
    lasso_model <- cv.glmnet(x = as.matrix(train_data), y = label[train_indices], alpha = 1,family = "binomial",nlambda=100)
    p1<-plot(lasso_model)
    p2<-plot(lasso_model$glmnet.fit,"lambda")
    best_lambda <- lasso_model$lambda.1se
    model <- glmnet(x = as.matrix(train_data), y = label[train_indices], alpha = 1, lambda = best_lambda, family = "binomial")
    non_zero_coeffs <- coef(model, s = best_lambda)
    non_zero_coeffs <- data.frame(Gene=non_zero_coeffs@Dimnames[[1]][non_zero_coeffs@i+1],
                                  coef=non_zero_coeffs@x)
    predicted_train_classes <- predict(model, newx = as.matrix(train_data), s = best_lambda, type = "response")
    predicted_train_classes <- ifelse(predicted_train_classes > 0.5, 1, 0)
    predicted_test_classes <- predict(model, newx = as.matrix(test_data), s = best_lambda, type = "response")
    predicted_test_classes <- ifelse(predicted_test_classes > 0.5, 1, 0)
    
    conf_train_matrix <- confusionMatrix(as.factor(as.numeric(as.factor(predicted_train_classes[,1]))), as.factor(as.numeric(as.factor(label[train_indices]))))
    conf_test_matrix <- confusionMatrix(as.factor(as.numeric(as.factor(predicted_test_classes[,1]))), as.factor(as.numeric(as.factor(label[-train_indices]))))
  }
  
  roc_train_curve <- roc(as.numeric(as.factor(label[train_indices])), as.numeric(as.factor(predicted_train_classes)))
  roc_test_curve <- roc(as.numeric(as.factor(label[-train_indices])), as.numeric(as.factor(predicted_test_classes)))
  roc_data <- data.frame(
    FPR = c(roc_train_curve$specificities,roc_test_curve$specificities),
    TPR = c(roc_train_curve$sensitivities,roc_test_curve$sensitivities),
    Dataset = c(rep("Train",length(roc_train_curve$specificities)),rep("Test",length(roc_test_curve$specificities)))
  )
  
  # 绘制ROC曲线并添加自定义样式
  roc_plot <- ggroc(list(Train=roc_train_curve,Test=roc_test_curve),linewidth=1.5) + 
    scale_color_manual(name="",values = c("#ff6f61","#004c6d"))+
    geom_abline(slope = 1, intercept = 1, linetype = "dashed", color = "gray") +
    labs(title = "ROC Curve",
         x = "False Positive Rate",
         y = "True Positive Rate") +
    theme_bw() +
    theme(
      plot.title = element_text(size = 16, hjust = 0.5),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12))+
    geom_text(x = -0.7, y = 0.3, label = paste("Train AUC =", round(auc(roc_train_curve), 2)), size = 5, color = "#ff6f61")+
    geom_text(x = -0.7, y = 0.2, label = paste("Test AUC =", round(auc(roc_test_curve), 2)), size = 5, color = "#004c6d")
  save.image(savepath)
  
}