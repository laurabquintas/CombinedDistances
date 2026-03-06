
library(dplyr)
library(parallel)
library(survival)
library(survminer)
library(glmnet)
library(grpreg)
library(penalized)
library(gglasso)
library(ggplot2)

################# Load Data ######################

df <- read.csv("/Users/Lau/Documents/GitHub/Combined-Distances/AML/GEN/Results/survival/t = 0.02/clinicaldata_train.csv")

df_encode <- read.csv("/Users/Lau/Documents/GitHub/Combined-Distances/AML/GEN/Results/survival/t = 0.02/clinicaldata_train_HE.csv")





######## USING GROUP LASSO ######

######### whole dataset with new
set.seed(42)

y <- Surv(df$Survival.in.days,df$VitalStatus)
group_vector <- c('age', 'Gender', 'maxCCF', 'TP53', 'IDH1', 'SRSF2', 'DNMT3A', 'FLT3', 'PTPN11', 'NRAS', 'IDH2', 'KRAS', 'ASXL1', 'TET2', 'NPM1', 'WT1', 'RUNX1', 'FLT3', 'DISC_union', 'DISC_union', 'DISC_intersection', 'DISC_intersection', 'DISC_intersection', 'DISC_intersection', 'DISC_intersection', 'DISC_intersection', 'DISC_intersection', 'DISC_intersection', 'DISC_intersection', 'DISC_intersection', 'DISC_intersection', 'DISC_intersection', 'CASet_intersection', 'CASet_intersection', 'CASet_intersection', 'CASet_union', 'CASet_union', 'BD', 'BD', '1BD', '1BD', '2BD', '2BD', 'AD', 'AD', 'CD', 'CD', 'PD', 'PD', 'PD', 'PD', 'PD', 'PD', 'PD', 'PCD', 'PCD', 'MLTD', 'MLTD', 'MP3_union', 'MP3_union', 'MP3_intersection', 'MP3_intersection')

group_factor <- factor(group_vector)


cvfit <- cv.grpsurv(df_encode, y, group_factor, penalty='grLasso', returnY = TRUE)
plot(cvfit)
summary(cvfit)

minlambda <- cvfit$lambda.min
plot(cvfit$fit)

plot(cvfit, type="rsq")

fit <- grpsurv(df_encode, y, group_factor, penalty = 'grLasso', returnY = TRUE)
pdf(file = "/Users/Lau/Documents/GitHub/Combined-Distances/AML/GEN/Results/survival/t = 0.02/GL_coef_lambda.pdf", width = 6, height = 6)

plot(fit)
dev.off()

df <- data.frame(lambda = fit$lambda, AIC = AIC(fit))

# Find the index of the lambda value in 'df' that matches 'cvfit$lambda.min'
index <- which(abs(df$lambda - cvfit$lambda.min) < .Machine$double.eps^0.5)

# Use this index to find the corresponding AIC value
aic <- df$AIC[index]

#best_lambda_index <- which.min(AIC(fit))
#best_lambda <- fit$lambda[best_lambda_index]
#AIC_min <- min(AIC(fit))
#print(AIC_min)

min_coef <- coef(fit,lambda=minlambda)

#print(AIC_coef)

p <- ggplot(df, aes(x = lambda, y = AIC)) +
  geom_point() +
  geom_vline(xintercept = minlambda, linetype = "dotted") +
  geom_text(aes(x = minlambda, y = aic, label = sprintf("AIC: %.2f", aic)), vjust = -1.9) +
  labs(x = "Lambda",
       y = "AIC") +
  theme_minimal()
ggsave(filename = "/Users/Lau/Documents/GitHub/Combined-Distances/AML/GEN/Results/survival/t = 0.02/GL_AIC_lambda.png", plot = p, width = 10, height = 6, dpi = 300)

# Display the plot
print(p)

### LOG LIKELIHOOD

df <- data.frame(lambda = fit$lambda, LL = logLik(fit))

# Find the index of the lambda value in 'df' that matches 'cvfit$lambda.min'
index <- which(abs(df$lambda - cvfit$lambda.min) < .Machine$double.eps^0.5)

# Use this index to find the corresponding AIC value
loglikelihood <- df$LL[index]

print(loglikelihood)


data_to_plot <- data.frame(Coefficients = min_coef, Names = attr(min_coef, "names"))
p <- ggplot(data_to_plot, aes(x=Names, y=Coefficients)) + 
  geom_point() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(filename = "/Users/Lau/Documents/GitHub/Combined-Distances/AML/GEN/Results/survival/t = 0.02/GL_coef_.png", plot = p, width = 20, height = 6, dpi = 300)

print(p)
# Save AIC coefficients as a CSV file
write.csv(data_to_plot, file = "/Users/Lau/Documents/GitHub/Combined-Distances/AML/GEN/Results/survival/t = 0.02/R_coef.csv", row.names = FALSE)



##### smt ####
cvfit <- cv.grpsurv(df_encode, y, group_factor, penalty='grLasso', returnY = TRUE)
plot(cvfit)
summary(cvfit)

minlambda <- cvfit$lambda.min
plot(cvfit$fit)

plot(cvfit, type="rsq")



fit <- grpsurv(df_encode, y, group_factor, penalty = 'grLasso', returnY = TRUE)
pdf(file = "/Users/Lau/Documents/GitHub/Combined-Distances/BREAST CANCER/GEN/Results/WARD/Survival/Simple LOESS/t = 0.02/GL_coef_lambda.pdf", width = 6, height = 6)

plot(fit)
dev.off()

######
folds <- createFolds(y, k = 5)

# Store results
results <- list()

# for each fold
for(i in 1:5){
  # Get the indices for the train and test sets
  train_indices <- setdiff(1:nrow(df_encode), folds[[i]])
  test_indices <- folds[[i]]
  
  # Subset the data
  train_data <- df_encode[train_indices, ]
  test_data <- df_encode[test_indices, ]
  
  # Subset the target variable
  train_y <- y[train_indices]
  test_y <- y[test_indices]
  
  # Fit the model on the training data
  fit <- grpsurv(train_data, train_y, group_factor, penalty = 'grLasso', returnY = TRUE)
  
  # Evaluate the model on the test data
  # (You'll need to replace this with the appropriate function to evaluate your model,
  # which might involve predicting the survival on the test data and comparing it to the actual survival.)
  # result <- evaluate_model(fit, test_data, test_y)
  
  min_coef <- coef(fit,lambda=minlambda)
  data_to_plot <- data.frame(Coefficients = min_coef, Names = attr(min_coef, "names"))
  p <- ggplot(data_to_plot, aes(x=Names, y=Coefficients)) + 
    geom_point() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  
  #ggsave(filename = "/Users/Lau/Documents/GitHub/Combined-Distances/BREAST CANCER/GEN/Results/WARD/Survival/Simple LOESS/t = 0.02/GL_AIC_lambda.png", plot = p, width = 10, height = 6, dpi = 300)
  
  # Display the plot
  print(p)
  
  # Store the result
  results[[i]] <- min_coef
}

# Aggregate results
# (This will depend on the form of your results, but might involve taking the mean of a performance metric across the 5 folds)
final_result <- aggregate_results(results)
######

df <- data.frame(lambda = fit$lambda, AIC = AIC(fit))

# Find the index of the lambda value in 'df' that matches 'cvfit$lambda.min'
index <- which(abs(df$lambda - cvfit$lambda.min) < .Machine$double.eps^0.5)

# Use this index to find the corresponding AIC value
aic <- df$AIC[index]

min_coef <- coef(fit,lambda=minlambda)

#print(AIC_coef)

p <- ggplot(df, aes(x = lambda, y = AIC)) +
  geom_point() +
  geom_vline(xintercept = minlambda, linetype = "dotted") +
  labs(x = "Lambda",
       y = "AIC") +
  theme_minimal()
ggsave(filename = "/Users/Lau/Documents/GitHub/Combined-Distances/BREAST CANCER/GEN/Results/WARD/Survival/Simple LOESS/t = 0.02/GL_AIC_lambda.png", plot = p, width = 10, height = 6, dpi = 300)

# Display the plot
print(p)



data_to_plot <- data.frame(Coefficients = min_coef, Names = attr(min_coef, "names"))
p <- ggplot(data_to_plot, aes(x=Names, y=Coefficients)) + 
  geom_point() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(filename = "/Users/Lau/Documents/GitHub/Combined-Distances/BREAST CANCER/GEN/Results/WARD/Survival/Simple LOESS/t = 0.02/GL_coef_AIC.png", plot = p, width = 20, height = 6, dpi = 300)

print(p)
# Save AIC coefficients as a CSV file
write.csv(data_to_plot, file = "/Users/Lau/Documents/GitHub/Combined-Distances/BREAST CANCER/GEN/Results/WARD/Survival/Simple LOESS/t = 0.02/R_coef.csv", row.names = FALSE)

# Filter to only include nonzero coefficients

data_to_plot_nonzero <- data_to_plot[data_to_plot$Coefficients != 0, ]

# Now, plot using ggplot
p <- ggplot(data_to_plot_nonzero, aes(x=Names, y=Coefficients)) + 
  geom_point() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(filename = "/Users/Lau/Documents/GitHub/Combined-Distances/AML/GEN/Results/survival/t = 0.02/GL_coef_AIC_filtered.png", plot = p, width = 8, height = 4, dpi = 300)

print(p)


######## whole dataset ########
set.seed(42)

y <- Surv(df$Survival.in.days,df$VitalStatus)
group_vector <- c('age', 'Gender', 'DISC_union', 'DISC_union', 'DISC_intersection', 'DISC_intersection', 'DISC_intersection', 'DISC_intersection', 'DISC_intersection', 'DISC_intersection', 'DISC_intersection', 'DISC_intersection', 'DISC_intersection', 'DISC_intersection', 'DISC_intersection', 'DISC_intersection', 'CASet_intersection', 'CASet_intersection', 'CASet_intersection', 'CASet_union', 'CASet_union', 'BD', 'BD', '1BD', '1BD', '2BD', '2BD', 'AD', 'AD', 'CD', 'CD', 'PD', 'PD', 'PD', 'PD', 'PD', 'PD', 'PD', 'PCD', 'PCD', 'MLTD', 'MLTD', 'MP3_union', 'MP3_union', 'MP3_sigmoid', 'MP3_sigmoid', 'MP3_intersection', 'MP3_intersection', 'MP3_geometric', 'MP3_geometric', 'CASet_geometric', 'CASet_geometric', 'CASet_geometric', 'CASet_geometric', 'CASet_geometric', 'CASet_geometric', 'CASet_sigmoid', 'CASet_sigmoid', 'DISC_geometric', 'DISC_geometric', 'DISC_geometric', 'DISC_geometric', 'DISC_geometric', 'DISC_geometric', 'DISC_geometric', 'DISC_geometric', 'DISC_geometric', 'DISC_geometric', 'DISC_geometric', 'DISC_geometric', 'DISC_sigmoid', 'DISC_sigmoid')

group_factor <- factor(group_vector)


cvfit <- cv.grpsurv(df_encode_2, y, group_factor, penalty='grLasso', returnY = TRUE)
plot(cvfit)
summary(cvfit)

fit <- grpsurv(df_encode_2, y, group_factor, penalty = 'grLasso', returnY = TRUE)
pdf(file = "/Users/Lau/Documents/GitHub/Combined-Distances/AML/Results/survival/t = 0.02/GL_coef_lambda_newmethods.pdf", width = 6, height = 6)

plot(fit)

dev.off()

df <- data.frame(lambda = fit$lambda, AIC = AIC(fit))

min_lambda <- df$lambda[which.min(df$AIC)]
min_aic <- min(df$AIC)

#best_lambda_index <- which.min(AIC(fit))
#best_lambda <- fit$lambda[best_lambda_index]
#AIC_min <- min(AIC(fit))
#print(AIC_min)

AIC_coef <- coef(fit,lambda=min_lambda)

print(AIC_coef)

p <- ggplot(df, aes(x = lambda, y = AIC)) +
  geom_point() +
  geom_vline(xintercept = min_lambda, linetype = "dotted") +
  geom_text(aes(x = min_lambda, y = min_aic, label = sprintf("Min AIC: %.2f", min_aic)), vjust = -1.9) +
  labs(x = "Lambda",
       y = "AIC") +
  theme_minimal()
ggsave(filename = "/Users/Lau/Documents/GitHub/Combined-Distances/AML/Results/survival/t = 0.02/GL_AIC_lambda_newmethods.png", plot = p, width = 10, height = 6, dpi = 300)

# Display the plot
print(p)


data_to_plot <- data.frame(Coefficients = AIC_coef, Names = attr(AIC_coef, "names"))
p <- ggplot(data_to_plot, aes(x=Names, y=Coefficients)) + 
  geom_point() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(filename = "/Users/Lau/Documents/GitHub/Combined-Distances/AML/Results/survival/t = 0.02/GL_coef_AIC_newmethods.png", plot = p, width = 20, height = 6, dpi = 300)

print(p)
# Save AIC coefficients as a CSV file
write.csv(data_to_plot, file = "/Users/Lau/Documents/GitHub/Combined-Distances/AML/Results/survival/t = 0.02/R_aic_coef_newmethods.csv", row.names = FALSE)
