# STAT 408 Final Project
# Term: Fall 2024 
# Data: Ion concentrations from simulated Streamlab experiment 
# Experiment was conducted at Loyola University Chicago
# R Script Author: MacKenzie Michaels 

# Libraries ==============
library(tidyverse)
library(reshape)
library(car)
library(ggplot2)
library(plyr)
library(agricolae)
library(lme4)
library(emmeans)
library(sjPlot)
library(lmtest)
library(regclass)

# Data Import ==========================

streamImport <- read.csv("Biochar_Stream_Lab_Time_Series_Data.csv", stringsAsFactors = T)

# Data Check =================================
str(streamImport)

#Transpose Blog ====================================
#https://www.r-statistics.com/tag/transpose/

streamTranspose <- melt(streamImport, id = c("block", "plot", "nutrient_levels", "salt_presence", "ID", "coloration"))
write.csv(streamTranspose, "Stream_Lab_Time_Series_Data_Analysis_Long.csv")

# Reimport Transposed Data ==============================
# Import data ===========
stream <- read.csv("Stream_Lab_Time_Series_Data_Analysis_Long.csv", stringsAsFactors = T)

# Data Check ==================================
str(stream)
summary(stream)
stream$salt_presence <- as.factor(stream$salt_presence)

# Clean data ==================================
streamK$block <- as.factor(streamK$block)
streamK$nutrient_levels <- as.factor(streamK$nutrient_levels)
stream$plot <- as.factor(stream$plot)
stream$week <- as.factor(stream$week)
stream$value[stream$value<0] <- 0

table(streamION$week)
levels(streamION$block)
levels(streamION$nutrient_levels)
table(streamION$block, streamION$plot)
vif(lm(value ~ nutrient_levels + block + plot, data = streamION, type = "predictor"))
alias(lm(value ~ nutrient_levels + block/plot, data = streamION))

# Subset by ION ===============================
streamION <- subset(stream, variable == "K")
streamION <- droplevels(streamION)

# data check ==============
summary(streamION)
str(streamION)

# Nutrient level ion ===============
# Visualize data =============
ggplot(data = streamION,aes(x=week,y=value,colour = nutrient_levels))+
  geom_point()

hist(streamION$value, main = "Histogram of PO4 ppm", xlab = "ppm", breaks = 30, col = "lightblue")

graphLine <- ggplot(data = streamION, aes(x = week, y = value, colour = block))
graphLine + 
  theme_classic() +                   
  geom_point(size = 2) +
  facet_grid(block ~ nutrient_levels) +
  labs(title = "PO4 concentration by week and treatment")
  geom_smooth(method = "lm", se = FALSE)

# model development ===============
modION <- lmer((value) ~ nutrient_levels + (1|week/block), data = streamION)
null_modION <- lm(value ~ 1, data = streamION)
  
summary(modION)
AIC(modION,null_modION)

# assess homoscedasticity ==========
plot(resid(modION) ~ fitted(modION),
     xlab = "Fitted Values",
     ylab = "Residuals",
     main = "Residuals vs. Fitted for Nitrate")
abline(h = 0, col = "red", lty = 2)

residuals_mod <- residuals(modION)
fitted_vals <- fitted(modION)
lm_bp <- lm(residuals_mod^2 ~ fitted_vals)
bptest(lm_bp)

# assess normality ==========
hist(resid(modION), breaks = 20, main = "Histogram of Residuals for Nitrate", xlab = "Residuals")

qqnorm(resid(modION))
qqline(resid(modION), col = "red")

res_raw <- residuals(modION)
res_sd <- sqrt(var(res_raw))
res_standard <- res_raw / res_sd
ks.test(res_standard, "pnorm")

# Effects and variances ============
# Estimate of the fixed effects
fixef(modION)

# Estimate of the random effects
ranef(modION)

# Covariance matrix of the fixed effects
vcov(modION)
object <- as.data.frame(print(VarCorr(modION),comp=c("Variance","Std.Dev.")))
sum(object$vcov)
object[1,4]/sum(object$vcov)*100 # variance % explained by block within week 
object[2,4]/sum(object$vcov)*100 # variance % explained by week alone 
object[3,4]/sum(object$vcov)*100 # unexplained variance % 

# compare predicted values and observed values 
plot(streamION$value, predict(modION),
     xlab = "Observed Values",
     ylab = "Predicted Values",
     main = "Nitrate (NO3+)")
abline(0, 1, col = "red", lty = 2)

# Weighted model of ion subsets ===========================
mod <- lmer(value ~ nutrient_levels + (1|week), 
            data = streamION)
# Calculate residuals and fitted values
res_mod <- residuals(mod)
fit_mod <- fitted(mod)

# Estimate the variance of residuals
fit_var <- var(fit_mod)
# Calculate weights as the inverse of residual variance
weights <- 1 / (fit_mod^2 + fit_var)

# Fit the WLS model using the calculated weights
modION_wls <- lmer(value ~ nutrient_levels + (1 | week), data = streamION, 
                weights = weights)

# Summarize the WLS model
summary(modION_wls)
# Plot to assess homoscedasticity 
plot(resid(modION_wls) ~ fitted(modION_wls),
     xlab = "Fitted Values",
     ylab = "Residuals",
     main = "Residuals vs. Fitted")
abline(h = 0, col = "red", lty = 2)

residuals_mod <- residuals(modION_wls)
fitted_vals <- fitted(modION)
lm_bp <- lm(residuals_mod^2 ~ fitted_vals)
bptest(lm_bp)

#Plots to assess normality 
hist(resid(modION_wls), breaks = 20, main = "Histogram of Residuals", xlab = "Residuals")
qqnorm(resid(modION_wls))
qqline(resid(modION_wls), col = "red")

res_raw <- residuals(modION_wls)
res_sd <- sqrt(var(res_raw))
res_standard <- res_raw / res_sd
ks.test(res_standard, "pnorm")

# compare predicted values and observed values 
plot(streamION$value, predict(modION_wls),
     xlab = "Observed Values",
     ylab = "Predicted Values",
     main = "Observed vs. Predicted Calcium Ion Concentrations")
abline(0, 1, col = "red", lty = 2)

# checking for influential observations subset Ion ==================
df <- 3
alpha <- 0.01
crit <- qchisq(1 - alpha, df)
which(abs(rstudent(modION_wls)) > crit) # no influetial values 

# Salt presence ions ==============
# Visualize data
ggplot(data = streamION,aes(x=week,y=value,colour = salt_presence))+
  geom_point()

hist(streamION$value, main = "Histogram of Ion ppm", xlab = "ppm", breaks = 30, col = "lightblue")

graphLine <- ggplot(data = streamION, aes(x = week, y = value, colour = block))
graphLine + 
  theme_classic() +                   
  geom_point(size = 2) +
  facet_grid(block ~ salt_presence) +
  labs(title = "Na concentration by week and treatment")
geom_smooth(method = "lm", se = FALSE)

# salt model development ===============
summary(streamION)
modION <- lmer(sqrt(value) ~ salt_presence + (1|week), data = streamION)
null_modION <- lm(value ~ 1, data = streamION)

summary(modION)
AIC(modION,null_modION)

# assess homoscedasticity in salt model ==========
plot(resid(modION) ~ fitted(modION),
     xlab = "Fitted Values",
     ylab = "Residuals",
     main = "Residuals vs. Fitted")
abline(h = 0, col = "red", lty = 2)

residuals_mod <- residuals(modION)
fitted_vals <- fitted(modION)
lm_bp <- lm(residuals_mod^2 ~ fitted_vals)
bptest(lm_bp)

# assess normality in salt model ==========
hist(resid(modION), breaks = 20, main = "Histogram of Residuals", xlab = "Residuals")

qqnorm(resid(modION))
qqline(resid(modION), col = "red")

res_raw <- residuals(modION)
res_sd <- sqrt(var(res_raw))
res_standard <- res_raw / res_sd
ks.test(res_standard, "pnorm")

# checking for influential observations salt model ==================
df <- 3
alpha <- 0.01
crit <- qchisq(1 - alpha, df)
which(abs(rstudent(modION)) > crit) # no influential values 

# salt model Effects and variances ============
# Estimate of the fixed effects
fixef(modION)

# Estimate of the random effects
ranef(modION)

# Covariance matrix of the fixed effects
vcov(modION)
object <- as.data.frame(print(VarCorr(modION),comp=c("Variance","Std.Dev.")))
sum(object$vcov)
object[1,4]/sum(object$vcov)*100 # variance explained by week  
object[2,4]/sum(object$vcov)*100 # unexplained variance 

# compare predicted values and observed values 
plot(streamION$value, predict(modION),
     xlab = "Observed Values",
     ylab = "Predicted Values",
     main = "Chloride (Cl-)")
abline(0, 1, col = "red", lty = 2)

# salt model Weighted model of ion subsets ===========================
mod <- lmer(value ~ salt_presence + (1|week), 
            data = streamION)
# Calculate residuals and fitted values
res_mod <- residuals(mod)
fit_mod <- fitted(mod)

# Estimate the variance of residuals
fit_var <- var(fit_mod)
# Calculate weights as the inverse of residual variance
weights <- 1 / (fit_mod^2 + fit_var)

# Fit the WLS model using the calculated weights
modION_wls <- lmer(value ~ salt_presence + (1 | week), data = streamION, 
                   weights = weights)

# Summarize the WLS model
summary(modION_wls)
# Plot to assess homoscedasticity 
plot(resid(modION_wls) ~ fitted(modION_wls),
     xlab = "Fitted Values",
     ylab = "Residuals",
     main = "Residuals vs. Fitted")
abline(h = 0, col = "red", lty = 2)

residuals_mod <- residuals(modION_wls)
fitted_vals <- fitted(modION)
lm_bp <- lm(residuals_mod^2 ~ fitted_vals)
bptest(lm_bp)

#Plots to assess normality 
hist(resid(modION_wls), breaks = 20, main = "Histogram of Residuals", xlab = "Residuals")
qqnorm(resid(modION_wls))
qqline(resid(modION_wls), col = "red")

res_raw <- residuals(modION_wls)
res_sd <- sqrt(var(res_raw))
res_standard <- res_raw / res_sd
ks.test(res_standard, "pnorm")

# compare predicted values and observed values 
plot(streamION$value, predict(modION_wls),
     xlab = "Observed Values",
     ylab = "Predicted Values",
     main = "Observed vs. Predicted Chloride Ion Concentrations")
abline(0, 1, col = "red", lty = 2)
