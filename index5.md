# Advanced Statistics

## GLMs

Generalized linear models are fit using the glm( ) function. The form of the glm function is

glm(formula, family=familytype(link=linkfunction), data=)

```
Family	Default Link Function
binomial	(link = "logit")
gaussian	(link = "identity")
Gamma	(link = "inverse")
inverse.gaussian	(link = "1/mu^2")
poisson	(link = "log")
quasi	(link = "identity", variance = "constant")
quasibinomial	(link = "logit")
quasipoisson	(link = "log")
```
See help(glm) for other modeling options. See help(family) for other allowable link functions for each family. Three subtypes of generalized linear models will be covered here: logistic regression, poisson regression, and survival analysis.

### Logistic Regression

Logistic regression is useful when you are predicting a binary outcome from a set of continuous predictor variables. It is frequently preferred over discriminant function analysis because of its less restrictive assumptions.

```
# Logistic Regression
# where F is a binary factor and 
# x1-x3 are continuous predictors 
fit <- glm(F~x1+x2+x3,data=mydata,family=binomial())
summary(fit) # display results
confint(fit) # 95% CI for the coefficients
exp(coef(fit)) # exponentiated coefficients
exp(confint(fit)) # 95% CI for exponentiated coefficients
predict(fit, type="response") # predicted values
residuals(fit, type="deviance") # residuals
```

You can use anova(fit1,fit2, test="Chisq") to compare nested models. Additionally, cdplot(F~x, data=mydata) will display the conditional density plot of the binary outcome F on the continuous x variable.

conditional density plot click to view

### Poisson Regression

Poisson regression is useful when predicting an outcome variable representing counts from a set of continuous predictor variables.

```
# Poisson Regression
# where count is a count and 
# x1-x3 are continuous predictors 
fit <- glm(count ~ x1+x2+x3, data=mydata, family=poisson())
summary(fit) display results
```

If you have overdispersion (see if residual deviance is much larger than degrees of freedom), you may want to use quasipoisson() instead of poisson().

### Survival Analysis

Survival analysis (also called event history analysis or reliability analysis) covers a set of techniques for modeling the time to an event. Data may be right censored - the event may not have occured by the end of the study or we may have incomplete information on an observation but know that up to a certain time the event had not occured (e.g. the participant dropped out of study in week 10 but was alive at that time).

While generalized linear models are typically analyzed using the glm( ) function, survival analyis is typically carried out using functions from the survival package . The survival package can handle one and two sample problems, parametric accelerated failure models, and the Cox proportional hazards model.

Data are typically entered in the format start time, stop time, and status (1=event occured, 0=event did not occur). Alternatively, the data may be in the format time to event and status (1=event occured, 0=event did not occur). A status=0 indicates that the observation is right cencored. Data are bundled into a Surv object via the Surv( ) function prior to further analyses.

survfit( ) is used to estimate a survival distribution for one or more groups.
survdiff( ) tests for differences in survival distributions between two or more groups. 
coxph( ) models the hazard function on a set of predictor variables.

```
# Mayo Clinic Lung Cancer Data
library(survival)

# learn about the dataset
help(lung)

# create a Surv object 
survobj <- with(lung, Surv(time,status))

# Plot survival distribution of the total sample
# Kaplan-Meier estimator 
fit0 <- survfit(survobj~1, data=lung)
summary(fit0)
plot(fit0, xlab="Survival Time in Days", 
   ylab="% Surviving", yscale=100,
   main="Survival Distribution (Overall)") 

# Compare the survival distributions of men and women 
fit1 <- survfit(survobj~sex,data=lung)

# plot the survival distributions by sex 
plot(fit1, xlab="Survival Time in Days", 
  ylab="% Surviving", yscale=100, col=c("red","blue"),
  main="Survival Distributions by Gender") 
  legend("topright", title="Gender", c("Male", "Female"),
  fill=c("red", "blue"))

# test for difference between male and female 
# survival curves (logrank test) 
survdiff(survobj~sex, data=lung) 

# predict male survival from age and medical scores 
MaleMod <- coxph(survobj~age+ph.ecog+ph.karno+pat.karno,
  data=lung, subset=sex==1)

# display results 
MaleMod

# evaluate the proportional hazards assumption 
cox.zph(MaleMod)
```

survival distribution for total group survival distribution by gender click to view

See Thomas Lumley's R news article on the survival package for more information. Other good sources include Mai Zhou's Use R Software to do Survival Analysis and Simulation and M. J. Crawley's chapter on Survival Analysis.

## Discriminant Function

The MASS package contains functions for performing linear and quadratic 
discriminant function analysis. Unless prior probabilities are specified, each assumes proportional prior probabilities (i.e., prior probabilities are based on sample sizes). In the examples below, lower case letters are numeric variables and upper case letters are categorical factors.

### Linear Discriminant Function
```
# Linear Discriminant Analysis with Jacknifed Prediction 
library(MASS)
fit <- lda(G ~ x1 + x2 + x3, data=mydata, 
   na.action="na.omit", CV=TRUE)
fit # show results
```

The code above performs an LDA, using listwise deletion of missing data. CV=TRUE generates jacknifed (i.e., leave one out) predictions. The code below assesses the accuracy of the prediction.

```
# Assess the accuracy of the prediction
# percent correct for each category of G
ct <- table(mydata$G, fit$class)
diag(prop.table(ct, 1))
# total percent correct
sum(diag(prop.table(ct)))
```

lda() prints discriminant functions based on centered (not standardized) variables. The "proportion of trace" that is printed is the proportion of between-class variance that is explained by successive discriminant functions. No significance tests are produced. Refer to the section on MANOVA for such tests.

### Quadratic Discriminant Function
To obtain a quadratic discriminant function use qda( ) instead of lda( ). Quadratic discriminant function does not assume homogeneity of variance-covariance matrices.

```
# Quadratic Discriminant Analysis with 3 groups applying 
# resubstitution prediction and equal prior probabilities. 
library(MASS)
fit <- qda(G ~ x1 + x2 + x3 + x4, data=na.omit(mydata),
  prior=c(1,1,1)/3))
```

Note the alternate way of specifying listwise deletion of missing data. Re-subsitution (using the same data to derive the functions and evaluate their prediction accuracy) is the default method unless CV=TRUE is specified. Re-substitution will be overly optimistic.

### Visualizing the Results
You can plot each observation in the space of the first 2 linear discriminant functions using the following code. Points are identified with the group ID.

```
# Scatter plot using the 1st two discriminant dimensions 
plot(fit) # fit from lda
```

linear discrimiant plot of points click to view

The following code displays histograms and density plots for the observations in each group on the first linear discriminant dimension. There is one panel for each group and they all appear lined up on the same graph.

```
# Panels of histograms and overlayed density plots
# for 1st discriminant function
plot(fit, dimen=1, type="both") # fit from lda
```

lda histogram/density plot click to view

The partimat( ) function in the klaR package can display the results of a linear or quadratic classifications 2 variables at a time.

```
# Exploratory Graph for LDA or QDA
library(klaR)
partimat(G~x1+x2+x3,data=mydata,method="lda")

partimat plot click to view
```

You can also produce a scatterplot matrix with color coding by group.

```
# Scatterplot for 3 Group Problem 
pairs(mydata[c("x1","x2","x3")], main="My Title ", pch=22, 
   bg=c("red", "yellow", "blue")[unclass(mydata$G)])
```
scatterplot matrix click to view

### Test Assumptions
See (M)ANOVA Assumptions for methods of evaluating multivariate normality and homogeneity of covariance matrices.

