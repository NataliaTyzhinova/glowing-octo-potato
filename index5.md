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

## Time Series

R has extensive facilities for analyzing time series data. This section describes the creation of a time series, seasonal decomposition, modeling with exponential and ARIMA models, and forecasting with the forecast package.

### Creating a time series
The ts() function will convert a numeric vector into an R time series object. The format is ts(vector, start=, end=, frequency=) where start and end are the times of the first and last observation and frequency is the number of observations per unit time (1=annual, 4=quartly, 12=monthly, etc.).

```
# save a numeric vector containing 72 monthly observations
# from Jan 2009 to Dec 2014 as a time series object
myts <- ts(myvector, start=c(2009, 1), end=c(2014, 12), frequency=12) 

# subset the time series (June 2014 to December 2014)
myts2 <- window(myts, start=c(2014, 6), end=c(2014, 12)) 

# plot series
plot(myts)
```

### Seasonal Decomposition
A time series with additive trend, seasonal, and irregular components can be decomposed using the stl() function. Note that a series with multiplicative effects can often by transformed into series with additive effects through a log transformation (i.e., newts <- log(myts)).

```
# Seasonal decomposition
fit <- stl(myts, s.window="period")
plot(fit)

# additional plots
monthplot(myts)
library(forecast)
seasonplot(myts)
```

### Exponential Models
Both the HoltWinters() function in the base installation, and the ets() function in the forecast package, can be used to fit exponential models.

```
# simple exponential - models level
fit <- HoltWinters(myts, beta=FALSE, gamma=FALSE)
# double exponential - models level and trend
fit <- HoltWinters(myts, gamma=FALSE)
# triple exponential - models level, trend, and seasonal components
fit <- HoltWinters(myts)

# predictive accuracy
library(forecast)
accuracy(fit)

# predict next three future values
library(forecast)
forecast(fit, 3)
plot(forecast(fit, 3))
```

### ARIMA Models
The arima() function can be used to fit an autoregressive integrated moving averages model. Other useful functions include:

```
lag(ts, k)	lagged version of time series, shifted back k observations
diff(ts, differences=d)	difference the time series d times
ndiffs(ts)	Number of differences required to achieve stationarity (from the forecast package)
acf(ts)	autocorrelation function
pacf(ts)	partial autocorrelation function
adf.test(ts)	Augemented Dickey-Fuller test. Rejecting the null hypothesis suggests that a time series is stationary (from the tseries package)
Box.test(x, type="Ljung-Box")	Pormanteau test that observations in vector or time series x are independent
Note that the forecast package has somewhat nicer versions of acf() and pacf() called Acf() and Pacf() respectively.
```

```
# fit an ARIMA model of order P, D, Q
fit <- arima(myts, order=c(p, d, q)

# predictive accuracy
library(forecast)
accuracy(fit)

# predict next 5 observations
library(forecast)
forecast(fit, 5)
plot(forecast(fit, 5))
```

### Automated Forecasting
The forecast package provides functions for the automatic selection of exponential and ARIMA models. The ets() function supports both additive and multiplicative models. The auto.arima() function can handle both seasonal and nonseasonal ARIMA models. Models are chosen to maximize one of several fit criteria.

```
library(forecast)
# Automated forecasting using an exponential model
fit <- ets(myts)

# Automated forecasting using an ARIMA model
fit <- auto.arima(myts)
```

### Going Further
There are many good online resources for learning time series analysis with R. These include A little book of R for time series by Avril Chohlan and DataCamp's manipulating time series in R course by Jeffrey Ryan.

## Factor Analysis 

This section covers principal components and factor analysis. The latter includes both exploratory and confirmatory methods.

### Principal Components
The princomp( ) function produces an unrotated principal component analysis.

```
# Pricipal Components Analysis
# entering raw data and extracting PCs 
# from the correlation matrix 
fit <- princomp(mydata, cor=TRUE)
summary(fit) # print variance accounted for 
loadings(fit) # pc loadings 
plot(fit,type="lines") # scree plot 
fit$scores # the principal components
biplot(fit)
```

scree plot biplot click to view

Use cor=FALSE to base the principal components on the covariance matrix. Use the covmat= option to enter a correlation or covariance matrix directly. If entering a covariance matrix, include the option n.obs=.

The principal( ) function in the psych package can be used to extract and rotate principal components.

```
# Varimax Rotated Principal Components
# retaining 5 components 
library(psych)
fit <- principal(mydata, nfactors=5, rotate="varimax")
fit # print results
```

mydata can be a raw data matrix or a covariance matrix. Pairwise deletion of missing data is used. rotate can "none", "varimax", "quatimax", "promax", "oblimin", "simplimax", or "cluster"

### Exploratory Factor Analysis
The factanal( ) function produces maximum likelihood factor analysis.

```
# Maximum Likelihood Factor Analysis
# entering raw data and extracting 3 factors, 
# with varimax rotation 
fit <- factanal(mydata, 3, rotation="varimax")
print(fit, digits=2, cutoff=.3, sort=TRUE)
# plot factor 1 by factor 2 
load <- fit$loadings[,1:2] 
plot(load,type="n") # set up plot 
text(load,labels=names(mydata),cex=.7) # add variable names
```

exploratory factor analysis click to view

The rotation= options include "varimax", "promax", and "none". Add the option scores="regression" or "Bartlett" to produce factor scores. Use the covmat= option to enter a correlation or covariance matrix directly. If entering a covariance matrix, include the option n.obs=.

The factor.pa( ) function in the psych package offers a number of factor analysis related functions, including principal axis factoring.

```
# Principal Axis Factor Analysis
library(psych)
fit <- factor.pa(mydata, nfactors=3, rotation="varimax")
fit # print results
```

mydata can be a raw data matrix or a covariance matrix. Pairwise deletion of missing data is used. Rotation can be "varimax" or "promax".

### Determining the Number of Factors to Extract
A crucial decision in exploratory factor analysis is how many factors to extract. The nFactors package offer a suite of functions to aid in this decision. Details on this methodology can be found in a PowerPoint presentation by Raiche, Riopel, and Blais. Of course, any factor solution must be interpretable to be useful.

```
# Determine Number of Factors to Extract
library(nFactors)
ev <- eigen(cor(mydata)) # get eigenvalues
ap <- parallel(subject=nrow(mydata),var=ncol(mydata),
  rep=100,cent=.05)
nS <- nScree(x=ev$values, aparallel=ap$eigen$qevpea)
plotnScree(nS)
```

number of factors click to view

### Going Further
The FactoMineR package offers a large number of additional functions for exploratory factor analysis. This includes the use of both quantitative and qualitative variables, as well as the inclusion of supplimentary variables and observations. Here is an example of the types of graphs that you can create with this package.

```
# PCA Variable Factor Map 
library(FactoMineR)
result <- PCA(mydata) # graphs generated automatically
```

factominer1 factominer2 click to view

Thye GPARotation package offers a wealth of rotation options beyond varimax and promax.

### Structual Equation Modeling
Confirmatory Factor Analysis (CFA) is a subset of the much wider Structural Equation Modeling (SEM) methodology. SEM is provided in R via the sem package. Models are entered via RAM specification (similar to PROC CALIS in SAS). While sem is a comprehensive package, my recommendation is that if you are doing significant SEM work, you spring for a copy of AMOS. It can be much more user-friendly and creates more attractive and publication ready output. Having said that, here is a CFA example using sem.

CFA Model

Assume that we have six observered variables (X1, X2, ..., X6). We hypothesize that there are two unobserved latent factors (F1, F2) that underly the observed variables as described in this diagram. X1, X2, and X3 load on F1 (with loadings lam1, lam2, and lam3). X4, X5, and X6 load on F2 (with loadings lam4, lam5, and lam6). The double headed arrow indicates the covariance between the two latent factors (F1F2). e1 thru e6 represent the residual variances (variance in the observed variables not accounted for by the two latent factors). We set the variances of F1 and F2 equal to one so that the parameters will have a scale. This will result in F1F2 representing the correlation between the two latent factors.

For sem, we need the covariance matrix of the observed variables - thus the cov( ) statement in the code below. The CFA model is specified using the specify.model( ) function. The format is arrow specification, parameter name, start value. Choosing a start value of NA tells the program to choose a start value rather than supplying one yourself. Note that the variance of F1 and F2 are fixed at 1 (NA in the second column). The blank line is required to end the RAM specification.

```
# Simple CFA Model
library(sem)
mydata.cov <- cov(mydata)
model.mydata <- specify.model() 
F1 ->  X1, lam1, NA
F1 ->  X2, lam2, NA 
F1 ->  X3, lam3, NA 
F2 ->  X4, lam4, NA 
F2 ->  X5, lam5, NA 
F2 ->  X6, lam6, NA 
X1 <-> X1, e1,   NA 
X2 <-> X2, e2,   NA 
X3 <-> X3, e3,   NA 
X4 <-> X4, e4,   NA 
X5 <-> X5, e5,   NA 
X6 <-> X6, e6,   NA 
F1 <-> F1, NA,    1 
F2 <-> F2, NA,    1 
F1 <-> F2, F1F2, NA

mydata.sem <- sem(model.mydata, mydata.cov, nrow(mydata))
# print results (fit indices, paramters, hypothesis tests) 
summary(mydata.sem)
# print standardized coefficients (loadings) 
std.coef(mydata.sem)
```

You can use the boot.sem( ) function to bootstrap the structual equation model. See help(boot.sem) for details. Additionally, the function mod.indices( ) will produce modification indices. Using modification indices to improve model fit by respecifying the parameters moves you from a confirmatory to an exploratory analysis.

For more information on sem, see Structural Equation Modeling with the sem Package in R, by John Fox.

## Correspondence Analysis

Correspondence analysis provides a graphic method of exploring the relationship between variables in a contingency table. There are many options for correspondence analysis in R. I recommend the ca package by Nenadic and Greenacre because it supports supplimentary points, subset analyses, and comprehensive graphics. You can obtain the package here.

Although ca can perform multiple correspondence analysis (more than two categorical variables), only simple correspondence analysis is covered here. See their article for details on multiple CA.

### Simple Correspondence Analysis
In the following example, A and B are categorical factors.

```
# Correspondence Analysis
library(ca)
mytable <- with(mydata, table(A,B)) # create a 2 way table
prop.table(mytable, 1) # row percentages
prop.table(mytable, 2) # column percentages
fit <- ca(mytable)
print(fit) # basic results 
summary(fit) # extended results 
plot(fit) # symmetric map
plot(fit, mass = TRUE, contrib = "absolute", map =
   "rowgreen", arrows = c(FALSE, TRUE)) # asymmetric map
```

The first graph is the standard symmetric representation of a simple correspondence analysis with rows and column represented by points.

correspondence analysis 1 click to view

Row points (column points) that are closer together have more similar column profiles (row profiles). Keep in mind that you can not interpret the distance between row and column points directly.

The second graph is asymmetric , with rows in the principal coordinates and columns in reconstructions of the standarized residuals. Additionally, mass is represented by points and columns are represented by arrows. Point intensity (shading) corresponds to the absolute contributions for the rows. This example is included to highlight some of the available options.

## Multidimensional Scaling

R provides functions for both classical and nonmetric multidimensional scaling. Assume that we have N objects measured on p numeric variables. We want to represent the distances among the objects in a parsimonious (and visual) way (i.e., a lower k-dimensional space).

### Classical MDS
You can perform a classical MDS using the cmdscale( ) function.

```
# Classical MDS
# N rows (objects) x p columns (variables)
# each row identified by a unique row name

d <- dist(mydata) # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
fit # view results

# plot solution 
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
  main="Metric MDS", type="n")
text(x, y, labels = row.names(mydata), cex=.7)
```

classical mds click to view

### Nonmetric MDS
Nonmetric MDS is performed using the isoMDS( ) function in the MASS package.

```
# Nonmetric MDS
# N rows (objects) x p columns (variables)
# each row identified by a unique row name

library(MASS)
d <- dist(mydata) # euclidean distances between the rows
fit <- isoMDS(d, k=2) # k is the number of dim
fit # view results

# plot solution 
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
  main="Nonmetric MDS", type="n")
text(x, y, labels = row.names(mydata), cex=.7)
```

nonmetric mds click to view

### Individual Difference Scaling
3-way or individual difference scaling can be completed using the indscal() function in the SensoMineR package. The smacof package offers a three way analysis of individual differences based on stress minimization of means of majorization.

## Cluster Analysis

R has an amazing variety of functions for cluster analysis. In this section, I will describe three of the many approaches: hierarchical agglomerative, partitioning, and model based. While there are no best solutions for the problem of determining the number of clusters to extract, several approaches are given below.

### Data Preparation
Prior to clustering data, you may want to remove or estimate missing data and rescale variables for comparability.

```
# Prepare Data
mydata <- na.omit(mydata) # listwise deletion of missing
mydata <- scale(mydata) # standardize variables
```

### Partitioning
K-means clustering is the most popular partitioning method. It requires the analyst to specify the number of clusters to extract. A plot of the within groups sum of squares by number of clusters extracted can help determine the appropriate number of clusters. The analyst looks for a bend in the plot similar to a scree test in factor analysis. See Everitt & Hothorn (pg. 251).

```
# Determine number of clusters
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(mydata, 
   centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
  ylab="Within groups sum of squares")

# K-Means Cluster Analysis
fit <- kmeans(mydata, 5) # 5 cluster solution
# get cluster means 
aggregate(mydata,by=list(fit$cluster),FUN=mean)
# append cluster assignment
mydata <- data.frame(mydata, fit$cluster)
```

A robust version of K-means based on mediods can be invoked by using pam( ) instead of kmeans( ). The function pamk( ) in the fpc package is a wrapper for pam that also prints the suggested number of clusters based on optimum average silhouette width.

### Hierarchical Agglomerative
There are a wide range of hierarchical clustering approaches. I have had good luck with Ward's method described below.

```
# Ward Hierarchical Clustering
d <- dist(mydata, method = "euclidean") # distance matrix
fit <- hclust(d, method="ward") 
plot(fit) # display dendogram
groups <- cutree(fit, k=5) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters 
rect.hclust(fit, k=5, border="red")
```

dendogram click to view

The pvclust( ) function in the pvclust package provides p-values for hierarchical clustering based on multiscale bootstrap resampling. Clusters that are highly supported by the data will have large p values. Interpretation details are provided Suzuki. Be aware that pvclust clusters columns, not rows. Transpose your data before using.

```
# Ward Hierarchical Clustering with Bootstrapped p values
library(pvclust)
fit <- pvclust(mydata, method.hclust="ward",
   method.dist="euclidean")
plot(fit) # dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit, alpha=.95)
```

clustering with p values click to view

### Model Based
Model based approaches assume a variety of data models and apply maximum likelihood estimation and Bayes criteria to identify the most likely model and number of clusters. Specifically, the Mclust( ) function in the mclust package selects the optimal model according to BIC for EM initialized by hierarchical clustering for parameterized Gaussian mixture models. (phew!). One chooses the model and number of clusters with the largest BIC. See help(mclustModelNames) to details on the model chosen as best.

```
# Model Based Clustering
library(mclust)
fit <- Mclust(mydata)
plot(fit) # plot results 
summary(fit) # display the best model
```

model based clustering cluster scatter plots click to view

### Plotting Cluster Solutions
It is always a good idea to look at the cluster results.

```
# K-Means Clustering with 5 clusters
fit <- kmeans(mydata, 5)

# Cluster Plot against 1st 2 principal components

# vary parameters for most readable graph
library(cluster) 
clusplot(mydata, fit$cluster, color=TRUE, shade=TRUE, 
   labels=2, lines=0)

# Centroid Plot against 1st 2 discriminant functions
library(fpc)
plotcluster(mydata, fit$cluster)
```
clusplot discriminant plot click to view

### Validating cluster solutions
The function cluster.stats() in the fpc package provides a mechanism for comparing the similarity of two cluster solutions using a variety of validation criteria (Hubert's gamma coefficient, the Dunn index and the corrected rand index)

```
# comparing 2 cluster solutions
library(fpc)
cluster.stats(d, fit1$cluster, fit2$cluster)
```

where d is a distance matrix among objects, and fit1$cluster and fit$cluster are integer vectors containing classification results from two different clusterings of the same data

## Tree-Based Models
Recursive partitioning is a fundamental tool in data mining. It helps us explore the stucture of a set of data, while developing easy to visualize decision rules for predicting a categorical (classification tree) or continuous (regression tree) outcome. This section briefly describes CART modeling, conditional inference trees, and random forests.

### CART Modeling via rpart
Classification and regression trees (as described by Brieman, Freidman, Olshen, and Stone) can be generated through the rpart package. Detailed information on rpart is available in An Introduction to Recursive Partitioning Using the RPART Routines. The general steps are provided below followed by two examples.

**1\. Grow the Tree**
To grow a tree, use
rpart(formula, data=, method=,control=) where

```
formula	is in the format 
outcome ~ predictor1+predictor2+predictor3+ect.
data=	specifies the data frame
method=	"class" for a classification tree 
"anova" for a regression tree
control=	optional parameters for controlling tree growth. For example, control=rpart.control(minsplit=30, cp=0.001) requires that the minimum number of observations in a node be 30 before attempting a split and that a split must decrease the overall lack of fit by a factor of 0.001 (cost complexity factor) before being attempted.
```

**2\. Examine the results**

The following functions help us to examine the results.

```
printcp(fit)	display cp table
plotcp(fit)	plot cross-validation results
rsq.rpart(fit)	plot approximate R-squared and relative error for different splits (2 plots). labels are only appropriate for the "anova" method.
print(fit)	print results
summary(fit)	detailed results including surrogate splits
plot(fit)	plot decision tree
text(fit)	label the decision tree plot
post(fit, file=)	create postscript plot of decision tree
In trees created by rpart( ), move to the LEFT branch when the stated condition is true (see the graphs below).
```

**3\. prune tree**
Prune back the tree to avoid overfitting the data. Typically, you will want to select a tree size that minimizes the cross-validated error, the xerror column printed by printcp( ).

Prune the tree to the desired size using
prune(fit, cp= )

Specifically, use printcp( ) to examine the cross-validated error results, select the complexity parameter associated with minimum error, and place it into the prune( ) function. Alternatively, you can use the code fragment

     fit$cptable[which.min(fit$cptable[,"xerror"]),"CP"]

to automatically select the complexity parameter associated with the smallest cross-validated error. Thanks to HSAUR for this idea.

### Classification Tree example
Let's use the data frame kyphosis to predict a type of deformation (kyphosis) after surgery, from age in months (Age), number of vertebrae involved (Number), and the highest vertebrae operated on (Start).

```
# Classification Tree with rpart
library(rpart)

# grow tree 
fit <- rpart(Kyphosis ~ Age + Number + Start,
   method="class", data=kyphosis)

printcp(fit) # display the results 
plotcp(fit) # visualize cross-validation results 
summary(fit) # detailed summary of splits

# plot tree 
plot(fit, uniform=TRUE, 
   main="Classification Tree for Kyphosis")
text(fit, use.n=TRUE, all=TRUE, cex=.8)

# create attractive postscript plot of tree 
post(fit, file = "c:/tree.ps", 
   title = "Classification Tree for Kyphosis")
```

cp Plot Classification Tree Classification Tree in Postscript click to view

```
# prune the tree 
pfit<- prune(fit, cp=   fit$cptable[which.min(fit$cptable[,"xerror"]),"CP"])

# plot the pruned tree 
plot(pfit, uniform=TRUE, 
   main="Pruned Classification Tree for Kyphosis")
text(pfit, use.n=TRUE, all=TRUE, cex=.8)
post(pfit, file = "c:/ptree.ps", 
   title = "Pruned Classification Tree for Kyphosis")
```

Pruned Classificaiton Tree Pruned Classification Tree in Postscript click to view

### Regression Tree example
In this example we will predict car mileage from price, country, reliability, and car type. The data frame is cu.summary.

```
# Regression Tree Example
library(rpart)

# grow tree 
fit <- rpart(Mileage~Price + Country + Reliability + Type, 
   method="anova", data=cu.summary)

printcp(fit) # display the results 
plotcp(fit) # visualize cross-validation results 
summary(fit) # detailed summary of splits

# create additional plots 
par(mfrow=c(1,2)) # two plots on one page 
rsq.rpart(fit) # visualize cross-validation results   

# plot tree 
plot(fit, uniform=TRUE, 
   main="Regression Tree for Mileage ")
text(fit, use.n=TRUE, all=TRUE, cex=.8)

# create attractive postcript plot of tree 
post(fit, file = "c:/tree2.ps", 
   title = "Regression Tree for Mileage ")
```

cp plot for regression tree rsquare plot for regression treeregression tree Regressio Tree in Post Script click to view

```
# prune the tree 
pfit<- prune(fit, cp=0.01160389) # from cptable   

# plot the pruned tree 
plot(pfit, uniform=TRUE, 
   main="Pruned Regression Tree for Mileage")
text(pfit, use.n=TRUE, all=TRUE, cex=.8)
post(pfit, file = "c:/ptree2.ps", 
   title = "Pruned Regression Tree for Mileage")
```

It turns out that this produces the same tree as the original.

### Conditional inference trees via party
The party package provides nonparametric regression trees for nominal, ordinal, numeric, censored, and multivariate responses. party: A laboratory for recursive partitioning, provides details.

You can create a regression or classification tree via the function

ctree(formula, data=)
The type of tree created will depend on the outcome variable (nominal factor, ordered factor, numeric, etc.). Tree growth is based on statistical stopping rules, so pruning should not be required.

The previous two examples are re-analyzed below.

```
# Conditional Inference Tree for Kyphosis
library(party)
fit <- ctree(Kyphosis ~ Age + Number + Start, 
   data=kyphosis)
plot(fit, main="Conditional Inference Tree for Kyphosis")
```

Condiitional Inference Tree for Kyphosis click to view

```
# Conditional Inference Tree for Mileage
library(party)
fit2 <- ctree(Mileage~Price + Country + Reliability + Type, 
   data=na.omit(cu.summary))
```

Conditional Inference Tree for Mileage click to view

### Random Forests
Random forests improve predictive accuracy by generating a large number of bootstrapped trees (based on random samples of variables), classifying a case using each tree in this new "forest", and deciding a final predicted outcome by combining the results across all of the trees (an average in regression, a majority vote in classification). Breiman and Cutler's random forest approach is implimented via the randomForest package.

Here is an example.

```
# Random Forest prediction of Kyphosis data
library(randomForest)
fit <- randomForest(Kyphosis ~ Age + Number + Start,   data=kyphosis)
print(fit) # view results 
importance(fit) # importance of each predictor
```

For more details see the comprehensive Random Forest website.

### Going Further
This section has only touched on the options available. To learn more, see the CRAN Task View on Machine & Statistical Learning.

