# Statistics

## Descriptive Statistics

R provides a wide range of functions for obtaining summary statistics. One method of obtaining descriptive statistics is to use the sapply( ) function with a specified summary statistic.

```
# get means for variables in data frame mydata
# excluding missing values 
sapply(mydata, mean, na.rm=TRUE)
```

Possible functions used in sapply include mean, sd, var, min, max, median, range, and quantile.

There are also numerous R functions designed to provide a range of descriptive statistics at once. For example

```
# mean,median,25th and 75th quartiles,min,max
summary(mydata)

# Tukey min,lower-hinge, median,upper-hinge,max
fivenum(x)
```

Using the Hmisc package

```
library(Hmisc)
describe(mydata) 
# n, nmiss, unique, mean, 5,10,25,50,75,90,95th percentiles 
# 5 lowest and 5 highest scores
```

Using the pastecs package

```
library(pastecs)
stat.desc(mydata) 
# nbr.val, nbr.null, nbr.na, min max, range, sum, 
# median, mean, SE.mean, CI.mean, var, std.dev, coef.var
```

Using the psych package

```
library(psych)
describe(mydata)
# item name ,item number, nvalid, mean, sd, 
# median, mad, min, max, skew, kurtosis, se
```

Summary Statistics by Group
A simple way of generating summary statistics by grouping variable is available in the psych package.

```
library(psych)
describe.by(mydata, group,...)
```

The doBy package provides much of the functionality of SAS PROC SUMMARY. It defines the desired table using a model formula and a function. Here is a simple example.

```
library(doBy)
summaryBy(mpg + wt ~ cyl + vs, data = mtcars, 
  FUN = function(x) { c(m = mean(x), s = sd(x)) } )
# produces mpg.m wt.m mpg.s wt.s for each 
# combination of the levels of cyl and vs
```

See also: aggregating data.

## Frequencies and Crosstabs

This section describes the creation of frequency and contingency tables from categorical variables, along with tests of independence, measures of association, and methods for graphically displaying results.

Generating Frequency Tables
R provides many methods for creating frequency and contingency tables. Three are described below. In the following examples, assume that A, B, and C represent categorical variables.

table
You can generate frequency tables using the table( ) function, tables of proportions using the prop.table( ) function, and marginal frequencies using margin.table( ).

```
# 2-Way Frequency Table 
attach(mydata)
mytable <- table(A,B) # A will be rows, B will be columns 
mytable # print table 

margin.table(mytable, 1) # A frequencies (summed over B) 
margin.table(mytable, 2) # B frequencies (summed over A)

prop.table(mytable) # cell percentages
prop.table(mytable, 1) # row percentages 
prop.table(mytable, 2) # column percentages
```

table( ) can also generate multidimensional tables based on 3 or more categorical variables. In this case, use the ftable( ) function to print the results more attractively.

```
# 3-Way Frequency Table 
mytable <- table(A, B, C) 
ftable(mytable)
```

Table ignores missing values. To include NA as a category in counts, include the table option exclude=NULL if the variable is a vector. If the variable is a factor you have to create a new factor using newfactor <- factor(oldfactor, exclude=NULL).

xtabs
The xtabs( ) function allows you to create crosstabulations using formula style input.

```
# 3-Way Frequency Table
mytable <- xtabs(~A+B+c, data=mydata)
ftable(mytable) # print table 
summary(mytable) # chi-square test of indepedence
```

If a variable is included on the left side of the formula, it is assumed to be a vector of frequencies (useful if the data have already been tabulated).

Crosstable
The CrossTable( ) function in the gmodels package produces crosstabulations modeled after PROC FREQ in SAS or CROSSTABS in SPSS. It has a wealth of options.

```
# 2-Way Cross Tabulation
library(gmodels)
CrossTable(mydata$myrowvar, mydata$mycolvar)
```

There are options to report percentages (row, column, cell), specify decimal places, produce Chi-square, Fisher, and McNemar tests of independence, report expected and residual values (pearson, standardized, adjusted standardized), include missing values as valid, annotate with row and column titles, and format as SAS or SPSS style output! 
See help(CrossTable) for details.

Tests of Independence
Chi-Square Test
For 2-way tables you can use chisq.test(mytable) to test independence of the row and column variable. By default, the p-value is calculated from the asymptotic chi-squared distribution of the test statistic. Optionally, the p-value can be derived via Monte Carlo simultation.

Fisher Exact Test
fisher.test(x) provides an exact test of independence. x is a two dimensional contingency table in matrix form.

Mantel-Haenszel test
Use the mantelhaen.test(x) function to perform a Cochran-Mantel-Haenszel chi-squared test of the null hypothesis that two nominal variables are conditionally independent in each stratum, assuming that there is no three-way interaction. x is a 3 dimensional contingency table, where the last dimension refers to the strata.

Loglinear Models
You can use the loglm( ) function in the MASS package to produce log-linear models. For example, let's assume we have a 3-way contingency table based on variables A, B, and C.

library(MASS)
mytable <- xtabs(~A+B+C, data=mydata)

We can perform the following tests:

Mutual Independence: A, B, and C are pairwise independent.
loglm(~A+B+C, mytable)

Partial Independence: A is partially independent of B and C (i.e., A is independent of the composite variable BC).
loglin(~A+B+C+B*C, mytable)

Conditional Independence: A is independent of B, given C.
loglm(~A+B+C+A*C+B*C, mytable)
No Three-Way Interaction
loglm(~A+B+C+A*B+A*C+B*C, mytable)

Martin Theus and Stephan Lauer have written an excellent article on Visualizing Loglinear Models, using mosaic plots.

Measures of Association
The assocstats(mytable) function in the vcd package calculates the phi coefficient, contingency coefficient, and Cramer's V for an rxc table. The kappa(mytable) function in the vcd package calculates Cohen's kappa and weighted kappa for a confusion matrix. See Richard Darlington's article on Measures of Association in Crosstab Tables for an excellent review of these statistics.

Visualizing results
Use bar and pie charts for visualizing frequencies in one dimension.

Use the vcd package for visualizing relationships among categorical data (e.g. mosaic and association plots).

Use the ca package for correspondence analysis (visually exploring relationships between rows and columns in contingency tables).

To practice making these charts, try the data visualization course at DataCamp.

Converting Frequency Tables to an "Original" Flat file
Finally, there may be times that you wil need the original "flat file" data frame rather than the frequency table. Marc Schwartz has provided code on the Rhelp mailing list for converting a table back into a data frame.

## Correlations

You can use the cor( ) function to produce correlations and the cov( ) function to produces covariances.

A simplified format is cor(x, use=, method= ) where

```
Option	Description
x	Matrix or data frame
use	Specifies the handling of missing data. Options are all.obs (assumes no missing data - missing data will produce an error), complete.obs (listwise deletion), and pairwise.complete.obs (pairwise deletion)
method	Specifies the type of correlation. Options are pearson, spearman or kendall.
```

```
# Correlations/covariances among numeric variables in 
# data frame mtcars. Use listwise deletion of missing data. 
cor(mtcars, use="complete.obs", method="kendall") 
cov(mtcars, use="complete.obs")
```

Unfortunately, neither cor( ) or cov( ) produce tests of significance, although you can use the cor.test( ) function to test a single correlation coefficient.

The rcorr( ) function in the Hmisc package produces correlations/covariances and significance levels for pearson and spearman correlations. However, input must be a matrix and pairwise deletion is used.

```
# Correlations with significance levels
library(Hmisc)
rcorr(x, type="pearson") # type can be pearson or spearman

#mtcars is a data frame 
rcorr(as.matrix(mtcars))
```

You can use the format cor(X, Y) or rcorr(X, Y) to generate correlations between the columns of X and the columns of Y. This similar to the VAR and WITH commands in SAS PROC CORR.

```
# Correlation matrix from mtcars
# with mpg, cyl, and disp as rows 
# and hp, drat, and wt as columns 
x <- mtcars[1:3]
y <- mtcars[4:6]
cor(x, y)
```

##Other Types of Correlations

```
# polychoric correlation
# x is a contingency table of counts
library(polycor)
polychor(x) 

# heterogeneous correlations in one matrix 
# pearson (numeric-numeric), 
# polyserial (numeric-ordinal), 
# and polychoric (ordinal-ordinal)
# x is a data frame with ordered factors 
# and numeric variables
library(polycor)
hetcor(x) 


# partial correlations
library(ggm)
data(mydata)
pcor(c("a", "b", "x", "y", "z"), var(mydata))
# partial corr between a and b controlling for x, y, z
```

##Visualizing Correlations
Use corrgram( ) to plot correlograms .

Use the pairs() or splom( ) to create scatterplot matrices.

## t-tests

The t.test( ) function produces a variety of t-tests. Unlike most statistical packages, the default assumes unequal variance and applies the Welsh df modification.

```
# independent 2-group t-test
t.test(y~x) # where y is numeric and x is a binary factor

# independent 2-group t-test
t.test(y1,y2) # where y1 and y2 are numeric

# paired t-test
t.test(y1,y2,paired=TRUE) # where y1 & y2 are numeric

# one sample t-test
t.test(y,mu=3) # Ho: mu=3
```

You can use the var.equal = TRUE option to specify equal variances and a pooled variance estimate. You can use the alternative="less" or alternative="greater" option to specify a one tailed test.

Nonparametric and resampling alternatives to t-tests are available.

Visualizing Results
Use box plots or density plots to visualize group differences.

## Nonparametric Statistics

R provides functions for carrying out Mann-Whitney U, Wilcoxon Signed Rank, Kruskal Wallis, and Friedman tests.


```
# independent 2-group Mann-Whitney U Test 
wilcox.test(y~A) 
# where y is numeric and A is A binary factor

# independent 2-group Mann-Whitney U Test
wilcox.test(y,x) # where y and x are numeric

# dependent 2-group Wilcoxon Signed Rank Test 
wilcox.test(y1,y2,paired=TRUE) # where y1 and y2 are numeric

# Kruskal Wallis Test One Way Anova by Ranks 
kruskal.test(y~A) # where y1 is numeric and A is a factor

# Randomized Block Design - Friedman Test 
friedman.test(y~A|B)
# where y are the data values, A is a grouping factor
# and B is a blocking factor
```

For the wilcox.test you can use the alternative="less" or alternative="greater" option to specify a one tailed test.

Parametric and resampling alternatives are available.

The package pgirmess provides nonparametric multiple comparisons. (Note: This package has been withdrawn but is still available in the CRAN archives.)

```
library(npmc)
npmc(x) 
# where x is a data frame containing variable 'var' 
# (response variable) and 'class' (grouping variable)
```

Visualizing Results
Use box plots or density plots to visual group differences.

## Multiple Regression
R provides comprehensive support for multiple linear regression. The topics below are provided in order of increasing complexity.

### Fitting the Model

```
# Multiple Linear Regression Example 
fit <- lm(y ~ x1 + x2 + x3, data=mydata)
summary(fit) # show results

# Other useful functions 
coefficients(fit) # model coefficients
confint(fit, level=0.95) # CIs for model parameters 
fitted(fit) # predicted values
residuals(fit) # residuals
anova(fit) # anova table 
vcov(fit) # covariance matrix for model parameters 
influence(fit) # regression diagnostics
```

### Diagnostic Plots
Diagnostic plots provide checks for heteroscedasticity, normality, and influential observerations.

```
# diagnostic plots 
layout(matrix(c(1,2,3,4),2,2)) # optional 4 graphs/page 
plot(fit)
```

regression diagnostic plots click to view

For a more comprehensive evaluation of model fit see regression diagnostics or the exercises in this interactive course on regression.

### Comparing Models
You can compare nested models with the anova( ) function. The following code provides a simultaneous test that x3 and x4 add to linear prediction above and beyond x1 and x2.

```
# compare models
fit1 <- lm(y ~ x1 + x2 + x3 + x4, data=mydata)
fit2 <- lm(y ~ x1 + x2)
anova(fit1, fit2)
```

### Cross Validation
You can do K-Fold cross-validation using the cv.lm( ) function in the DAAG package.

```
# K-fold cross-validation
library(DAAG)
cv.lm(df=mydata, fit, m=3) # 3 fold cross-validation
```

Sum the MSE for each fold, divide by the number of observations, and take the square root to get the cross-validated standard error of estimate.

You can assess R2 shrinkage via K-fold cross-validation. Using the crossval() function from the bootstrap package, do the following:

```
# Assessing R2 shrinkage using 10-Fold Cross-Validation 

fit <- lm(y~x1+x2+x3,data=mydata) 

library(bootstrap)
# define functions 
theta.fit <- function(x,y){lsfit(x,y)}
theta.predict <- function(fit,x){cbind(1,x)%*%fit$coef} 

# matrix of predictors
X <- as.matrix(mydata[c("x1","x2","x3")])
# vector of predicted values
y <- as.matrix(mydata[c("y")]) 

results <- crossval(X,y,theta.fit,theta.predict,ngroup=10)
cor(y, fit$fitted.values)**2 # raw R2 
cor(y,results$cv.fit)**2 # cross-validated R2
```

### Variable Selection
Selecting a subset of predictor variables from a larger set (e.g., stepwise selection) is a controversial topic. You can perform stepwise selection (forward, backward, both) using the stepAIC( ) function from the MASS package. stepAIC( ) performs stepwise model selection by exact AIC.

```
# Stepwise Regression
library(MASS)
fit <- lm(y~x1+x2+x3,data=mydata)
step <- stepAIC(fit, direction="both")
step$anova # display results
```

Alternatively, you can perform all-subsets regression using the leaps( ) function from the leaps package. In the following code nbest indicates the number of subsets of each size to report. Here, the ten best models will be reported for each subset size (1 predictor, 2 predictors, etc.).

```
# All Subsets Regression
library(leaps)
attach(mydata)
leaps<-regsubsets(y~x1+x2+x3+x4,data=mydata,nbest=10)
# view results 
summary(leaps)
# plot a table of models showing variables in each model.
# models are ordered by the selection statistic.
plot(leaps,scale="r2")
# plot statistic by subset size 
library(car)
subsets(leaps, statistic="rsq")
```

all subsets regression 1 all subsets regression 2 click to view

Other options for plot( ) are bic, Cp, and adjr2. Other options for plotting with 
subset( ) are bic, cp, adjr2, and rss.

### Relative Importance
The relaimpo package provides measures of relative importance for each of the predictors in the model. See help(calc.relimp) for details on the four measures of relative importance provided.

```
# Calculate Relative Importance for Each Predictor
library(relaimpo)
calc.relimp(fit,type=c("lmg","last","first","pratt"),
   rela=TRUE)

# Bootstrap Measures of Relative Importance (1000 samples) 
boot <- boot.relimp(fit, b = 1000, type = c("lmg", 
  "last", "first", "pratt"), rank = TRUE, 
  diff = TRUE, rela = TRUE)
booteval.relimp(boot) # print result
plot(booteval.relimp(boot,sort=TRUE)) # plot result
```

bootstrapped relative importance click to view

### Graphic Enhancements
The car package offers a wide variety of plots for regression, including added variable plots, and enhanced diagnostic and Scatterplots.

### Going Further
#### Nonlinear Regression
The nls package provides functions for nonlinear regression. See John Fox's Nonlinear Regression and Nonlinear Least Squares for an overview. Huet and colleagues' Statistical Tools for Nonlinear Regression: A Practical Guide with S-PLUS and R Examples is a valuable reference book.

#### Robust Regression
There are many functions in R to aid with robust regression. For example, you can perform robust regression with the rlm( ) function in the MASS package. John Fox's (who else?) Robust Regression provides a good starting overview. The UCLA Statistical Computing website has Robust Regression Examples.

The robust package provides a comprehensive library of robust methods, including regression. The robustbase package also provides basic robust statistics including model selection methods. And David Olive has provided an detailed online review of Applied Robust Statistics with sample R code.

