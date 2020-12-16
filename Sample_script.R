
# This program will run an analysis that has 1 penalized function of time. The
# other script "Sample_script_Pcov.R" runs a model with a penalized function
# for multiple groups.


#Load in the programs.  Requires the following packages ggplot2, nlme, MASS, splines, Matrix

source("Programs/Programs.R")


####### Reading in the necessary data ########

all_data <- read.csv("Data/Example_data.csv")
### all_data has the following columns:
#    - country: country level identifier (must be a factor).
#    - year
#    - Y: prevalence on (0,1) scale.
#    - SE_var: Estimated sampling standard error of Y
#    - Region: the region of the country  (must be a factor).
#    - source_cov: a covariate used to denote that a the prevalence came from a
#              source that may have systematic bias.
#    - source_cov2: a covariate used to denote that a the prevalence came from a
#              different source which may also have systematic bias.
#    - cov1: Example Covariate 1
#    - cov2: Example Covariate 2


# We'll now transform the data from (0,1) to (-Inf,Inf) using a logit transformation.
# The standard errors are transformed using the delta method.
all_data$SE_var <- all_data$SE_var*(1/all_data$Y + 1/(1-all_data$Y))
all_data$Y <- log(all_data$Y/(1-all_data$Y))


# Define all items needed for analysis:
# Main data file
data <- all_data
# Formula for the fixed effect
formula <- Y ~ source_cov + source_cov2 + cov1 + cov2 + Region
# Formula for the random effects
re_formula <- ~ 1 + year
# Names of area, time and standard error variables.
area_var <- "country"
time_var <- "year"
SE_var <- "SE_var"
# Names of variables to set to zero for prediction 
zero_covs <- c("source_cov", "source_cov2")



# Location of penalized splines
DF_P <- seq(1993,2020,1)
# Location of random splines (note: the last 2 random splines are removed to 
# avoid edge effects)
DF_R <- quantile( data$year, probs = c(0.5))


# Specifying the covariance matrix of the random splines (the covariance matrix
# of the other random effects is always unstructured)  
#   - cov_mat="CS" (default) -> compound symmetric covariance matrix
#   - cov_mat="UN"  -> unstructured symmetric covariance matrix
#   - cov_mat="VC"  -> Diagonal only covariance matrix
cov_mat <- "CS"


Estimation <- cmnpe(formula, re_formula, data = data, area_var = area_var, 
                    time_var = time_var, SE_var = SE_var, DF_P = DF_P, DF_R =  
                    DF_R, cov_mat = cov_mat, zero_covs = zero_covs)

# Estimated degrees of freedom
Estimation$df
# AICc model fit criteria
2*Estimation$df -2*Estimation$model$logLik

####Summary of linear model####
summary(Estimation$model)

##### Open and run Plot_script.R to output results and create figures ######



