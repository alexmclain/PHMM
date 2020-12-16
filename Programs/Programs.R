require(ggplot2)
require(nlme)
require(MASS)
require(splines)
require(Matrix)



cmnpe <- function(formula, re_formula, data, area_var, time_var, SE_var, DF_P, DF_R, 
                  Pen_var = NULL, Pen_type = "Single", cov_mat="CS", 
                  zero_covs = NULL){
  
  if(missing(formula)){formula <- Y ~ 1}
  if(missing(re_formula)){re_formula <- ~ 1}
  
  ###### Order the data by country and year (required) #######
  cnames <- colnames(data)
  data <- data[order(data[,cnames==area_var],data[,cnames==time_var]),]
  
  ##### Run estimation. #####  
  Estimation <- model_fit(formula, re_formula, data, area_var, time_var, 
                          SE_var, DF_P, DF_R, cov_mat, Pen_var, Pen_type, 
                          zero_covs)
  
  return(Estimation)
}



model_fit <- function(formula, re_formula, data, area_var, time_var, SE_var, DF_P, 
                      DF_R, cov_mat, Pen_var, Pen_type, zero_covs){
  
  cnames <- colnames(data)
  # Making the penalized splines
  year <- data[ , cnames == time_var]
  B.knots <- c(min(year)-1,max(year)+2)
  # Order of penalization
  q.order <-  2
  bsplines <- bs( year, knots = DF_P, Boundary.knots = B.knots, intercept= F)
  B.matrix <- as.matrix(bsplines)
  b.col <- ncol(B.matrix)
  b.row <- nrow(B.matrix)
  D.matrix <- diff(diag(b.col),differences = q.order)
  P.matrix <- crossprod(D.matrix)
  Pi.matrix.svd <- svd(P.matrix)
  Ui.matrix <-(Pi.matrix.svd$u)[,1:(b.col-q.order)]  		# matrix of eigenvectors
  eigen.vec <-(Pi.matrix.svd$d)[1:(b.col-q.order)]   		# vector of eigenvalues
  Sigmai.inv <- diag(1/sqrt(eigen.vec))  	  				    # diagonal matrix of eigenvalues
  Zi.matrix <<- B.matrix%*%Ui.matrix%*%Sigmai.inv        # instead of Z=BU
  
  # Checking to see if the time-var is included in the formula, if not adding it
  mod_factors <- colnames(attr(terms(formula),"factors"))
  if(! any( mod_factors == time_var ) ){
    char_form <- paste( as.character( formula ) )
    formula   <- as.formula( paste( char_form[2], char_form[1], char_form[3], 
                                    "+", time_var ) )
    mod_factors <- c(mod_factors,time_var)
  }
  
  # Making random splines
  rand_spline <- as.matrix( bs( year, knots = DF_R, Boundary.knots = B.knots, 
                                intercept = F) )
  rand_spline <- rand_spline[,-c( (ncol(rand_spline)-1) : (ncol(rand_spline))) ]
  
  
  
  colnames(rand_spline) <- paste( "Rand", 1:(ncol(rand_spline)), sep = "")
  colnames(Zi.matrix) <- paste( "Zi.matrix", 1:(ncol(Zi.matrix)), sep = "")
  
  data <- data.frame(data, Zi.matrix, rand_spline, intercept = 1)
  cnames <- colnames(data)
  data$area_var2 <- data[ , cnames == area_var]
  
  d.names <<- cnames
  
  # Setting the weight and RE formulas
  data[ , cnames == SE_var] <- (data[ , cnames == SE_var])^2
  weight_form = varSum(form = as.formula(paste("~",SE_var)))
  
  ctrl <- lmeControl(opt = c("nlminb"),maxIter = 1000, msMaxIter = 1000, 
                     niterEM = 500,  msMaxEval = 2000,tolerance = 1e-6)
  
  Zi.formula <- as.formula("~0 + Zi.matrix")
  Zi.structure = pdIdent(Zi.formula)
  if( !is.null(Pen_var)){
    if( !is.factor( data[ , cnames == Pen_var])){
      data[ , cnames == Pen_var] <- as.factor( data[ , cnames == Pen_var])
    }
    
    if(Pen_type == "Single"){
      Pcov_data <<- model.matrix( as.formula( paste("~",Pen_var)) , data = data)[ , -1]
      Zi.formula <- as.formula("~0 + Pcov_data:Zi.matrix")
      Zi.structure = pdIdent(Zi.formula)
    }else{
      Pcov_data <<- model.matrix( as.formula( paste("~0+",Pen_var)) , data = data)
      Pen_vals <- levels( data[ , cnames == Pen_var])
      L <- length( Pen_vals)
      Zi.formula <- list()
      for(j in 1:L){
        assign( paste("int",j,sep = ""), Pcov_data[,j]*Zi.matrix, envir = .GlobalEnv)
        Zi.formula[[(j)]] <- pdIdent( as.formula( paste("~0+int",j,sep = "")))
      }
      Zi.structure = pdBlocked(Zi.formula)
    }
  }
  
  
  rand_formula <- "~0"
  for(i in 1:ncol(rand_spline)){rand_formula <- paste(rand_formula,"+",colnames(rand_spline)[i])}
  rand_formula <- as.formula( rand_formula )
  
  
  random_list <- list(temp1 = Zi.structure, temp2 = pdSymm(re_formula),
                      temp3 = pdCompSymm(rand_formula))
  if(cov_mat == "UN"){
    random_list <- list(temp1 = Zi.structure, temp2 = pdSymm(re_formula),
                        temp3 = pdSymm(rand_formula))
  }
  if(cov_mat == "VC"){
    random_list <- list(temp1 = Zi.structure, temp2 = pdSymm(re_formula),
                        temp3 = pdDiag(rand_formula))
  }
  names(random_list) <- c("intercept", area_var, "area_var2")
  
  # Estimating the mixed model
  fitted_model <- lme( formula, data = data, random = random_list, weights = 
                         weight_form, control = ctrl , na.action = na.omit)
  
  na_default <- options("na.action")
  
  # Performing the prediction
  pred_results <- prediction( data, formula, re_formula, Zi.formula, rand_formula, 
                              fitted_model, SE_var, area_var, time_var, zero_covs,
                              Pen_type, Pen_var)
  options(na.action = na_default$na.action)
  
  mis_ordered <- pred_results$newdata
  SIGMA_T_mis_ordered <- pred_results$SIGMA_T
  SIGMA_Y_mis_ordered <- pred_results$SIGMA_Y
  
  ### Reordering the prediction data and covariance matrices
  mis_ordered[ , colnames(mis_ordered) == area_var ] <- 
    as.character( mis_ordered[ , colnames(mis_ordered) == area_var ])
  order_list <- order( mis_ordered[ , colnames(mis_ordered) == area_var ], 
                       mis_ordered[ , colnames(mis_ordered) == time_var ])
  SIGMA_T <- SIGMA_T_mis_ordered[ order_list, ]
  SIGMA_T <- SIGMA_T[ , order_list ]
  SIGMA_Y <- SIGMA_Y_mis_ordered[ order_list, ]
  SIGMA_Y <- SIGMA_Y[ , order_list ]
  ordered <- mis_ordered[ order_list,]
  
  # Calculate sd's to get CI's and PI's
  sigma_T_est <- sqrt( ordered$sigma_all )
  sigma_Y_est <- sqrt( ordered$sigma_Y_all )
  lower_CI <- ordered$pred - 1.96*sigma_T_est
  upper_CI <- ordered$pred + 1.96*sigma_T_est
  lower_PI <- ordered$pred - 1.96*sigma_Y_est
  upper_PI <- ordered$pred + 1.96*sigma_Y_est
  
  # Combining results to export.
  data <- data[order( as.character( data[ , colnames(data) == area_var]), 
                      data[ , colnames(data) == time_var]),]
  pred_data <- data.frame(data$area_var2, 
                          data[ , colnames(data) == time_var], 
                          data[ , colnames(data) == as.character(formula)[2]], 
                          data[ , colnames(data) == SE_var], 
                          pred = ordered$pred, lower_CI = lower_CI, upper_CI = upper_CI, 
                          lower_PI = lower_PI, upper_PI = upper_PI, 
                          sigma_T_est = sigma_T_est, sigma_Y_est = sigma_Y_est,
                          resid = ordered$resid, pred_fixed = ordered$pred_all_fixed, 
                          pred_fixpen = ordered$pred_all_fixpen)
  colnames(pred_data)[1:4] <- c(area_var, time_var, as.character(formula)[2], SE_var)
  
  remove(Zi.matrix, pos = ".GlobalEnv")
  if( !is.null(Pen_var)){
    remove(Pcov_data, pos = ".GlobalEnv")
    remove(list = ls(pattern = "int"), pos = ".GlobalEnv")
  }
  
  fin_res <- list(model = fitted_model, pred_data = pred_data, full_pars = 
                    pred_results$b, COV_beta_b = pred_results$COV_beta_b, df = 
                    pred_results$df, SIGMA_T = SIGMA_T, SIGMA_Y = SIGMA_Y, DF_P = 
                    DF_P, DF_R = DF_R, B.knots = B.knots)
  
  return(fin_res)
}


prediction <- function(data, formula, re_formula, Zi.formula, rand_formula, 
                       fitted_model, SE_var, area_var, time_var, zero_covs, 
                       Pen_type, Pen_var = NULL){
  
  options(na.action='na.pass')
  # Getting the X and Z matrices for observations with data
  cnames <- colnames(data)
  Y             <- data[ , cnames == as.character(formula)[2] ]
  SE_vec        <- data[ , cnames == SE_var]
  country       <- data$area_var2
  X_matrix <- model.matrix(formula,data=data)
  Z_matrix <- Z.matrix( Zi.formula, re_formula, rand_formula, data, NA_ind = NULL,
                        Pen_type)
  NA_ind <- complete.cases(Y, SE_vec, country, X_matrix, Z_matrix)
  
  ### is there a way to do the NA check with the mixed model results?
  
  data_obs <- data[ NA_ind, ]
  Y_obs             <- data_obs[ , cnames == as.character(formula)[2] ]
  SE_vec_obs        <- data_obs[ , cnames == SE_var]
  X_matrix_obs <- model.matrix( formula, data = data_obs)
  factor_list <- unlist( lapply( data_obs[1, ], is.factor) )
  if( any(factor_list)){ 
    if( sum(factor_list) == 1){ 
      data_obs[ , factor_list] <- as.character( data_obs[ , factor_list ])
    }else{ 
      data_obs[ , factor_list] <- apply( data_obs[ , factor_list ], 2, as.character)
    }
  }
  Z_matrix_obs <- Z.matrix( Zi.formula, re_formula, rand_formula, data_obs, 
                            NA_ind, Pen_type)
  # Checking Pen_var to see if it's approproate. 
  if( !is.null(Pen_var)){
    t_N <- length( data_obs[ , cnames == Pen_var] )
    tab_vals <- table( data_obs[ , cnames == Pen_var])
    if( min(tab_vals) < max( c( 0.01*t_N,10))){stop(paste( Pen_var, "has at least 1 group that has too few observations to use as a penalized variable."))}
  }
  
  # Estimating all fixed and random effects with covariance
  beta_results <- beta.fun( Y_obs, X_matrix_obs, Z_matrix_obs, 
                            SE_vec_obs, fitted_model)
  full_b <- beta_results$b
  S <- beta_results$S
  
  # Predicting for areas with data
  data[ , cnames %in% zero_covs] <- 0
  data_obs_c <- data[ data$area_var2 %in% data_obs$area_var2, ]
  SE_vec_obs_c        <- median( data_obs_c[ , cnames == SE_var], na.rm = TRUE)
  X_matrix_obs_c <- model.matrix( formula, data = data_obs_c)
  factor_list <- unlist( lapply( data_obs_c[1,], is.factor))
  if( any(factor_list)){ 
    if( sum(factor_list) == 1){ 
      data_obs_c[ , factor_list] <- as.character( data_obs_c[ , factor_list ])
    }else{
      data_obs_c[ , factor_list] <- apply( data_obs_c[ , factor_list ], 2, as.character)
    }
  }
  Z_matrix_obs_c <- Z.matrix( Zi.formula, re_formula, rand_formula, data_obs_c, 
                              data$area_var2 %in% data_obs$area_var2, Pen_type)
  pred_mat_new <- pred( X_matrix_obs_c, Z_matrix_obs_c, SE_vec_obs_c, 
                        fitted_model, full_b, S)
  
  # Predicting for areas with no data
  data_noobs_c <- data[ !( data$area_var2 %in% data_obs$area_var2 ), ]
  SE_vec_noobs_c <- median( SE_vec, na.rm = TRUE)
  N2             <- length( unique( data_noobs_c$area_var2))
  X_matrix_noobs_c <- X_matrix[ !( data$area_var2 %in% data_obs$area_var2 ), ]
  X_matrix_noobs_c[ , colnames(X_matrix_noobs_c) %in% zero_covs] <- 0
  factor_list <- unlist( lapply( data_noobs_c[ 1, ], is.factor))
  if( any(factor_list)){ 
    if( sum(factor_list) == 1){ 
      data_noobs_c[ , factor_list] <- as.character( data_noobs_c[ , factor_list ])
    }else{
      data_noobs_c[ , factor_list] <- apply( data_noobs_c[ , factor_list ], 2, as.character)
    }
  }
  Z_matrix_noobs_c <- Z.matrix(Zi.formula, re_formula, rand_formula, data_noobs_c, 
                               NA_ind = !( data$area_var2 %in% data_obs$area_var2),
                               Pen_type)
  num_pen <- dim(as.matrix(fitted_model$modelStruct$reStruct$intercept))[2]
  pred_mat_nob <- pred.noRand(X_matrix_noobs_c, Z_matrix_noobs_c, SE_vec_noobs_c, 
                              fitted_model, N2, full_b, S, num_pen)
  
  # Merging prediction results together
  data_w_out_obs_c <- data.frame(data_obs_c$area_var2, 
                                 data_obs_c[ , cnames == time_var], 
                                 data_obs_c[ , cnames == as.character(formula)[2]], 
                                 data_obs_c[ , cnames == SE_var])
  data_w_out_noobs_c <- data.frame(data_noobs_c$area_var2, 
                                   data_noobs_c[ , cnames == time_var], 
                                   data_noobs_c[ , cnames == as.character(formula)[2]], 
                                   data_noobs_c[ , cnames == SE_var])
  colnames(data_w_out_obs_c) <- colnames(data_w_out_noobs_c) <- 
    c(area_var, time_var, as.character(formula)[2], SE_var)
  data_w_out <- rbind(data_w_out_obs_c,data_w_out_noobs_c)
  pred_all  <- c(pred_mat_new$mu_T,pred_mat_nob$mu_T)
  pred_all_fixed   <- c(pred_mat_new$mu_F,pred_mat_nob$mu_F)
  pred_all_fixpen  <- c(pred_mat_new$mu_RF,pred_mat_nob$mu_T)
  SIGMA_T=blockMatrixDiagonal(pred_mat_new$sigma_T,pred_mat_nob$sigma_T)
  SIGMA_Y=blockMatrixDiagonal(pred_mat_new$sigma_Y,pred_mat_nob$sigma_Y)
  sigma_all <- diag(SIGMA_T)
  sigma_Y_all <- diag(SIGMA_Y)
  resid_final <- data_w_out$Y - pred_all
  cov_w_preds <- data.frame(data_w_out, pred = pred_all, sigma_all = sigma_all, 
                            sigma_Y_all = sigma_Y_all, pred_all_fixed = pred_all_fixed,
                            pred_all_fixpen = pred_all_fixpen, resid = resid_final)
  
  # Exporting all of the prediction results.
  pred_results <- list(newdata = cov_w_preds, b = full_b, D = beta_results$D, COV_beta_b = 
                         beta_results$S, SIGMA_T = SIGMA_T, SIGMA_Y = SIGMA_Y, df = 
                         beta_results$df, resid = resid_final)
  return(pred_results)
}


pred <- function(X,Z,SE_vec,fitted_model,full_b,S){  
  
  gamma <- exp(coef(fitted_model$modelStruct$varStruct))
  sigma2 <- fitted_model$sigma^2
  
  C <- as.matrix(cbind(X,Z))
  R <- sigma2*(1 + gamma*SE_vec)*diag(dim(Z)[1])
  
  mu_F <- as.matrix(X)%*%full_b[1:(dim(X)[2])]
  num <- (dim(X)[2]) + dim(as.matrix(fitted_model$modelStruct$reStruct$intercept))[2]
  mu_RF <- C[,1:num]%*%full_b[1:num]
  mu_T <- C%*%full_b
  sigma_T <- C%*%S%*%t(C)
  sigma_Y <- R+sigma_T 
  
  # Export the mean, covariance matrix, and random effects 
  preds <- list(mu_T=mu_T,mu_F=mu_F,mu_RF=mu_RF,sigma_T=sigma_T,sigma_Y=sigma_Y)
  return(preds)
}

pred.noRand <- function(X,Z,SE_vec,fitted_model,N,full_b,S,num_pen){  
  
  gamma <- exp(coef(fitted_model$modelStruct$varStruct))
  sigma2 <- fitted_model$sigma^2
  num_fix <- (dim(X)[2]) 
  
  Z_pen <- Z[,(1:num_pen)]
  Z_ran <- Z[,-(1:num_pen)]
  C <- as.matrix(cbind(X,Z_pen))
  R <- sigma2*(1 + gamma*SE_vec^2)*diag(dim(Z)[1])
  
  b_var <- as.matrix(fitted_model$modelStruct$reStruct[[2]])*sigma2
  l_var <- as.matrix(fitted_model$modelStruct$reStruct[[1]])*sigma2
  bl_var <- blockMatrixDiagonal(b_var,l_var)
  # Construct the block diagonal covariance matrix of the random effects.  
  Dpr1 <- do.call(blockMatrixDiagonal,replicate(N, bl_var, simplify=FALSE))
  l_bl <- ncol(bl_var)
  # Now reorder to match Z
  colnames(Dpr1) <- row.names(Dpr1) <- paste( rep( paste( "C", 1:l_bl), N), 
                                              rep( 1:N, each = l_bl))
  ord_colnames   <- paste( "C", rep( 1:l_bl, each = N), rep( 1:N, l_bl))
  D <- Dpr1[ , ord_colnames]
  D <- D[ ord_colnames, ]
  
  
  mu_F <- as.matrix( X)%*%full_b[ 1:(num_fix)]
  mu_T <- C%*%full_b[ 1:(num_fix + num_pen)]
  fixed_S <- S[ 1:(num_fix + num_pen), 1:(num_fix + num_pen)]
  sigma_T <- Z_ran%*%D%*%t(Z_ran)+C%*%fixed_S%*%t(C)
  sigma_Y <- R+sigma_T 
  
  # Export the mean, covariance matrix, and random effects 
  preds <- list(mu_T=mu_T,mu_F=mu_F,sigma_T=sigma_T,sigma_Y=sigma_Y)
  return(preds)
}

beta.fun <- function(Y,X,Z,SE_vec,fitted_model){
  
  gamma <- exp(coef(fitted_model$modelStruct$varStruct))
  sigma2 <- fitted_model$sigma^2
  
  C <- as.matrix(cbind(X,Z))
  R <- sigma2*(1 + gamma*SE_vec)*diag(dim(Z)[1])
  df<- NULL
  N <- length(levels(fitted_model$groups$area_var2))
  u_var <- as.matrix(fitted_model$modelStruct$reStruct[[3]])*sigma2
  b_var <- as.matrix(fitted_model$modelStruct$reStruct[[2]])*sigma2
  l_var <- as.matrix(fitted_model$modelStruct$reStruct[[1]])*sigma2
  bl_var <- blockMatrixDiagonal(b_var,l_var)
  # Construct the block diagonal covariance matrix of the random effects.  
  Dpr1 <- do.call(blockMatrixDiagonal,replicate(N, bl_var, simplify=FALSE))
  l_bl <- ncol(bl_var)
  # Now reorder to match Z
  colnames(Dpr1) <- row.names(Dpr1) <- paste( rep( paste( "C", 1:l_bl), N), 
                                              rep( 1:N, each = l_bl))
  ord_colnames   <- paste( "C", rep( 1:l_bl, each = N), rep( 1:N, l_bl))
  Dpr2 <- Dpr1[ , ord_colnames]
  Dpr2 <- Dpr2[ ord_colnames, ]
  
  D <- blockMatrixDiagonal(u_var,Dpr2)
  
  t1 <- try(R_inv <- solve(R),silent = TRUE)
  if(!is.null(attr(t1,"class"))){
    cat("Warning: Generalized Inverse of R required.\n")
    R_inv <- ginv(R)}
  t1 <- try(D_inv <- solve(D),silent = TRUE)
  if(!is.null(attr(t1,"class"))){
    if(any(diag(D)<(10e-13))){cat("Warning: Some variance components are zero (with machine precision), setting them to a \n non-zero small value to increase computational stability. See lme output to see which\n variance components are essentially zero and consider revising model. \n")}
    diag(u_var)[diag(u_var)<(10e-13)] <- 10e-13
    diag(Dpr2)[diag(Dpr2)<(10e-13)] <- 10e-13
    D <- blockMatrixDiagonal(u_var,Dpr2)
    t1 <- try(D_inv <- solve(D),silent = TRUE)
    if(!is.null(attr(t1,"class"))){
      cat("Warning: Generalized Inverse of D required.\n")
      D_inv <- ginv(D)}
  }
  
  
  if(!is.null(attr(t1,"class"))){
    diag(u_var)[diag(u_var)<(10e-13)] <- 10e-13
    D <- blockMatrixDiagonal(u_var,Dpr1,Dpr2)
    
    t1 <- try(D_inv <- solve(D),silent = TRUE)
    if(!is.null(attr(t1,"class"))){
      cat("Warning: Generalized Inverse of D required.\n")
      D_inv <- ginv(D)}
  }
  
  G <- t(C)%*%R_inv%*%C
  B <- blockMatrixDiagonal(0*diag(dim(X)[2]),D_inv)
  
  t1 <- try(S <- solve(G + B),silent = TRUE)
  if(!is.null(attr(t1,"class"))){
    cat("Warning: Approximate Inverse of S required.\n")
    FHU <- as.matrix(nearPD(G + B)$mat)
    S <- solve(FHU)
  }
  df <- sum(diag(S%*%G))
  
  
  
  Sigma <- Z%*%D%*%t(Z) + R
  t1 <- try(Sigma_inv <- solve(Sigma),silent = TRUE)
  if(!is.null(attr(t1,"class"))){
    cat("Warning: Generalized Inverse of Sigma required.\n")
    t1 <- try(Sigma_inv <- ginv(Sigma),silent = TRUE)}
  if(!is.null(attr(t1,"class"))){stop("Inverse of covariance matrix failed.")}
  
  t1 <- try(beta_hat <- solve(t(X)%*%Sigma_inv%*%X)%*%t(X)%*%Sigma_inv%*%Y,silent = TRUE)
  if(!is.null(attr(t1,"class"))){
    cat("Warning: Generalized Inverse of X'Sigma X required.\n")
    t1 <- try(beta_hat <- ginv(t(X)%*%Sigma_inv%*%X)%*%t(X)%*%Sigma_inv%*%Y,silent = TRUE)}
  full_b <- c(beta_hat,D%*%t(Z)%*%Sigma_inv%*%(Y - X%*%beta_hat))
  
  
  # Export the coefficients and covariance matrices 
  preds <- list(b=full_b,D=D,S=S,df=df)
  return(preds)
}

##### Making the Z matrix.
Z.matrix <- function(Zi.formula, re_formula, rand_formula, data, NA_ind = NULL, 
                     Pen_type){
  if( Pen_type == "Multi"){
    Zi.formula <- ~0 + Zi.matrix:Pcov_data
  }
  Pen_mat <- model.matrix( Zi.formula, data = data)
  if( !is.null( NA_ind)){ Pen_mat <- Pen_mat[NA_ind,] }
  
  re_form_area <- paste("~0+",area_var)
  if( length(attr(terms(re_formula),"factors")) > 0 ){
    
    for(i in 1:(length(colnames(attr(terms(re_formula),"factors"))))){
      re_form_area <- paste(re_form_area,"+",colnames(attr(terms(re_formula),"factors"))[i],":",area_var)
    }
  }
  re_form_area <- as.formula(re_form_area)
  
  re_mat <- model.matrix( re_form_area, data = data)
  
  rand_form_area <- "~0"
  for(i in 1:(length(colnames(attr(terms(rand_formula),"factors"))))){
    rand_form_area <- paste(rand_form_area,"+",colnames(attr(terms(rand_formula),"factors"))[i],":",area_var)
  }
  rand_form_area <- as.formula(rand_form_area)
  
  rand_mat <- model.matrix(rand_form_area, data = data)
  
  full_Z_matrix <- cbind(Pen_mat,re_mat,rand_mat)
  return(full_Z_matrix) 
}








varSum <-
  function(value = numeric(0), form = ~ fitted(.), fixed = NULL)
  {
    value <- unlist(value)		# may be given as a list
    fixed <- attr(value, "fixed") <- unlist(fixed)
    attr(value, "formula") <- form <- asOneSidedFormula(form)
    if (length(all.vars(getCovariateFormula(form))) == 0) {
      stop("\"form\" must have a covariate")
    }
    if (!is.null(getGroupsFormula(form))) {
      if (is.null(grpNames <- names(value)) && (length(value) > 1)) {
        stop("Initial values must have groups names in varPower")
      }
      if (!is.null(fixed)) {
        if (is.null(names(fixed))) {
          stop("Fixed parameters must have groups names in varPower")
        }
      }
      attr(value, "groupNames") <- c(grpNames, names(fixed))
    } else {
      attr(value, "whichFix") <- !is.null(fixed)
    }
    class(value) <- c("varSum", "varFunc")
    value
  }

###*# Methods for standard generics

coef.varSum <-
  function(object, unconstrained = TRUE, allCoef = FALSE, ...)
  {
    if (((length(object) == 0) &&
         (!allCoef || is.null(attr(object, "fixed")))) ||
        is.null( wPar <- attr(object, "whichFix"))) {
      return(numeric(0))
    }
    val <- double(length(wPar))
    if (any(wPar)) {
      val[wPar] <- attr(object, "fixed")
    }
    if (any(!wPar)) {
      val[!wPar] <- as.vector(object)
    }
    if (!is.null(getGroupsFormula(object))) {
      ##different values per group
      names(val) <- attr(object, "groupNames")
    } else {
      names(val) <- "covarCoef"
    }
    if (!allCoef) {
      val <- val[!wPar]
    }
    val
  }

"coef<-.varSum" <-
  function(object, ..., value)
  {
    if (length(object) > 0) {		# varying parameters
      value <- as.numeric(value)
      if (length(value) != length(object)) {
        stop(paste("Cannot change the length of the varStruct",
                   "parameter after initialization"))
      }
      object[] <- value
      aux <- coef(object, FALSE, allCoef = TRUE)
      if (!is.null(grps <- getGroups(object))) {
        aux <- aux[grps]
      }
      attr(object, "logLik") <-
        sum(log(attr(object, "weights") <- 1/sqrt(1+exp(aux)* getCovariate(object))))
    } else {
      stop(paste("Cannot change coefficients before initialization or",
                 "when all parameters are fixed"))
    }
    object
  }

Initialize.varSum <-
  function(object, data, ...)
  {
    form <- formula(object)
    if (all(!is.na(match(all.vars(getCovariateFormula(form)), names(data))))) {
      ## can evaluate covariate on data
      attr(object, "needUpdate") <- FALSE
      attr(object, "covariate") <- getCovariate(data, form)
    } else {
      attr(object, "needUpdate") <- TRUE
    }
    if (!is.null(grpForm <- getGroupsFormula(form))) {
      strat <- as.character(getGroups(data, form,
                                      level = length(splitFormula(grpForm, sep = "*")),
                                      sep = "*"))
      uStrat <- unique(strat)
      if (length(uStrat) > 1) {		# multi-groups
        attr(object, "groups") <- strat
        if (!is.null(attr(object, "fixed"))) {
          fixNames <- names(attr(object, "fixed"))
          if (is.null(fixNames)) {
            stop("Fixed parameters must have group names")
          }
          if (any(is.na(match(fixNames, uStrat)))) {
            stop("Mismatch between group names and fixed values names")
          }
        } else {
          fixNames <- NULL
        }
        uStratVar <- uStrat[is.na(match(uStrat, fixNames))]
        nStratVar <- length(uStratVar)
        attr(object, "whichFix") <- !is.na(match(uStrat, fixNames))
        if (nStratVar > 0) {
          if (length(object) <= 1) {
            ## repeat for all groups
            names(object) <- NULL
            oldAttr <- attributes(object)
            if (length(object) > 0) {
              object <- rep(as.vector(object), nStratVar)
            } else {
              object <- rep(0, nStratVar)
            }
            attributes(object) <- oldAttr
            attr(object, "groupNames") <- uStrat
            names(object) <- uStratVar
          } else {
            if (length(as.vector(object)) != nStratVar) {
              stop(paste("Initial value for \"varSum\" should be of length",
                         nStratVar))
            }
            stN <- attr(object, "groupNames") #must have names
            if ((length(stN) != length(uStrat)) ||
                any(sort(stN) != sort(uStrat))) {
              stop("Nonexistent groups names for initial values in varSum")
            }
          }
        } else {
          if (all(attr(object, "fixed") == 0)) {
            ## equal variances structure
            return(Initialize(varIdent(), data))
          } else {
            oldAttr <- attributes(object)
            object <- numeric(0)
            attributes(object) <- oldAttr
            attr(object, "groupNames") <- uStrat
          }
        }
      } else {                            # single stratum
        attr(object, "formula") <- getCovariateFormula(formula(object))
        attr(object, "whichFix") <- !is.null(attr(object, "fixed"))
      }
    }
    if (is.null(getGroupsFormula(object))) {
      ## single stratum
      if (attr(object, "whichFix")) {
        if (!attr(object, "fixed")) {
          ## equal variances structure
          return(Initialize(varIdent(), data))
        } else {
          oldAttr <- attributes(object)
          object <- numeric(0)
          attributes(object) <- oldAttr
        }
      } else {
        len <- length(as.vector(object))
        if (len == 0) {			# uninitialized
          oldAttr <- attributes(object)
          object <- 0
          attributes(object) <- oldAttr
        } else if (len > 1) {
          stop("Initial value for \"varSum\" should be of length 1.")
        }
      }
    }
    if (!is.null(covar <- getCovariate(object))) {
      natPar <- coef(object, allCoef = TRUE)
      if (!is.null(grps <- getGroups(object))) {
        natPar <- natPar[grps]
      }
      attr(object, "logLik") <-
        sum(log(attr(object, "weights") <-  1/sqrt(1+exp(natPar)* getCovariate(object))))
      object
    } else {
      NextMethod()
    }
  }


summary.varSum <-
  function(object, structName = "Sum of one and constant times variance covariate", ...)
  {
    if (!is.null(getGroupsFormula(object))) {
      structName <- paste(structName, " different strata", sep = ",")
    }
    summary.varFunc(object, structName)
  }

update.varSum <-
  function(object, data, ...)
  {
    val <- NextMethod()
    if (length(val) == 0) {		# chance to update weights
      aux <- coef(val, allCoef = TRUE)
      if (!is.null(grps <- getGroups(val))) {
        aux <- aux[grps]
      }
      attr(val, "logLik") <-
        sum(log(attr(val, "weights") <- 1/sqrt(1+exp(aux)* getCovariate(object))))
    }
    val
  }

summary.varFunc <-
  function(object, structName = class(object)[1], ...)
  {
    attr(object, "structName") <- structName
    attr(object, "oClass") <- class(object)
    class(object) <- "summary.varFunc"
    object
  }

blockMatrixDiagonal<-function(...){  
  ## This is a function for making block diagonal matrices.
  matrixList<-list(...)
  if(is.list(matrixList[[1]])) matrixList<-matrixList[[1]]
  
  dimensions<-sapply(matrixList,FUN=function(x) dim(x)[1])
  finalDimension<-sum(dimensions)
  finalMatrix<-matrix(0,nrow=finalDimension,ncol=finalDimension)
  index<-1
  for(k in 1:length(dimensions)){
    finalMatrix[index:(index+dimensions[k]-1),index:(index+dimensions[k]-1)]<-matrixList[[k]]
    index<-index+dimensions[k]
  }
  finalMatrix
}





