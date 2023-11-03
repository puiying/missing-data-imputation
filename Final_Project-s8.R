
# load library
library(MASS)
library(mice)
library(dplyr)
library(tableone)
library(labelled)
library(knitr)
library(rms)
library(forestplot)
library(optmatch)
library(Matching)
library(reshape2)
library(survey)
library(boot)
library(ggplot2)
library(MatchIt)
library(locfit)
library(hash)
library(data.table)
library(tidyverse)


####### Data Generation Mechanism
# set seed
set.seed(0010)
RR <- 1 # or 2

### create a function to generate datasets based on parameters
simulate_dataset <- function(p_missing, n, p, RR){
  # define theta c to filter relative risk ratio
  if (RR == 1){
    theta_c <- 0 # or 1.221 or 1.289
  }
  
  ### generate 3 covariates x1, x2, x3
  # define mean and variance-covariance
  cov_mu <- c(0, 0, 0)                                   
  cov_sigma <- matrix(c(1, p, p, p, 1, p, p, p, 1), ncol = 3)
  
  # create multivariate normal distribution
  X <- mvrnorm(n = n, mu = cov_mu, Sigma = cov_sigma)
  # dichotomize to a threshold of 0 for binary variable, x3
  X[, 3] <- (X[, 3] >= 0) * 1.0
  # we need approx ~ 0.5 prevalence
  print(paste("The prevalence of binary covariate, X3 is", mean(X[, 3]))) 
  prev_X3 <- mean(X[, 3])
  
  ### treatment assignment
  Z <- rbinom(n, size=1, prob=expit(-1.15 + 0.7*X[,1] + 0.6*X[,2] + 0.6*X[,3]) )
  # we need approx ~ 0.3 treatment prevalence
  print(paste("The expected value of treatment assignment is", mean(Z)))
  treatment_mean <- mean(Z)
  
  ### binary outcome
  Y <- rbinom(n, size=1, prob=expit(-1.5 + 0.5*X[,1] + 0.5*X[,2] + 0.3*X[,3] + theta_c*Z) )
  # we need approx ~ 1 relative risk ratio
  print(paste("The nominal relative risk ratio is", RR, ', and the unadjusted RR estimated from dataset is', mean(Y[Z==1])/mean(Y[Z==0])))
  RR_unadj <- mean(Y[Z==1])/mean(Y[Z==0])
  
  
  ### missingness mechanism
  if (p_missing == 0.1){
    pm = -1.6
  }
  else if (p_missing == 0.3){
    pm = 0
  }
  else if (p_missing == 0.6){
    pm = 1.6
  }
  M_ind <- matrix(0, n, 3)
  M_ind[, 1] <- rbinom(n, size=1, prob=expit(-1.5 + pm + Z + X[,2]) )
  M_ind[, 3] <- rbinom(n, size=1, prob=expit(-1.5 + pm + Z + X[,2]) )
  # we need approx ~ 0.3 missing data
  print("The missing indicator on average are")
  M_ind_ave <- colMeans(M_ind)
  print(colMeans(M_ind))
  
  ### Full data 
  df <- data.frame(X1_orig = X[,1],
                   X2 = X[,2], 
                   X3_orig = X[,3], 
                   X1 = ifelse(M_ind[,1]==1, NA, X[,1]),
                   X3 = ifelse(M_ind[,3]==1, NA, X[,3]),
                   M1 = M_ind[,1],
                   M3 = M_ind[,3],
                   Z = Z,
                   Y = Y,
                   sample_size = n,
                   correlation = p,
                   missingness = p_missing,
                   prev_X3 = prev_X3,
                   treatment_mean = treatment_mean,
                   RR_nominal = RR,
                   RR_unadj = RR_unadj,
                   M_ind_ave_X1 = M_ind_ave[1],
                   M_ind_ave_X3 = M_ind_ave[3]
  )
  
  return(df)
}


####### Missing Data Imputations
# create a function to do data imputation/approaches
impute_data <- function(imp_type, df, M_imp=M){
  
  # define missing sequence from least to most missingness
  if (sum(is.na(data[["S1"]]$X1)) > sum(is.na(data[["S1"]]$X3))){
    missing_seq <- c('Y', 'Z', 'X2', 'X3', 'X1')
  }
  else{
    missing_seq <- c('Y', 'Z', 'X2', 'X1', 'X3')
  }
  
  ### Full data 
  if (imp_type == "full"){
    # define final df
    df.out <- df
    # define propensity score model
    f.ps <- quote(Z ~ X1_orig + X2 + X3_orig)
  }
  ### Complete Case Analysis
  else if (imp_type == "cc"){
    # define final df
    df.out <- na.omit(df)
    # define propensity score model
    f.ps <- quote(Z ~ X1 + X2 + X3)
  }
  ### Missing Pattern (with missing indicator)
  else if (imp_type == "mp"){
    # define final df
    df.out <- df
    df.out$X1[is.na(df.out$X1)] <- 0
    df.out$X3[is.na(df.out$X3)] <- 0
    # define propensity score model
    f.ps <- quote(Z ~ X1 + X2 + X3 + M1 + M3)
  }
  
  ### MIte
  else if (imp_type == "mite"){
    # rearrange columns based on amount of missingness
    df <- df[missing_seq]
    
    # impute using MICE algorithm
    df.out <- mice(df, 
                   method=c('', '', '', 'pmm', 'pmm'), 
                   maxit=25,
                   m=M_imp, printFlag =FALSE)
    # define propensity score model
    f.ps <- quote(Z ~ X1 + X2 + X3)
  }
  ### MIps
  else if (imp_type == "mips"){
    # rearrange columns based on amount of missingness
    df <- df[missing_seq]
    
    # impute using MICE algorithm
    imputed <- mice(df, 
                    method=c('', '', '', 'pmm', 'pmm'), 
                    maxit=25,
                    m=M_imp, printFlag =FALSE)
    
    # get ps model
    ps.m <- with(imputed, glm(Z ~ X1 + X2 + X3, family="binomial"))
    
    # get an id value
    id <- tibble(id = seq(from = 1, to = nrow(df), length.out = nrow(df)))
    
    # predict ps
    ps.pred <-
      tibble(.imputed = 1:M) %>% 
      mutate(p = map(.imputed, ~ predict(ps.m$analyses[[.]], type="response", 
                                         se.fit = TRUE) %>%
                       data.frame())
      )
    
    # unnest the results and then augment a little
    ps.pred <- ps.pred %>% 
      unnest(p) %>% # add in the nd predictor data
      bind_cols(
        bind_rows(replicate(M, id, simplify = FALSE))
      ) 
    
    # average out ps
    df["ps"] <- ps.pred %>% 
      group_by(id) %>% 
      summarise(fit_bar = mean(fit))
    
    # final df
    df.out <- df
    
    # define propensity score model
    f.ps <- quote(Z ~ ps)
  }
  # return a list with dataframe and formula
  return (list(df.out, f.ps))
}


####### Propensity Score Matching
# create a function that will perform ps matching
ps_matching <- function(data, match_type, f_exp){
  # one-to-one Optimal Matching
  if (match_type == "optimal"){
    mod_match <- with(data, 
                      matchit(as.formula(f_exp), 
                              method="optimal", distance="logit", ratio=1))
  }
  # one-to-one NNM matching: no replacement, no caliper
  else if (match_type == "nnm-no-no"){
    mod_match <- with(data,
                      matchit(as.formula(f_exp), method="nearest", replace=FALSE, 
                              m.order="random", distance="logit", ratio=1))
  }
  # one-to-one NNM matching: with replacement, no caliper
  else if (match_type == "nnm-yes-no"){
    mod_match <- with(data, 
                      matchit(as.formula(f_exp), method="nearest", replace=TRUE, 
                              m.order="random", distance="logit", ratio=1))
  }
  # one-to-one NNM matching: no replacement, with caliper
  else if (match_type == "nnm-no-yes"){
    mod_match <- with(data, 
                      matchit(as.formula(f_exp), method="nearest", replace=FALSE,
                              caliper=1, std.caliper=TRUE, 
                              m.order="random", distance="logit", ratio=1))
  }
  # one-to-one NNM matching: with replacement, with caliper
  else if (match_type == "nnm-yes-yes"){
    mod_match <- with(data,
                      matchit(as.formula(f_exp), method="nearest", replace=TRUE,
                              caliper=1, std.caliper=TRUE,
                              m.order="random", distance="logit", ratio=1))
  }
  return (mod_match)
}


####### Performance Analysis
# define formula that calculates att 
cal_att <- function(scenario, imputation, matching, datasets, match_matrix, repetition){
  
  # Get mean per group
  p1 <- mean(datasets$Y[datasets$Z == 1], na.rm = TRUE)
  p0 <- mean(datasets$Y[as.numeric(match_matrix)], na.rm = TRUE)
  # Get sd per group
  s1 <- sd(datasets$Y[datasets$Z == 1], na.rm = TRUE)
  s0 <- sd(datasets$Y[as.numeric(match_matrix)], na.rm = TRUE)
  # Get number of observations per group
  n1 <- nrow(datasets[datasets$Z == 1, ])
  n0 <- nrow(datasets[as.numeric(match_matrix), ])
  
  # Get unadjusted average treatment effects
  ATT <- p1 - p0
  
  # Get difference pooled standard error for discrete
  SE <- sqrt(p1*(1-p1)/n1 + p0*(1-p0)/n0)
  
  # Get upper and lower 95% confidence interval
  CI <- ATT + c(-1, 1)*qnorm(0.975)*SE
  
  # Get all information in a df and output
  df <- data.frame(scenario = scenario, 
                   imputation = imputation, 
                   matching = matching, 
                   repetition = repetition,
                   ATT = ATT, 
                   SE = SE, 
                   CI.Lower = CI[1], 
                   CI.Upper = CI[2])
  df$ATT <- as.numeric(df$ATT)
  df$SE <- as.numeric(df$SE)
  df$CI.Lower <- as.numeric(df$CI.Lower)
  df$CI.Upper <- as.numeric(df$CI.Upper)
  return (df)
}


########## Putting them together
# define imputation
M <- 50

imputation.methods <- c("cc", "full", "mp", "mips", "mite") # 
matching.methods <- c("optimal", "nnm-no-no", "nnm-yes-no", "nnm-no-yes", "nnm-yes-yes")
var <- c("X1", "X2", "X3")

# create empty placeholders
performance <- data.frame()
tabmatches <- new.env(hash = TRUE, parent = emptyenv(), size = 100L)

# define sample size n
n_set <- c(500, 2000)
# define correlation between covariates, p
p_set <- c(0.3, 0.6)
# define missingness rate, p_missing
p_missing_set <- c(0.1, 0.6)

# # simulate different datasets based on different scenarios
# data_simulation_table <- tibble(`Senario` = character(),  `Cor Coef`=numeric(), 
#                                 `Missing Rate (Nominal)`=numeric(), 
#                                 `Missing Rate X1 (Observed)` = numeric(),
#                                 `Missing Rate X3 (Observed)` = numeric(),
#                                 `Sample Size`=numeric(), `Relative Risk (Nominal)` = numeric(),
#                                 `Prevalance of X3` = numeric(), `Treatment Assignment Rate` = numeric(),
#                                 `Relative Risk (Obeserved, unadjusted)` = numeric())

counter <- 8
for (p in c(0.6)){
  for (n in c(2000)){
    for (p_missing in c(0.6)){
      datasets <- data.frame()
      data_match <- data.frame()
      
      for (r in 1:50){  
        # set a name for simulation set
        print(paste("Scenario", counter, "(S", counter, "- Iteration", r, "):"))
        print(paste("n =", n, ", missing =", p_missing, ", correlation =", p))
        sim_name <- paste0("S", counter)
        
        # simulate dataset
        data <- simulate_dataset(p_missing, n, p, RR)
        datasets <- data.frame(rbindlist(list(datasets, data)))
        
        #tabmatches[[sim_name]] <- list()
        for (i in imputation.methods){
          
          # impute dataset
          out <- impute_data(i, data)
          
          #tabmatches[[sim_name]][[i]] <- list()
          for (mm in matching.methods){
            
            # ps matching 1:1
            mod_match <- ps_matching(data=out[[1]], match_type=mm, f_exp=out[[2]])
            
            # calculate target performance metrics
            # for mite only
            if (i == "mite"){
              imp_list_df <- map(1:M, function(x) complete(out[[1]], x))
              temp <- data.frame()
              for (m in 1:M){
                # prepare for att
                temp <- rbind(temp, cal_att(scenario=sim_name,
                                            imputation=i,
                                            matching=mm,
                                            datasets=imp_list_df[[m]],
                                            match_matrix=mod_match$analyses[[m]]$match.matrix,
                                            repetition=r))
                
                # prepare for smd plot
                dmatch <- data.frame(rbind(imp_list_df[[m]][imp_list_df[[m]]$Z == 1, ],
                                           imp_list_df[[m]][as.numeric(mod_match$analyses[[m]]$match.matrix), ]))
                
                dmatch$scenario <- sim_name
                dmatch$imputation <- i
                dmatch$matching <- mm
                dmatch$repetition <- r
                data_match <- data.frame(rbindlist(list(data_match, dmatch[c("Y", "Z", "X1", "X2", "X3", "scenario", "imputation", "matching", "repetition")])))
                
              }
              # get average ATT
              att <- as.numeric(colMeans(temp[5]))
              # get within variance from Rubin's Rule
              V <- as.numeric(colMeans(temp[6]**2))
              # get between variance from Rubin's Rule
              B <- as.numeric(var(temp[5]))
              # Rubin's rule sd
              se <- V + (1 + 1/M) * B
              # calculate CI
              CI <- att + c(-1, 1)*qnorm(0.975)*se
              # average of att
              performance <- rbind(performance, c(sim_name, i, mm, r, att, se, CI[1], CI[2]))
              # # average of tableone
              # tabmatches[[sim_name]][[i]][[mm]] <- CreateTableOne(vars = var, strata = "Z", data = data_match, test = FALSE)
              
            }
            # for others
            else {
              # prepare for smd plot
              dmatch <- data.frame(rbind(out[[1]][out[[1]]$Z == 1, ],
                                         out[[1]][as.numeric(mod_match$match.matrix), ]))
              dmatch$scenario <- sim_name
              dmatch$imputation <- i
              dmatch$matching <- mm
              dmatch$repetition <- r
              data_match <- data.frame(rbindlist(list(data_match, dmatch[c("Y", "Z", "X1", "X2", "X3", "scenario", "imputation", "matching", "repetition")])))
              # tabmatches[[sim_name]][[i]][[mm]] <- CreateTableOne(vars = var, strata = "Z", data = data_match, test = FALSE)
              
              # calculate att
              performance <- rbind(performance, cal_att(scenario=sim_name,
                                                        imputation=i,
                                                        matching=mm,
                                                        repetition=r,
                                                        datasets=out[[1]],
                                                        match_matrix=mod_match$match.matrix))
            }
          }
        }
      }
      tabmatches[[sim_name]] <- list()
      # get unweighted SMD results
      tabUnmatched <- CreateTableOne(vars = var, strata = "Z", datasets, test = FALSE)
      
      for (i in imputation.methods){
        tabmatches[[sim_name]][[i]] <- list()
        for (mm in matching.methods){
          tabmatches[[sim_name]][[i]][[mm]] <- CreateTableOne(vars = var, strata = "Z", data = data_match[(data_match$scenario==sim_name) & (data_match$imputation==i) & (data_match$matching==mm), ], test = FALSE)
        }
      }
      # counter <- counter + 1
    }
  }
}

performance['bias'] <- abs(0 - as.numeric(performance$ATT))
performance$ATT <- as.numeric(performance$ATT)
performance$SE <- as.numeric(performance$SE)
performance$CI.Lower <- as.numeric(performance$CI.Lower)
performance$CI.Upper <- as.numeric(performance$CI.Upper)
write.csv(performance,"performance.csv", row.names = FALSE)
save.image(file = "final_project.RData")