############################################
## 0. Package loading and parallel setup
############################################

library(FDRreg)
library(locfdr)
library(foreach)
library(doMC)
library(parallel)
library(mvtnorm)
library(faraway)
library(glmnet)
library(biglasso)
library(writexl)

set.seed(6)

tstart = proc.time()
registerDoMC(20)
options(error = NULL)


############################################
## 1. Global simulation configuration
############################################

LD_simulation = FALSE 
extract_some_snps = TRUE
no_of_snps_from_map = 100000

if (LD_simulation==TRUE && extract_some_snps)   {nobserve = no_of_snps_from_map}
if (LD_simulation==TRUE && !extract_some_snps)  {nobserve = nrow(map)} 
if (LD_simulation==F)                           {nobserve = 1000}


nobserve = 1000  # Can change the sample size here

nsims = 100 
fdr_level = 0.1
total_no_covar = 200
nfolds_enet = 10


############################################
## 2. Signal probability functions
############################################

flist_no_covar = c(2, 3, 4, 5, 5, 5, 5, 2, 3, 4, 5, 5, 5, 5, 3, 3, 3, 5, 5, 5)

flist = list(
  # Group 1: Exploring low number of covariates (medium coefficients, medium intercept)
  function(x1, x2) -3.0 + 0.8*x1 + 1.0*x2,
  function(x1, x2, x3) -3.0 + 0.8*x1 + 1.0*x2 + 1.2*x3,
  function(x1, x2, x3, x4) -3.0 + 0.8*x1 + 1.0*x2 + 1.2*x3 + 1.4*x4,
  function(x1, x2, x3, x4, x5) -3.0 + 0.8*x1 + 1.0*x2 + 1.2*x3 + 1.4*x4 + 1.6*x5,
  
  # Group 2: Exploring low coefficients (medium number of covariates, medium intercept)
  function(x1, x2, x3, x4, x5) -3.0 + 0.4*x1 + 0.6*x2 + 0.8*x3 + 1.0*x4 + 1.2*x5,
  function(x1, x2, x3, x4, x5) -3.0 + 0.2*x1 + 0.4*x2 + 0.6*x3 + 0.8*x4 + 1.0*x5,
  function(x1, x2, x3, x4, x5) -3.0 + 0.08*x1 + 0.16*x2 + 0.24*x3 + 0.32*x4 + 0.4*x5,
  
  # Group 3: Combination of low number of covariates and low coefficients
  function(x1, x2) -3.0 + 0.2*x1 + 0.4*x2,
  function(x1, x2, x3) -3.0 + 0.2*x1 + 0.4*x2 + 0.6*x3,
  function(x1, x2, x3, x4) -3.0 + 0.2*x1 + 0.4*x2 + 0.6*x3 + 0.8*x4,
  function(x1, x2, x3, x4, x5) -3.0 + 0.2*x1 + 0.4*x2 + 0.6*x3 + 0.8*x4 + 1.0*x5,
  
  # Group 4: Exploring intercept effect with low coefficients
  function(x1, x2, x3, x4, x5) -3.0 + 0.08*x1 + 0.16*x2 + 0.24*x3 + 0.32*x4 + 0.4*x5,
  function(x1, x2, x3, x4, x5) -2.0 + 0.08*x1 + 0.16*x2 + 0.24*x3 + 0.32*x4 + 0.4*x5,
  function(x1, x2, x3, x4, x5) -1.0 + 0.08*x1 + 0.16*x2 + 0.24*x3 + 0.32*x4 + 0.4*x5,
  
  # Group 5: Exploring intercept effect with low number of covariates
  function(x1, x2, x3) -3.0 + 0.8*x1 + 1.0*x2 + 1.2*x3,
  function(x1, x2, x3) -2.0 + 0.8*x1 + 1.0*x2 + 1.2*x3,
  function(x1, x2, x3) -1.0 + 0.8*x1 + 1.0*x2 + 1.2*x3,
  
  # Group 6: Exploring intercept effect with combination of low covariates and low coefficients
  function(x1, x2, x3, x4, x5) -3.0 + 0.2*x1 + 0.4*x2 + 0.6*x3 + 0.8*x4 + 1.0*x5,
  function(x1, x2, x3, x4, x5) -2.0 + 0.2*x1 + 0.4*x2 + 0.6*x3 + 0.8*x4 + 1.0*x5,
  function(x1, x2, x3, x4, x5) -1.0 + 0.2*x1 + 0.4*x2 + 0.6*x3 + 0.8*x4 + 1.0*x5
)



############################################
## 3. Effect size mixture model
############################################

# Parameter settings for the mixture model
parlist = list(
  list(weights = c(0.48, 0.04, 0.48), mu = c(-2, 0, 2), sigma2 = c(1, 5, 1)),
  list(weights = c(0.4, 0.2, 0.4), mu = c(-1.25, 0, 1.25), sigma2 = c(2, 4, 2))
)


############################################
## 4. Utility functions
############################################

# Function to calculate error rates
GetErrorRates = function(is_signal, is_finding) {
  FDP_actual = sum(is_signal[is_finding == 1] == 0) / sum(is_finding == 1)
  power_actual = sum(is_finding[is_signal == 1] == 1) / sum(is_signal == 1)
  return(c(FDP_actual, power_actual))
}


rnormix = function(n, weights, mu, sigma2) {
  k = length(weights)
  component = sample(1:k, size = n, replace = TRUE, prob = weights)
  rnorm(n, mean = mu[component], sd = sqrt(sigma2[component]))
}


# Function for marginal association test (Sure independence screening)
marginal_association_test = function(X, y, threshold = 0.05) {
  p_values = apply(X, 2, function(x) cor.test(x, y)$p.value)
  return(p_values <= threshold)
}



############################################
## 5. Result containers
############################################


results_list = list()
counter = 1


result_power_FDR = list()
result_time = list()
result_NArate = list()
result_correct_selection_rates = list()
result_correct_selection_rates2 = list()
result_correct_selection_rates3 = list()
save_signal_findings = list()
save_selected_var = list()
rate  = matrix(NA,nsims,length(flist))



############################################
## 6. Choose signal density scenario
############################################
jj = 1 # jj = 1 or 2
cat("\tSignal Density", jj, "\n")
mypars = parlist[[jj]]
############################################



############################################
## 7. Simulation 
############################################


for (i in seq_along(flist)) {
  print(c("first loop = ",i))

  # Elastic Net screening 
  ENN = 0 
  # Lasso
  BLL = 0 
  # MA test 
  MAA = 0
  
  no_assoc_covar = flist_no_covar[i]
  cat("Function", i, "\n")
  myxfunction = flist[[i]]
  
  results  = foreach(j = 1:nsims) %dopar% {

    print(c("second loop = ",j))

    # Generate data
    X = matrix(abs(rnorm(nobserve * no_assoc_covar)), nrow = nobserve, ncol = no_assoc_covar)
    log_odds_expr = paste0("log_odds = myxfunction(", paste0("X[,", 1:no_assoc_covar, "]", collapse = ", "), ")")
    eval(parse(text = log_odds_expr))
    success_prob = faraway::ilogit(log_odds)
    

    is_signal = rbinom(nobserve, 1, success_prob) # Whether the snp is significant or not
    signal_true = list(is_signal)

    mu = rnormix(nobserve, mypars$weights, mypars$mu, mypars$sigma2)

     # Sets the effect size to exactly 0 for all observations that are not true signals
    mu[is_signal == 0] = 0  
    y = mu + rnorm(nobserve, 0, 1)
    
    # Generate irrelevant covariates
    num_notRelevant_covar = total_no_covar - no_assoc_covar
    X_other = matrix(abs(rnorm(nobserve * num_notRelevant_covar)), nrow = nobserve, ncol = num_notRelevant_covar)
    X_final = cbind(X, X_other)
    
    ########################################
    ## Method 1: Original FDRreg
    ########################################
    t_1 = proc.time() 
    out3 = FDRreg(y, covars= X_final, nulltype = 'theoretical')
    t_origFDRreg = proc.time() - t_1
  

    FDRreg_result = out3$FDR 
    FDRreg_signal = list(FDRreg_result)
    findings = which(out3$FDR <= fdr_level)
    
    is_finding = rep(0, nobserve)
    is_finding[findings] = 1
    rates_EBayesFDRreg = GetErrorRates(is_signal, is_finding)
    
    #######################################################################################################################################################################
    
    
    ########################################
    ## Method 2: Elastic Net screening
    ########################################

    lambda_seq = exp(seq(log(0.0001), log(1), length.out = 100))
    cvfit = cv.glmnet(x = X_final, y = abs(y), nfolds = nfolds_enet, alpha = 0.5,lambda = lambda_seq)
  
    coef_list = as.vector(coef(cvfit, s = "lambda.min"))[-1]
    
    t_2 = proc.time()
    
    if (sum(coef_list != 0) == 1) {
      mat = matrix(X_final[, coef_list != 0], nobserve , 1)
      out3_enet = FDRreg(y, mat, nulltype = 'theoretical' )
      
      EN_result = out3_enet$FDR 
      EN_signal = list(EN_result)  
      EN_selected_var = list(which(coef_list != 0))

      t_enet = proc.time() - t_2
      
      findings_enet = which(out3_enet$FDR <= fdr_level)
      is_finding_enet = rep(0, nobserve)
      is_finding_enet[findings_enet] = 1
      rates_EBayesFDRreg_enet = GetErrorRates(is_signal, is_finding_enet)
      
      # Calculate true positive rate for variable selection
      enet_true_selected = 1:no_assoc_covar  # Original variables
      enet_selected_vars = which(coef_list != 0)
      enet_correct_rate = sum(enet_selected_vars %in% enet_true_selected) / no_assoc_covar
      enet_correct_rate2 = sum(enet_selected_vars %in% enet_true_selected) / length(enet_selected_vars)
      enet_correct_rate3 = 2 * (enet_correct_rate * enet_correct_rate2) / (enet_correct_rate + enet_correct_rate2)
    }
    
    if (sum(coef_list != 0) > 1) {
      out3_enet = FDRreg(y, X_final[, coef_list != 0], nulltype = 'theoretical' )
      
      EN_result = out3_enet$FDR 
      EN_signal = list(EN_result)
      EN_selected_var = list(which(coef_list != 0))
      ##################################################################

      t_enet = proc.time() - t_2
      
      findings_enet = which(out3_enet$FDR <= fdr_level)
      is_finding_enet = rep(0, nobserve)
      is_finding_enet[findings_enet] = 1
      rates_EBayesFDRreg_enet = GetErrorRates(is_signal, is_finding_enet)
      
      # Calculate true positive rate for variable selection
      enet_true_selected = 1:no_assoc_covar  # Original variables
      enet_selected_vars = which(coef_list != 0)
      enet_correct_rate = sum(enet_selected_vars %in% enet_true_selected) / no_assoc_covar
      enet_correct_rate2 = sum(enet_selected_vars %in% enet_true_selected) / length(enet_selected_vars)
      enet_correct_rate3 = 2 * (enet_correct_rate * enet_correct_rate2) / (enet_correct_rate + enet_correct_rate2)
      ##################################################################
      
      
    } else{

      # If no variables were selected, set results to NA or some default value
      t_enet = NA
      findings_enet = integer(0)
      is_finding_marginal = rep(0, nobserve)
      rates_EBayesFDRreg_enet = c(NA, NA)  # or c(0, 0) depending on how you want to handle this case
      
      
      enet_correct_rate = NA
      enet_correct_rate2 = NA
      enet_correct_rate3 = NA
      
      EN_signal = list()
      EN_selected_var = list()
      
      ENN = ENN + 1
    }
    ###############################################################################
    
    
    
    
    
    ########################################
    ## Method 3: Lasso screening
    ########################################
    
    lambda_seq = exp(seq(log(0.0001), log(10), length.out = 100))
    cvfit2 = cv.biglasso(X = as.big.matrix(X_final, backingfile = ""), 
                          y = abs(y), 
                          lambda = lambda_seq,
                          ncores = 1,  
                          nfolds = nfolds_enet, 
                          screen = "SSR",
                          trace = FALSE)

    coef_list2 = as.vector(coef(cvfit2, s = "lambda.min"))[-1]
    
    t_3 = proc.time()
    
    if (sum(coef_list2 != 0) == 1) {
      mat = matrix(X_final[, coef_list2 != 0], nobserve , 1)
      out3_biglasso = FDRreg(y, mat, nulltype = 'theoretical' )
      
      LASSO_result = out3_biglasso$FDR 
      LASSO_signal = list(LASSO_result)
      LASSO_selected_var = list(which(coef_list2 != 0))

      t_biglasso = proc.time() - t_3
      
      findings_biglasso = which(out3_biglasso$FDR <= fdr_level)
      is_finding_biglasso = rep(0, nobserve)
      is_finding_biglasso[findings_biglasso] = 1
      rates_EBayesFDRreg_biglasso = GetErrorRates(is_signal, is_finding_biglasso)
      
      # Calculate true positive rate for variable selection
      lasso_true_selected = 1:no_assoc_covar  # Original variables
      lasso_selected_vars = which(coef_list2 != 0)
      lasso_correct_rate = sum(lasso_selected_vars %in% lasso_true_selected) / no_assoc_covar
      lasso_correct_rate2 = sum(lasso_selected_vars %in% lasso_true_selected) / length(lasso_selected_vars)
      lasso_correct_rate3 = 2 * (lasso_correct_rate * lasso_correct_rate2) / (lasso_correct_rate + lasso_correct_rate2)
      #############################################################################################
    }
    
    
    
    
    if (sum(coef_list2 != 0) > 1) {
      out3_biglasso = FDRreg(y, X_final[, coef_list2 != 0], nulltype = 'theoretical' )
      
      LASSO_result = out3_biglasso$FDR 
      LASSO_signal = list(LASSO_result)
      LASSO_selected_var = list(which(coef_list2 != 0))
           
      t_biglasso = proc.time() - t_3
      
      findings_biglasso = which(out3_biglasso$FDR <= fdr_level)
      is_finding_biglasso = rep(0, nobserve)
      is_finding_biglasso[findings_biglasso] = 1
      rates_EBayesFDRreg_biglasso = GetErrorRates(is_signal, is_finding_biglasso)
      
      # Calculate true positive rate for variable selection
      lasso_true_selected = 1:no_assoc_covar  # Original variables
      lasso_selected_vars = which(coef_list2 != 0)
      lasso_correct_rate = sum(lasso_selected_vars %in% lasso_true_selected) / no_assoc_covar
      lasso_correct_rate2 = sum(lasso_selected_vars %in% lasso_true_selected) / length(lasso_selected_vars)
      lasso_correct_rate3 = 2 * (lasso_correct_rate * lasso_correct_rate2) / (lasso_correct_rate + lasso_correct_rate2)
      #############################################################################################
      
      
    }else{

      # If no variables were selected, set results to NA or some default value
      t_biglasso = NA
      findings_biglasso = integer(0)
      is_finding_biglasso = rep(0, nobserve)
      rates_EBayesFDRreg_biglasso = c(NA, NA) 
      
      
      lasso_correct_rate = NA
      lasso_correct_rate2 = NA
      lasso_correct_rate3 = NA
      
      LASSO_signal = list()
      LASSO_selected_var = list()
      
      
      
      BLL = BLL + 1
      
    }
    
    ########################################
    ## Method 4: Marginal association test
    ########################################
    
    selected_marginal = marginal_association_test(X=X_final, y=abs(y))
    
    t_4 = proc.time()
    if (sum(selected_marginal != 0) == 1) {
      mat = matrix(X_final[, selected_marginal != 0], nobserve , 1)
      out3_marginal = FDRreg(y, mat, nulltype = 'theoretical')
      
 
      MA_result = out3_marginal$FDR 
      MA_signal = list(MA_result)
      MA_selected_var = list(which(selected_marginal != 0))

      
      
      t_marginal = proc.time() - t_4
      
      findings_marginal = which(out3_marginal$FDR <= fdr_level)
      is_finding_marginal = rep(0, nobserve)
      is_finding_marginal[findings_marginal] = 1
      rates_EBayesFDRreg_marginal = GetErrorRates(is_signal, is_finding_marginal)
      

      # Calculate true positive rate for variable selection
      MA_true_selected = 1:no_assoc_covar  # Original variables
      MA_selected_vars = which(selected_marginal)
      MA_correct_rate = sum(MA_selected_vars %in% MA_true_selected) / no_assoc_covar
      MA_correct_rate2 = sum(MA_selected_vars %in% MA_true_selected) / length(MA_selected_vars)
      MA_correct_rate3 = 2 * (MA_correct_rate * MA_correct_rate2) / (MA_correct_rate + MA_correct_rate2)

    }
    
    
    if (sum(selected_marginal) > 1) {
      out3_marginal = FDRreg(y, X_final[, selected_marginal], nulltype = 'theoretical')
      
      MA_result = out3_marginal$FDR 
      MA_signal = list(MA_result)
      MA_selected_var = list(which(selected_marginal != 0))
      
      t_marginal = proc.time() - t_4
      
      findings_marginal = which(out3_marginal$FDR <= fdr_level)
      is_finding_marginal = rep(0, nobserve)
      is_finding_marginal[findings_marginal] = 1
      rates_EBayesFDRreg_marginal = GetErrorRates(is_signal, is_finding_marginal)

      # Calculate true positive rate for variable selection
      MA_true_selected = 1:no_assoc_covar  # Original variables
      MA_selected_vars = which(selected_marginal)
      MA_correct_rate = sum(MA_selected_vars %in% MA_true_selected) / no_assoc_covar
      MA_correct_rate2 = sum(MA_selected_vars %in% MA_true_selected) / length(MA_selected_vars)
      MA_correct_rate3 = 2 * (MA_correct_rate * MA_correct_rate2) / (MA_correct_rate + MA_correct_rate2)
      ##################################################################
      
    } else {
      # If no variables were selected, set results to NA or some default value
      t_marginal = NA
      findings_marginal = integer(0)
      is_finding_marginal = rep(0, nobserve)
      rates_EBayesFDRreg_marginal = c(NA, NA)  # or c(0, 0) depending on how you want to handle this case
      
      MA_correct_rate = NA
      MA_correct_rate2 = NA
      MA_correct_rate3 = NA
      
      MA_signal = list()
      MA_selected_var = list()
      
      
      MAA = MAA + 1
    }
    
    
    
    FP_matrix =c(rates_EBayesFDRreg,rates_EBayesFDRreg_enet,rates_EBayesFDRreg_biglasso,
                 rates_EBayesFDRreg_marginal)
    
    Time_matrix = c(t_origFDRreg[3],t_enet[3],t_biglasso[3],t_marginal[3])
    NA_matrix = c(0, ENN, BLL, MAA)
    correct_selection_rates = c(1, enet_correct_rate, lasso_correct_rate, MA_correct_rate)  ### recall 
    correct_selection_rates2 = c(1, enet_correct_rate2, lasso_correct_rate2, MA_correct_rate2)  ### precision 
    correct_selection_rates3 = c(1, enet_correct_rate3, lasso_correct_rate3, MA_correct_rate3) ### f1 score
    
    
    # Return result list
    list(
      save_signal_findings_new = list(signal = signal_true, FDRreg = FDRreg_signal, EN = EN_signal, LASSO = LASSO_signal, MA = MA_signal),
      selected_var_new = list(EN_var = EN_selected_var , LASSO_var = LASSO_selected_var, MA_var = MA_selected_var),
      FP_new = FP_matrix, 
      Time_new = Time_matrix, 
      NA_new = NA_matrix,     
      cr_new = correct_selection_rates,    
      cr2_new = correct_selection_rates2,   
      cr3_new = correct_selection_rates3     
    )
    
  }
  
  ########################################
  ## Aggregate results across simulations
  ########################################

  result_power_FDR[[i]] = do.call(rbind, lapply(results, `[[`, "FP_new"))
  result_time[[i]] = do.call(rbind, lapply(results, `[[`, "Time_new"))
  result_NArate[[i]] = do.call(rbind, lapply(results, `[[`, "NA_new"))
  
  result_correct_selection_rates[[i]] = do.call(rbind, lapply(results, `[[`, "cr_new"))  ### recall 
  result_correct_selection_rates2[[i]] = do.call(rbind, lapply(results, `[[`, "cr2_new"))  ### precision 
  result_correct_selection_rates3[[i]] = do.call(rbind, lapply(results, `[[`, "cr3_new")) ### f1 score
  
  
  save_signal_findings[[i]] = results$save_signal_findings_new
  save_selected_var[[i]] = results$selected_var_new
  
  # General result
  results_list[[counter]] = list(func=i, result= results)
  counter = counter+1
  
}


r1 = result_power_FDR

r2 = result_time

r3 = result_NArate

r4 = result_correct_selection_rates

r5 = result_correct_selection_rates2

r6 = result_correct_selection_rates3

r7 = save_signal_findings

r8 = save_selected_var



############################################
## 8. Post-processing and output
############################################

colMeans_na = function(x) {
  colMeans(x, na.rm = TRUE)
}

save.image("26_01.RData")


result_meanfdr_power = lapply(r1, colMeans_na)
result_meantime = lapply(r2, colMeans_na)

result_NA_rate = lapply(r3, colSums)
result_select_rate = lapply(r4, colMeans_na)
result_select_rate2 = lapply(r5, colMeans_na)
result_select_rate3 = lapply(r6, colMeans_na)


rr = list()
for(ii in seq_along(flist)){
  m1_1 = matrix(result_meanfdr_power[[ii]],2,4)
  m1_2 = matrix(result_meantime[[ii]],1,4)
  m1_3 = matrix(result_NA_rate[[ii]]/nsims,1,4)
  m1_4 = matrix(result_select_rate[[ii]],1,4)
  m1_5 = matrix(result_select_rate2[[ii]],1,4)
  m1_6 = matrix(result_select_rate3[[ii]],1,4)
  
  mm = rbind(m1_1,m1_2,m1_3,m1_4,m1_5,m1_6)
  
  rr[[ii]] = mm
}

final_result = do.call(rbind, rr)


final_result = round(final_result,4)

my_df = as.data.frame(final_result)


# Save in excel
df_name = paste0("N_",nobserve,"_paragroup_", jj,"_sim_",nsims,".xlsx")
write_xlsx(my_df, df_name)

