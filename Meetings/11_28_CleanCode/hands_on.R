library(MASS)
library('bigstep')
###########
#zad 1
zad_1 = 1
n = 1000
p = 950
X = matrix(rnorm(n*p, 0, 1/sqrt(1000)), n, p)
beta = rep(0, p)
beta[1:5] = rep(3, 5)
Y =  X%*%beta
Y1 = Y + rnorm(n)
k = c(5, 10, 50, 100, 500, 950)
RSS = matrix(0, ncol = 6, nrow= 100)
PE = matrix(0, ncol = 6, nrow= 100)
PE_RSS_sigma_z0a = matrix(0, ncol = 6, nrow= 100)
PE_RSS_sigma_niez0a = matrix(0, ncol = 6, nrow= 100)
CV = matrix(0, ncol = 6, nrow = 100)
AIC_sigma_z0a = matrix(0, ncol = 6, nrow = 100)
AIC_sigma_niez0a = matrix(0, ncol = 6, nrow = 100)
for (j in 1:100){
  g = rnorm(n)
  Yg = Y+g
  for(i in 1:length(k)){
    model = lm(Yg~X[,1:k[i]]-1)
    b_est = model$coeff
    #b_est = lsfit(X[,1:k[i]], Yg, intercept = FALSE)$coeff
    #b_est_full = lsfit(X, Y+g, intercept = FALSE)$coeff
    RSS[j,i] = sum(((Y1)-X[,1:k[i]]%*%b_est)^2)
    #RSS_full = sum(((Y+g)-X%*%b_est_full)^2)
    M = X[,1:k[i]]%*%solve(t(X[,1:k[i]])%*%X[,1:k[i]])%*%t(X[,1:k[i]])
    tr = sum(diag(M))
    PE[j, i] = sum((X[,1:k[i]]%*%beta[1:k[i]] - X[,1:k[i]]%*%b_est[1:k[i]] + g)^2)
    PE_RSS_sigma_z0a[j,i] = RSS[j,i] + 2*1*tr
    PE_RSS_sigma_niez0a[j,i] = RSS[j,i] +2*(RSS[j,i]/((n-k[i])*k[i]))
    CV[j,i] = sum((model$residuals/(1-diag(M)))^2)
    AIC_sigma_z0a[j,i] = RSS[j,i] + 2* 1 * k[i]+1000*log(2*pi)
    AIC_sigma_niez0a[j,i] = 2*k[i]+n*log(RSS[j,i]/n)+1000*log(2*pi)
    cat(zad_1, " : ",j, "--", i, "\n")
  }
}



##############
#zad 2 (AIC, BIC)
zad_2 = 2
k= c(20,100,500,950)
prawdziwe_odkrycia_AIC = matrix(0, nrow=100, ncol = length(k))
falszywe_odkrycia_AIC = matrix(0, nrow=100, ncol = length(k))
MSP_AIC = matrix(0, nrow=100, ncol = length(k))
moc_AIC = matrix(0, nrow=100, ncol = length(k))
fdp_AIC = matrix(0, nrow=100, ncol = length(k))

prawdziwe_odkrycia_BIC = matrix(0, nrow=100, ncol = length(k))
falszywe_odkrycia_BIC = matrix(0, nrow=100, ncol = length(k))
MSP_BIC = matrix(0, nrow=100, ncol = length(k))
moc_BIC = matrix(0, nrow=100, ncol = length(k))
fdp_BIC = matrix(0, nrow=100, ncol = length(k))

for (j in 1:100){
  g = rnorm(1000)
  Yg = Y+g
  for(i in 1:length(k)){
    data = prepare_data(Y+g, X[,1:k[i]])
    wynik_aic = fast_forward(data, crit = aic)
    wynik_bic = fast_forward(data, crit = bic)
    #
    b_est_aic = as.integer(wynik_aic$model)
    if(length(b_est_aic) != 0) {
      if(length(b_est_aic) > 1){
        model = lm(Yg~X[,b_est_aic]-1)
        estymatory_bet_aic = model$coeff
        prawdziwe_odkrycia_AIC[j,i] = sum(b_est_aic < 6)
        falszywe_odkrycia_AIC[j,i] = sum(b_est_aic > 5)
        MSP_AIC[j,i] = sum((Y - X[,b_est_aic]%*%estymatory_bet_aic)^2)
        moc_AIC[j,i] = sum(b_est_aic < 6)/5
        fdp_AIC[j,i] = sum(b_est_aic > 5)/max(length(b_est_aic), 1)
      }else{
        model = lm(Yg~X[,b_est_aic]-1)
        estymatory_bet_aic = model$coeff
        prawdziwe_odkrycia_AIC[j,i] = sum(b_est_aic < 6)
        falszywe_odkrycia_AIC[j,i] = sum(b_est_aic > 5)
        MSP_AIC[j,i] = sum((Y - X[,b_est_aic]%*%estymatory_bet_aic)^2)
        moc_AIC[j,i] = sum(b_est_aic < 6)/5
        fdp_AIC[j,i] = sum(b_est_aic > 5)/max(length(b_est_aic), 1)
      }
    }else{
      prawdziwe_odkrycia_AIC[j,i] = 0
      falszywe_odkrycia_AIC[j,i] = 0
      MSP_AIC[j,i] = sum(Yg^2)
      moc_AIC[j,i] = 0
      fdp_AIC[j,i] = 0
    }
    
    #
    b_est_bic = as.integer(wynik_bic$model)
    if(length(b_est_bic) != 0){
      if(length(b_est_bic) > 1){
        model = lm(Yg~X[,b_est_bic]-1)
        estymatory_bet_bic = model$coeff
        prawdziwe_odkrycia_BIC[j,i] = sum(b_est_bic < 6)
        falszywe_odkrycia_BIC[j,i] = sum(b_est_bic > 5)
        MSP_BIC[j,i] = sum((Y - X[,b_est_bic]%*%estymatory_bet_bic)^2)
        moc_BIC[j,i] = sum(b_est_bic < 6)/5
        fdp_BIC[j,i] = sum(b_est_bic > 5)/max(length(b_est_bic), 1)
      }else{
        model = lm(Yg~X[,b_est_bic]-1)
        estymatory_bet_bic = model$coeff
        prawdziwe_odkrycia_BIC[j,i] = sum(b_est_bic < 6)
        falszywe_odkrycia_BIC[j,i] = sum(b_est_bic > 5)
        MSP_BIC[j,i] = sum((Y - X[,b_est_bic]*estymatory_bet_bic)^2)
        moc_BIC[j,i] = sum(b_est_bic < 6)/5
        fdp_BIC[j,i] = sum(b_est_bic > 5)/max(length(b_est_bic), 1)
      }
      
    }else{
      prawdziwe_odkrycia_BIC[j,i] = 0
      falszywe_odkrycia_BIC[j,i] = 0
      MSP_BIC[j,i] = sum(Yg^2)
      moc_BIC[j,i] = 0
      fdp_BIC[j,i] = 0
    }
    cat(zad_2, " : ",j, "--", i, "\n")
  }
}


############
#zad 2 (RIC, mBic i mBic2)
k= c(20,100,500,950)
prawdziwe_odkrycia_RIC = matrix(0, nrow=100, ncol = length(k))
falszywe_odkrycia_RIC = matrix(0, nrow=100, ncol = length(k))
MSP_RIC = matrix(0, nrow=100, ncol = length(k))
moc_RIC = matrix(0, nrow=100, ncol = length(k))
fdp_RIC = matrix(0, nrow=100, ncol = length(k))

prawdziwe_odkrycia_mBIC = matrix(0, nrow=100, ncol = length(k))
falszywe_odkrycia_mBIC = matrix(0, nrow=100, ncol = length(k))
MSP_mBIC = matrix(0, nrow=100, ncol = length(k))
moc_mBIC = matrix(0, nrow=100, ncol = length(k))
fdp_mBIC = matrix(0, nrow=100, ncol = length(k))

prawdziwe_odkrycia_mBIC2 = matrix(0, nrow=100, ncol = length(k))
falszywe_odkrycia_mBIC2 = matrix(0, nrow=100, ncol = length(k))
MSP_mBIC2 = matrix(0, nrow=100, ncol = length(k))
moc_mBIC2 = matrix(0, nrow=100, ncol = length(k))
fdp_mBIC2 = matrix(0, nrow=100, ncol = length(k))

beta[1:5] = rep(4, 5)
Y =  X%*%beta
Y1 = Y + rnorm(n)

for (j in 1:100){
  g=rnorm(1000)
  for(i in 1:length(k)){
    beta[1:5] = rep(sqrt(2*log(p/k[i])), 5)
    Y =  X%*%beta
    Y1 = Y + rnorm(n)
    data = prepare_data(Y+g, X[,1:k[i]])
    my_crit <- function(loglik,k,p=950) {
      -2*loglik+2*k*log(p)
    }
    wynik_ric = stepwise(data, crit = my_crit)
    wynik_mbic = stepwise(data, crit = mbic)
    wynik_mbic2 = stepwise(data, crit = mbic2)
    #
    b_est_ric = as.integer(wynik_ric$model)
    if(length(b_est_ric) != 0){
      if(length(b_est_ric) > 1){
        model = lm(Y1~X[,b_est_ric]-1)
        estymatory_bet_ric = model$coeff
        prawdziwe_odkrycia_RIC[j,i] = sum(b_est_ric < 6)
        falszywe_odkrycia_RIC[j,i] = sum(b_est_ric > 5)
        MSP_RIC[j,i] = sum((Y - X[,b_est_ric]%*%estymatory_bet_ric)^2)
        moc_RIC[j,i] = sum(b_est_ric < 6)/5
        fdp_RIC[j,i] = sum(b_est_ric > 5)/max(length(b_est_ric), 1)
      }else{
        model = lm(Y1~X[,b_est_ric]-1)
        estymatory_bet_ric = model$coeff
        prawdziwe_odkrycia_RIC[j,i] = sum(b_est_ric < 6)
        falszywe_odkrycia_RIC[j,i] = sum(b_est_ric > 5)
        MSP_RIC[j,i] = sum((Y - X[,b_est_ric]*estymatory_bet_ric)^2)
        moc_RIC[j,i] = sum(b_est_ric < 6)/5
        fdp_RIC[j,i] = sum(b_est_ric > 5)/max(length(b_est_ric), 1)
      }
      
    }else{
      prawdziwe_odkrycia_RIC[j,i] = 0
      falszywe_odkrycia_RIC[j,i] = 0
      MSP_RIC[j,i] = sum(Y1^2)
      moc_RIC[j,i] = 0
      fdp_RIC[j,i] = 0
    }
    #
    b_est_mbic = as.integer(wynik_mbic$model)
    if(length(b_est_mbic) != 0){
      if(length(b_est_mbic) > 1){
        model = lm(Y1~X[,b_est_mbic]-1)
        estymatory_bet_mbic = model$coeff
        prawdziwe_odkrycia_mBIC[j,i] = sum(b_est_mbic < 6)
        falszywe_odkrycia_mBIC[j,i] = sum(b_est_mbic > 5)
        MSP_mBIC[j,i] = sum((Y - X[,b_est_mbic]%*%estymatory_bet_mbic)^2)
        moc_mBIC[j,i] = sum(b_est_mbic < 6)/5
        fdp_mBIC[j,i] = sum(b_est_mbic > 5)/max(length(b_est_mbic), 1)
      }else{
        model = lm(Y1~X[,b_est_mbic]-1)
        estymatory_bet_mbic = model$coeff
        prawdziwe_odkrycia_mBIC[j,i] = sum(b_est_mbic < 6)
        falszywe_odkrycia_mBIC[j,i] = sum(b_est_mbic > 5)
        MSP_mBIC[j,i] = sum((Y - X[,b_est_mbic]*estymatory_bet_mbic)^2)
        moc_mBIC[j,i] = sum(b_est_mbic < 6)/5
        fdp_mBIC[j,i] = sum(b_est_mbic > 5)/max(length(b_est_mbic), 1)
      }
    }else{
      prawdziwe_odkrycia_mBIC[j,i] = 0
      falszywe_odkrycia_mBIC[j,i] = 0
      MSP_mBIC[j,i] = sum(Y1^2)
      moc_mBIC[j,i] = 0
      fdp_mBIC[j,i] = 0
    }
    #
    b_est_mbic2 = as.integer(wynik_mbic2$model)
    if(length(b_est_mbic2) != 0){
      
      if(length(b_est_mbic2) > 1){
        model = lm(Y1~X[,b_est_mbic2]-1)
        estymatory_bet_mbic2 = model$coeff
        prawdziwe_odkrycia_mBIC2[j,i] = sum(b_est_mbic2 < 6)
        falszywe_odkrycia_mBIC2[j,i] = sum(b_est_mbic2 > 5)
        MSP_mBIC2[j,i] = sum((Y - X[,b_est_mbic2]%*%estymatory_bet_mbic2)^2)
        moc_mBIC2[j,i] = sum(b_est_mbic2 < 6)/5
        fdp_mBIC2[j,i] = sum(b_est_mbic2 > 5)/max(length(b_est_mbic2), 1)
      }else{
        model = lm(Y1~X[,b_est_mbic2]-1)
        estymatory_bet_mbic2 = model$coeff
        prawdziwe_odkrycia_mBIC2[j,i] = sum(b_est_mbic2 < 6)
        falszywe_odkrycia_mBIC2[j,i] = sum(b_est_mbic2 > 5)
        MSP_mBIC2[j,i] = sum((Y - X[,b_est_mbic2]*estymatory_bet_mbic2)^2)
        moc_mBIC2[j,i] = sum(b_est_mbic2 < 6)/5
        fdp_mBIC2[j,i] = sum(b_est_mbic2 > 5)/max(length(b_est_mbic2), 1)
      }
    }else{
      prawdziwe_odkrycia_mBIC2[j,i] = 0
      falszywe_odkrycia_mBIC2[j,i] = 0
      MSP_mBIC2[j,i] = sum(Y1^2)
      moc_mBIC2[j,i] = 0
      fdp_mBIC2[j,i] = 0
    }
    cat(zad_2, " : ",j, "--", i, "\n")
  }
}
colMeans(prawdziwe_odkrycia_mBIC2)
colMeans(falszywe_odkrycia_RIC)

#############
# zad 3
zad_3 = 3
k=c(950)
prawdziwe_odkrycia_RIC3 = matrix(0, nrow=100, ncol = length(k))
falszywe_odkrycia_RIC3 = matrix(0, nrow=100, ncol = length(k))
MSP_RIC3 = matrix(0, nrow=100, ncol = length(k))
moc_RIC3 = matrix(0, nrow=100, ncol = length(k))
fdp_RIC3 = matrix(0, nrow=100, ncol = length(k))

prawdziwe_odkrycia_mBIC3 = matrix(0, nrow=100, ncol = length(k))
falszywe_odkrycia_mBIC3 = matrix(0, nrow=100, ncol = length(k))
MSP_mBIC3 = matrix(0, nrow=100, ncol = length(k))
moc_mBIC3 = matrix(0, nrow=100, ncol = length(k))
fdp_mBIC3 = matrix(0, nrow=100, ncol = length(k))

prawdziwe_odkrycia_mBIC23 = matrix(0, nrow=100, ncol = length(k))
falszywe_odkrycia_mBIC23 = matrix(0, nrow=100, ncol = length(k))
MSP_mBIC23 = matrix(0, nrow=100, ncol = length(k))
moc_mBIC23 = matrix(0, nrow=100, ncol = length(k))
fdp_mBIC23 = matrix(0, nrow=100, ncol = length(k))
beta[1:50] = rep(4, 50)
Y =  X%*%beta
Y1 = Y + rnorm(n)
#
for (j in 1:100){
  g=rnorm(1000)
  for(i in 1:length(k)){
    beta[1:50] = rep(3, 50)
    Y =  X%*%beta
    Y1 = Y + rnorm(n)
    data = prepare_data(Y+g, X[,1:k[i]])
    my_crit <- function(loglik,k,p=950) {
      -2*loglik+2*k*log(p)
    }
    wynik_ric = stepwise(data, crit = my_crit)
    wynik_mbic = stepwise(data, crit = mbic)
    wynik_mbic2 = stepwise(data, crit = mbic2)
    #
    b_est_ric = as.integer(wynik_ric$model)
    if(length(b_est_ric) != 0){
      
      if(length(b_est_ric) > 1){
        model = lm(Y1~X[,b_est_ric]-1)
        estymatory_bet_ric = model$coeff
        prawdziwe_odkrycia_RIC3[j,i] = sum(b_est_ric < 51)
        falszywe_odkrycia_RIC3[j,i] = sum(b_est_ric > 50)
        MSP_RIC3[j,i] = sum((Y - X[,b_est_ric]%*%estymatory_bet_ric)^2)
        moc_RIC3[j,i] = sum(b_est_ric < 51)/50
        fdp_RIC3[j,i] = sum(b_est_ric > 50)/max(length(b_est_ric), 1)
      }else{
        model = lm(Y1~X[,b_est_ric]-1)
        estymatory_bet_ric = model$coeff
        prawdziwe_odkrycia_RIC3[j,i] = sum(b_est_ric < 51)
        falszywe_odkrycia_RIC3[j,i] = sum(b_est_ric > 50)
        MSP_RIC3[j,i] = sum((Y - X[,b_est_ric]*estymatory_bet_ric)^2)
        moc_RIC3[j,i] = sum(b_est_ric < 51)/50
        fdp_RIC3[j,i] = sum(b_est_ric > 50)/max(length(b_est_ric), 1)
      }
      
    }else{
      prawdziwe_odkrycia_RIC3[j,i] = 0
      falszywe_odkrycia_RIC3[j,i] = 0
      MSP_RIC3[j,i] = sum(Y1^2)
      moc_RIC3[j,i] = 0
      fdp_RIC3[j,i] = 0
    }
    #
    b_est_mbic = as.integer(wynik_mbic$model)
    if(length(b_est_mbic) != 0){
      if(length(b_est_mbic) > 1){
        model = lm(Y1~X[,b_est_mbic]-1)
        estymatory_bet_mbic = model$coeff
        prawdziwe_odkrycia_mBIC3[j,i] = sum(b_est_mbic < 51)
        falszywe_odkrycia_mBIC3[j,i] = sum(b_est_mbic > 50)
        MSP_mBIC3[j,i] = sum((Y - X[,b_est_mbic]%*%estymatory_bet_mbic)^2)
        moc_mBIC3[j,i] = sum(b_est_mbic < 51)/50
        fdp_mBIC3[j,i] = sum(b_est_mbic > 50)/max(length(b_est_mbic), 1)
      }else{
        model = lm(Y1~X[,b_est_mbic]-1)
        estymatory_bet_mbic = model$coeff
        prawdziwe_odkrycia_mBIC3[j,i] = sum(b_est_mbic < 51)
        falszywe_odkrycia_mBIC3[j,i] = sum(b_est_mbic > 50)
        MSP_mBIC3[j,i] = sum((Y - X[,b_est_mbic]*estymatory_bet_mbic)^2)
        moc_mBIC3[j,i] = sum(b_est_mbic < 51)/50
        fdp_mBIC3[j,i] = sum(b_est_mbic > 50)/max(length(b_est_mbic), 1)
      }
    }else{
      prawdziwe_odkrycia_mBIC3[j,i] = 0
      falszywe_odkrycia_mBIC3[j,i] = 0
      MSP_mBIC3[j,i] = sum(Y1^2)
      moc_mBIC3[j,i] = 0
      fdp_mBIC3[j,i] = 0
    }
    #
    b_est_mbic2 = as.integer(wynik_mbic2$model)
    if(length(b_est_mbic2) != 0){
      if(length(b_est_mbic2) > 1){
        model = lm(Y1~X[,b_est_mbic2]-1)
        estymatory_bet_mbic2 = model$coeff
        prawdziwe_odkrycia_mBIC23[j,i] = sum(b_est_mbic2 < 51)
        falszywe_odkrycia_mBIC23[j,i] = sum(b_est_mbic2 > 50)
        MSP_mBIC23[j,i] = sum((Y - X[,b_est_mbic2]%*%estymatory_bet_mbic2)^2)
        moc_mBIC23[j,i] = sum(b_est_mbic2 < 51)/50
        fdp_mBIC23[j,i] = sum(b_est_mbic2 > 50)/max(length(b_est_mbic2), 1)
      }else{
        model = lm(Y1~X[,b_est_mbic2]-1)
        estymatory_bet_mbic23 = model$coeff
        prawdziwe_odkrycia_mBIC23[j,i] = sum(b_est_mbic2 < 51)
        falszywe_odkrycia_mBIC23[j,i] = sum(b_est_mbic2 > 50)
        MSP_mBIC23[j,i] = sum((Y - X[,b_est_mbic2]*estymatory_bet_mbic2)^2)
        moc_mBIC23[j,i] = sum(b_est_mbic2 < 51)/50
        fdp_mBIC23[j,i] = sum(b_est_mbic2 > 50)/max(length(b_est_mbic2), 1)
      }
    }else{
      prawdziwe_odkrycia_mBIC23[j,i] = 0
      falszywe_odkrycia_mBIC23[j,i] = 0
      MSP_mBIC23[j,i] = sum(Y1^2)
      moc_mBIC23[j,i] = 0
      fdp_mBIC23[j,i] = 0
    }
    cat(zad_3, " : ",j, "--", i, "\n")
  }
}


########
# zad 4
zad_4 = 4
n = 1000
p = 950
X = matrix(rnorm(n*p, 0, 1/sqrt(1000)), n, p)
beta= c(rep(10,30), rep(0,920))
Y = X%*%beta
k = c(950)
prawdziwe_odkrycia_mBIC_exp = matrix(0, nrow=100, ncol = length(k))
falszywe_odkrycia_mBIC_exp  = matrix(0, nrow=100, ncol = length(k))
MSP_mBIC_exp  = matrix(0, nrow=100, ncol = length(k))
moc_mBIC_exp  = matrix(0, nrow=100, ncol = length(k))
fdp_mBIC_exp  = matrix(0, nrow=100, ncol = length(k))

prawdziwe_odkrycia_mBIC2_exp = matrix(0, nrow=100, ncol = length(k))
falszywe_odkrycia_mBIC2_exp  = matrix(0, nrow=100, ncol = length(k))
MSP_mBIC2_exp  = matrix(0, nrow=100, ncol = length(k))
moc_mBIC2_exp  = matrix(0, nrow=100, ncol = length(k))
fdp_mBIC2_exp  = matrix(0, nrow=100, ncol = length(k))

prawdziwe_odkrycia_rmBIC_exp = matrix(0, nrow=100, ncol = length(k))
falszywe_odkrycia_rmBIC_exp  = matrix(0, nrow=100, ncol = length(k))
MSP_rmBIC_exp  = matrix(0, nrow=100, ncol = length(k))
moc_rmBIC_exp  = matrix(0, nrow=100, ncol = length(k))
fdp_rmBIC_exp  = matrix(0, nrow=100, ncol = length(k))

prawdziwe_odkrycia_rmBIC2_exp = matrix(0, nrow=100, ncol = length(k))
falszywe_odkrycia_rmBIC2_exp  = matrix(0, nrow=100, ncol = length(k))
MSP_rmBIC2_exp  = matrix(0, nrow=100, ncol = length(k))
moc_rmBIC2_exp  = matrix(0, nrow=100, ncol = length(k))
fdp_rmBIC2_exp  = matrix(0, nrow=100, ncol = length(k))
MSE_exp = matrix(0, nrow = 100, ncol = length(k))

prawdziwe_odkrycia_mBIC_cauchy = matrix(0, nrow=100, ncol = length(k))
falszywe_odkrycia_mBIC_cauchy  = matrix(0, nrow=100, ncol = length(k))
MSP_mBIC_cauchy  = matrix(0, nrow=100, ncol = length(k))
moc_mBIC_cauchy  = matrix(0, nrow=100, ncol = length(k))
fdp_mBIC_cauchy  = matrix(0, nrow=100, ncol = length(k))

prawdziwe_odkrycia_mBIC2_cauchy = matrix(0, nrow=100, ncol = length(k))
falszywe_odkrycia_mBIC2_cauchy  = matrix(0, nrow=100, ncol = length(k))
MSP_mBIC2_cauchy  = matrix(0, nrow=100, ncol = length(k))
moc_mBIC2_cauchy  = matrix(0, nrow=100, ncol = length(k))
fdp_mBIC2_cauchy  = matrix(0, nrow=100, ncol = length(k))

prawdziwe_odkrycia_rmBIC_cauchy = matrix(0, nrow=100, ncol = length(k))
falszywe_odkrycia_rmBIC_cauchy  = matrix(0, nrow=100, ncol = length(k))
MSP_rmBIC_cauchy  = matrix(0, nrow=100, ncol = length(k))
moc_rmBIC_cauchy  = matrix(0, nrow=100, ncol = length(k))
fdp_rmBIC_cauchy  = matrix(0, nrow=100, ncol = length(k))

prawdziwe_odkrycia_rmBIC2_cauchy = matrix(0, nrow=100, ncol = length(k))
falszywe_odkrycia_rmBIC2_cauchy  = matrix(0, nrow=100, ncol = length(k))
MSP_rmBIC2_cauchy  = matrix(0, nrow=100, ncol = length(k))
moc_rmBIC2_cauchy  = matrix(0, nrow=100, ncol = length(k))
fdp_rmBIC2_cauchy  = matrix(0, nrow=100, ncol = length(k))
MSE_cauchy = matrix(0, nrow = 100, ncol = length(k))
beta= c(rep(10,30), rep(0,920))
Y = X%*%beta
i = 1
for (j in 1:100){
  Y_exp = Y + rexp(1000,1)
  data_1_exp = prepare_data(Y_exp, X)
  wynik_mbic = stepwise(data_1_exp, crit = mbic)
  wynik_mbic2 = stepwise(data_1_exp, crit = mbic2)
  Y_exp_rank = rank(Y_exp)
  data_1_exp_rank = prepare_data(Y_exp_rank, X)
  wynik_rmbic = stepwise(data_1_exp_rank, crit = mbic)
  wynik_rmbic2 = stepwise(data_1_exp_rank, crit = mbic2)
  #
  b_est_mbic = as.integer(wynik_mbic$model)
  if(length(b_est_mbic) != 0){
    if(length(b_est_mbic) > 1){
      model = lm(Y_exp~X[,b_est_mbic]-1)
      estymatory_bet_mbic = model$coeff
      prawdziwe_odkrycia_mBIC_exp[j,i] = sum(b_est_mbic < 31)
      falszywe_odkrycia_mBIC_exp[j,i] = sum(b_est_mbic > 30)
      MSP_mBIC_exp[j,i] = sum((Y - X[,b_est_mbic]%*%estymatory_bet_mbic)^2)
      moc_mBIC_exp[j,i] = sum(b_est_mbic < 31)/30
      fdp_mBIC_exp[j,i] = sum(b_est_mbic > 30)/max(length(b_est_mbic), 1)
    }else{
      model = lm(Y_exp~X[,b_est_mbic]-1)
      estymatory_bet_mbic = model$coeff
      prawdziwe_odkrycia_mBIC_exp[j,i] = sum(b_est_mbic < 31)
      falszywe_odkrycia_mBIC_exp[j,i] = sum(b_est_mbic > 30)
      MSP_mBIC_exp[j,i] = sum((Y - X[,b_est_mbic]*estymatory_bet_mbic)^2)
      moc_mBIC_exp[j,i] = sum(b_est_mbic < 31)/30
      fdp_mBIC_exp[j,i] = sum(b_est_mbic > 30)/max(length(b_est_mbic), 1)
    }
  }else{
    prawdziwe_odkrycia_mBIC_exp[j,i] = 0
    falszywe_odkrycia_mBIC_exp[j,i] = 0
    MSP_mBIC_exp[j,i] = sum(Y_exp^2)
    moc_mBIC_exp[j,i] = 0
    fdp_mBIC_exp[j,i] = 0
  }
  #
  b_est_mbic2 = as.integer(wynik_mbic2$model)
  if(length(b_est_mbic2) != 0){
    if(length(b_est_mbic2) > 1){
      model = lm(Y_exp~X[,b_est_mbic2]-1)
      estymatory_bet_mbic2 = model = lm(Y1~X[,b_est_mbic2]-1)$coeff
      prawdziwe_odkrycia_mBIC2_exp[j,i] = sum(b_est_mbic2 < 31)
      falszywe_odkrycia_mBIC2_exp[j,i] = sum(b_est_mbic2 > 30)
      MSP_mBIC2_exp[j,i] = sum((Y - X[,b_est_mbic2]%*%estymatory_bet_mbic2)^2)
      moc_mBIC2_exp[j,i] = sum(b_est_mbic2 < 31)/30
      fdp_mBIC2_exp[j,i] = sum(b_est_mbic2 > 30)/max(length(b_est_mbic2), 1)
    }else{
      model = lm(Y_exp~X[,b_est_mbic2]-1)
      estymatory_bet_mbic2 = model$coeff
      prawdziwe_odkrycia_mBIC2_exp[j,i] = sum(b_est_mbic2 < 31)
      falszywe_odkrycia_mBIC2_exp[j,i] = sum(b_est_mbic2 > 30)
      MSP_mBIC2_exp[j,i] = sum((Y - X[,b_est_mbic2]*estymatory_bet_mbic2)^2)
      moc_mBIC2_exp[j,i] = sum(b_est_mbic2 < 31)/30
      fdp_mBIC2_exp[j,i] = sum(b_est_mbic2 > 30)/max(length(b_est_mbic2), 1)
    }
  }else{
    prawdziwe_odkrycia_mBIC2_exp[j,i] = 0
    falszywe_odkrycia_mBIC2_exp[j,i] = 0
    MSP_mBIC2_exp[j,i] = sum(Y_exp^2)
    moc_mBIC2_exp[j,i] = 0
    fdp_mBIC2_exp[j,i] = 0
  }
  #
  b_est_rmbic = as.integer(wynik_rmbic$model)
  if(length(b_est_rmbic) != 0){
    if(length(b_est_rmbic) > 1){
      model = lm(Y_exp~X[,b_est_rmbic]-1)
      estymatory_bet_rmbic = model$coeff
      prawdziwe_odkrycia_rmBIC_exp[j,i] = sum(b_est_rmbic < 31)
      falszywe_odkrycia_rmBIC_exp[j,i] = sum(b_est_rmbic > 30)
      MSP_rmBIC_exp[j,i] = sum((Y - X[,b_est_rmbic]%*%estymatory_bet_rmbic)^2)
      moc_rmBIC_exp[j,i] = sum(b_est_rmbic < 31)/30
      fdp_rmBIC_exp[j,i] = sum(b_est_rmbic > 30)/max(length(b_est_rmbic), 1)
    }else{
      model = lm(Y_exp~X[,b_est_rmbic]-1)
      estymatory_bet_rmbic = model$coeff
      prawdziwe_odkrycia_rmBIC_exp[j,i] = sum(b_est_rmbic < 31)
      falszywe_odkrycia_rmBIC_exp[j,i] = sum(b_est_rmbic > 30)
      MSP_rmBIC_exp[j,i] = sum((Y - X[,b_est_rmbic]*estymatory_bet_rmbic)^2)
      moc_rmBIC_exp[j,i] = sum(b_est_rmbic < 31)/30
      fdp_rmBIC_exp[j,i] = sum(b_est_rmbic > 30)/max(length(b_est_rmbic), 1)
    }
  }else{
    prawdziwe_odkrycia_rmBIC_exp[j,i] = 0
    falszywe_odkrycia_rmBIC_exp[j,i] = 0
    MSP_rmBIC_exp[j,i] = sum(Y_exp^2)
    moc_rmBIC_exp[j,i] = 0
    fdp_rmBIC_exp[j,i] = 0
  }
  b_est_rmbic2 = as.integer(wynik_rmbic2$model)
  if(length(b_est_rmbic2) != 0){
    if(length(b_est_rmbic2) > 1){
      model = lm(Y_exp~X[,b_est_rmbic2]-1)
      estymatory_bet_rmbic2 = model$coeff
      prawdziwe_odkrycia_rmBIC2_exp[j,i] = sum(b_est_rmbic2 < 31)
      falszywe_odkrycia_rmBIC2_exp[j,i] = sum(b_est_rmbic2 > 30)
      MSP_rmBIC2_exp[j,i] = sum((Y - X[,b_est_rmbic2]%*%estymatory_bet_rmbic2)^2)
      moc_rmBIC2_exp[j,i] = sum(b_est_rmbic2 < 31)/30
      fdp_rmBIC2_exp[j,i] = sum(b_est_rmbic2 > 30)/max(length(b_est_rmbic2), 1)
    }else{
      model = lm(Y_exp~X[,b_est_rmbic2]-1)
      estymatory_bet_rmbic2 = model$coeff
      prawdziwe_odkrycia_rmBIC2_exp[j,i] = sum(b_est_rmbic2 < 31)
      falszywe_odkrycia_rmBIC2_exp[j,i] = sum(b_est_rmbic2 > 30)
      MSP_rmBIC2_exp[j,i] = sum((Y - X[,b_est_rmbic2]*estymatory_bet_rmbic2)^2)
      moc_rmBIC2_exp[j,i] = sum(b_est_rmbic2 < 31)/30
      fdp_rmBIC2_exp[j,i] = sum(b_est_rmbic2 > 30)/max(length(b_est_rmbic2), 1)
    }
  }else{
    prawdziwe_odkrycia_rmBIC2_exp[j,i] = 0
    falszywe_odkrycia_rmBIC2_exp[j,i] = 0
    MSP_rmBIC2_exp[j,i] = sum(Y_exp^2)
    moc_rmBIC2_exp[j,i] = 0
    fdp_rmBIC2_exp[j,i] = 0
  }
  coeff = rlm(Y~X[,b_est_rmbic2]-1)$coeff
  MSE_exp[j,i] = sum((beta[b_est_rmbic2]-coeff)^2)
  
  Y_cauchy = Y + rcauchy(1000)
  data_1_cauchy = prepare_data(Y_cauchy, X)
  wynik_mbic = stepwise(data_1_cauchy, crit = mbic)
  wynik_mbic2 = stepwise(data_1_cauchy, crit = mbic2)
  Y_cauchy_rank = rank(Y_cauchy)
  data_1_cauchy_rank = prepare_data(Y_cauchy_rank, X)
  wynik_rmbic = stepwise(data_1_cauchy_rank, crit = mbic)
  wynik_rmbic2 = stepwise(data_1_cauchy_rank, crit = mbic2)
  
  b_est_mbic = as.integer(wynik_mbic$model)
  if(length(b_est_mbic) != 0){
    if(length(b_est_mbic) > 1){
      model = lm(Y_cauchy~X[,b_est_mbic]-1)
      estymatory_bet_mbic = model$coeff
      prawdziwe_odkrycia_mBIC_cauchy[j,i] = sum(b_est_mbic < 31)
      falszywe_odkrycia_mBIC_cauchy[j,i] = sum(b_est_mbic > 30)
      MSP_mBIC_cauchy[j,i] = sum((Y - X[,b_est_mbic]%*%estymatory_bet_mbic)^2)
      moc_mBIC_cauchy[j,i] = sum(b_est_mbic < 31)/30
      fdp_mBIC_cauchy[j,i] = sum(b_est_mbic > 30)/max(length(b_est_mbic), 1)
    }else{
      model = lm(Y_cauchy~X[,b_est_mbic]-1)
      estymatory_bet_mbic = model$coeff
      prawdziwe_odkrycia_mBIC_cauchy[j,i] = sum(b_est_mbic < 31)
      falszywe_odkrycia_mBIC_cauchy[j,i] = sum(b_est_mbic > 30)
      MSP_mBIC_cauchy[j,i] = sum((Y - X[,b_est_mbic]*estymatory_bet_mbic)^2)
      moc_mBIC_cauchy[j,i] = sum(b_est_mbic < 31)/30
      fdp_mBIC_cauchy[j,i] = sum(b_est_mbic > 30)/max(length(b_est_mbic), 1)
    }
  }else{
    prawdziwe_odkrycia_mBIC_cauchy[j,i] = 0
    falszywe_odkrycia_mBIC_cauchy[j,i] = 0
    MSP_mBIC_cauchy[j,i] = sum(Y_cauchy^2)
    moc_mBIC_cauchy[j,i] = 0
    fdp_mBIC_cauchy[j,i] = 0
  }
  #
  b_est_mbic2 = as.integer(wynik_mbic2$model)
  if(length(b_est_mbic2) != 0){
    if(length(b_est_mbic2) > 1){
      model = lm(Y_cauchy~X[,b_est_mbic2]-1)
      estymatory_bet_mbic2 = model$coeff
      prawdziwe_odkrycia_mBIC2_cauchy[j,i] = sum(b_est_mbic2 < 31)
      falszywe_odkrycia_mBIC2_cauchy[j,i] = sum(b_est_mbic2 > 30)
      MSP_mBIC2_cauchy[j,i] = sum((Y - X[,b_est_mbic2]%*%estymatory_bet_mbic2)^2)
      moc_mBIC2_cauchy[j,i] = sum(b_est_mbic2 < 31)/30
      fdp_mBIC2_cauchy[j,i] = sum(b_est_mbic2 > 30)/max(length(b_est_mbic2), 1)
    }else{
      model = lm(Y_cauchy~X[,b_est_mbic2]-1)
      estymatory_bet_mbic2 = model$coeff
      prawdziwe_odkrycia_mBIC2_cauchy[j,i] = sum(b_est_mbic2 < 31)
      falszywe_odkrycia_mBIC2_cauchy[j,i] = sum(b_est_mbic2 > 30)
      MSP_mBIC2_cauchy[j,i] = sum((Y - X[,b_est_mbic2]*estymatory_bet_mbic2)^2)
      moc_mBIC2_cauchy[j,i] = sum(b_est_mbic2 < 31)/30
      fdp_mBIC2_cauchy[j,i] = sum(b_est_mbic2 > 30)/max(length(b_est_mbic2), 1)
    }
  }else{
    prawdziwe_odkrycia_mBIC2_cauchy[j,i] = 0
    falszywe_odkrycia_mBIC2_cauchy[j,i] = 0
    MSP_mBIC2_cauchy[j,i] = sum(Y_cauchy^2)
    moc_mBIC2_cauchy[j,i] = 0
    fdp_mBIC2_cauchy[j,i] = 0
  }
  #
  b_est_rmbic = as.integer(wynik_rmbic$model)
  if(length(b_est_rmbic) != 0){
    if(length(b_est_rmbic) > 1){
      model = lm(Y_cauchy~X[,b_est_rmbic]-1)
      estymatory_bet_rmbic = model$coeff
      prawdziwe_odkrycia_rmBIC_cauchy[j,i] = sum(b_est_rmbic < 31)
      falszywe_odkrycia_rmBIC_cauchy[j,i] = sum(b_est_rmbic > 30)
      MSP_rmBIC_cauchy[j,i] = sum((Y - X[,b_est_rmbic]%*%estymatory_bet_rmbic)^2)
      moc_rmBIC_cauchy[j,i] = sum(b_est_rmbic < 31)/30
      fdp_rmBIC_cauchy[j,i] = sum(b_est_rmbic > 30)/max(length(b_est_rmbic), 1)
    }else{
      model = lm(Y_cauchy~X[,b_est_rmbic]-1)
      estymatory_bet_rmbic = model$coeff
      prawdziwe_odkrycia_rmBIC_cauchy[j,i] = sum(b_est_rmbic < 31)
      falszywe_odkrycia_rmBIC_cauchy[j,i] = sum(b_est_rmbic > 30)
      MSP_rmBIC_cauchy[j,i] = sum((Y - X[,b_est_rmbic]*estymatory_bet_rmbic)^2)
      moc_rmBIC_cauchy[j,i] = sum(b_est_rmbic < 31)/30
      fdp_rmBIC_cauchy[j,i] = sum(b_est_rmbic > 30)/max(length(b_est_rmbic), 1)
    }
  }else{
    prawdziwe_odkrycia_rmBIC_cauchy[j,i] = 0
    falszywe_odkrycia_rmBIC_cauchy[j,i] = 0
    MSP_rmBIC_cauchy[j,i] = sum(Y_cauchy^2)
    moc_rmBIC_cauchy[j,i] = 0
    fdp_rmBIC_cauchy[j,i] = 0
  }
  b_est_rmbic2 = as.integer(wynik_rmbic2$model)
  if(length(b_est_rmbic2) != 0){
    if(length(b_est_rmbic2) > 1){
      model = lm(Y_cauchy~X[,b_est_rmbic2]-1)
      estymatory_bet_rmbic2 = model$coeff
      prawdziwe_odkrycia_rmBIC2_cauchy[j,i] = sum(b_est_rmbic2 < 31)
      falszywe_odkrycia_rmBIC2_cauchy[j,i] = sum(b_est_rmbic2 > 30)
      MSP_rmBIC2_cauchy[j,i] = sum((Y - X[,b_est_rmbic2]%*%estymatory_bet_rmbic2)^2)
      moc_rmBIC2_cauchy[j,i] = sum(b_est_rmbic2 < 31)/30
      fdp_rmBIC2_cauchy[j,i] = sum(b_est_rmbic2 > 30)/max(length(b_est_rmbic2), 1)
    }else{
      model = lm(Y_cauchy~X[,b_est_rmbic2]-1)
      estymatory_bet_rmbic2 = model$coeff
      prawdziwe_odkrycia_rmBIC2_cauchy[j,i] = sum(b_est_rmbic2 < 31)
      falszywe_odkrycia_rmBIC2_cauchy[j,i] = sum(b_est_rmbic2 > 30)
      MSP_rmBIC2_cauchy[j,i] = sum((Y - X[,b_est_rmbic2]*estymatory_bet_rmbic2)^2)
      moc_rmBIC2_cauchy[j,i] = sum(b_est_rmbic2 < 31)/30
      fdp_rmBIC2_cauchy[j,i] = sum(b_est_rmbic2 > 30)/max(length(b_est_rmbic2), 1)
    }
  }else{
    prawdziwe_odkrycia_rmBIC2_cauchy[j,i] = 0
    falszywe_odkrycia_rmBIC2_cauchy[j,i] = 0
    MSP_rmBIC2_cauchy[j,i] = sum(Y_cauchy^2)
    moc_rmBIC2_cauchy[j,i] = 0
    fdp_rmBIC2_cauchy[j,i] = 0
  }
  coeff = rlm(Y~X[,b_est_rmbic2]-1)$coeff
  MSE_cauchy[j,i] = sum((beta[b_est_rmbic2]-coeff)^2)
  cat(zad_4, " : --", j, "\n")
}


