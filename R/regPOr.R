regPOr <- function(time,event,X,order = 3,equal_space = T,nknot, myknots,conv_cri = 1e-6,c_initial = F,initial){

  library(nleqslv)
  library(Matrix)
  library(ggplot2)

  #function for inverse when the information is singular
  myinversev<- function(v, b, g, tol=1e-7)
  { rg=1:length(g); ind=(g<tol); index=rg[ind]+length(b); vv=v[-index,-index]
  var.b=diag(solve(vv))[1:length(b)];
  se.b=sqrt(var.b);
  return(se.b)
  }

  #iv is the X's
  iv <- as.matrix(X)
  #the number of coefficeints of beta
  P <- ncol(iv)
  #del is delta/event 0=alive, 1=dead
  del <- event
  #total number of observations
  n <- length(time)

  if(equal_space == T){
    #the number of interior knots
    ik <- nknot-2
    #equal space knots
    #knots at the right end + 0.01 because all the basis equals to 0 at the end points if not
    mx <-  max(time) + 0.01
    knots <- seq(0,mx,length.out = ik + 2)
  }else{
    knots <- myknots
    ik <- length(knots) - 2
  }

  #Degree for Ispline
  dgr <- order
  K <- dgr + ik

  M <- Mspline(time,dgr,knots)
  I <- Ispline(time,dgr,knots)

  #set an initial value
  if(c_initial == F){
    bt <- rep(1,ncol(iv))
    gama <- rep(1,K)
  }else{
    bt <- initial[1:P]
    gama <- initial[(P+1):(P+K)]
  }

  df <- rep(1,length(bt))
  ite <- 0

  while( t(df)%*%df > conv_cri | t(df)%*%df == 0){
    #Ephi and epsi are 200*1 vectors
    Ephi <- (t(I)%*%gama*exp(iv%*%bt) + 1)^(-1)

    Epsi <- (t(I)%*%gama*exp(iv%*%bt) + 1)^(-1)*del
    #Eu is a 200*K matrix, the ^del is used to get rid of the 0 denomator in right censored data
    Eu <- ( t(M*gama)/as.vector(t(M)%*%gama)^del)*del
    #set the equations needed to be solved
    btf <- function(x){
      y <- numeric(length(bt))
      for(h in 1:length(bt)){
        y[h] <- del%*%iv[,h] - sum( ( t(del%*%Eu)*(I%*%(exp(iv%*%x)*(Ephi + Epsi)*iv[,h])) )  / ( I%*%(exp(iv%*%x)*(Ephi + Epsi)) ) )
      }
      y
    }

    #btstart <- rnorm(length(bt),0,1)
    btstart <- bt
    #solve and get the updated bt
    sol <- nleqslv(btstart,btf,method="Newton")
    btnew <- sol$x

    gamanew <-  t(Eu)%*%del/( I%*%(exp(iv%*%btnew)*(Ephi + Epsi)) )

    ite <- ite + 1
    df <- btnew -bt
    if(t(df)%*%df == 0){bt <- rnorm(ncol(iv),1,0)
    gama <- rexp(K)} else {bt <- as.vector(btnew)
    gama <- as.vector(gamanew)}

  }
  ##########################calculate covariance matrix for the selected model####################
  Ephi <- (t(I)%*%gama*exp(iv%*%bt) + 1)^(-1)
  Epsi[del == 1] <- (t(I[,del == 1])%*%gama*exp(as.matrix(iv[del == 1,])%*%bt) + 1)^(-1)
  #Eu is a 200*K matrix, the ^del is used to get rid of the 0 denomator in right censored data
  Eu[del ==1,] <- t(M[,del == 1]*gama)/as.vector(t(M[,del == 1])%*%gama)
  #variance of phi + psi*del is a n*1 vector
  V <- Ephi^2 + Epsi^2
  ####part 1
  Q <- matrix(0,(K+P),(K+P))
  #elements in xx^t as a vector, for entries as a group
  e <- as.vector(t(do.call(cbind,replicate(P,iv,simplify = FALSE)))) * rep(as.vector(t(iv)),each = P)
  #coefficients for each matrix, repeat each elememt p times so it match the above vector
  co <-  rep(t(I)%*%gama*exp(iv%*%bt)*(Ephi + Epsi), each = P^2)
  #fill in the corresponding place in Q
  Q[seq(1,P),seq(1,P)] <- - matrix(apply(matrix(e*co,P^2,n),1,sum),P,P)

  Q[seq(P+1,K+P),seq(1,P)] <- - t(t(I)*as.vector(exp(iv%*%bt)*(Ephi + Epsi))) %*% iv

  Q[seq(1,P),seq(P+1,P+K)] <- t( Q[seq(P+1,K+P),seq(1,P)])

  diag(Q[seq(1+P,P+K),seq(1+P,P+K)]) <- - gama^(-2)* apply(Eu,2,sum)

  ####part 2
  vc <- matrix(0,(K+P),(K+P))

  covc <-  rep((t(I)%*%gama*exp(iv%*%bt))^2*V, each = P^2)
  vc[seq(1,P),seq(1,P)] <- matrix(apply(matrix(e*covc,P^2,n),1,sum),P,P)

  vc[seq(1,P),seq(P+1,P+K)] <- t(iv)%*%(t(I)*as.vector(t(I)%*%gama*exp(iv%*%bt*2)*V))
  vc[seq(P+1,K+P),seq(1,P)] <- t(vc[seq(1,P),seq(P+1,P+K)])

  for(l in 1:K){
    for(m in 1:K){
      if(l == m){
        vc[P+l,P+m] <- sum( Eu[,l]*(1 - Eu[,l])*del/(gama[l]^2) + I[l,]*I[m,]*exp(iv%*%bt*2)*V )
      } else{
        vc[P+l,P+m] <- sum( -M[l,]*M[m,]*del/((t(M)%*%gama)^(del*2)) + I[l,]*I[m,]*exp(iv%*%bt*2)*V )
      }
    }
  }

  v <- -(Q + vc)
  tol <- 1e-7
  ####part3
  if(sum(gama < tol) == 0){
    vv <- v
  } else {
    rg=1:K
    index=rg[gama < tol]+P
    vv=v[-index,-index]
  }

  if( rcond(vv) > .Machine$double.eps ){
    se_theta = sqrt(diag(solve(vv))[1:P])
  } else {
    se_theta = sqrt(diag(solve(vv + diag(1e-4,nrow(vv),ncol(vv))))[1:P])
  }
  ###############CI###########################
  CI_lower <- bt - 1.96*se_theta
  CI_upper <- bt + 1.96*se_theta
  ci <- cbind(CI_lower,CI_upper)

  #define function log??likelihood AND calculate AIC and BIC
  llhd <- sum(log(((t(M)%*%gama) * exp(iv%*%bt))^del * (t(I)%*%gama*exp(iv%*%bt) + 1)^(-del-1)))
  AIC <- 2*(P+K) - 2*llhd
  BIC <- (P+K)*log(n) - 2*llhd
  #########################################################
  # plots:odds, survival,hazard
  ########################################################
  #grids of time points
  tgrids <- seq(0,mx,0.1)
  #calculate baseline odds
  b <- Ispline(tgrids,dgr,knots)
  m <- Mspline(tgrids,dgr,knots)
  odds <- t(as.matrix(gama)) %*% b
  ff <- t(as.matrix(gama))%*%m

  #check baseline hazard rate
  hzd <- ff/(1+odds)
  bsl_hz = ggplot() + geom_line(data = as.data.frame( cbind(tgrids,t(hzd) ) ),aes(tgrids,t(hzd) )) + labs(x="t",y="h(t)")
  #check baseline odds
  bsl_odds = ggplot() + geom_line(data = as.data.frame( cbind(tgrids,t(odds) ) ),aes(tgrids,t(odds) )) + labs(x="t",y="Baseline odds")
  #check baseline survival
  sur <- 1/(1+odds)
  bsl_surv = ggplot() + geom_line(data = as.data.frame( cbind(tgrids,t(sur) ) ),aes(tgrids,t(sur) )) + labs(x="t",y="S(t)")

  #Return the result
  list(beta = bt,beta_se = se_theta, CI = ci, spline_coef = gama, knots = knots, AIC = AIC, BIC = BIC, Baseline_Surv = bsl_surv, Baseline_hazard = bsl_hz, Baseline_odds = bsl_odds )

}



