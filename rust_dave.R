########
# rust #
########

library(R.matlab)
library(ggplot2)
library(Hmisc)

data <- readMat('../rust.mat')

i <- data$it
x <- data$xt

dim(x)

#### transition probability

# j = 0,1,2
zero <- vector()
one <- vector()
two <- vector()

for (k in 1:100){

        print(k)
        
        # for a given bus k
        xk <- x[,k]
        ik <- i[,k]
        jk <- xk[-1]-xk[-1000]
        
        # for periods that i=0
        j <- jk[ik==0]

        # count the number of 0,1,2
        zero[k] <- length(j[j==0])
        one[k] <- length(j[j==1])
        two[k] <- length(j[j==2])
        
        }

# estimated transition probabilities (independent of x_t)
theta_30 = sum(zero)/(sum(zero)+sum(one)+sum(two))
theta_31 = sum(one)/(sum(zero)+sum(one)+sum(two))
theta_32 = sum(two)/(sum(zero)+sum(one)+sum(two))

#### value function iteration at true value of theta

theta_1 = .05

# cost function
cost <- function(x){
        return (theta_1*x)}

theta_30 = .3
theta_31 = .5
theta_32 = .2

beta =.99

RC = 10

# flow utility (without structural errors)
u <- function(x,i){
        -RC*i - cost(x*(1-i))
}

# initial values for EV(x,i), rows represent x=0,1,...,30,31,32
# need extra rows because we need EV(31,i), EV(32,i) to compute EV(30,i)
  EV <- matrix(100,33,2)
  EV2 <- matrix(0,33,2)
  
  value.Iterate <- function(EV){
          for (x in 1:31){ # indexing for x value 0:30
                  # for i=0
          EV2[x,1] <- log(exp(u(x,0)+beta*EV[x,1])+exp(u(x,1)+beta*EV[x,2]))*theta_30 # j = 0
                  + log(exp(u(x=1,0)+beta*EV[x+1,1])+exp(u(x+1,1)+beta*EV[x+1,2]))*theta_31 # j=1
                  + log(exp(u(x+2,0)+beta*EV[x+2,1])+exp(u(x+2,1)+beta*EV[x+2,2]))*theta_32 # j=2
                  
                  # for i=1
          EV2[x,2] <- log(exp(u(1,0)+beta*EV[1,1])+exp(u(1,1)+beta*EV[1,2]))
                  
          }
          return(EV2)
  }
  
  # critical value 
  cri <- 10^(-8)
  
  while(max(abs(EV-EV2))>cri){
    print(EV[1,2])
    EV <- EV2
    EV2 <- value.Iterate(EV)
  }

# get EV(x,i) for x=0,1,2,..,30
EV <- EV[1:31,]
EV

# plot 
df <- data.frame('x'=1:30,'EV0'=EV[2:31,1],'EV1'=EV[2:31,2])

ggplot(aes(x,EV0),data=df) + theme_minimal()+geom_point()+geom_point(aes(x,EV1),color='red')+
        xlab('mileage')+ylab('EV')

#plot(seq(1,201),data$EV[,1])

#### nested fixed point with outer loop (over parameters)

# estimate choice probability from EV

p_i <- matrix(0,31,2)

choice.prob.Estimate <- function(p_i){
        
        for (x in 1:31){
                # for i=0
                p_i[x,1] <- exp(u(x-1,0)+beta*EV[x,1])/(exp(u(x-1,0)+beta*EV[x,1])+exp(u(x-1,1)+beta*EV[x,2]))
                p_i[x,2] <- 1- p_i[x,1]
        }
        
        return(p_i)
}

p_i <- choice.prob.Estimate(p_i)

# log liklihood function

log_choice_prob <- 0
log_transition_prob <- 0
total <- 0

log.likelihood.Compute <- function(){
        for (bus in 1:100){
                for (t in 1:999){
                        if (x[t,bus] <= 30){
                                log_choice_prob <- log(p_i[x[t,bus]+1,i[t,bus]+1]) + log_choice_prob 
                        } else {
                                log_choice_prob <- log(p_i[31,i[t,bus]+1]) + log_choice_prob
                        }
                        
                        
                        if (x[t+1,bus]-x[t,bus]==0){
                                log_transition_prob <- log(theta_30)+log_transition_prob
                        } else if (x[t+1,bus]-x[t,bus]==1) {
                                log_transition_prob <- log(theta_31)+log_transition_prob
                        } else if (x[t+1,bus]-x[t,bus]==2) {
                                log_transition_prob <- log(theta_32)+log_transition_prob
                        }
                }
                
                total <- (log_choice_prob+log_transition_prob) + total
        }
        return (total)
}

log.likelihood.Compute()

# putting things together

# 1) compute theta_30,31,32 directly from the data, this step is not in the loop

# 2) create a grid for theta_1, beta and RC
theta_1_range <- seq(.01,.10,.01)
beta_range <- seq(.90,.99,.01)
RC_range <- seq(6,15,1)

# 3) compute nested fixed point
likelihood <- data.frame('theta_1'=rep(0),'beta'=rep(0),'RC'=rep(0),'log.likelihood'=rep(0))
#likelihood <- data.frame('theta_1'='','beta'='','RC'='','log.likelihood'='')

for (theta_1 in theta_1_range){
        print (theta_1) # just to check the progress
        for (beta in beta_range){
                for (RC in RC_range){
                        # fix the value of parameters
                        
                        # value function iteration 
                        EV <- matrix(100,33,2)
                        EV2 <- matrix(0,33,2)
                        
                        while(max(abs(EV-EV2))>cri){
                                EV <- EV2
                                EV2 <- value.Iterate(EV)
                        }
                        
                        EV <- EV[1:31,]
                        
                        # choice probability 
                        
                        p_i <- matrix(0,31,2)
                        p_i <- choice.prob.Estimate(p_i)
                        
                        # likelihood
                        
                        log_choice_prob <- 0
                        log_transition_prob <- 0
                        total <- 0
                        
                        likelihood <- rbind(likelihood,c(theta_1,beta,RC,log.likelihood.Compute()))
                        
                        
                }
        }
}

likelihood <- likelihood[-1,]
likelihood[which.max(likelihood[,4]),]

#### Hotz-Miller approach 


zero <- vector()
one <- vector()
two <- vector()
d_1_prob <- rep(0, max(x))
x_prob <- rep(0, max(x))

for (k in 1:100){
  
  print(k)
  
  # for a given bus k
  xk <- x[,k]
  ik <- i[,k]
  jk <- xk[-1]-xk[-1000]
  
  # for periods that i=0
  j <- jk[ik==0]
  
  # count the number of 0,1,2
  zero[k] <- length(j[j==0])
  one[k] <- length(j[j==1])
  two[k] <- length(j[j==2])
  
  print(table((xk+1)*ik))
  # Get probabilities for Hotz + Miller
  for (j in names(table((xk+1)*ik))) {
    if (j != '0') { 
      print(j)
      print(as.numeric(j))
      table(xk*ik)[j]
    d_1_prob[as.numeric(j)] = d_1_prob[as.numeric(j)] + table((xk+1)*ik)[j]
    }
  }
  for (j in names(table(xk+1))) {
    x_prob[as.numeric(j)] = x_prob[as.numeric(j)] + table(xk+1)[j]
  }
  
}

# estimated transition probabilities (independent of x_t)
theta_30 = sum(zero)/(sum(zero)+sum(one)+sum(two))
theta_31 = sum(one)/(sum(zero)+sum(one)+sum(two))
theta_32 = sum(two)/(sum(zero)+sum(one)+sum(two))

P_d_1_given_x = d_1_prob/x_prob
P_d_1_given_x = P_d_1_given_x[1:length(P_d_1_given_x)-1]
eulers_constant = .577

theta_1_range <- seq(.01,.10,.01)
S = 100
theta_best <- 1000
diff_best <- 100000
for (theta_1 in theta_1_range){
  V = matrix(0, nrow = length(P_d_1_given_x), ncol = 2)
  for (x in 1:length(P_d_1_given_x)) {
    for (d in 0:1) {
      V_tilde_components = rep(0, S)
      for (s in (1:S)){
        action_sequence = matrix(0, 10, 3)
        action_sequence[1, 1] = x
        action_sequence[1, 2] = d
        if (action_sequence[1, 2] == 1) {
          action_sequence[1, 3] = P_d_1_given_x[action_sequence[1, 1]]
        } else {
          action_sequence[1, 3] = (1-P_d_1_given_x[action_sequence[1, 1]])
        }
        for (j in (2:10)) {
          if (action_sequence[j-1, 2] == 0) { 
            action_sequence[j, 1] = action_sequence[j-1, 1] + multinomial_draw(theta_30, theta_31, theta_32)
            increment = rbinom(1, 1, P_d_1_given_x[action_sequence[j, 1]])
            action_sequence[j, 2] = ifelse(is.na(increment), 1, increment)
            if (action_sequence[j, 2] == 1) {
              action_sequence[j, 3] = P_d_1_given_x[action_sequence[j, 1]]
            } else {
              action_sequence[j, 3] = (1 - P_d_1_given_x[action_sequence[j, 1]])
            }
          } else if (action_sequence[j-1, 2] == 1) { 
            action_sequence[j, 1] = 1
            increment = rbinom(1, 1, P_d_1_given_x[action_sequence[j, 1]])
            action_sequence[j, 2] = ifelse(is.na(increment), 1, increment)
            if (action_sequence[j, 2] == 1) {
              action_sequence[j, 3] = P_d_1_given_x[action_sequence[j, 1]]
            } else {
              action_sequence[j, 3] = (1 - P_d_1_given_x[action_sequence[j, 1]])
            }
          } else {
            action_sequence[j, 1] = 0
            action_sequence[j, 2] = 0
          }
        }
        V_tilde_components[s] = u(action_sequence[1,1], action_sequence[1,2]) + sum((u(action_sequence[2:10,1], 
                                                                                       action_sequence[2:10,2]) + 
                                                                                       eulers_constant - 
                                                                                       log(action_sequence[2:10,3]))*(beta^seq(1,9,1)))
      }
      V[x,d+1] = mean(V_tilde_components, na.rm=TRUE)
      }
  }
  prob_array <- exp(V[,1])/(exp(V[,1]) + exp(V[,2]))
  diff <- sum(prob_array - P_d_1_given_x)^2
  print(diff)
  print(diff_best)
  if (diff < diff_best) { 
    theta_best <- theta_1
    diff_best <- diff
    V_best <- V
  }
}

multinomial_draw <- function(theta_30, theta_31, theta_32) {
  rn <- runif(1, 0, 1)
  if (rn <= theta_30) {
    return(0)
  } else if (rn <= (theta_30 + theta_31)) {
    return(1)
  } else {
    return(2)
  }
}
  