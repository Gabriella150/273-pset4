########
# rust #
########

library(R.matlab)
library(ggplot2)

data <- readMat('rust.mat')

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
        EV2[x,1] <- log(exp(u(x-1,0)+beta*EV[x,1])+exp(u(x-1,1)+beta*EV[x,2]))*theta_30 # j = 0
                + log(exp(u(x,0)+beta*EV[x+1,1])+exp(u(x,1)+beta*EV[x+1,2]))*theta_31 # j=1
                + log(exp(u(x+1,0)+beta*EV[x+2,1])+exp(u(x+1,1)+beta*EV[x+2,2]))*theta_32 # j=2
                
                # for i=1
        EV2[x,2] <- log(exp(u(0,0)+beta*EV[1,1])+exp(u(0,1)+beta*EV[1,2]))
                
        }
        return(EV2)
}

# critical value 
cri <- 10^(-8)

while(max(abs(EV-EV2))>cri){
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

plot(seq(1,201),data$EV[,1])

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

# value matrix 
V <- matrix(0,33,2)

# estiamted choice probabiity
p_ix <- matrix(0,33,2)

ones <- vector()
total <- vector()

for (state in 0:32){
        
        a <- 0
        b <- 0
        
        for (bus in 1:100){
                
                a <- sum(i[which(x[,bus]==state),bus]) + a 
                b <- length(i[which(x[,bus]==state),bus]) + b
        }
        
        ones[state+1] <- a 
        total[state+1] <- b
}

ones
total

p_ix[,1] <- 1-ones/total
p_ix[,2] <- ones/total

# simulation

S = 1000
gamma = 0

# when i=0
for (state in 0:30){
        a = 0
        for (s in 1:S){
                x_prime = state + sample(c(0,1,2),1,replace = T, prob = c(theta_30,theta_31,theta_32))
                i_prime = sample(c(0,1),1,replace=T,prob = c(p_ix[x_prime+1,1],p_ix[x_prime+1,2]))
                a = (u(state,0) + beta*(u(x_prime,i_prime)+gamma-log(p_ix[x_prime+1,i_prime+1]))) + a
        }
        
        V[state+1,1] = a/S
}

# when i=1
for (state in 0:30){
        b=0
        for (s in 1:S){
                # x_prime = 0 
                i_prime = sample(c(0,1),1,replace=T,prob = c(p_ix[1,1],p_ix[1,2]))
                b = (u(state,1) + beta*(u(0,i_prime)+gamma-log(p_ix[1,i_prime+1]))) + b
        }
        
        V[state+1,2] = b/S
}

V
p_ix

# estimated choice probability
p_ix_hat <- matrix(0,33,2)

for (state in 0:30){
        p_ix_hat[state+1,1] <- exp(V[state+1,1])/(exp(V[state+1,1])+exp(V[state+1,2]))
        p_ix_hat[state+1,2] <- 1- p_ix_hat[state+1,1]
}
p_ix_hat

p_ix <- p_ix[1:31,]
p_ix_hat <- p_ix_hat[1:31,]

p_ix
p_ix_hat

max(p_ix-p_ix_hat)

# add in the outer loop
theta_1_range <- seq(.01,.10,.01)
beta_range <- seq(.90,.99,.01)
RC_range <- seq(6,15,1)

difference <- data.frame('theta_1'=rep(0),'beta'=rep(0),'RC'=rep(0),'difference'=rep(0))

for (theta_1 in theta_1_range){
        print (theta_1)
        for (beta in beta_range){
                for (RC in RC_range){
                        
                        V <- matrix(0,33,2)
                        
                        # when i=0
for (state in 0:30){
        a = 0
        for (s in 1:S){
                x_prime = state + sample(c(0,1,2),1,replace = T, prob = c(theta_30,theta_31,theta_32))
                i_prime = sample(c(0,1),1,replace=T,prob = c(p_ix[x_prime+1,1],p_ix[x_prime+1,2]))
                a = (u(state,0) + beta*(u(x_prime,i_prime)+gamma-log(p_ix[x_prime+1,i_prime+1]))) + a
        }
        
        V[state+1,1] = a/S
}

# when i=1
for (state in 0:30){
        b=0
        for (s in 1:S){
                # x_prime = 0 
                i_prime = sample(c(0,1),1,replace=T,prob = c(p_ix[1,1],p_ix[1,2]))
                b = (u(state,1) + beta*(u(0,i_prime)+gamma-log(p_ix[1,i_prime+1]))) + b
        }
        
        V[state+1,2] = b/S
}

                p_ix_hat <- matrix(0,33,2)

for (state in 0:30){
        p_ix_hat[state+1,1] <- exp(V[state+1,1])/(exp(V[state+1,1])+exp(V[state+1,2]))
        p_ix_hat[state+1,2] <- 1- p_ix_hat[state+1,1]
}

difference <- rbind(difference,c(theta_1,beta,RC,max(abs(p_ix[1:31,]-p_ix_hat[1:31,]))))
                        
                        }
                        
        }
}

head(difference)
difference <- difference[-1,]
difference[which.min(difference[,4]),]
