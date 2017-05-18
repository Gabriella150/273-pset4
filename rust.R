#### This code uses the methodologies of both Rust (1987) and Hotz and Miller (1993) to estimate the parameters of a single
#### agent dynamic problem where an agent (Harold Zurcher) must choose when to have the engines replaced in a fleet of buses.

## Import libraries
library(R.matlab)
library(ggplot2)

## Read in data
setwd('~/Dropbox (MIT)/MIT/Spring_2017/14.273/HW4/273-pset4/')
data <- readMat('../rust.mat')

## Extract bus replacement events
i <- data$it
## Extract bus mileage counts (in increments of 5,000 miles)
x <- data$xt


## Buses transition from different mileage states, and can jump forward zero, one, or two 5,000 mile buckets. This block
## of code estimates the transition probabilities empirically from the data.

### Initialize empty vectors to hold a count of how many times each jump happens
zero <- vector()
one <- vector()
two <- vector()

### Loop through the 100 buses in the dataset
for (k in 1:100) {
        
        ### Given a bus k, grab the mileage counts
        xk <- x[,k]
        ### Also grab the engine replacement events
        ik <- i[,k]
        ### Get a modified array which gives the change in mileage buckets from period j to period j+1
        jk <- xk[-1]-xk[-1000]
        
        ### We only care about periods where i=0 for transition probabilities, since i=1 will always send
        ### x back to 0. This selects out only time periods for this bus where i = 0
        j <- jk[ik==0]

        ### This counts up how many times the mileage bucket counter, x, moves up by 0, 1, or 2 when i=0
        zero[k] <- length(j[j==0])
        one[k] <- length(j[j==1])
        two[k] <- length(j[j==2])
}

## Estimate the x_t-independent transition probabilities by dividing the number of times for each transition by the  
## total number of transitions
theta_30 = sum(zero)/(sum(zero)+sum(one)+sum(two))
theta_31 = sum(one)/(sum(zero)+sum(one)+sum(two))
theta_32 = sum(two)/(sum(zero)+sum(one)+sum(two))

################
# Question 2.3 #
################

#### We'll now take the true values of the parameter values as given, and use the method described in Rust (1987) to iteratively 
#### estimate the value function (or in this case, the EV function). 

## Initialize parameters to their true values
theta_1 = .05
theta_30 = .3
theta_31 = .5
theta_32 = .2
beta =.99
RC = 10

### Define the linear cost function. If an engine is not replaced, the bus incurs cost theta_1*x, so cost
### increases linearly as a bus gets older.
cost <- function(x){
        return (theta_1*x)}

### Define the utility function at mileage x from action i. If the agent chooses to replace the engine in a bus,
### it costs RC. If they choose _not_ to replace the engine, they incur the cost of running the bus at mileage x.
u <- function(x,i){
        -RC*i - cost(x*(1-i))
}


### The value function can be estimated through an iteration procedure. We start with some initial guess for EV, 
### calculate EV with an expression that includes our initial guess of EV, and continue iterating until the difference
### between subsequent EV estimates becomes small.

#### value.Iterate is a function to iteratively update the value function according to the methodology in Rust. The function
#### takes as an argument a current estimate of EV, and returns an updated estimate of EV. EV is an x by d matrix - we want the
#### EV values for each decision d at every possible current mileage value x.
value.Iterate <- function(EV){
  
  ### First iterate through each of the 30 x states
  for (x in 1:31){
    ## Update the EV value corresponding to not replacing the engine. There are three contributions here - one from the 
    ## j = 0 case, one from the j = 1 case, and one from the j = 2 case. Note the indexing here. When x =1, the state is 
    ## equal to 0 (this is the x that needs to be passed into u()), but we want to grab the EV corresponding to the 1st entry.
    EV2[x,1] <- log(exp(u(x-1,0)+beta*EV[x,1])+exp(u(x-1,1)+beta*EV[x,2]))*theta_30 
    + log(exp(u(x,0)+beta*EV[x+1,1])+exp(u(x,1)+beta*EV[x+1,2]))*theta_31
    + log(exp(u(x+1,0)+beta*EV[x+2,1])+exp(u(x+1,1)+beta*EV[x+2,2]))*theta_32 
    
    ## Update the EV value corresponding to replacing the engine. When the engine is replaced, x at the next period will
    ## deterministically reset to x = 0.
    EV2[x,2] <- log(exp(u(0,0)+beta*EV[1,1])+exp(u(0,1)+beta*EV[1,2]))
  }
  
  ## Return the updated EV values.
  return(EV2)
}

### Set a critical value for to measure the deviation between iterative updates of EV. The distance between the two EV matrices
### is the infinity norm of the difference
cri <- 10^(-8)

### Set an initial value for the EV matrix (all 0s, EV), and another EV object to hold the updated estimates, EV2.
EV <- matrix(100,33,2)
EV2 <- matrix(0,33,2)

## While the infinity norm is less than the threshold, iterate
while(max(abs(EV-EV2))>cri){
  
  ### Set the current EV to the previous updated EV
  EV <- EV2
  ### Compute a new updated EV by iterating on the current EV
  EV2 <- value.Iterate(EV)
}

### Do one last update to set EV equal to the last EV2
EV <- EV2

# get EV(x,i) for x=0,1,2,..,30
### EV contains extra states, which we needed to compute the above computation. Throw them away.
EV <- EV[1:31,]

### Plot the EV of both replacing the engine (i = 1) and not replacing the engine (i = 0) at every x 
### between 1 and 30
df <- data.frame('x'=c(1:30, 1:30),'EV'=c(EV[2:31,1], EV[2:31,2]),'Action' = c(rep('i = 0', 30), rep('i = 1', 30)))

### Generate a plot that compares the EV of replacing the engine and not replacing the engine
ev_plot <- ggplot(df, aes(x=x, y=EV, color=Action)) + geom_point() + xlab('Mileage') + ylab('EV') + 
  ggtitle('EV as a function of mileage and action \n at x between 1 and 30')

### This is a plot to see the EV data in the attached rust matlab file. The state space is different than ours (200 states),
### so its hard to compare. Our's is linear (seems wrong), whereas the provided data is not. However, the first 30 states 
### _do_ look approximately linear, so maybe we're not so far off. 

df_rust <- data.frame('x'=c(seq(1,201), seq(1, 201)), 'EV'=c(data$EV[,1], data$EV[,2]), 'Action' = c(rep('i = 0', 201),
                                                                                                     rep('i = 1', 201)))
ev_plot_rust <- ggplot(df_rust, aes(x=x, y=EV, color=Action)) + geom_point() + xlab('Mileage') + ylab('EV') + 
  ggtitle('Rust dataset EV as a function of mileage and action \n at x between 1 and 201')

################
# Question 2.4 #
################



################
# Question 3.1 #
################

### This code will estimate the parameters beta, theta_1, and RC using the nested fixed-point algorithm
### described in Rust. 

### A function to compute the probability of Zurcher's choices using the EVs we calculate using the EV calculation
### framework above. The probability of choosing different actions basically acts like a multichoice logit function.
choice.prob.Estimate <- function(){
  ### Iterate through all of the possible mileage states, x.
  for (x in 1:31){
    ## For each mileage state x, calculate the probability that Zurcher will choose i = 0.
    p_i[x,1] <- exp(u(x-1,0)+beta*EV[x,1])/(exp(u(x-1,0)+beta*EV[x,1])+exp(u(x-1,1)+beta*EV[x,2]))
    ## Calculate the probability of choosing i = 1 at each state, which is just 1 - P(i = 0).
    p_i[x,2] <- 1- p_i[x,1]
  }
  
  ### Return the updated p_i object
  return(p_i)
}

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

### First we initialize a p_i matrix of all zeros. This is where we will store our estimates at each step.
p_i <- matrix(0,31,2)
log_choice_prob <- 0
log_transition_prob <- 0
total <- 0

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

#p_ix <- p_ix[1:31,]
#p_ix_hat <- p_ix_hat[1:31,]

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
