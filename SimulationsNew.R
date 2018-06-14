#######################################
#    SETUP
#######################################
  # Libraries
    library(lme4)
    library(lmerTest)

  # Variables 
    # variance/sd epsilon
    seps <- 1  

    # Number of simulations
    asim <- 2500

    # Effect size
    delta <- 1

    # Number of participants
    n <- 30
    n.1 <- n/2      # in the first group
    n.2 <- n/2      # in the second group

    # Level of statistical significance 
    alpha <- 0.05


#######################################
# SIMULATIONS
#######################################
  # instead of specifying the correlation between first and second measurement, i
    # would set seps equal to the standard deviation of the true score and let the
    # error variance vary. effect size is defined based on seps. larger error 
    # variance will then imply smaller correlation but also smaller effect size and hence lower power.
  # basically, we then use the following measurement model: 
    # observed score x = true score t + error variance e
  # the predicted correlation between 2 parallel measurements is then equal to var(t)/(var(t)+var(e))


# Define rho, correlation between two measurements
  # to have reasonable range of correlation, define rho and then infer error variance
  # to test, should be smaller steps - i artifically added two rho's
  # we now expect that the maximum possible power (given es and sample size) will
  # be achieved when correlation is perfect and is equal to using a single 
  # measurement with sample size n
    rho <- c(seq(0.01,0.99,0.1),0.95,0.99)
    errvar<-(1-rho)*seps^2/rho
    errsd<-sqrt(errvar)


# Prepare variables to store results
  pow.mean1<-vector("numeric",length(errsd))
  pow.mean3<-vector("numeric",length(errsd))
  pow.mean2<-vector("numeric",length(errsd))
  cori<-vector("numeric",length(errsd))



# Loop over preset correlations between measure 1 and measure 2
for(i in 1:length(errsd)){
  # Create objects to store power in for every simulations
  print(i)
  pow.1<-vector("numeric",asim)
  pow.3<-vector("numeric",asim)
  pow.2<-vector("numeric",asim)
  corobs<-vector("numeric",asim)

  for(k in 1:asim){
    # Scenario 1
      # compare 2 groups of 15 subjects (n=30). 
      # assume that measurements are normally distributed within each group 
      # (sd = 1) and a mean difference of 1 between groups. 
      # In addition, we assume that measurements are perfectly reliable, 
      # i.e. a second measurement for each subject would be identical to the first measurement.
      
      x<-c(rep(1,n.1),rep(0,n.2))
    
      #Vector with observations in the set of participants
      y<-rnorm(n,0,seps)       
      #Add an effect size to the first group
      y[1:n.1]<-y[1:n.1]+delta    
      # Boolean of whether an effect is detected, this is later used to compute the power
      pow.1[k]<-summary(lm(y~x))$coef[2,4]<alpha
    
    # Scenario 2: 
      # additional measurement error is added. This means that measurements 
      # are no longer perfectly reliable. The higher the measurement error, 
      # the smaller the reliability (and hence the smaller the correlation 
      # between consecutive measurements).
      
      # Generating true score t according to measurement model
        y3T<-rnorm(n,0,seps)
        
      # Add effect size
        y3T[1:n.1]<-y3T[1:n.1]+delta

      # Add measurement error
        y3<-y3T+rnorm(n,0,errsd[i])

      # Boolean of whether an effect is detected, this is later used to compute the power
        pow.2[k]<-summary(lm(y3~x))$coef[2,4]<alpha
    
    # Scenario 3
      # we reconsider scenario 2 but now add a second measurement (scan) for each subject. 
      # The more reliable measurements are, the less information will be added by 
      # using this second scan. This scenario thus also converges to the first scenario 
      # when reliability increases.

        x3<-c(rep(1,n.1),rep(0,n.2),rep(1,n.1),rep(0,n.2))
      
      # Vector with observations in the set of participants, correlated with previous vector
        y3.2<-y3T+rnorm(n,0,errsd[i])
        y3o<-c(y3,y3.2)

      # Check correlations
        corobs[k]<-cor(y3,y3.2)

      # Boolean of whether an effect is detected, this is later used to compute the power
        # Define subject numbers
        subject<-rep(1:n,2)
        # Construct mixed model
        mm<-lmer(y3o ~ x3 + (1 | subject))   
        # Boolean of whether an effect is detected, this is later used to compute the power
        pow.3[k]<-summary(mm)$coef[2,5]<alpha
  }
  pow.mean1[i]<-mean(pow.1)
  pow.mean3[i]<-mean(pow.3)
  pow.mean2[i]<-mean(pow.2)
  cori[i]<-mean(corobs)
}


#######################################
# PLOTS
#######################################
  # save results
    write.table(pow.mean1, file = "pow1.RData")
    write.table(pow.mean2, file = "pow2.RData")
    write.table(pow.mean3, file = "pow3.RData")
    write.table(cori, file = "cori.RData")

  # load data if necessary
    pow.mean1 <- read.csv("pow1.RData")
    pow.mean2 <- read.csv("pow2.RData")
    pow.mean3 <- read.csv("pow3.RData")
    cori <- read.csv("cori.RData")  

  # Check correlation between measurement 1 and 2
    plot(rho,cori,xlab="Predicted correlation",ylab="Observed correlation")
    abline(0,1)

plot(cori,pow.mean1,col=1,ylim=c(0,1), type = "l", xlab = '', ylab = '')
par(new = TRUE)
plot(cori,pow.mean2,col=2,ylim=c(0,1), type = "l", xlab = '', ylab = '')
par(new = TRUE)
plot(cori,pow.mean3,xlab="Correlation",ylab="Power",ylim=c(0,1), type = "l", main = "Power in function of correlation", col = 3)
legend(0.7, 0.4, legend = c("Scenario 1", "Scenario 2", "Scenario 3"), col = c(1, 2, 3), lwd = 3)



plot(errsd,pow.mean3,xlab="Measurement error",ylab="Power",ylim=c(0,1), type = "l")
par(new = TRUE)
plot(errsd,pow.mean1,col=2,ylim=c(0,1), type = "l", xlab = '', ylab = '')
par(new = TRUE)
plot(errsd,pow.mean2,col=3,ylim=c(0,1), type = "l", xlab = '', ylab = '')
