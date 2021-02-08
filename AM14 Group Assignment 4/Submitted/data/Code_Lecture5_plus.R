#delete all the existing variables, just for completeness
rm(list = ls())

#make sure these libraries are installed, if they are missing, first run the code like 
# install.packages("tseries")
library(tseries)
library(data.table)
library(zoo)
library(ggplot2)
library(reshape2)
library(lars)

#portfolio returns are from Ken French data library: https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html
#don't forget to keep data files in the same source folder as the code you're running
Returns <- read.csv("industry49.csv", header=T)
Date1 <- "1969-07-01"
Date2 <- "2019-06-01"

#just in case, also save the risk-free rate and the dates
rf <- read.csv("rf.csv", header=T)
rf$Date<-as.Date(rf$Date)
dates<-rf$Date[rf$Date >= Date1 & rf$Date <= Date2]
rf<-rf[rf$Date >= Date1 & rf$Date <= Date2,2]
remove(Date1,Date2)

#number of time series observations
t<-as.numeric(dim(Returns)[1])
#number of portfolios
nport<-as.numeric(dim(Returns)[2])

flag.rolling<-0 #set to 0, if you consider non-overlapping periods of, say, 5 years

T.est<-180 #number of months in training sample, to estimate the mean and variance
T.test<-60 #number of months out of sample, to study investment outcomes

phi<-0.75 #sample value of shrinkage to the mean: the larger it is, larger is shrinkage 
    #of 49 returns to their cross-sectional average of 49
delta<-0.001 #value of the ridge penalty: the larger it is, the higher is the shrinkage of the variance-covariance
  # matrix towards the diagonal 

#number of periods for out-of-sample performance evaluation
if (flag.rolling==0){
  k<-floor((t-T.test-T.est)/T.test)+1
  T.oos<-k*T.test
} else{
  k<-(t-T.test-T.est)+1 
  T.oos<-t-T.est
}

test.periods<-c(1:k)

#creating the variables to store the information about Sharpe ratio and out-of-sample returns for 
  #different strategies
SR<-data.frame("n1"= array(0,c(k,1))) # 1/n portfolio
SR$ridge<-array(0,c(k,1))
SR$lasso<-array(0,c(k,1))
SR$tangency<-array(0,c(k,1))
SR$minvar<-array(0,c(k,1))
SR$shrink<-array(0,c(k,1)) #shrinkage towards the cross-sectional mean

oosret<-data.frame("n1"= array(0,c(T.oos,1)))
oosret$ridge<-array(0,c(T.oos,1))
oosret$lasso<-array(0,c(T.oos,1))
oosret$tangency<-array(0,c(T.oos,1))
oosret$minvar<-array(0,c(T.oos,1))
oosret$shrink<-array(0,c(T.oos,1))


#the loop for OOS evaluation starts
for (i in 1:k)
{
  if (flag.rolling==0){
    est.start<-(i-1)*(T.test)+1 #find the month when parameters estimation starts
  } else{
    est.start<-i 
  }
  est.end<-T.est+est.start-1 #find the month when parameter estimation ends
  oos.start<- est.end+1 #start of the oos investment
  oos.end <-oos.start+T.test-1 #end of the oos investment for this iteration
  
  #this prints out the estimation step
  cat("i=", i, "estimation period:", est.start, " to ", est.end, "evaluation: ", oos.start, " to ", oos.end, "\n")
  
  #extracts the vector of returns for estimation (a matrix of 49 timee series)
  est.returns<-as.matrix(Returns[est.start:est.end,])
  sample.vcov<-cov(est.returns) #estimate variances and covariance of 49 portfolios, all put in a 49 by 49 matrix 
  sample.expret <- colMeans(est.returns) #estimate sample average returns, a 49 by 1 vector
  
  w.tangency<-solve(sample.vcov) %*% sample.expret #standard tangency portfolio that should give in theory the highest SR
  w.tangency<-w.tangency/sum(w.tangency) #normalize the weights to sum up to 1
  w.n1<-array(c(1/nport), c(nport,1)) #1/n strategy
  w.minvar<-solve(sample.vcov) %*% array(1, c(nport,1))  #minimum variance portfolio (in theory)
  w.minvar<-w.minvar/sum(w.minvar) #normalize the weights
  w.shrink<-phi*w.minvar + (1-phi)*w.tangency #weights of the portfolio that shrinks average returns towards the cross-sectional mean
  w.ridge<-solve(sample.vcov+delta * diag(nport)) %*% sample.expret #ridge portfolio weights
  w.ridge<-w.ridge/sum(w.ridge) #normalize to sum up to 1
  
  #lasso bit - reformulating the mean-variance approach by introducing weird y and X variables
  #use the SVD decomposition to build a matrix that is Sigma^(1/2), a square root of matrix Sigma
  vcov.svd<-svd(sample.vcov)
  vcov.sqrt<-vcov.svd$v %*% diag(sqrt(vcov.svd$d)) %*% solve(vcov.svd$v)
  
  #formulate the weird X and y variables, as in the lecture notes
  y<-vcov.sqrt%*%sample.expret
  X<-vcov.sqrt
  
  #compute the lasso path for every possible tuning value of a grid of 100 lambda points
  test.lasso<-lars(X, y, normalize = FALSE, intercept=FALSE, max.steps=100) 
  
  #use a 5-fold cross-validation to find the CV loss function for each tuning value
  test.lasso.cv<-cv.lars(X,y, intercept=FALSE, plot.it = FALSE, K = 5, mode="step")
  
  # s.min constains the value of lambda that leads to the smallest CV value
  s.min <- which.min(test.lasso.cv$cv)
  
  #find 1 standard error of the cross-validating criterion and compute lambda within 1 s.e. or s.min
  cv1se <- min(test.lasso.cv$cv) + test.lasso.cv$cv.error[s.min]
  s.1se <- test.lasso.cv$index[min(which(test.lasso.cv$cv < cv1se))]
  
  #get portfolio weights that correspond to the chosen value of lambda, in this particular case, lambda = s.min 
  coef.lasso <- predict(test.lasso, s=s.min, type="coefficient", mode="step")
  w.lasso<-coef.lasso$coefficients/sum(coef.lasso$coefficients) #normalize the coefficients to sum up to 1
  
  #extract the time series of 49 portfolio returns "out-of-sample" for this iteration
  oos.returns<-as.matrix(Returns[oos.start:oos.end,])
  oos.rf<-rf[oos.start:oos.end]
  
  #find the time series of the 1/n portfolio returns out-of-sample, and its SR
  ret.n1<-oos.returns%*%w.n1
  oosret$n1[(oos.start-T.est):(oos.end-T.est)]<- ret.n1
  SR$n1[i]<- sqrt(12)*mean(ret.n1-oos.rf)/sd(ret.n1)
  
  #returns on the tangency portfolio and its SR
  ret.tangency<-oos.returns%*%w.tangency
  oosret$tangency[(oos.start-T.est):(oos.end-T.est)]<-ret.tangency
  SR$tangency[i]<-sqrt(12)*mean(ret.tangency-oos.rf)/sd(ret.tangency)
  
  #returns on the minimum variance portfolio and its SR
  ret.minvar<-oos.returns%*%w.minvar
  oosret$minvar[(oos.start-T.est):(oos.end-T.est)]<-ret.minvar
  SR$minvar[i]<-sqrt(12)*mean(ret.minvar-oos.rf)/sd(ret.minvar)
  
  #returns on the portfolio constructed with shrinkage towards a cross-sectional mean
  ret.shrink<-oos.returns%*%w.shrink
  oosret$shrink[(oos.start-T.est):(oos.end-T.est)]<-ret.shrink
  SR$shrink[i]<- sqrt(12)*mean(ret.shrink-oos.rf)/sd(ret.shrink)
  
  #rfeturns on the ridge portfolio and its SR
  ret.ridge<-oos.returns%*%w.ridge
  oosret$ridge[(oos.start-T.est):(oos.end-T.est)]<-ret.ridge
  SR$ridge[i]<- sqrt(12)*mean(ret.ridge-oos.rf)/sd(ret.ridge)
  
  #return on the portfolio constructed with lasso, and its SR
  ret.lasso<-oos.returns%*%w.lasso
  oosret$lasso[(oos.start-T.est):(oos.end-T.est)]<-ret.lasso
  SR$lasso[i]<- sqrt(12)*mean(ret.lasso-oos.rf)/sd(ret.lasso)
}

data<- reshape2::melt(SR) #put all the SR info (time series, for each out-of-sample period) together

#plotting the empirical distribution of the SR for all the strategies
graph<-ggplot(data,aes(x=value, fill=variable)) 
graph+
  geom_density(alpha=0.25)+
  labs(x="Annualized SR", y = "Density") +
  scale_fill_discrete(name="", 
                      breaks=c("n1",  "tangency","shrink", "ridge", "lasso", "minvar"),
                      labels=c("1/n",   "Sample tangency", paste("Mean shrinkage", 100*phi, "%"), "Ridge", "Lasso", "Minimum variance")) +
  theme(text = element_text(size=15), legend.position="right", legend.title = element_blank(), legend.text=element_text(size=15))+
  scale_x_continuous(labels = scales::percent)


#computing cumulated total return, starting from 1 dollar for each of the strategies
cumret<-data.frame("n1"=array(c(1, cumprod(1+oosret$n1)), c((T.oos+1),1)))
cumret$ridge<-array(c(1, cumprod(1+oosret$ridge)), c((T.oos+1),1))
cumret$lasso<- array(c(1, cumprod(1+oosret$lasso)), c((T.oos+1),1))
cumret$tangency<-array(c(1, cumprod(1+oosret$tangency)), c((T.oos+1),1))
cumret$shrink<- array(c(1, cumprod(1+oosret$shrink)), c((T.oos+1),1))
cumret$minvar<- array(c(1, cumprod(1+oosret$minvar)), c((T.oos+1),1))
#adding dates info
date.oos.start<-dates[T.est]
date.oos.end<-dates[T.est+T.oos]
dates.oos<-dates[dates >= (date.oos.start-1) & dates <= date.oos.end]
cumret<-cbind(cumret, date=dates.oos)

#plotting cumulative returns for all the strategies
data2<-  reshape2::melt(cumret, id="date")
graph2<-ggplot(data2,aes(x=date, y=value, color=variable)) 
graph2+geom_line()+
  labs(y="Portfolio value, USD", x = "time") +
  scale_color_discrete(name="", 
                       breaks = c("n1", "ridge", "lasso", "tangency", "shrink", "minvar"),
                       labels = c("1/n", paste("Ridge", delta), "Lasso (with minCV)", "Sample tangency", paste("Mean shrinkage", 100*phi, "%"), "Minimum variance"))+
  theme(text = element_text(size=15), legend.position="right", legend.title = element_blank(), legend.text=element_text(size=15))





