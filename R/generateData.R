
#' Generate Data
#'
#' @param size Size of dataset you want to draw from when doing plasmode simulations
#'
#' @return Survival dataset of desired size with approximately 20% censored data with 10 variables
#'
#' @importFrom stats rnbinom
#' @import survival

#' @export generateData
#' @examples
#' generateData(100)
#' generateData(1000)
#'
#'
generateData <- function(size){
  b1= 0.69
  treat <-rbinom(size,1,0.6)
  sex <- rbinom(size,1,0.4)
  hosp <- rbinom(size,1,0.2)
  smok <- rbinom(size,1,0.25)
  age <- rnorm(size, 0.6)
  bp <- rnorm(size, mean(age)+0.3*sex+0.9*treat)

  corr <- 0.65
  sigma <- matrix(c(1,corr,corr,corr,corr,1,corr,corr,corr,corr,1,corr,corr,corr,corr,1),nrow=4)
  mu <- as.vector(c(0.3*mean(bp)+0.9*mean(age),0.3*mean(bp)+0.9*mean(age),0.3*mean(bp)+0.9*mean(age),0.3*mean(bp)+0.9*mean(age)))
  mvn <- mvrnorm(size,mu=mu,Sigma= sigma)
  cov1 <- mvn[,1]
  cov2 <- mvn[,2]
  cov3 <- mvn[,3]
  cov4 <- mvn[,4]

  x <- cbind(treat,sex,age,bp,hosp,smok,cov1,cov2,cov3,cov4)

  entryt1 <- runif(size, 0,12)
  lambda1 = exp(0.08+b1*1 -0.1*x[,1]-0.3*x[,2]+0.5*x[,3]+0.15*x[,4]-0.01*x[,5]-0.2*x[,6]+0.3*x[,7]+0.3*x[,8]+0.8*x[,9]+0.1*x[,10])
  t1 = rexp(size,lambda1)
  c1 = rep(10-entryt1,1)
  Data <- as.data.frame(cbind(t1,c1))
  Data$time <- apply(Data,1,min)
  Data$Event <- ifelse(t1<c1,1,0)
  patientid <- seq(1,size)
  sdat.c <- as.data.frame(cbind(patientid,Data$time,Data$Event,x))
  colnames(sdat.c)<- c("patientid","time1","event1",colnames(x))
  return(sdat.c)
}
