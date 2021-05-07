#' Title plotResults
#'
#' @param mech Missing data mechanism of interest (e.g. MCAR,MAR,MNAR)
#' @param lowerlim Lower limit for plotting
#' @param upperlim Upper limit for plotting
#' @param size Size of dataset you are comparing imputation methods on
#' @param wd Specify working directory (should contain results from RestructureResults)
#'
#' @import ggplot2
#' @import gridExtra
#'@import ggpubr
#'@importFrom utils read.csv write.csv
#' @importFrom grDevices dev.off pdf
#' @return plots results comparing imputation methods and saves pdfs
#' @export plotResults
#'
#' @examples
#'  workingdir <-  "C:/Users/kgetz1/Documents/Year2/ComputingProject/GitHub/CompareMI"
#' plotResults("MCAR",-4,4,size=100,wd=workingdir)
#'
plotResults <- function(mech,lowerlim,upperlim,size,wd){
  setwd(wd)
  title <- mech


  if(mech=="MCAR"){
    MARdata <- read.csv(paste0("BiasMCAR",size,".csv"))
    MARdataExtra <- read.csv(paste0("MCARextraFull",size,".csv"))
  }

  if(mech=="MAR"){
    MARdata <- read.csv(paste0("BiasMAR",size,".csv"))
    MARdataExtra <- read.csv(paste0("MARextraFull",size,".csv"))
  }
  if(mech=="MNAR"){
    MARdata <- read.csv(paste0("BiasMNAR",size,".csv"))
    MARdataExtra <- read.csv(paste0("MNARextraFull",size,".csv"))
  }




  MARdatatreat <- MARdata[c(1,3,5,7),]
  MARdatatreat[,1] <- c("MiceTreat","OracleTreat","ccTreat","daeTreat")
  MARdataExtraTreat <- MARdataExtra[c(1,2,5,6,9,10,13,14),]


  dfMAR0.1 <- data.frame(means =MARdatatreat[,2],
                         lower = MARdataExtraTreat[c(1,3,5,7),2],
                         upper = MARdataExtraTreat[c(2,4,6,8),2],
                         feats =MARdatatreat[,1],
                         group = rep(c("MICE","Oracle", "cc","dae"), each = 1))

  plot1 <- ggplot(dfMAR0.1, aes(x = feats, color = group)) + ggtitle(paste("Treat", mech, "0.1")) +  scale_x_discrete(breaks=NULL) +
    xlab("Methods") + ylab("Bias")+ylim(lowerlim,upperlim)+
    geom_errorbar(aes(ymax =upper, ymin = lower),
                  position = "dodge",size = 1.5)+ geom_hline(yintercept=0, linetype="dashed", color = "black")

  dfMAR0.3 <- data.frame(means =MARdatatreat[,3],
                         lower = MARdataExtraTreat[c(1,3,5,7),3],
                         upper = MARdataExtraTreat[c(2,4,6,8),3],
                         feats =MARdatatreat[,1],
                         group = rep(c("MICE","Oracle", "cc","dae"), each = 1))

  plot2 <- ggplot(dfMAR0.3, aes(x = feats, color = group)) +ggtitle(paste("Treat", mech, "0.3")) +scale_x_discrete(breaks=NULL) +
    xlab("Methods") + ylab("Bias")+ylim(lowerlim,upperlim)+
    geom_errorbar(aes(ymax = upper, ymin = lower),
                  position = "dodge",size = 1.5)+geom_hline(yintercept=0, linetype="dashed", color = "black")

  dfMAR0.5 <- data.frame(means =MARdatatreat[,4],
                         lower = MARdataExtraTreat[c(1,3,5,7),4],
                         upper = MARdataExtraTreat[c(2,4,6,8),4],
                         feats =MARdatatreat[,1],
                         group = rep(c("MICE","Oracle", "cc","dae"), each = 1))

  plot3 <- ggplot(dfMAR0.5, aes(x = feats, color = group)) + ggtitle(paste("Treat", mech, "0.5")) + scale_x_discrete(breaks=NULL) +
    xlab("Methods") + ylab("Bias")+ylim(lowerlim,upperlim)+
    geom_errorbar(aes(ymax = upper, ymin =lower),
                  position = "dodge",size = 1.5)+geom_hline(yintercept=0, linetype="dashed", color = "black")

  pdf(paste0(mech,"PlotTreat",size,".pdf"))

  print(ggarrange(plot1, plot2, plot3, ncol=3, common.legend = TRUE, legend="bottom"))
  dev.off()

  MARdatatreat <- MARdata[c(2,4,6,8),]
  MARdatatreat[,1] <- c("MiceBp","OracleBp","ccBp","daeBp")
  MARdataExtraBp <- MARdataExtra[c(3,4,7,8,11,12,15,16),]

  dfMAR0.1 <- data.frame(means =MARdatatreat[,2],
                         lower = MARdataExtraBp[c(1,3,5,7),2],
                         upper = MARdataExtraBp[c(2,4,6,8),2],
                         feats =MARdatatreat[,1],
                         group = c("MICE","Oracle", "cc","dae"))

  plot1 <- ggplot(dfMAR0.1, aes(x = feats, color = group)) +ggtitle(paste("Bp", mech, "0.1")) + scale_x_discrete(breaks=NULL) +
    xlab("Methods") + ylab("Bias")+ ylim(lowerlim,upperlim)+
    geom_errorbar(aes(ymax = upper, ymin =lower),
                  position = "dodge",size = 1.5)+ geom_hline(yintercept=0, linetype="dashed", color = "black")

  dfMAR0.3 <- data.frame(means =MARdatatreat[,3],
                         lower = MARdataExtraBp[c(1,3,5,7),3],
                         upper = MARdataExtraBp[c(2,4,6,8),3],
                         feats =MARdatatreat[,1],
                         group = c("MICE","Oracle", "cc","dae"))

  plot2 <- ggplot(dfMAR0.3, aes(x = feats, color = group)) +ggtitle(paste("Bp", mech, "0.3")) +scale_x_discrete(breaks=NULL) +
    xlab("Methods") + ylab("Bias")+ ylim(lowerlim,upperlim)+
    geom_errorbar(aes(ymax =upper, ymin = lower),
                  position = "dodge",size = 1.5)+geom_hline(yintercept=0, linetype="dashed", color = "black")

  dfMAR0.5 <- data.frame(means =MARdatatreat[,4],
                         lower = MARdataExtraBp[c(1,3,5,7),4],
                         upper = MARdataExtraBp[c(2,4,6,8),4],
                         feats =MARdatatreat[,1],
                         group = c("MICE","Oracle", "cc","dae"))

  plot3 <- ggplot(dfMAR0.5, aes(x = feats, color = group)) +ggtitle(paste("Bp", mech, "0.5")) + scale_x_discrete(breaks=NULL) +
    xlab("Methods") + ylab("Bias")+ylim(lowerlim,upperlim)+
    geom_errorbar(aes(ymax = upper, ymin = lower),
                  position = "dodge",size = 1.5)+geom_hline(yintercept=0, linetype="dashed", color = "black")
  pdf(paste0(mech,"PlotBp",size,".pdf"))
  print(ggarrange(plot1, plot2, plot3, ncol=3, common.legend = TRUE, legend="bottom"))
  dev.off()
}
