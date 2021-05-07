#' RestructureResults
#'
#' @param mech  Missing data mechanism of interest (e.g. MCAR,MAR,MNAR)
#' @param size Size of dataset you want to compare imputation methods on
#' @param wd Specify working directory (should contain AEDatsets folder and python script)
#'
#' @importFrom utils read.csv write.csv
#' @return Writes csv files with reorganized structure to prepare for plotting
#' @export RestructureResults
#'
#' @examples
#'  workingdir <- "C:/Users/kgetz1/Documents/Year2/ComputingProject/GitHub/CompareMI"
#' RestructureResults(mech="MCAR",size=100,wd=workingdir)
#'
RestructureResults<- function(mech,size,wd){
  setwd(wd)
  if(mech=="MCAR"){
    MCAR0.1 <- read.csv("MCAR0.1.csv")
    rownames(MCAR0.1) <- MCAR0.1[,1]
    MCAR0.1 <- t(MCAR0.1[,-1])
    colnames(MCAR0.1) <- paste0(colnames(MCAR0.1),"_","MCAR","_",round(MCAR0.1[1,2],2))
    MCAR0.1 <- MCAR0.1[,-2]
    MCAR0.3 <- read.csv("MCAR0.3.csv")
    rownames(MCAR0.3) <- MCAR0.3[,1]
    MCAR0.3 <- t(MCAR0.3[,-1])
    colnames(MCAR0.3) <- paste0(colnames(MCAR0.3),"_","MCAR","_",round(MCAR0.3[1,2],2))
    MCAR0.3 <- MCAR0.3[,-2]
    MCAR0.5 <- read.csv("MCAR0.5.csv")
    rownames(MCAR0.5) <- MCAR0.5[,1]
    MCAR0.5 <- t(MCAR0.5[,-1])
    colnames(MCAR0.5) <- paste0(colnames(MCAR0.5),"_","MCAR","_",round(MCAR0.5[1,2],2))
    MCAR0.5 <- MCAR0.5[,-2]

    MCAR0.1Extra <- read.csv("MCARExtra0.1.csv")
    MCAR0.3Extra <- read.csv("MCARExtra0.3.csv")
    MCAR0.5Extra <- read.csv("MCARExtra0.5.csv")

    MCARextra <- as.data.frame(rbind(MCAR0.1Extra,MCAR0.3Extra,MCAR0.5Extra))
    MCARextra <- t(MCARextra[,-1])
    colnames(MCARextra) <- c("0.1","0.3","0.5")
    write.csv(MCARextra,paste0("MCARextraFull",size,".csv"))
    MCAR <- as.data.frame(cbind(MCAR0.1,MCAR0.3,MCAR0.5))
    write.csv(MCAR,paste0("BiasMCAR",size,".csv"))

  }

  if(mech=="MAR"){
    MAR0.1 <- read.csv("MAR0.1.csv")
    rownames(MAR0.1) <- MAR0.1[,1]
    MAR0.1 <- t(MAR0.1[,-1])
    MAR0.1 <- MAR0.1[,-2]
    MAR0.3 <- read.csv("MAR0.3.csv")
    rownames(MAR0.3) <- MAR0.3[,1]
    MAR0.3 <- t(MAR0.3[,-1])
    MAR0.3 <- MAR0.3[,-2]
    MAR0.5 <- read.csv("MAR0.5.csv")
    rownames(MAR0.5) <- MAR0.5[,1]
    MAR0.5 <- t(MAR0.5[,-1])
    MAR0.5 <- MAR0.5[,-2]

    MAR0.1Extra <- read.csv("MARExtra0.1.csv")
    MAR0.3Extra <- read.csv("MARExtra0.3.csv")
    MAR0.5Extra <- read.csv("MARExtra0.5.csv")

    MARextra <- as.data.frame(rbind(MAR0.1Extra,MAR0.3Extra,MAR0.5Extra))
    MARextra <- t(MARextra[,-1])
    colnames(MARextra) <- c("0.1","0.3","0.5")
    write.csv(MARextra,paste0("MARextraFull",size,".csv"))
    MAR <- as.data.frame(cbind(MAR0.1,MAR0.3,MAR0.5))
    write.csv(MAR,paste0("BiasMAR",size,".csv"))

  }

  if(mech=="MNAR"){
    MNAR0.1 <- read.csv("MNAR0.1.csv")
    rownames(MNAR0.1) <- MNAR0.1[,1]
    MNAR0.1 <- t(MNAR0.1[,-1])
    colnames(MNAR0.1) <- paste0(colnames(MNAR0.1),"_","MNAR","_",round(MNAR0.1[1,2],2))
    MNAR0.1 <- MNAR0.1[,-2]
    MNAR0.3 <- read.csv("MNAR0.3.csv")
    rownames(MNAR0.3) <- MNAR0.3[,1]
    MNAR0.3 <- t(MNAR0.3[,-1])
    colnames(MNAR0.3) <- paste0(colnames(MNAR0.3),"_","MNAR","_",round(MNAR0.3[1,2],2))
    MNAR0.3 <- MNAR0.3[,-2]
    MNAR0.5 <- read.csv("MNAR0.5.csv")
    rownames(MNAR0.5) <- MNAR0.5[,1]
    MNAR0.5 <- t(MNAR0.5[,-1])
    colnames(MNAR0.5) <- paste0(colnames(MNAR0.5),"_","MNAR","_",round(MNAR0.5[1,2],2))
    MNAR0.5 <- MNAR0.5[,-2]

    MNAR0.1Extra <- read.csv("MNARExtra0.1.csv")
    MNAR0.3Extra <- read.csv("MNARExtra0.3.csv")
    MNAR0.5Extra <- read.csv("MNARExtra0.5.csv")
    MNARextra <- as.data.frame(rbind(MNAR0.1Extra,MNAR0.3Extra,MNAR0.5Extra))
    MNARextra <- t(MNARextra[,-1])
    colnames(MNARextra) <- c("0.1","0.3","0.5")
    write.csv(MNARextra,paste0("MNARextraFull",size,".csv"))
    MNAR <- as.data.frame(cbind(MNAR0.1,MNAR0.3,MNAR0.5))
    write.csv(MNAR,paste0("BiasMNAR",size,".csv"))

  }
}
