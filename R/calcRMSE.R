
#' Title calcRMSE
#'
#' @param mech  Missing data mechanism of interest (e.g. MCAR,MAR,MNAR)
#' @param size  Size of dataset you want to compare imputation methods on
#' @param sdat.c Dataset to sample from for plasmode simulations
#' @param pythonPath specify where on computer python is (desired version <3.8)
#' @param wd set working directory (should contain AEDatsets folder and python script)
#' @return Writes csv with calculated RMSE for DAE and MICE
#'#' @import MASS
#' @import Metrics
#' @import Plasmode
#' @import caret
#' @import data.table
#' @import doParallel
#' @import dplyr
#' @import foreach
#' @import mice
#' @import parallel
#' @import reticulate
#' @import survival
#' @import tibble
#' @import tidyr
#' @import tidyverse
#' @import testthat
#' @importFrom stats complete.cases quantile rbinom rexp rnorm runif sd time var
#' @importFrom utils read.csv write.csv
#'@export calcRMSE
#'
#' @examples
#' Data <- generateData(100)
#'  workingdir <-  "C:/Users/kgetz1/Documents/Year2/ComputingProject/GitHub/CompareMI"
#' path <-  "C:/Users/kgetz1/AppData/Local/Programs/Python/Python38/python.exe"
#' calcRMSE("MCAR",size=100,sdat.c=Data,pythonPath=path,wd=workingdir)
#'
calcRMSE <- function(mech,size,sdat.c,pythonPath,wd){
  setwd(wd)
  if(mech=="MCAR"){
    do.call(file.remove, list(list.files("AEDatasets", full.names = TRUE)))

    cons <- 1.2
    propName <-0.2
    timeEvent <- sdat.c[,c(1,2,3)]
    scaledCov <- as.data.frame(scale(sdat.c[,-c(1,2,3)]))
    sdat.c <- cbind(timeEvent,scaledCov)

    os1 <- coxph(Surv(sdat.c$time1 ,sdat.c$event1) ~ treat+sex+age+bp+hosp+smok+cov1+cov2+cov3+cov4, data=sdat.c, x=T)
    oc1 <- coxph(Surv(sdat.c$time1 ,sdat.c$event1) ~ treat+sex+age+bp+hosp+smok+cov1+cov2+cov3+cov4, data=sdat.c, x=T)

    Iter <- 1:1

    final <- foreach(iterNum = Iter) %dopar% {

      use_python(pythonPath,required=T)
      py_config()
      i <- iterNum
      # Run simulation #1: effect OR = 1.0
      sor <- PlasmodeSur(
        objectOut = os1,
        objectCen = oc1,
        idVar = sdat.c$patientid,
        effectOR =0.5,
        nsim=1,
        size=size)

      ### Generate missingness in simulated data #####################################
      sdat.m <- sdat.c

      cur_sim <- sor
      s.list <- list() # initialize empty storage list
      sim.index <- 1

      # Pull simulation i data
      # Pull simulation i data
      wdat <- cur_sim[["Sim_Data"]]
      names(wdat) <- c("patientid", "event", "time")
      # Attach full data to simulation i
      id <- wdat$patientid
      wdat <- dplyr::left_join(wdat, sdat.m, by="patientid")
      wdat <- wdat[,-1]
      wdat$event <- ifelse(wdat$event == T, 1, 0)
      wdat<-wdat[,-c(3,4)]
      haz <- nelsonaalen(wdat, time, event)
      wdat <- add_column(wdat, hazard = haz, .after = 1)
      CONS <- cons
      miss.coef <- c(rep(log(1),12))
      X <- as.matrix(wdat[,-2])
      p <- exp(CONS+X%*%miss.coef)/(1+exp(CONS+X%*%miss.coef))

      miss.ind <- rbinom(nrow(X),1,p)
      wdatNA <- wdat
      bpCol <- which(colnames(wdatNA)=="bp")
      wdatNA[miss.ind==0,bpCol] <- NA
      s.list[[paste("MissingSim", sim.index, sep = "")]] <- cbind(id,wdatNA)
      s.list[[paste("FullSim", sim.index, sep = "")]] <- cbind(id,wdat)

      prop <- sum(is.na(wdatNA))/nrow(wdatNA)
      trueTreat <- sor$TrueOutBeta[1]
      trueBp <- sor$TrueOutBeta[4]

      mData <-s.list$MissingSim1
      cData <-s.list$FullSim1

      imp <- mice(mData,maxit=0)
      predM = imp$predictorMatrix


      # Setting values of variables I'd like to leave out to 0 in the predictor matrix
      predM[, c("id")]=0
      predM[,c("time")]=0

      imp2 <- mice(mData, m=5,maxit = 5,
                   predictorMatrix = predM,
                   method = "pmm", print =  FALSE)

      RMSE <- c()
      for(q in 1:5){
        anesimp_long <- mice::complete(imp2, action=q, include = F)
        NArows <- which(is.na(mData$bp))
        reconBp <- anesimp_long[NArows,8]
        trueBp <- cData[NArows,8]
        RMSE[q] <- sqrt(sum((reconBp-trueBp)^2)/length(trueBp))
      }
      meanRMSEMICE<- mean(RMSE)


      train <- mData[,-1]

      write.csv(train, paste0("tempTrain",1,".csv"),row.names=F)
      #input_dropout_ratio = 0.5 means half features set tao missing for each row
      daeData<-c()
      daeTreatList <- c()
      daeBpList <- c()
      VarTreatdaeList <- c()
      VarBpdaeList <- c()
      dataReturned <- c()

      source_python("AEDescProj.Py")
      runDAEMIDA(train,i)
      fileList <-list.files(paste0("AEDatasets"))
      while(length(fileList[grepl(paste0("ae_Data_iter",1,"_"),fileList)])<5){
      }
      RMSE <- c()
      for (j in 1:5){
        data <- read.csv(paste0("AEDatasets/dae_Data_iter",1,"_num_",j,".csv"),row.names=1) #read in instead
        NArows <- which(is.na(mData$bp))
        reconBp <- anesimp_long[NArows,8]
        trueBp <- cData[NArows,8]
        RMSE[j] <- sqrt(sum((reconBp-trueBp)^2)/length(trueBp))
      }
      meanRMSEDAE <-  mean(RMSE)

      returnData <- cbind(meanRMSEMICE,meanRMSEDAE)
      return(returnData)
    }

    finalFrame <- as.data.frame(matrix(unlist(final), ncol = 2, byrow = TRUE))
    colnames(finalFrame)<- c('MICE','DAE')


    write.csv(finalFrame,paste0("MCARRMSE",propName,".csv"))
  }

  if(mech=="MAR"){
    consList <- list(1.5,0,-0.5)

    for(w in 1:3){
      cons <- 1.5
      propName <-0.2
      timeEvent <- sdat.c[,c(1,2,3)]
      scaledCov <- as.data.frame(scale(sdat.c[,-c(1,2,3)]))
      sdat.c <- cbind(timeEvent,scaledCov)


      os1 <- coxph(Surv(sdat.c$time1 ,sdat.c$event1) ~ treat+sex+age+bp+hosp+smok+cov1+cov2+cov3+cov4, data=sdat.c, x=T)
      # Censoring hazard
      oc1 <- coxph(Surv(sdat.c$time1 ,sdat.c$event1) ~ treat+sex+age+bp+hosp+smok+cov1+cov2+cov3+cov4, data=sdat.c, x=T)

      Iter <- 1:1

      final <- foreach(iterNum = Iter) %dopar% {
        use_python(pythonPath,required=T)
        py_config()
        i <- iterNum
        # Run simulation #1: effect OR = 1.0
        sor <- PlasmodeSur(
          objectOut = os1,
          objectCen = oc1,
          idVar = sdat.c$patientid,
          effectOR =0.5,
          nsim=1,
          size=size)

        ### Generate missingness in simulated data #####################################
        sdat.m <- sdat.c

        cur_sim <- sor
        s.list <- list() # initialize empty storage list
        sim.index <- 1

        # Pull simulation i data
        wdat <- cur_sim[["Sim_Data"]]
        names(wdat) <- c("patientid", "event", "time")
        # Attach full data to simulation i
        id <- wdat$patientid
        wdat <- dplyr::left_join(wdat, sdat.m, by="patientid")
        wdat <- wdat[,-1]
        wdat$event <- ifelse(wdat$event == T, 1, 0)
        wdat<-wdat[,-c(3,4)]
        haz <- nelsonaalen(wdat, time, event)
        wdat <- add_column(wdat, hazard = haz, .after = 1)
        CONS <- cons
        miss.coef <- c(rep(log(1.1),5),0,rep(log(1.1),6))
        X <- as.matrix(wdat[,-2])
        p <- exp(CONS+X%*%miss.coef)/(1+exp(CONS+X%*%miss.coef))
        miss.ind <- rbinom(nrow(X),1,p)
        wdatNA <- wdat
        bpCol <- which(colnames(wdatNA)=="bp")
        wdatNA[miss.ind==0,bpCol] <- NA
        s.list[[paste("MissingSim", sim.index, sep = "")]] <- cbind(id,wdatNA)
        s.list[[paste("FullSim", sim.index, sep = "")]] <- cbind(id,wdat)
        prop <- sum(is.na(wdatNA))/nrow(wdatNA)
        trueTreat <- sor$TrueOutBeta[1]
        trueBp <- sor$TrueOutBeta[4]


        #MAR data
        mData <-s.list$MissingSim1
        cData <-s.list$FullSim1


        imp <- mice(mData,maxit=0)
        predM = imp$predictorMatrix


        # Setting values of variables I'd like to leave out to 0 in the predictor matrix
        predM[, c("id")]=0
        predM[,c("time")]=0

        imp2 <- mice(mData, m=5,maxit = 5,
                     predictorMatrix = predM,
                     method = "pmm", print =  FALSE)

        RMSE <- c()
        for(q in 1:5){
          anesimp_long <- mice::complete(imp2, action=q, include = F)
          NArows <- which(is.na(mData$bp))
          reconBp <- anesimp_long[NArows,8]
          trueBp <- cData[NArows,8]
          RMSE[q] <-sqrt(sum((reconBp-trueBp)^2)/length(trueBp))
        }
        meanRMSEMICE<- sqrt(mean(RMSE))

        train <- mData[,-1]

        write.csv(train, paste0("tempTrain",i,".csv"),row.names=F)
        #input_dropout_ratio = 0.5 means half features set tao missing for each row
        daeData<-c()
        daeTreatList <- c()
        daeBpList <- c()
        VarTreatdaeList <- c()
        VarBpdaeList <- c()
        dataReturned <- c()

        source_python("AEDescProj.Py")
        runDAEMIDA(train,1)
        fileList <-list.files(paste0("AEDatasets"))
        while(length(fileList[grepl(paste0("ae_Data_iter",1,"_"),fileList)])<5){
        }
        RMSE <- c()
        for (j in 1:5){
          data <- read.csv(paste0("AEDatasets/dae_Data_iter",1,"_num_",j,".csv"),row.names=1) #read in instead
          NArows <- which(is.na(mData$bp))
          reconBp <- anesimp_long[NArows,8]
          trueBp <- cData[NArows,8]
          RMSE[j] <-sqrt(sum((reconBp-trueBp)^2)/length(trueBp))
        }
        meanRMSEDAE <-  mean(RMSE)

        returnData <- cbind(meanRMSEMICE,meanRMSEDAE)
        return(returnData)
      }

      finalFrame <- as.data.frame(matrix(unlist(final), ncol = 2, byrow = TRUE))
      colnames(finalFrame)<- c('MICE','DAE')

      write.csv(finalFrame,paste0("MARRMSE",propName,".csv"))
    }
  }




  if(mech=="MNAR"){
    cons <- 3
    propName <-0.2
    timeEvent <- sdat.c[,c(1,2,3)]
    scaledCov <- as.data.frame(scale(sdat.c[,-c(1,2,3)]))
    sdat.c <- cbind(timeEvent,scaledCov)

    sdat.c$zeta <- rnorm(n=nrow(sdat.c),mean=sdat.c$bp)

    os1 <- coxph(Surv(sdat.c$time1 ,sdat.c$event1) ~ treat+sex+age+bp+hosp+smok+cov1+cov2+cov3+cov4+zeta, data=sdat.c, x=T)
    # Censoring hazard
    oc1 <- coxph(Surv(sdat.c$time1 ,sdat.c$event1) ~ treat+sex+age+bp+hosp+smok+cov1+cov2+cov3+cov4+zeta, data=sdat.c, x=T)


    Iter <- 1:1

    final <- foreach(iterNum = Iter) %dopar% {
      iterStart <-Sys.time()
      use_python(pythonPath,required=T)
      py_config()
      i <- iterNum

      # Run simulation #1: effect OR = 1.0
      sor <- PlasmodeSur(
        objectOut = os1,
        objectCen = oc1,
        idVar = sdat.c$patientid,
        effectOR =0.5,
        nsim=1,
        size=size)
      sdat.m <- sdat.c

      cur_sim <- sor
      s.list <- list() # initialize empty storage list

      sim.index <- 1

      # Pull simulation i data
      wdat <- cur_sim[["Sim_Data"]]
      names(wdat) <- c("patientid", "event", "time")
      # Attach full data to simulation i
      id <- wdat$patientid
      wdat <- dplyr::left_join(wdat, sdat.m, by="patientid")
      wdat <- wdat[,-1]
      wdat$event <- ifelse(wdat$event == T, 1, 0)
      wdat<-wdat[,-c(3,4)]
      haz <- nelsonaalen(wdat, time, event)
      wdat <- add_column(wdat, hazard = haz, .after = 1)
      CONS <- cons
      miss.coef <- c(rep(log(1.1),5),log(5),rep(log(1.1),6),log(5))
      X <- as.matrix(wdat[,-2])
      p <- exp(CONS+X%*%miss.coef)/(1+exp(CONS+X%*%miss.coef))
      miss.ind <- rbinom(nrow(X),1,p)
      wdatNA <- wdat
      bpCol <- which(colnames(wdatNA)=="bp")
      wdatNA[miss.ind==0,bpCol] <- NA
      s.list[[paste("MissingSim", sim.index, sep = "")]] <- cbind(id,wdatNA)
      s.list[[paste("FullSim", sim.index, sep = "")]] <- cbind(id,wdat)
      prop <- sum(is.na(wdatNA))/nrow(wdatNA)

      mData <-s.list$MissingSim1
      cData <-s.list$FullSim1
      imp <- mice(mData,maxit=0)
      predM = imp$predictorMatrix

      # Setting values of variables I'd like to leave out to 0 in the predictor matrix
      predM[, c("id")]=0
      predM[,c("time")]=0
      predM[, c("zeta")]=0

      imp2 <- mice(mData, m=5,maxit = 5,
                   predictorMatrix = predM,
                   method = "pmm", print =  FALSE)

      RMSE <- c()
      for(q in 1:5){
        anesimp_long <- mice::complete(imp2, action=q, include = F)
        NArows <- which(is.na(mData$bp))
        reconBp <- anesimp_long[NArows,8]
        trueBp <- cData[NArows,8]
        RMSE[q] <- sqrt(sum((reconBp-trueBp)^2)/length(trueBp))
      }
      meanRMSEMICE<- mean(RMSE)
      zetaCol <- which(colnames(mData)=="zeta")
      train <- mData[,-c(1,15)]

      write.csv(train, paste0("tempTrain",i,".csv"),row.names=F)
      #input_dropout_ratio = 0.5 means half features set tao missing for each row
      daeData<-c()


      daeTreatList <- c()
      daeBpList <- c()
      VarTreatdaeList <- c()
      VarBpdaeList <- c()
      dataReturned <- c()

      source_python("AEDescProj.Py")
      runDAEMIDA(train,i)
      fileList <-list.files(paste0("AEDatasets"))
      while(length(fileList[grepl(paste0("ae_Data_iter",1,"_"),fileList)])<5){
      }
      RMSE <- c()
      for (j in 1:5){
        data <- read.csv(paste0("AEDatasets/dae_Data_iter",1,"_num_",j,".csv"),row.names=1) #read in instead
        NArows <- which(is.na(mData$bp))
        reconBp <- anesimp_long[NArows,8]
        trueBp <- cData[NArows,8]
        RMSE[j] <- sqrt(sum((reconBp-trueBp)^2)/length(trueBp))
      }
      meanRMSEDAE <-  mean(RMSE)

      returnData <- cbind(meanRMSEMICE,meanRMSEDAE)
      return(returnData)
    }

    finalFrame <- as.data.frame(matrix(unlist(final), ncol = 2, byrow = TRUE))
    colnames(finalFrame)<- c('MICE','DAE')

    write.csv(finalFrame,paste0("MNARRMSE",propName,".csv"))
  }
}
