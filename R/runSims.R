#' runSims
#'
#' @param mech Missing data mechanism of interest (e.g. MCAR,MAR,MNAR)
#' @param size Size of dataset you want to compare imputation methods on
#' @param iterations Number of iterations for simulation
#' @param sdat.c Dataset to sample from for plasmode simulations
#' @param parallelOption Set equal to T if you want to parallelize
#' @param coreNum Specify number of core if parallelizing
#' @param pythonPath specify where on computer python is (desired version <3.8)
#' @param wd set working directory (should contain AEDatsets folder and python script)
#' @import MASS
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


#' @return Writes files with information about bias and percentiles for each imputation method
#' @export runSims
#'
#' @examples
#' Data <- generateData(100)
#' path <-  "C:/Users/kgetz1/AppData/Local/Programs/Python/Python38/python.exe"
#' workingdir <-  "C:/Users/kgetz1/Documents/Year2/ComputingProject/GitHub/CompareMI"
#' runSims("MCAR",size=100,iterations=2,sdat.c=Data,parallelOption=FALSE,pythonPath=path,wd=workingdir)


runSims <- function(mech,size,iterations,sdat.c,parallelOption,coreNum=1,pythonPath,wd){
  setwd(wd)
  startWholetime <- Sys.time()
  tensorflow <- NULL
  math <- NULL
  numpy <- NULL
  random <- NULL
  sklearn <- NULL
  sys <- NULL
  os <- NULL
  .onLoad <- function(libname, pkgname) {
    tensorflow <<- import("tensorflow", delay_load = TRUE)
  }
  .onLoad <- function(libname, pkgname) {
    math <<- import("math", delay_load = TRUE)
  }
  .onLoad <- function(libname, pkgname) {
    numpy <<-  import("numpy", delay_load = TRUE)
  }
  .onLoad <- function(libname, pkgname) {
    random <<- import("random", delay_load = TRUE)
  }
  .onLoad <- function(libname, pkgname) {
    sklearn <<- import("sklearn", delay_load = TRUE)
  }
  .onLoad <- function(libname, pkgname) {
    sys <<- import("sys", delay_load = TRUE)
  }
  .onLoad <- function(libname, pkgname) {
    os <<- import("os", delay_load = TRUE)
  }
  skip_if_no_tensorflow <- function() {
    have_tensorflow  <- py_module_available("tensorflow")
    if (!have_tensorflow)
      skip("tensorflow not available for testing")
  }
  test_that("Things work as expected", {
    skip_if_no_tensorflow()
    # test code here...
  })
  skip_if_no_math <- function() {
    have_math  <- py_module_available("math")
    if (!have_math)
      skip("math not available for testing")
  }
  test_that("Things work as expected", {
    skip_if_no_math()
    # test code here...
  })
  skip_if_no_numpy <- function() {
    have_numpy  <- py_module_available("numpy")
    if (!have_numpy)
      skip("numpy not available for testing")
  }
  test_that("Things work as expected", {
    skip_if_no_numpy()
    # test code here...
  })
  skip_if_no_random <- function() {
    have_random  <- py_module_available("random")
    if (!have_random)
      skip("random not available for testing")
  }
  test_that("Things work as expected", {
    skip_if_no_random()
    # test code here...
  })

  skip_if_no_sklearn <- function() {
    have_sklearn <- py_module_available("sklearn")
    if (!have_sklearn)
      skip("sklearn not available for testing")
  }
  test_that("Things work as expected", {
    skip_if_no_sklearn()
    # test code here...
  })
  skip_if_no_sys <- function() {
    have_sys  <- py_module_available("sys")
    if (!have_sys)
      skip("sys not available for testing")
  }
  test_that("Things work as expected", {
    skip_if_no_sys()
    # test code here...
  })
  skip_if_no_os <- function() {
    have_os  <- py_module_available("os")
    if (!have_os)
      skip("os not available for testing")
  }
  test_that("Things work as expected", {
    skip_if_no_sys()
    # test code here...
  })


  if(mech=="MCAR"){
    registerDoSEQ()

    for(w in 1:3){
      do.call(file.remove, list(list.files("AEDatasets", full.names = TRUE)))

      consList <- list(2,1,0)
      propList <- list(0.1,0.3,0.5)
      cons <- consList[[w]]
      propName <-propList[[w]]


      os1 <- coxph(Surv(sdat.c$time1 ,sdat.c$event1) ~ treat+sex+age+bp+hosp+smok+cov1+cov2+cov3+cov4, data=sdat.c, x=T)
      # Censoring hazard
      oc1 <- coxph(Surv(sdat.c$time1 ,sdat.c$event1) ~ treat+sex+age+bp+hosp+smok+cov1+cov2+cov3+cov4, data=sdat.c, x=T)


      if(parallelOption==T){
       # doParallel::stopImplicitCluster()
      #  registerDoParallel(coreNum)
        cl <- makeCluster(coreNum,type="SOCK")
        registerDoParallel(cl)
      }

      Iter <- 1:iterations

      final <- foreach(iterNum = Iter) %dopar% {
        library(Plasmode)
        library(dplyr)
        library(survival)
        library(MASS)
        library(tidyverse)
        library(caret)
        library(mice)
        library(tidyr)
        library(Metrics)
        library(data.table)
        library(parallel)
        library(reticulate)
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


        anesimp_long <- mice::complete(imp2, action="long", include = TRUE)
        anesimp_long_mids<-as.mids(anesimp_long)

        fitimp1 <- with(anesimp_long_mids,
                        coxph(Surv(time, event) ~ treat+sex + age +
                                bp + hosp + smok + cov1 + cov2 + cov3 +
                                cov4))


        pool <- summary(pool(fitimp1))

        TreatRow <-which(pool$term =="treat")
        bpRow <-which(pool$term=="bp")

        MICEtreat <- pool$estimate[TreatRow]
        MICEBp <-  pool$estimate[bpRow]

        biasTreatMICE <- as.numeric(MICEtreat-trueTreat)
        biasBpMICE <- as.numeric(MICEBp-trueBp)


        #compare complete dataset to OS1 estimates
        CompleteFit <-  coxph(Surv(time, event) ~ treat+sex + age +
                                bp + hosp + smok + cov1 + cov2 + cov3 +
                                cov4,cData)


        Completetreat <- CompleteFit$coefficients[TreatRow]
        CompleteBp <-  CompleteFit$coefficients[bpRow]

        biasTreatComplete <- as.numeric(Completetreat-trueTreat)
        biasBpComplete <- as.numeric(CompleteBp-trueBp)


        ExcludeData <- mData[complete.cases(mData),]


        ExcludeFit <-  coxph(Surv(time, event) ~treat+sex + age + bp + hosp + smok + cov1 + cov2 + cov3 +cov4,data=ExcludeData)

        Excludetreat <- ExcludeFit$coefficients[TreatRow]
        ExcludeBp <-  ExcludeFit$coefficients[bpRow]

        biasTreatExclude <- as.numeric(Excludetreat-trueTreat)
        biasBpExclude <- as.numeric(ExcludeBp-trueBp)


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
        runDAEMIDA(train,i)
        fileList <-list.files(paste0("AEDatasets"))
        while(length(fileList[grepl(paste0("ae_Data_iter",i,"_"),fileList)])<5){
        }
        for (j in 1:5){
          data <- read.csv(paste0("AEDatasets/dae_Data_iter",i,"_num_",j,".csv"),row.names=1) #read in instead
          daeFit <- coxph(Surv(time, event) ~treat+sex + age + bp + hosp + smok + cov1 + cov2 + cov3 +cov4,data=data)
          daeTreatList[j] <- daeFit$coefficients[1]
          daeBpList[j] <- daeFit$coefficients[7]
          VarTreatdaeList[j] <- (summary(daeFit)$coefficients[TreatRow,3])^2
          VarBpdaeList[j] <- (summary(daeFit)$coefficients[bpRow,3])^2
        }
        #now try to pool them
        daetreat <-mean(daeTreatList)
        daeBp <-mean(daeBpList)

        biasTreatdae <- as.numeric(daetreat-trueTreat)
        biasBpdae <- as.numeric(daeBp-trueBp)

        UbarTreat <- mean(VarTreatdaeList)
        UbarBp <- mean(VarBpdaeList)

        BTreat <- var(daeTreatList)
        BBp <- var(daeBpList)

        Ttreat <- UbarTreat + (1+(1/5))*BTreat #variance
        TBp <- UbarBp + (1+(1/5))*BBp

        returnData <- cbind(biasTreatMICE,biasBpMICE,biasTreatComplete,biasBpComplete,biasTreatExclude,biasBpExclude,biasTreatdae,biasBpdae)
        return(returnData)
      }
      if(parallelOption==T){
      #  doParallel::stopImplicitCluster()
        stopCluster(cl)
      }


      finalFrame <- as.data.frame(matrix(unlist(final), ncol = 8, byrow = TRUE))
      colnames(finalFrame)<- c('biasTreatMICE','biasBpMICE','biasTreatComplete','biasBpComplete','biasTreatExclude','biasBpExclude','biasTreatdae','biasBpdae')


      AvgbiasTreatMICE <- mean(finalFrame$biasTreatMICE)
      AvgbiasBpMICE <- mean(finalFrame$biasBpMICE)
      AvgbiasTreatComplete <- mean(finalFrame$biasTreatComplete)
      AvgbiasBpComplete <- mean(finalFrame$biasBpComplete)
      AvgbiasTreatExclude <- mean(finalFrame$biasTreatExclude)
      AvgbiasBpExclude <- mean(finalFrame$biasBpExclude)
      AvgbiasTreatdae <- mean(finalFrame$biasTreatdae)
      AvgbiasBpdae <- mean(finalFrame$biasBpdae)

      PercentileTreatMICE <-  quantile(finalFrame$biasTreatMICE, c(0.025, 0.975))
      PercentileBpMICE <- quantile(finalFrame$biasBpMICE, c(0.025, 0.975))
      PercentileTreatComplete <- quantile(finalFrame$biasTreatComplete, c(0.025, 0.975))
      PercentileBpComplete <-quantile(finalFrame$biasBpComplete, c(0.025, 0.975))
      PercentileTreatExclude <-quantile(finalFrame$biasTreatExclude, c(0.025, 0.975))
      PercentileBpExclude <-quantile(finalFrame$biasBpExclude, c(0.025, 0.975))
      PercentileTreatDAE <-quantile(finalFrame$biasTreatdae, c(0.025, 0.975))
      PercentileBpDAE <-quantile(finalFrame$biasBpdae, c(0.025, 0.975))



      names <- c("MiceTreat","MiceBp","oracleTreat","oracleBp","complete-caseTreat","complete-caseBp","daetreat","daeBp")
      bias <- c(AvgbiasTreatMICE,AvgbiasBpMICE,AvgbiasTreatComplete,AvgbiasBpComplete,AvgbiasTreatExclude,AvgbiasBpExclude,AvgbiasTreatdae,AvgbiasBpdae)
      Percentile <- c(PercentileTreatMICE,PercentileBpMICE,PercentileTreatComplete,PercentileBpComplete,PercentileTreatExclude,PercentileBpExclude,PercentileTreatDAE,PercentileBpDAE)
      Perc <- t(as.data.frame(unlist(Percentile)))
      namesPerc <- c("MiceTreatLower","MiceTreatUpper","MiceBpLower","MiceBpUpper","oracleTreatLower","oracleTreatUpper","oracleBpLower","oracleBpUpper","complete-caseTreatLower","complete-caseTreatUpper","complete-caseBpLower","complete-caseBpUpper","daetreatLower","daetreatUpper","daeBpLower","daeBpUpper")
      colnames(Perc) <- namesPerc
      rownames(Perc) <- c()
      Extra <- cbind(Perc)
      results <- data.frame(rbind(bias,propName))
      colnames(results) <- names

      write.csv(results,paste0("MCAR",propName,".csv"))
      write.csv(Extra,paste0("MCARExtra",propName,".csv"))

    }

    wholeTime <-Sys.time() -startWholetime
    print(wholeTime)
  }








  if(mech=="MAR"){
    registerDoSEQ()

    for(w in 1:3){
      do.call(file.remove, list(list.files("AEDatasets", full.names = TRUE)))

      consList <- list(1.5,0,-0.5)
      propList <- list(0.1,0.3,0.5)
      cons <- consList[[w]]
      propName <-propList[[w]]
      os1 <- coxph(Surv(sdat.c$time1 ,sdat.c$event1) ~ treat+sex+age+bp+hosp+smok+cov1+cov2+cov3+cov4, data=sdat.c, x=T)
      # Censoring hazard
      oc1 <- coxph(Surv(sdat.c$time1 ,sdat.c$event1) ~ treat+sex+age+bp+hosp+smok+cov1+cov2+cov3+cov4, data=sdat.c, x=T)

      if(parallelOption==T){
        doParallel::stopImplicitCluster()
        registerDoParallel(coreNum)
      }


      Iter <- 1:iterations

      final <- foreach(iterNum = Iter) %dopar% {
        library(Plasmode)
        library(dplyr)
        library(survival)
        library(MASS)
        library(tidyverse)
        library(caret)
        library(mice)
        library(tidyr)
        library(Metrics)
        library(data.table)
        library(parallel)
        library(reticulate)
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


        anesimp_long <- mice::complete(imp2, action="long", include = TRUE)
        anesimp_long_mids<-as.mids(anesimp_long)

        fitimp1 <- with(anesimp_long_mids,
                        coxph(Surv(time, event) ~ treat+sex + age +
                                bp + hosp + smok + cov1 + cov2 + cov3 +
                                cov4))



        pool <- summary(pool(fitimp1))

        TreatRow <-which(pool$term =="treat")
        bpRow <-which(pool$term=="bp")

        MICEtreat <- pool$estimate[TreatRow]
        MICEBp <-  pool$estimate[bpRow]

        biasTreatMICE <- as.numeric(MICEtreat-trueTreat)
        biasBpMICE <- as.numeric(MICEBp-trueBp)


        #compare complete dataset to OS1 estimates
        CompleteFit <-  coxph(Surv(time, event) ~ treat+sex + age +
                                bp + hosp + smok + cov1 + cov2 + cov3 +
                                cov4,cData)


        Completetreat <- CompleteFit$coefficients[TreatRow]
        CompleteBp <-  CompleteFit$coefficients[bpRow]

        biasTreatComplete <- as.numeric(Completetreat-trueTreat)
        biasBpComplete <- as.numeric(CompleteBp-trueBp)


        ExcludeData <- mData[complete.cases(mData),]


        ExcludeFit <-  coxph(Surv(time, event) ~treat+sex + age + bp + hosp + smok + cov1 + cov2 + cov3 +cov4,data=ExcludeData)

        Excludetreat <- ExcludeFit$coefficients[TreatRow]
        ExcludeBp <-  ExcludeFit$coefficients[bpRow]

        biasTreatExclude <- as.numeric(Excludetreat-trueTreat)
        biasBpExclude <- as.numeric(ExcludeBp-trueBp)

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
        runDAEMIDA(train,i)
        fileList <-list.files(paste0("AEDatasets"))
        while(length(fileList[grepl(paste0("ae_Data_iter",i,"_"),fileList)])<5){
        }
        for (j in 1:5){
          data <- read.csv(paste0("AEDatasets/dae_Data_iter",i,"_num_",j,".csv"),row.names=1) #read in instead
          daeFit <- coxph(Surv(time, event) ~treat+sex + age + bp + hosp + smok + cov1 + cov2 + cov3 +cov4,data=data)
          daeTreatList[j] <- daeFit$coefficients[1]
          daeBpList[j] <- daeFit$coefficients[7]
          VarTreatdaeList[j] <- (summary(daeFit)$coefficients[TreatRow,3])^2
          VarBpdaeList[j] <- (summary(daeFit)$coefficients[bpRow,3])^2
        }
        #now try to pool them
        daetreat <-mean(daeTreatList)
        daeBp <-mean(daeBpList)

        biasTreatdae <- as.numeric(daetreat-trueTreat)
        biasBpdae <- as.numeric(daeBp-trueBp)

        UbarTreat <- mean(VarTreatdaeList)
        UbarBp <- mean(VarBpdaeList)

        BTreat <- var(daeTreatList)
        BBp <- var(daeBpList)

        Ttreat <- UbarTreat + (1+(1/5))*BTreat #variance
        TBp <- UbarBp + (1+(1/5))*BBp

        returnData <- cbind(biasTreatMICE,biasBpMICE,biasTreatComplete,biasBpComplete,biasTreatExclude,biasBpExclude,biasTreatdae,biasBpdae)
        return(returnData)
      }
      if(parallelOption==T){
        doParallel::stopImplicitCluster()
      }

      finalFrame <- as.data.frame(matrix(unlist(final), ncol = 8, byrow = TRUE))
      colnames(finalFrame)<- c('biasTreatMICE','biasBpMICE','biasTreatComplete','biasBpComplete','biasTreatExclude','biasBpExclude','biasTreatdae','biasBpdae')


      AvgbiasTreatMICE <- mean(finalFrame$biasTreatMICE)
      AvgbiasBpMICE <- mean(finalFrame$biasBpMICE)
      AvgbiasTreatComplete <- mean(finalFrame$biasTreatComplete)
      AvgbiasBpComplete <- mean(finalFrame$biasBpComplete)
      AvgbiasTreatExclude <- mean(finalFrame$biasTreatExclude)
      AvgbiasBpExclude <- mean(finalFrame$biasBpExclude)
      AvgbiasTreatdae <- mean(finalFrame$biasTreatdae)
      AvgbiasBpdae <- mean(finalFrame$biasBpdae)

      PercentileTreatMICE <-  quantile(finalFrame$biasTreatMICE, c(0.025, 0.975))
      PercentileBpMICE <- quantile(finalFrame$biasBpMICE, c(0.025, 0.975))
      PercentileTreatComplete <- quantile(finalFrame$biasTreatComplete, c(0.025, 0.975))
      PercentileBpComplete <-quantile(finalFrame$biasBpComplete, c(0.025, 0.975))
      PercentileTreatExclude <-quantile(finalFrame$biasTreatExclude, c(0.025, 0.975))
      PercentileBpExclude <-quantile(finalFrame$biasBpExclude, c(0.025, 0.975))
      PercentileTreatDAE <-quantile(finalFrame$biasTreatdae, c(0.025, 0.975))
      PercentileBpDAE <-quantile(finalFrame$biasBpdae, c(0.025, 0.975))



      names <- c("MiceTreat","MiceBp","oracleTreat","oracleBp","complete-caseTreat","complete-caseBp","daetreat","daeBp")
      bias <- c(AvgbiasTreatMICE,AvgbiasBpMICE,AvgbiasTreatComplete,AvgbiasBpComplete,AvgbiasTreatExclude,AvgbiasBpExclude,AvgbiasTreatdae,AvgbiasBpdae)
      Percentile <- c(PercentileTreatMICE,PercentileBpMICE,PercentileTreatComplete,PercentileBpComplete,PercentileTreatExclude,PercentileBpExclude,PercentileTreatDAE,PercentileBpDAE)
      Perc <- t(as.data.frame(unlist(Percentile)))
      namesPerc <- c("MiceTreatLower","MiceTreatUpper","MiceBpLower","MiceBpUpper","oracleTreatLower","oracleTreatUpper","oracleBpLower","oracleBpUpper","complete-caseTreatLower","complete-caseTreatUpper","complete-caseBpLower","complete-caseBpUpper","daetreatLower","daetreatUpper","daeBpLower","daeBpUpper")
      colnames(Perc) <- namesPerc
      rownames(Perc) <- c()
      Extra <- cbind(Perc)
      results <- data.frame(rbind(bias,propName))
      colnames(results) <- names

      write.csv(results,paste0("MAR",propName,".csv"))
      write.csv(Extra,paste0("MARExtra",propName,".csv"))

    }

    wholeTime <-Sys.time() -startWholetime
    print(wholeTime)
  }

  if(mech=="MNAR"){
    registerDoSEQ()

    for(w in 1:3){
      do.call(file.remove, list(list.files("AEDatasets", full.names = TRUE)))

      consList <- list(0,-3,-5) #0.03,0.1,0.12
      propList <- list(0.1,0.3,0.5)
      cons <- consList[[w]]
      propName <-propList[[w]]
      sdat.c$zeta <- rnorm(n=nrow(sdat.c),mean=sdat.c$bp)

      os1 <- coxph(Surv(sdat.c$time1 ,sdat.c$event1) ~ treat+sex+age+bp+hosp+smok+cov1+cov2+cov3+cov4+zeta, data=sdat.c, x=T)
      # Censoring hazard
      oc1 <- coxph(Surv(sdat.c$time1 ,sdat.c$event1) ~ treat+sex+age+bp+hosp+smok+cov1+cov2+cov3+cov4+zeta, data=sdat.c, x=T)
      # Get true estimates
      os1True <- coxph(Surv(sdat.c$time1 ,sdat.c$event1) ~ treat+sex+age+bp+hosp+smok+cov1+cov2+cov3+cov4, data=sdat.c, x=T)

      oc1True <- coxph(Surv(sdat.c$time1 ,sdat.c$event1) ~ treat+sex+age+bp+hosp+smok+cov1+cov2+cov3+cov4,data=sdat.c, x=T)

      sor <- PlasmodeSur(
        objectOut = os1True,
        objectCen = oc1True,
        idVar = sdat.c$patientid,
        effectOR = 0.5,
        nsim=1,
        size=1000000)

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
      s.list[[paste("FullSim", sim.index, sep = "")]] <- cbind(id,wdat)

      cData <-s.list$FullSim1
      TrueFit <- coxph(Surv(time,event) ~ treat+sex+age+bp+hosp+smok+cov1+cov2+cov3+cov4, data=cData, x=T)


      TreatRow <-which(names(TrueFit$coefficients) =="treat")
      BpRow <-which(names(TrueFit$coefficients)=="bp")
      trueTreat <-TrueFit$coefficients[TreatRow]
      trueBp <-  TrueFit$coefficients[BpRow]

      if(parallelOption==T){
        doParallel::stopImplicitCluster()
        registerDoParallel(coreNum)
      }


      Iter <- 1:iterations

      final <- foreach(iterNum = Iter) %dopar% {
        library(Plasmode)
        library(dplyr)
        library(survival)
        library(MASS)
        library(tidyverse)
        library(caret)
        library(mice)
        library(tidyr)
        library(Metrics)
        library(data.table)
        library(parallel)
        library(reticulate)
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


        anesimp_long <- mice::complete(imp2, action="long", include = TRUE)
        anesimp_long_mids<-as.mids(anesimp_long)


        fitimp1 <- with(anesimp_long_mids,
                        coxph(Surv(time, event) ~ treat+sex + age +
                                bp + hosp + smok + cov1 + cov2 + cov3 +
                                cov4))




        pool <- summary(pool(fitimp1))

        TreatRow <-which(pool$term =="treat")
        bpRow <-which(pool$term=="bp")

        MICEtreat <- pool$estimate[TreatRow]
        MICEBp <-  pool$estimate[bpRow]

        biasTreatMICE <- as.numeric(MICEtreat-trueTreat)
        biasBpMICE<- as.numeric(MICEBp-trueBp)


        #compare complete dataset to OS1 estimates
        CompleteFit <-  coxph(Surv(time, event) ~ treat+sex+age+bp+hosp+smok+cov1+cov2+cov3+cov4,cData)

        Completetreat <- CompleteFit$coefficients[TreatRow]
        CompleteBp <-  CompleteFit$coefficients[bpRow]

        biasTreatComplete <- as.numeric(Completetreat-trueTreat)
        biasBpComplete <- as.numeric(CompleteBp-trueBp)


        ExcludeData <- mData[complete.cases(mData),]


        ExcludeFit <-  coxph(Surv(time, event) ~treat+sex+age+bp+hosp+smok+cov1+cov2+cov3+cov4,data=ExcludeData)

        Excludetreat <- ExcludeFit$coefficients[TreatRow]
        ExcludeBp <-  ExcludeFit$coefficients[bpRow]

        biasTreatExclude <- as.numeric(Excludetreat-trueTreat)
        biasBpExclude <- as.numeric(ExcludeBp-trueBp)

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
        while(length(fileList[grepl(paste0("ae_Data_iter",i,"_"),fileList)])<5){
        }
        for (j in 1:5){
          data <- read.csv(paste0("AEDatasets/dae_Data_iter",i,"_num_",j,".csv"),row.names=1) #read in instead
          daeFit <- coxph(Surv(time, event) ~treat+sex + age + bp + hosp + smok + cov1 + cov2 + cov3 +cov4,data=data)
          daeTreatList[j] <- daeFit$coefficients[1]
          daeBpList[j] <- daeFit$coefficients[7]
          VarTreatdaeList[j] <- (summary(daeFit)$coefficients[TreatRow,3])^2
          VarBpdaeList[j] <- (summary(daeFit)$coefficients[bpRow,3])^2
        }
        #now try to pool them
        daetreat <-mean(daeTreatList)
        daeBp <-mean(daeBpList)

        biasTreatdae <- as.numeric(daetreat-trueTreat)
        biasBpdae <- as.numeric(daeBp-trueBp)

        UbarTreat <- mean(VarTreatdaeList)
        UbarBp <- mean(VarBpdaeList)

        BTreat <- var(daeTreatList)
        BBp <- var(daeBpList)

        Ttreat <- UbarTreat + (1+(1/5))*BTreat #variance
        TBp <- UbarBp + (1+(1/5))*BBp

        iterTime <- Sys.time() -iterStart
        returnData <- cbind(biasTreatMICE,biasBpMICE,biasTreatComplete,biasBpComplete,biasTreatExclude,biasBpExclude,biasTreatdae,biasBpdae,iterTime)
        return(returnData)
      }
      if(parallelOption==T){
        doParallel::stopImplicitCluster()
      }

      finalFrame <- as.data.frame(matrix(unlist(final), ncol = 9, byrow = TRUE))
      colnames(finalFrame)<- c('biasTreatMICE','biasBpMICE','biasTreatComplete','biasBpComplete','biasTreatExclude','biasBpExclude','biasTreatdae','biasBpdae','iterTime')

      AvgbiasTreatMICE <- mean(finalFrame$biasTreatMICE)
      AvgbiasBpMICE <- mean(finalFrame$biasBpMICE)
      AvgbiasTreatComplete <- mean(finalFrame$biasTreatComplete)
      AvgbiasBpComplete <- mean(finalFrame$biasBpComplete)
      AvgbiasTreatExclude <- mean(finalFrame$biasTreatExclude)
      AvgbiasBpExclude <- mean(finalFrame$biasBpExclude)
      AvgbiasTreatdae <- mean(finalFrame$biasTreatdae)
      AvgbiasBpdae <- mean(finalFrame$biasBpdae)

      PercentileTreatMICE <-  quantile(finalFrame$biasTreatMICE, c(0.025, 0.975))
      PercentileBpMICE <- quantile(finalFrame$biasBpMICE, c(0.025, 0.975))
      PercentileTreatComplete <- quantile(finalFrame$biasTreatComplete, c(0.025, 0.975))
      PercentileBpComplete <-quantile(finalFrame$biasBpComplete, c(0.025, 0.975))
      PercentileTreatExclude <-quantile(finalFrame$biasTreatExclude, c(0.025, 0.975))
      PercentileBpExclude <-quantile(finalFrame$biasBpExclude, c(0.025, 0.975))
      PercentileTreatDAE <-quantile(finalFrame$biasTreatdae, c(0.025, 0.975))
      PercentileBpDAE <-quantile(finalFrame$biasBpdae, c(0.025, 0.975))

      AvgIterTime <- mean(finalFrame$iterTime)
      SDIterTime <- sd(finalFrame$iterTime)
      IterTimeData <- cbind(AvgIterTime,SDIterTime)
      colnames(IterTimeData) <- c("Mean","SD")
      write.csv(IterTimeData,paste0("MNARIterTime",propName,".csv"))

      names <- c("MiceTreat","MiceBp","oracleTreat","oracleBp","complete-caseTreat","complete-caseBp","daetreat","daeBp")
      bias <- c(AvgbiasTreatMICE,AvgbiasBpMICE,AvgbiasTreatComplete,AvgbiasBpComplete,AvgbiasTreatExclude,AvgbiasBpExclude,AvgbiasTreatdae,AvgbiasBpdae)
      Percentile <- c(PercentileTreatMICE,PercentileBpMICE,PercentileTreatComplete,PercentileBpComplete,PercentileTreatExclude,PercentileBpExclude,PercentileTreatDAE,PercentileBpDAE)
      Perc <- t(as.data.frame(unlist(Percentile)))
      namesPerc <- c("MiceTreatLower","MiceTreatUpper","MiceBpLower","MiceBpUpper","oracleTreatLower","oracleTreatUpper","oracleBpLower","oracleBpUpper","complete-caseTreatLower","complete-caseTreatUpper","complete-caseBpLower","complete-caseBpUpper","daetreatLower","daetreatUpper","daeBpLower","daeBpUpper")
      colnames(Perc) <- namesPerc
      rownames(Perc) <- c()
      Extra <- cbind(Perc)
      results <- data.frame(rbind(bias,propName))
      colnames(results) <- names
      test_that("SimResultsGood", {
        expect_equal(sum(is.na(results)),0)
      })
      test_that("SimExtraGood", {
        expect_equal(sum(is.na(Extra)),0)
      })

      write.csv(results,paste0("MNAR",propName,".csv"))
      write.csv(Extra,paste0("MNARExtra",propName,".csv"))

    }

    wholeTime <-Sys.time() -startWholetime
    print(wholeTime)
  }
}
