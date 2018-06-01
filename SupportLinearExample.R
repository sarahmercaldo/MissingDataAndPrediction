library(plyr)
library(rms)
library(mice)
library(cvTools)
library(caret)

#################################################
# This is the support dataset example
###############################################

source('/Users/SFM/Dropbox/GraduateSchool/Dissertation/MissingDataWithPredictionModels/Functions/PMKSfunctions.R')
dat <- read.csv(file = "/Users/SFM/Dropbox/GraduateSchool/Dissertation/MissingDataWithPredictionModels/support2/support2.csv")
dd <- datadist(dat)
options(datadist='dd')

#either include dzclass or dzgroup
dat <- subset(dat, select = c(sps,
                              pafi,
                              meanbp,
                              wblc,
                              alb,
                              resp,
                              temp,
                              hrt,
                              bili,
                              crea,
                              sod))


#remove the one person with a missing outcome
dat <- dat[-which(is.na(dat$sps)),]

#remove the few people with only one missing 
dat <- dat[-which(is.na(dat$meanbp)),]
#dat <- dat[-which(is.na(dat$resp)),]
#dat <- dat[-which(is.na(dat$temp)),]
#dat <- dat[-which(is.na(dat$hrt)),]
#dat <- dat[-which(is.na(dat$sod)),]


dat.I <- create.miss.ind(dat) 
dd <- datadist(dat.I)
options(datadist='dd')

################################
#Fit Model -- look at indicators to see if any are significant
################################

imp.alldat <- aregImpute( ~ sps  +
                       pafi + 
                       meanbp + 
                       wblc + 
                       alb +
                       resp + 
                       temp + 
                       hrt + 
                       bili + 
                       crea + 
                       sod , data=dat.I, n.impute=10,pr=FALSE)
mod.MIMI.alldat <- fit.mult.impute(sps ~  
                                     pafi + 
                                     meanbp + 
                                     wblc + 
                                     alb +
                                     resp + 
                                     temp + 
                                     hrt + 
                                     bili + 
                                     crea + 
                                     sod +
                                     m.pafi + 
                                     m.wblc +
                                     m.alb +
                                     m.bili +
                                     m.crea, ols,imp.alldat,data=dat.I,pr=FALSE)


summary(mod.MIMI.alldat)

mod.MIMI.alldat <- with(data = imps, glm(sps ~  
                                    pafi + 
                                    meanbp + 
                                    wblc + 
                                    alb +
                                    resp + 
                                    temp + 
                                    hrt + 
                                    bili + 
                                    crea + 
                                    sod +
                                    m.pafi + 
                                    m.wblc +
                                    m.alb +
                                    m.bili +
                                    m.crea,
                                  family='gaussian'))

pooled.mod.MIMI.alldat <- mice::pool(mod.MIMI.alldat)

summary(pooled.mod.MIMI.alldat)



################################
#Cross Validate prediction model
################################


nimp = 10  #number of imputations  
k=10
results <- data.frame(brier.MI.y      = rep(NA, k),
                      brier.MI.noy    = rep(NA, k),
                      brier.MIMI.y    = rep(NA, k),
                      brier.MIMI.noy  = rep(NA, k)   
)

#results.pattern <- lapply(vector('list',k),function(ZZ) vector('list',3))
results.pattern <- vector('list',length = k)


flds <- createFolds(dat$sps, k=k, list=TRUE, returnTrain=FALSE)

for( i in 1:k){
  trainingset <- dat.I[-flds[[i]],]
  testset <- dat.I[flds[[i]],]
  
  
  imps  <- mice(trainingset, m=nimp, printFlag=FALSE,
                seed = 1, diagnostics = FALSE)
  
  mod <- with(data = imps,
              exp = glm(sps ~
                          pafi + 
                          meanbp + 
                          wblc + 
                          alb +
                          resp + 
                          temp + 
                          hrt + 
                          bili + 
                          crea + 
                          sod,
                        family='gaussian'))
  
  pooled.mod <- mice::pool(mod)
  
  mod.MIMI <- with(data = imps, glm(sps ~  
                                      pafi + 
                                      meanbp + 
                                      wblc + 
                                      alb +
                                      resp + 
                                      temp + 
                                      hrt + 
                                      bili + 
                                      crea + 
                                      sod +
                                      m.pafi + 
                                      m.wblc +
                                      m.alb +
                                      m.bili +
                                      m.crea,
                                    family='gaussian'))
  
  pooled.mod.MIMI <- mice::pool(mod.MIMI)
  
  
  ############################
  # Start Using the test data
  ############################
  imputed.datasets.withy <- replicate(nimp, testset, simplify=FALSE)
  imputed.datasets.noy   <- replicate(nimp, testset, simplify=FALSE)
  
  for(cc in which(!complete.cases(testset))){
    addNewPatient <- rbind.match.columns(trainingset, testset[cc,])
    Validateimps_y  <- mice(with(addNewPatient,data.frame(sps,
                                                          pafi,
                                                          meanbp,
                                                          wblc,
                                                          alb,
                                                          resp,
                                                          temp,
                                                          hrt,
                                                          bili,
                                                          crea,
                                                          sod)),method = 'norm.predict', m=nimp, maxit = 0)
    
    
    #Impute Validation Set without Y
    Validateimps_noy  <- mice(with(addNewPatient,data.frame( pafi,
                                                             meanbp,
                                                             wblc,
                                                             alb,
                                                             resp,
                                                             temp,
                                                             hrt,
                                                             bili,
                                                             crea,
                                                             sod)), m=nimp,method = 'norm.predict',
                              maxit = 0)
    
    
    for (j in 1:nimp){
      
      tmp.imp.y <- mice::complete(Validateimps_y)[nrow(addNewPatient),]
      
      imputed.datasets.withy[[j]][cc,names(tmp.imp.y) ] <- tmp.imp.y
      
      tmp.imp.noy <- mice::complete(Validateimps_noy)[nrow(addNewPatient),]
      
      imputed.datasets.noy[[j]][cc,names(tmp.imp.noy)] <-   tmp.imp.noy				
    }	
  }
  
  
  pred.withy      <- lapply(imputed.datasets.withy, 
                            function(z) predict.oos.mice(pooled.mod, NEWDATA=z))
  
  pred.noy        <- lapply(imputed.datasets.noy, 
                            function(z) predict.oos.mice(pooled.mod, NEWDATA=z))   
  
  pred.MIMI.withy <- lapply(imputed.datasets.withy, 
                            function(z) predict.oos.mice(pooled.mod.MIMI, NEWDATA=z))
  
  pred.MIMI.noy   <- lapply(imputed.datasets.noy, 
                            function(z) predict.oos.mice(pooled.mod.MIMI, NEWDATA=z))  
  
  
  #Combine prediction estimates			
  pred.withy.avg   	  <- Reduce("+", pred.withy) / length(pred.withy)
  pred.noy.avg	      <- Reduce("+", pred.noy) / length(pred.noy)
  pred.MIMI.withy.avg <- Reduce("+", pred.MIMI.withy) / length(pred.MIMI.withy)
  pred.MIMI.noy.avg	  <- Reduce("+", pred.MIMI.noy) / length(pred.MIMI.noy)
  
  
  #######################
  # Get Brier Score Data
  #######################
  
  results[i,'brier.MI.y']      <-  brier.score(pred.withy.avg,
                                               testset$sps)
  results[i,'brier.MI.noy']    <-  brier.score(pred.noy.avg,testset$sps)
  results[i,'brier.MIMI.y']    <-  brier.score(pred.MIMI.withy.avg,
                                               testset$sps)
  results[i,'brier.MIMI.noy']  <-  brier.score(pred.MIMI.noy.avg,
                                               testset$sps)
  
  ex.fit <- "sps ~ pafi + meanbp + wblc + alb + resp +  temp + 
  hrt +  bili + crea + sod"
  pattern <- which.pattern(testset, model=ex.fit)
  
  ##pred.withy.avg
  MI.withy.avg.df <- data.frame(pred.withy.avg, testset$sps, pattern)
  metric.MI.withy <- get_metric_linear(DATA=MI.withy.avg.df)
  
  ##pred.noy.avg
  MI.noy.avg.df <- data.frame(pred.noy.avg, testset$sps, pattern)
  metric.MI.noy <- get_metric_linear(DATA=MI.noy.avg.df)
  
  ##pred.noy.avg
  MIMI.withy.avg.df <- data.frame(pred.MIMI.withy.avg, testset$sps, pattern)
  metric.MIMI.withy <- get_metric_linear(DATA=MIMI.withy.avg.df)
  
  ##pred.noy.avg
  MIMI.noy.avg.df <- data.frame(pred.MIMI.noy.avg, testset$sps, pattern)
  metric.MIMI.noy <- get_metric_linear(DATA=MIMI.noy.avg.df)
  
  
  results.pattern[[i]] <- list(
    metric.MI.withy = metric.MI.withy,
    metric.MI.noy = metric.MI.noy,
    metric.MIMI.withy = metric.MIMI.withy,
    metric.MIMI.noy = metric.MIMI.noy)
  
}

save(results, file = "/Users/SFM/Desktop/supportResultsLinearSM.Rda")
save(results.pattern, file = "/Users/SFM/Desktop/supportResultsPatternLinearSM.Rda")

#####################################################
# Cross Validate within each missing data pattern
####################################################
set.seed(3)
k=10


results.pmks <- vector('list',length = k)
results.ccsm <- vector('list',length = k)



ex.fit <- "sps ~ pafi + meanbp + wblc + alb + resp +  temp + hrt +  bili + crea + sod"

flds <- createFolds(dat$sps, k=k, list=TRUE, returnTrain=FALSE)

for( i in 1:k){
  
  trainingset <- dat[-flds[[i]],]
  testset <- dat[flds[[i]],]
  
  ex.ccsm <- ccsm(DATA=trainingset, model=ex.fit, logistic=FALSE)
  ex.pmks <- pmks(DATA=trainingset, model=ex.fit, logistic=FALSE)
  
  pred.pmks <- tryCatch(predict.sm(prediction.data=testset, model=ex.fit, pmks.object=ex.pmks, logistic = FALSE),
                        error = function(e) {NA})
  if(is.na(pred.pmks)) pred.pmks.res <- NA
  if(!is.na(pred.pmks)) pred.pmks.res <- do.call(rbind,lapply(pred.pmks, function(z) data.frame(pattern = z$numeric.pattern, brier = z$strat.brier, prop.pattern = (length(z$truth)/nrow(testset)))))
  
  pred.ccsm <- tryCatch(predict.sm(prediction.data=testset, model=ex.fit, pmks.object=ex.ccsm, logistic = FALSE),
                        error = function(e) {NA})
  if(is.na(pred.ccsm)) pred.ccsm.res <- NA
  if(!is.na(pred.ccsm)) pred.ccsm.res <- do.call(rbind,lapply(pred.ccsm, function(z) data.frame(pattern = z$numeric.pattern,brier = z$strat.brier,  prop.pattern = (length(z$truth)/nrow(testset)))))
  
  results.pmks[[i]] <- pred.pmks.res
  results.ccsm[[i]] <- pred.ccsm.res
  
}

save(results.pmks, file = "/Users/SFM/Desktop/supportResultsPMKSSM.Rda")
save(results.ccsm, file = "/Users/SFM/Desktop/supportResultsPatternCCSMSM.Rda")


##############
# Plot
###############
#Average results of non pattern metrics
results.avg <- colMeans(results)

#Average results of pattern metrics
results.pattern.tmp <- lapply(results.pattern, function(zz) lapply(zz, function(q) data.frame(q,pattern = rownames(q))))
results.pattern.MI.withy <- NULL
results.pattern.MI.noy <- NULL
results.pattern.MIMI.withy <- NULL
results.pattern.MIMI.noy <- NULL

for(i in 1:k){
  results.pattern.MI.withy <- rbind(results.pattern.MI.withy,as.data.frame(results.pattern.tmp[[i]]['metric.MI.withy']))
  results.pattern.MI.noy <- rbind(results.pattern.MI.noy,as.data.frame(results.pattern.tmp[[i]]['metric.MI.noy']))
  results.pattern.MIMI.withy <- rbind(results.pattern.MIMI.withy,as.data.frame(results.pattern.tmp[[i]]['metric.MIMI.withy']))
  results.pattern.MIMI.noy <- rbind(results.pattern.MIMI.noy,as.data.frame(results.pattern.tmp[[i]]['metric.MIMI.noy']))
}


##MI no y
tmps <- split(results.pattern.MI.noy[,1:2], results.pattern.MI.noy$metric.MI.noy.pattern)
results.pattern.MI.noy <- do.call(rbind, lapply(tmps, function(ZZ) {
  zz <- ZZ
  zz[!apply(zz,2,is.finite)] <- NA
  if(nrow(zz)>1) { 
    zz <- colMeans(zz,na.rm=TRUE)
    zz[!is.finite(zz)] <- NA
  }
  zz}))

##MI with y
tmps <- split(results.pattern.MI.withy[,1:2], results.pattern.MI.withy$metric.MI.withy.pattern)
results.pattern.MI.withy <- do.call(rbind, lapply(tmps, function(ZZ) {
  zz <- ZZ
  zz[!apply(zz,2,is.finite)] <- NA
  if(nrow(zz)>1) { 
    zz <- colMeans(zz,na.rm=TRUE)
    zz[!is.finite(zz)] <- NA
  }
  zz}))

##MIMI no y
tmps <- split(results.pattern.MIMI.noy[,1:2], results.pattern.MIMI.noy$metric.MIMI.noy.pattern)
results.pattern.MIMI.noy <- do.call(rbind, lapply(tmps, function(ZZ) {
  zz <- ZZ
  zz[!apply(zz,2,is.finite)] <- NA
  if(nrow(zz)>1) { 
    zz <- colMeans(zz,na.rm=TRUE)
    zz[!is.finite(zz)] <- NA
  }
  zz}))


##MIMI wtihy
tmps <- split(results.pattern.MIMI.withy[,1:2], results.pattern.MIMI.withy$metric.MIMI.withy.pattern)
results.pattern.MIMI.withy <- do.call(rbind, lapply(tmps, function(ZZ) {
  zz <- ZZ
  zz[!apply(zz,2,is.finite)] <- NA
  if(nrow(zz)>1) { 
    zz <- colMeans(zz,na.rm=TRUE)
    zz[!is.finite(zz)] <- NA
  }
  zz}))



##Average the results of the PMKS
tmp  <- do.call(rbind, results.pmks)
tmps <- split(tmp, tmp$pattern)
results.pmks.avg <- do.call(rbind, lapply(tmps, function(ZZ) {
  zz <- ZZ[,-1]
  zz[!apply(zz,2,is.finite)] <- NA
  if(nrow(zz)>1) { 
    zz <- colMeans(zz,na.rm=TRUE)
    zz[!is.finite(zz)] <- NA
  }
  zz}))



##Average the results of the CCSM
tmp  <- do.call(rbind, results.ccsm)
tmps <- split(tmp, tmp$pattern)
results.ccsm.avg <- data.frame(do.call(rbind, lapply(tmps, function(ZZ) {
  zz <- ZZ[,-1]
  zz[!apply(zz,2,is.finite)] <- NA
  if(nrow(zz)>1) { 
    zz <- colMeans(zz,na.rm=TRUE)
    zz[!is.finite(zz)] <- NA
  }
  zz})))

emp.pattern <- as.data.frame(table(which.pattern(DATA = dat, ex.fit)))
names(emp.pattern) <- c('Row.names', 'prop.pattern')
emp.pattern[,2] <- emp.pattern[,2]/nrow(dat)
results.pmks.avg$pmks.weight.avg = results.pmks.avg[,1]*results.pmks.avg[,2]
results.ccsm.avg$ccsm.weight.avg = results.ccsm.avg[,1]*results.ccsm.avg[,2]


results.pattern.MI.noy$MI.noy.weight.avg = results.pattern.MI.noy[,1]*results.pattern.MI.noy[,2]
results.pattern.MI.withy$MI.withy.weight.avg = results.pattern.MI.withy[,1]*results.pattern.MI.withy[,2]
results.pattern.MIMI.withy$MIMI.withy.weight.avg = results.pattern.MIMI.withy[,1]*results.pattern.MIMI.withy[,2]
results.pattern.MIMI.noy$MIMI.noy.weight.avg = results.pattern.MIMI.noy[,1]*results.pattern.MIMI.noy[,2]



m1 <- merge(results.pmks.avg,results.ccsm.avg, by = 0)
m2 <- merge(results.pattern.MIMI.noy, results.pattern.MI.noy, by = 0)
m3 <- merge(results.pattern.MI.withy,results.pattern.MIMI.withy, by = 0 )
all.results <- merge(m1,m2, by = 'Row.names')
all.results <- merge(all.results,m3, by = 'Row.names')
all.results <- merge(all.results, emp.pattern, by = 'Row.names')


plot(seq(1, length(all.results$Row.names), by=1),all.results$brier.x*all.results$prop.pattern,
     pch = NA, xlab = "Pattern", ylab = "Weighted MSE", xaxt='n', ylim = c(-2,60),
     main = 'sps ~ pafi + meanbp + wblc + alb + resp + temp + hrt + bili + crea + sod', cex.main = 1) #, xlim = c(0,5), ylim = c(0,5))


points(jitter(seq(1, length(all.results$Row.names), by=1)), all.results$metric.MI.noy.brier*all.results$prop.pattern, pch = 20, col = 'red')
points(jitter(seq(1, length(all.results$Row.names), by=1)), all.results$metric.MIMI.noy.brier*all.results$prop.pattern, pch = 20, col = 'blue')
points(jitter(seq(1, length(all.results$Row.names), by=1)), all.results$brier.y*all.results$prop.pattern, pch = 20, col = 'darkgreen')
points(jitter(seq(1, length(all.results$Row.names), by=1)), all.results$brier.x*all.results$prop.pattern, pch = 20, col = 'black')


text(x = seq(1, length(all.results$Row.names), by=1), par("usr")[3] - 2.5, 
     labels = all.results$Row.names
     , srt = 90, pos = 1, xpd = TRUE, cex = 0.5)
abline(v = c(2:length(all.results$Row.names))-0.5, lty = 2, lwd=0.5, col='gray')

legend('topright', legend = c("PMKS", "MI", "MIMI", "CCSM",round(sum(all.results$brier.x*all.results$prop.pattern),1),
                              round(sum(all.results$metric.MI.noy.brier*all.results$prop.pattern),1),
                              round(sum(all.results$metric.MIMI.noy.brier*all.results$prop.pattern),1),
                              round(sum(all.results$brier.y*all.results$prop.pattern),1)), ncol = 2, bty='n',
       title = "Total MSE", col = c('black','red','blue','darkgreen'), pch = c(20,20,20,20,NA,NA,NA,NA))

text(x = seq(1, length(all.results$Row.names), by=1), par("usr")[3] + 2.8, 
     labels = round(all.results$prop.pattern,2)
     ,pos = 1, xpd = TRUE, cex = 0.5)
text(x =-1,  par("usr")[3] + 2.8, labels = "Prop/Pattern",pos = 1, xpd = TRUE, cex = 0.5)



save(all.results, file = '/Users/SFM/Desktop/supportSM')
















#################################################################
#Exagerate the effects by inducing a greater amount of missingness
#
###############################

m = rbinom(nrow(dat), 1, 
           expit(miss.int(miss.mech=dat[,c('meanbp','sps')], betaM=2,0.9)
                 + 2*(((dat[,'sps']/sd(dat[,'sps'])) + dat[,'meanbp'])/as.numeric(sqrt(2*(1+cor(dat[,'sps'],dat[,'meanbp'])))))))

dat[m==1,'meanbp'] <- NA

dat.I <- create.miss.ind(dat) 
dd <- datadist(dat.I)
options(datadist='dd')

################################
#Fit Model -- look at indicators to see if any are significant
################################

imp.alldat <- aregImpute( ~ sps  +
                            pafi + 
                            meanbp + 
                            wblc + 
                            alb +
                            resp + 
                            temp + 
                            hrt + 
                            bili + 
                            crea + 
                            sod , data=dat.I, n.impute=10,pr=FALSE)
mod.MIMI.alldat <- fit.mult.impute(sps ~  
                                     pafi + 
                                     meanbp + 
                                     wblc + 
                                     alb +
                                     resp + 
                                     temp + 
                                     hrt + 
                                     bili + 
                                     crea + 
                                     sod +
                                     m.pafi + 
                                     m.wblc +
                                     m.alb +
                                     m.bili +
                                     m.crea + m.meanbp
                                   , ols,imp.alldat,data=dat.I,pr=FALSE)


summary(mod.MIMI.alldat)

# mod.MIMI.alldat <- with(data = imps, glm(sps ~  
#                                            pafi + 
#                                            meanbp + 
#                                            wblc + 
#                                            alb +
#                                            resp + 
#                                            temp + 
#                                            hrt + 
#                                            bili + 
#                                            crea + 
#                                            sod +
#                                            m.pafi + 
#                                            m.wblc +
#                                            m.alb +
#                                            m.bili +
#                                            m.crea,
#                                          family='gaussian'))
# 
# pooled.mod.MIMI.alldat <- mice::pool(mod.MIMI.alldat)
# 
# summary(pooled.mod.MIMI.alldat)
# 
# 

################################
#Cross Validate prediction model
################################


nimp = 10  #number of imputations  
k=10
results <- data.frame(brier.MI.y      = rep(NA, k),
                      brier.MI.noy    = rep(NA, k),
                      brier.MIMI.y    = rep(NA, k),
                      brier.MIMI.noy  = rep(NA, k)   
)

#results.pattern <- lapply(vector('list',k),function(ZZ) vector('list',3))
results.pattern <- vector('list',length = k)


flds <- createFolds(dat$sps, k=k, list=TRUE, returnTrain=FALSE)

for( i in 1:k){
  trainingset <- dat.I[-flds[[i]],]
  testset <- dat.I[flds[[i]],]
  
  
  imps  <- mice(trainingset, m=nimp, printFlag=FALSE,
                seed = 1, diagnostics = FALSE)
  
  mod <- with(data = imps,
              exp = glm(sps ~
                          pafi + 
                          meanbp + 
                          wblc + 
                          alb +
                          resp + 
                          temp + 
                          hrt + 
                          bili + 
                          crea + 
                          sod,
                        family='gaussian'))
  
  pooled.mod <- mice::pool(mod)
  
  mod.MIMI <- with(data = imps, glm(sps ~  
                                      pafi + 
                                      meanbp + 
                                      wblc + 
                                      alb +
                                      resp + 
                                      temp + 
                                      hrt + 
                                      bili + 
                                      crea + 
                                      sod +
                                      m.pafi + 
                                      m.wblc +
                                      m.alb +
                                      m.bili +
                                      m.crea + m.meanbp,
                                    family='gaussian'))
  
  pooled.mod.MIMI <- mice::pool(mod.MIMI)
  
  
  ############################
  # Start Using the test data
  ############################
  imputed.datasets.withy <- replicate(nimp, testset, simplify=FALSE)
  imputed.datasets.noy   <- replicate(nimp, testset, simplify=FALSE)
  
  for(cc in which(!complete.cases(testset))){
    addNewPatient <- rbind.match.columns(trainingset, testset[cc,])
    Validateimps_y  <- mice(with(addNewPatient,data.frame(sps,
                                                          pafi,
                                                          meanbp,
                                                          wblc,
                                                          alb,
                                                          resp,
                                                          temp,
                                                          hrt,
                                                          bili,
                                                          crea,
                                                          sod)),method = 'norm.predict', m=nimp, maxit = 0)
    
    
    #Impute Validation Set without Y
    Validateimps_noy  <- mice(with(addNewPatient,data.frame( pafi,
                                                             meanbp,
                                                             wblc,
                                                             alb,
                                                             resp,
                                                             temp,
                                                             hrt,
                                                             bili,
                                                             crea,
                                                             sod)), m=nimp,method = 'norm.predict',
                              maxit = 0)
    
    
    for (j in 1:nimp){
      
      tmp.imp.y <- mice::complete(Validateimps_y)[nrow(addNewPatient),]
      
      imputed.datasets.withy[[j]][cc,names(tmp.imp.y) ] <- tmp.imp.y
      
      tmp.imp.noy <- mice::complete(Validateimps_noy)[nrow(addNewPatient),]
      
      imputed.datasets.noy[[j]][cc,names(tmp.imp.noy)] <-   tmp.imp.noy				
    }	
  }
  
  
  pred.withy      <- lapply(imputed.datasets.withy, 
                            function(z) predict.oos.mice(pooled.mod, NEWDATA=z))
  
  pred.noy        <- lapply(imputed.datasets.noy, 
                            function(z) predict.oos.mice(pooled.mod, NEWDATA=z))   
  
  pred.MIMI.withy <- lapply(imputed.datasets.withy, 
                            function(z) predict.oos.mice(pooled.mod.MIMI, NEWDATA=z))
  
  pred.MIMI.noy   <- lapply(imputed.datasets.noy, 
                            function(z) predict.oos.mice(pooled.mod.MIMI, NEWDATA=z))  
  
  
  #Combine prediction estimates			
  pred.withy.avg   	  <- Reduce("+", pred.withy) / length(pred.withy)
  pred.noy.avg	      <- Reduce("+", pred.noy) / length(pred.noy)
  pred.MIMI.withy.avg <- Reduce("+", pred.MIMI.withy) / length(pred.MIMI.withy)
  pred.MIMI.noy.avg	  <- Reduce("+", pred.MIMI.noy) / length(pred.MIMI.noy)
  
  
  #######################
  # Get Brier Score Data
  #######################
  
  results[i,'brier.MI.y']      <-  brier.score(pred.withy.avg,
                                               testset$sps)
  results[i,'brier.MI.noy']    <-  brier.score(pred.noy.avg,testset$sps)
  results[i,'brier.MIMI.y']    <-  brier.score(pred.MIMI.withy.avg,
                                               testset$sps)
  results[i,'brier.MIMI.noy']  <-  brier.score(pred.MIMI.noy.avg,
                                               testset$sps)
  
  ex.fit <- "sps ~ pafi + meanbp + wblc + alb + resp +  temp + 
  hrt +  bili + crea + sod"
  pattern <- which.pattern(testset, model=ex.fit)
  
  ##pred.withy.avg
  MI.withy.avg.df <- data.frame(pred.withy.avg, testset$sps, pattern)
  metric.MI.withy <- get_metric_linear(DATA=MI.withy.avg.df)
  
  ##pred.noy.avg
  MI.noy.avg.df <- data.frame(pred.noy.avg, testset$sps, pattern)
  metric.MI.noy <- get_metric_linear(DATA=MI.noy.avg.df)
  
  ##pred.noy.avg
  MIMI.withy.avg.df <- data.frame(pred.MIMI.withy.avg, testset$sps, pattern)
  metric.MIMI.withy <- get_metric_linear(DATA=MIMI.withy.avg.df)
  
  ##pred.noy.avg
  MIMI.noy.avg.df <- data.frame(pred.MIMI.noy.avg, testset$sps, pattern)
  metric.MIMI.noy <- get_metric_linear(DATA=MIMI.noy.avg.df)
  
  
  results.pattern[[i]] <- list(
    metric.MI.withy = metric.MI.withy,
    metric.MI.noy = metric.MI.noy,
    metric.MIMI.withy = metric.MIMI.withy,
    metric.MIMI.noy = metric.MIMI.noy)
  
}

#####################################################
# Cross Validate within each missing data pattern
####################################################
set.seed(3)
k=20


results.pmks <- vector('list',length = k)
results.ccsm <- vector('list',length = k)



ex.fit <- "sps ~ pafi + meanbp + wblc + alb + resp +  temp + hrt +  bili + crea + sod"

flds <- createFolds(dat$sps, k=k, list=TRUE, returnTrain=FALSE)

for( i in 1:k){
  
  trainingset <- dat[-flds[[i]],]
  testset <- dat[flds[[i]],]
  
  ex.ccsm <- ccsm(DATA=trainingset, model=ex.fit, logistic=FALSE)
  ex.pmks <- pmks(DATA=trainingset, model=ex.fit, logistic=FALSE)
  
  pred.pmks <- tryCatch(predict.sm(prediction.data=testset, model=ex.fit, pmks.object=ex.pmks, logistic = FALSE),
                        error = function(e) {NA})
  if(is.na(pred.pmks)) pred.pmks.res <- NA
  if(!is.na(pred.pmks)) pred.pmks.res <- do.call(rbind,lapply(pred.pmks, function(z) data.frame(pattern = z$numeric.pattern, brier = z$strat.brier, prop.pattern = (length(z$truth)/nrow(testset)))))
  
  pred.ccsm <- tryCatch(predict.sm(prediction.data=testset, model=ex.fit, pmks.object=ex.ccsm, logistic = FALSE),
                        error = function(e) {NA})
  if(is.na(pred.ccsm)) pred.ccsm.res <- NA
  if(!is.na(pred.ccsm)) pred.ccsm.res <- do.call(rbind,lapply(pred.ccsm, function(z) data.frame(pattern = z$numeric.pattern,brier = z$strat.brier,  prop.pattern = (length(z$truth)/nrow(testset)))))
  
  results.pmks[[i]] <- pred.pmks.res
  results.ccsm[[i]] <- pred.ccsm.res
  
}

##############
# Plot
###############

#Average results of non pattern metrics
results.avg <- colMeans(results)

#Average results of pattern metrics


results.pattern.tmp <- lapply(results.pattern, function(zz) lapply(zz, function(q) data.frame(q,pattern = rownames(q))))


results.pattern.MI.withy <- NULL
results.pattern.MI.noy <- NULL
results.pattern.MIMI.withy <- NULL
results.pattern.MIMI.noy <- NULL

for(i in 1:k){
  results.pattern.MI.withy <- rbind(results.pattern.MI.withy,as.data.frame(results.pattern.tmp[[i]]['metric.MI.withy']))
  results.pattern.MI.noy <- rbind(results.pattern.MI.noy,as.data.frame(results.pattern.tmp[[i]]['metric.MI.noy']))
  results.pattern.MIMI.withy <- rbind(results.pattern.MIMI.withy,as.data.frame(results.pattern.tmp[[i]]['metric.MIMI.withy']))
  results.pattern.MIMI.noy <- rbind(results.pattern.MIMI.noy,as.data.frame(results.pattern.tmp[[i]]['metric.MIMI.noy']))
}


##MI no y
tmps <- split(results.pattern.MI.noy[,1:2], results.pattern.MI.noy$metric.MI.noy.pattern)
results.pattern.MI.noy <- do.call(rbind, lapply(tmps, function(ZZ) {
  zz <- ZZ
  zz[!apply(zz,2,is.finite)] <- NA
  if(nrow(zz)>1) { 
    zz <- colMeans(zz,na.rm=TRUE)
    zz[!is.finite(zz)] <- NA
  }
  zz}))

##MI with y
tmps <- split(results.pattern.MI.withy[,1:2], results.pattern.MI.withy$metric.MI.withy.pattern)
results.pattern.MI.withy <- do.call(rbind, lapply(tmps, function(ZZ) {
  zz <- ZZ
  zz[!apply(zz,2,is.finite)] <- NA
  if(nrow(zz)>1) { 
    zz <- colMeans(zz,na.rm=TRUE)
    zz[!is.finite(zz)] <- NA
  }
  zz}))

##MIMI no y
tmps <- split(results.pattern.MIMI.noy[,1:2], results.pattern.MIMI.noy$metric.MIMI.noy.pattern)
results.pattern.MIMI.noy <- do.call(rbind, lapply(tmps, function(ZZ) {
  zz <- ZZ
  zz[!apply(zz,2,is.finite)] <- NA
  if(nrow(zz)>1) { 
    zz <- colMeans(zz,na.rm=TRUE)
    zz[!is.finite(zz)] <- NA
  }
  zz}))


##MIMI wtihy
tmps <- split(results.pattern.MIMI.withy[,1:2], results.pattern.MIMI.withy$metric.MIMI.withy.pattern)
results.pattern.MIMI.withy <- do.call(rbind, lapply(tmps, function(ZZ) {
  zz <- ZZ
  zz[!apply(zz,2,is.finite)] <- NA
  if(nrow(zz)>1) { 
    zz <- colMeans(zz,na.rm=TRUE)
    zz[!is.finite(zz)] <- NA
  }
  zz}))



##Average the results of the PMKS
tmp  <- do.call(rbind, results.pmks)
tmps <- split(tmp, tmp$pattern)
results.pmks.avg <- do.call(rbind, lapply(tmps, function(ZZ) {
  zz <- ZZ[,-1]
  zz[!apply(zz,2,is.finite)] <- NA
  if(nrow(zz)>1) { 
    zz <- colMeans(zz,na.rm=TRUE)
    zz[!is.finite(zz)] <- NA
  }
  zz}))



##Average the results of the CCSM
tmp  <- do.call(rbind, results.ccsm)
tmps <- split(tmp, tmp$pattern)
results.ccsm.avg <- data.frame(do.call(rbind, lapply(tmps, function(ZZ) {
  zz <- ZZ[,-1]
  zz[!apply(zz,2,is.finite)] <- NA
  if(nrow(zz)>1) { 
    zz <- colMeans(zz,na.rm=TRUE)
    zz[!is.finite(zz)] <- NA
  }
  zz})))

emp.pattern <- as.data.frame(table(which.pattern(DATA = dat, ex.fit)))
names(emp.pattern) <- c('Row.names', 'prop.pattern')
emp.pattern[,2] <- emp.pattern[,2]/nrow(dat)
results.pmks.avg$pmks.weight.avg = results.pmks.avg[,1]*results.pmks.avg[,2]
results.ccsm.avg$ccsm.weight.avg = results.ccsm.avg[,1]*results.ccsm.avg[,2]


results.pattern.MI.noy$MI.noy.weight.avg = results.pattern.MI.noy[,1]*results.pattern.MI.noy[,2]
results.pattern.MI.withy$MI.withy.weight.avg = results.pattern.MI.withy[,1]*results.pattern.MI.withy[,2]
results.pattern.MIMI.withy$MIMI.withy.weight.avg = results.pattern.MIMI.withy[,1]*results.pattern.MIMI.withy[,2]
results.pattern.MIMI.noy$MIMI.noy.weight.avg = results.pattern.MIMI.noy[,1]*results.pattern.MIMI.noy[,2]



m1 <- merge(results.pmks.avg,results.ccsm.avg, by = 0)
m2 <- merge(results.pattern.MIMI.noy, results.pattern.MI.noy, by = 0)
m3 <- merge(results.pattern.MI.withy,results.pattern.MIMI.withy, by = 0 )
all.results <- merge(m1,m2, by = 'Row.names')
all.results <- merge(all.results,m3, by = 'Row.names')
all.results <- merge(all.results, emp.pattern, by = 'Row.names')


plot(all.results$brier.x*all.results$prop.pattern, all.results$metric.MI.noy.brier*all.results$prop.pattern,
     pch = 20, xlab = "PMKS Weighted MSE", ylab = "MI Weighted MSE") #, xlim = c(0,5), ylim = c(0,5))
abline(0,1)
legend('bottomright', legend = c("PMKS", "MI", "MIMI", "CCSM",round(sum(all.results$brier.x*all.results$prop.pattern),2),
                                 round(sum(all.results$metric.MI.noy.brier*all.results$prop.pattern),2),
                                 round(sum(all.results$metric.MIMI.noy.brier*all.results$prop.pattern),2),
                                 round(sum(all.results$brier.y*all.results$prop.pattern),2)), ncol = 2, bty='n',
       title = "Weighted MSE", col = c('black','red','blue','darkgreen'), pch = c(1,1,1,1,NA,NA,NA,NA))
points(all.results$brier.x*all.results$prop.pattern, all.results$metric.MI.noy.brier*all.results$prop.pattern, pch = 1, col = 'red')
points(all.results$brier.x*all.results$prop.pattern, all.results$metric.MIMI.noy.brier*all.results$prop.pattern, pch = 1, col = 'blue')
points(all.results$brier.x*all.results$prop.pattern, all.results$brier.y*all.results$prop.pattern, pch = 1, col = 'darkgreen')


plot(seq(1, length(all.results$Row.names), by=1),all.results$brier.x*all.results$prop.pattern,
     pch = 1, xlab = "Pattern", ylab = "Weighted MSE", xaxt='n') #, xlim = c(0,5), ylim = c(0,5))


points(seq(1, length(all.results$Row.names), by=1), all.results$metric.MI.noy.brier*all.results$prop.pattern, pch = 1, col = 'red')
points(seq(1, length(all.results$Row.names), by=1), all.results$metric.MIMI.noy.brier*all.results$prop.pattern, pch = 1, col = 'blue')
points(seq(1, length(all.results$Row.names), by=1), all.results$brier.y*all.results$prop.pattern, pch = 1, col = 'darkgreen')
text(x = seq(1, length(all.results$Row.names), by=1), par("usr")[3] - 0.6, 
     labels = all.results$Row.names
     , srt = 90, pos = 1, xpd = TRUE, cex = 0.5)
abline(v = c(2:length(all.results$Row.names))-0.5, lty = 2, col='gray')

legend('topright', legend = c("PMKS", "MI", "MIMI", "CCSM",round(sum(all.results$brier.x*all.results$prop.pattern),1),
                              round(sum(all.results$metric.MI.noy.brier*all.results$prop.pattern),1),
                              round(sum(all.results$metric.MIMI.noy.brier*all.results$prop.pattern),1),
                              round(sum(all.results$brier.y*all.results$prop.pattern),1)), ncol = 2, bty='n',
       title = "Weighted MSE", col = c('black','red','blue','darkgreen'), pch = c(1,1,1,1,NA,NA,NA,NA))


















################################################################
#Exagerate the effects by inducing a greater amount of missingness
# ADD 25 to the intercept
###############################

dat <- read.csv(file = "/Users/SFM/Dropbox/GraduateSchool/Dissertation/Paper1-MissingDataWithPredictionModels/support2/support2.csv")
dd <- datadist(dat)
options(datadist='dd')

#either include dzclass or dzgroup
dat <- subset(dat, select = c(sps,
                              pafi,
                              meanbp,
                              wblc,
                              alb,
                              resp,
                              temp,
                              hrt,
                              bili,
                              crea,
                              sod))


#remove the one person with a missing outcome
dat <- dat[-which(is.na(dat$sps)),]

#remove the few people with only one missing 
dat <- dat[-which(is.na(dat$meanbp)),]
#dat <- dat[-which(is.na(dat$resp)),]
#dat <- dat[-which(is.na(dat$temp)),]
#dat <- dat[-which(is.na(dat$hrt)),]
#dat <- dat[-which(is.na(dat$sod)),]


dat[which(is.na(dat[,'pafi'])),'sps'] <- dat[which(is.na(dat[,'pafi'])),'sps'] + 25


dat.I <- create.miss.ind(dat) 
dd <- datadist(dat.I)
options(datadist='dd')

################################
#Fit Model -- look at indicators to see if any are significant
################################

imp.alldat <- aregImpute( ~ sps  +
                            pafi + 
                            meanbp + 
                            wblc + 
                            alb +
                            resp + 
                            temp + 
                            hrt + 
                            bili + 
                            crea + 
                            sod , data=dat.I, n.impute=10,pr=FALSE)
mod.MIMI.alldat <- fit.mult.impute(sps ~  
                                     pafi + 
                                     meanbp + 
                                     wblc + 
                                     alb +
                                     resp + 
                                     temp + 
                                     hrt + 
                                     bili + 
                                     crea + 
                                     sod +
                                     m.pafi + 
                                     m.wblc +
                                     m.alb +
                                     m.bili +
                                     m.crea
                                   , ols,imp.alldat,data=dat.I,pr=FALSE)


summary(mod.MIMI.alldat)

# mod.MIMI.alldat <- with(data = imps, glm(sps ~  
#                                            pafi + 
#                                            meanbp + 
#                                            wblc + 
#                                            alb +
#                                            resp + 
#                                            temp + 
#                                            hrt + 
#                                            bili + 
#                                            crea + 
#                                            sod +
#                                            m.pafi + 
#                                            m.wblc +
#                                            m.alb +
#                                            m.bili +
#                                            m.crea,
#                                          family='gaussian'))
# 
# pooled.mod.MIMI.alldat <- mice::pool(mod.MIMI.alldat)
# 
# summary(pooled.mod.MIMI.alldat)
# 
# 

################################
#Cross Validate prediction model
################################


nimp = 5  #number of imputations  
k=10
results <- data.frame(brier.MI.y      = rep(NA, k),
                      brier.MI.noy    = rep(NA, k),
                      brier.MIMI.y    = rep(NA, k),
                      brier.MIMI.noy  = rep(NA, k)   
)

#results.pattern <- lapply(vector('list',k),function(ZZ) vector('list',3))
results.pattern <- vector('list',length = k)


flds <- createFolds(dat$sps, k=k, list=TRUE, returnTrain=FALSE)

for( i in 1:k){
  trainingset <- dat.I[-flds[[i]],]
  testset <- dat.I[flds[[i]],]
  
  
  imps  <- mice(trainingset, m=nimp, printFlag=FALSE,
                seed = 1, diagnostics = FALSE)
  
  mod <- with(data = imps,
              exp = glm(sps ~
                          pafi + 
                          meanbp + 
                          wblc + 
                          alb +
                          resp + 
                          temp + 
                          hrt + 
                          bili + 
                          crea + 
                          sod,
                        family='gaussian'))
  
  pooled.mod <- mice::pool(mod)
  
  mod.MIMI <- with(data = imps, glm(sps ~  
                                      pafi + 
                                      meanbp + 
                                      wblc + 
                                      alb +
                                      resp + 
                                      temp + 
                                      hrt + 
                                      bili + 
                                      crea + 
                                      sod +
                                      m.pafi + 
                                      m.wblc +
                                      m.alb +
                                      m.bili +
                                      m.crea ,
                                    family='gaussian'))
  
  pooled.mod.MIMI <- mice::pool(mod.MIMI)
  
  
  ############################
  # Start Using the test data
  ############################
  n.imp.oos <- 1
  imputed.datasets.withy <- replicate( n.imp.oos, testset, simplify=FALSE)
  imputed.datasets.noy   <- replicate( n.imp.oos, testset, simplify=FALSE)
  
  for(cc in which(!complete.cases(testset))){
    addNewPatient <- rbind.match.columns(trainingset, testset[cc,])
    Validateimps_y  <- mice(with(addNewPatient,data.frame(sps,
                                                          pafi,
                                                          meanbp,
                                                          wblc,
                                                          alb,
                                                          resp,
                                                          temp,
                                                          hrt,
                                                          bili,
                                                          crea,
                                                          sod)),method = 'norm.predict', 
                            m= n.imp.oos, maxit = 0)
    
    
    #Impute Validation Set without Y
    Validateimps_noy  <- mice(with(addNewPatient,data.frame( pafi,
                                                             meanbp,
                                                             wblc,
                                                             alb,
                                                             resp,
                                                             temp,
                                                             hrt,
                                                             bili,
                                                             crea,
                                                             sod)), m= n.imp.oos,method = 'norm.predict',
                              maxit = 0)
    
    
    for (j in 1: n.imp.oos){
      
      tmp.imp.y <- mice::complete(Validateimps_y)[nrow(addNewPatient),]
      
      imputed.datasets.withy[[j]][cc,names(tmp.imp.y) ] <- tmp.imp.y
      
      tmp.imp.noy <- mice::complete(Validateimps_noy)[nrow(addNewPatient),]
      
      imputed.datasets.noy[[j]][cc,names(tmp.imp.noy)] <-   tmp.imp.noy				
    }	
  }
  
  
  pred.withy      <- lapply(imputed.datasets.withy, 
                            function(z) predict.oos.mice(pooled.mod, NEWDATA=z))
  
  pred.noy        <- lapply(imputed.datasets.noy, 
                            function(z) predict.oos.mice(pooled.mod, NEWDATA=z))   
  
  pred.MIMI.withy <- lapply(imputed.datasets.withy, 
                            function(z) predict.oos.mice(pooled.mod.MIMI, NEWDATA=z))
  
  pred.MIMI.noy   <- lapply(imputed.datasets.noy, 
                            function(z) predict.oos.mice(pooled.mod.MIMI, NEWDATA=z))  
  
  
  #Combine prediction estimates			
  pred.withy.avg   	  <- Reduce("+", pred.withy) / length(pred.withy)
  pred.noy.avg	      <- Reduce("+", pred.noy) / length(pred.noy)
  pred.MIMI.withy.avg <- Reduce("+", pred.MIMI.withy) / length(pred.MIMI.withy)
  pred.MIMI.noy.avg	  <- Reduce("+", pred.MIMI.noy) / length(pred.MIMI.noy)
  
  
  #######################
  # Get Brier Score Data
  #######################
  
  results[i,'brier.MI.y']      <-  brier.score(pred.withy.avg,
                                               testset$sps)
  results[i,'brier.MI.noy']    <-  brier.score(pred.noy.avg,testset$sps)
  results[i,'brier.MIMI.y']    <-  brier.score(pred.MIMI.withy.avg,
                                               testset$sps)
  results[i,'brier.MIMI.noy']  <-  brier.score(pred.MIMI.noy.avg,
                                               testset$sps)
  
  ex.fit <- "sps ~ pafi + meanbp + wblc + alb + resp +  temp + 
  hrt +  bili + crea + sod"
  pattern <- which.pattern(testset, model=ex.fit)
  
  ##pred.withy.avg
  MI.withy.avg.df <- data.frame(pred.withy.avg, testset$sps, pattern)
  metric.MI.withy <- get_metric_linear(DATA=MI.withy.avg.df)
  
  ##pred.noy.avg
  MI.noy.avg.df <- data.frame(pred.noy.avg, testset$sps, pattern)
  metric.MI.noy <- get_metric_linear(DATA=MI.noy.avg.df)
  
  ##pred.noy.avg
  MIMI.withy.avg.df <- data.frame(pred.MIMI.withy.avg, testset$sps, pattern)
  metric.MIMI.withy <- get_metric_linear(DATA=MIMI.withy.avg.df)
  
  ##pred.noy.avg
  MIMI.noy.avg.df <- data.frame(pred.MIMI.noy.avg, testset$sps, pattern)
  metric.MIMI.noy <- get_metric_linear(DATA=MIMI.noy.avg.df)
  
  
  results.pattern[[i]] <- list(
    metric.MI.withy = metric.MI.withy,
    metric.MI.noy = metric.MI.noy,
    metric.MIMI.withy = metric.MIMI.withy,
    metric.MIMI.noy = metric.MIMI.noy)
  
}

#####################################################
# Cross Validate within each missing data pattern
####################################################
set.seed(5)
k=10


results.pmks <- vector('list',length = k)
results.ccsm <- vector('list',length = k)



ex.fit <- "sps ~ pafi + meanbp + wblc + alb + resp +  temp + hrt +  bili + crea + sod"

flds <- createFolds(dat$sps, k=k, list=TRUE, returnTrain=FALSE)

for( i in 1:k){
  
  trainingset <- dat[-flds[[i]],]
  testset <- dat[flds[[i]],]
  
  ex.ccsm <- ccsm(DATA=trainingset, model=ex.fit, logistic=FALSE)
  ex.pmks <- pmks(DATA=trainingset, model=ex.fit, logistic=FALSE)
  
  pred.pmks <- tryCatch(predict.sm(prediction.data=testset, model=ex.fit, pmks.object=ex.pmks, logistic = FALSE),
                        error = function(e) {NA})
  if(is.na(pred.pmks)) pred.pmks.res <- NA
  if(!is.na(pred.pmks)) pred.pmks.res <- do.call(rbind,lapply(pred.pmks, function(z) data.frame(pattern = z$numeric.pattern, brier = z$strat.brier, prop.pattern = (length(z$truth)/nrow(testset)))))
  
  pred.ccsm <- tryCatch(predict.sm(prediction.data=testset, model=ex.fit, pmks.object=ex.ccsm, logistic = FALSE),
                        error = function(e) {NA})
  if(is.na(pred.ccsm)) pred.ccsm.res <- NA
  if(!is.na(pred.ccsm)) pred.ccsm.res <- do.call(rbind,lapply(pred.ccsm, function(z) data.frame(pattern = z$numeric.pattern,brier = z$strat.brier,  prop.pattern = (length(z$truth)/nrow(testset)))))
  
  results.pmks[[i]] <- pred.pmks.res
  results.ccsm[[i]] <- pred.ccsm.res
  
}

##############
# Plot
###############

#Average results of non pattern metrics
results.avg <- colMeans(results)

#Average results of pattern metrics


results.pattern.tmp <- lapply(results.pattern, function(zz) lapply(zz, function(q) data.frame(q,pattern = rownames(q))))


results.pattern.MI.withy <- NULL
results.pattern.MI.noy <- NULL
results.pattern.MIMI.withy <- NULL
results.pattern.MIMI.noy <- NULL

for(i in 1:k){
  results.pattern.MI.withy <- rbind(results.pattern.MI.withy,as.data.frame(results.pattern.tmp[[i]]['metric.MI.withy']))
  results.pattern.MI.noy <- rbind(results.pattern.MI.noy,as.data.frame(results.pattern.tmp[[i]]['metric.MI.noy']))
  results.pattern.MIMI.withy <- rbind(results.pattern.MIMI.withy,as.data.frame(results.pattern.tmp[[i]]['metric.MIMI.withy']))
  results.pattern.MIMI.noy <- rbind(results.pattern.MIMI.noy,as.data.frame(results.pattern.tmp[[i]]['metric.MIMI.noy']))
}


##MI no y
tmps <- split(results.pattern.MI.noy[,1:2], results.pattern.MI.noy$metric.MI.noy.pattern)
results.pattern.MI.noy <- do.call(rbind, lapply(tmps, function(ZZ) {
  zz <- ZZ
  zz[!apply(zz,2,is.finite)] <- NA
  if(nrow(zz)>1) { 
    zz <- colMeans(zz,na.rm=TRUE)
    zz[!is.finite(zz)] <- NA
  }
  zz}))

##MI with y
tmps <- split(results.pattern.MI.withy[,1:2], results.pattern.MI.withy$metric.MI.withy.pattern)
results.pattern.MI.withy <- do.call(rbind, lapply(tmps, function(ZZ) {
  zz <- ZZ
  zz[!apply(zz,2,is.finite)] <- NA
  if(nrow(zz)>1) { 
    zz <- colMeans(zz,na.rm=TRUE)
    zz[!is.finite(zz)] <- NA
  }
  zz}))

##MIMI no y
tmps <- split(results.pattern.MIMI.noy[,1:2], results.pattern.MIMI.noy$metric.MIMI.noy.pattern)
results.pattern.MIMI.noy <- do.call(rbind, lapply(tmps, function(ZZ) {
  zz <- ZZ
  zz[!apply(zz,2,is.finite)] <- NA
  if(nrow(zz)>1) { 
    zz <- colMeans(zz,na.rm=TRUE)
    zz[!is.finite(zz)] <- NA
  }
  zz}))


##MIMI wtihy
tmps <- split(results.pattern.MIMI.withy[,1:2], results.pattern.MIMI.withy$metric.MIMI.withy.pattern)
results.pattern.MIMI.withy <- do.call(rbind, lapply(tmps, function(ZZ) {
  zz <- ZZ
  zz[!apply(zz,2,is.finite)] <- NA
  if(nrow(zz)>1) { 
    zz <- colMeans(zz,na.rm=TRUE)
    zz[!is.finite(zz)] <- NA
  }
  zz}))



##Average the results of the PMKS
tmp  <- do.call(rbind, results.pmks)
tmps <- split(tmp, tmp$pattern)
results.pmks.avg <- do.call(rbind, lapply(tmps, function(ZZ) {
  zz <- ZZ[,-1]
  zz[!apply(zz,2,is.finite)] <- NA
  if(nrow(zz)>1) { 
    zz <- colMeans(zz,na.rm=TRUE)
    zz[!is.finite(zz)] <- NA
  }
  zz}))



##Average the results of the CCSM
tmp  <- do.call(rbind, results.ccsm)
tmps <- split(tmp, tmp$pattern)
results.ccsm.avg <- data.frame(do.call(rbind, lapply(tmps, function(ZZ) {
  zz <- ZZ[,-1]
  zz[!apply(zz,2,is.finite)] <- NA
  if(nrow(zz)>1) { 
    zz <- colMeans(zz,na.rm=TRUE)
    zz[!is.finite(zz)] <- NA
  }
  zz})))

emp.pattern <- as.data.frame(table(which.pattern(DATA = dat, ex.fit)))
names(emp.pattern) <- c('Row.names', 'prop.pattern')
emp.pattern[,2] <- emp.pattern[,2]/nrow(dat)
results.pmks.avg$pmks.weight.avg = results.pmks.avg[,1]*results.pmks.avg[,2]
results.ccsm.avg$ccsm.weight.avg = results.ccsm.avg[,1]*results.ccsm.avg[,2]


results.pattern.MI.noy$MI.noy.weight.avg = results.pattern.MI.noy[,1]*results.pattern.MI.noy[,2]
results.pattern.MI.withy$MI.withy.weight.avg = results.pattern.MI.withy[,1]*results.pattern.MI.withy[,2]
results.pattern.MIMI.withy$MIMI.withy.weight.avg = results.pattern.MIMI.withy[,1]*results.pattern.MIMI.withy[,2]
results.pattern.MIMI.noy$MIMI.noy.weight.avg = results.pattern.MIMI.noy[,1]*results.pattern.MIMI.noy[,2]



m1 <- merge(results.pmks.avg,results.ccsm.avg, by = 0)
m2 <- merge(results.pattern.MIMI.noy, results.pattern.MI.noy, by = 0)
m3 <- merge(results.pattern.MI.withy,results.pattern.MIMI.withy, by = 0 )
all.results <- merge(m1,m2, by = 'Row.names')
all.results <- merge(all.results,m3, by = 'Row.names')
all.results <- merge(all.results, emp.pattern, by = 'Row.names')


plot(all.results$brier.x*all.results$prop.pattern, all.results$metric.MI.noy.brier*all.results$prop.pattern,
     pch = 20, xlab = "PMKS Weighted MSE", ylab = "MI Weighted MSE") #, xlim = c(0,5), ylim = c(0,5))
abline(0,1)
legend('bottomright', legend = c("PMKS", "MI", "MIMI", "CCSM",round(sum(all.results$brier.x*all.results$prop.pattern),2),
                                 round(sum(all.results$metric.MI.noy.brier*all.results$prop.pattern),2),
                                 round(sum(all.results$metric.MIMI.noy.brier*all.results$prop.pattern),2),
                                 round(sum(all.results$brier.y*all.results$prop.pattern),2)), ncol = 2, bty='n',
       title = "Weighted MSE", col = c('black','red','blue','darkgreen'), pch = c(1,1,1,1,NA,NA,NA,NA))
points(all.results$brier.x*all.results$prop.pattern, all.results$metric.MI.noy.brier*all.results$prop.pattern, pch = 1, col = 'red')
points(all.results$brier.x*all.results$prop.pattern, all.results$metric.MIMI.noy.brier*all.results$prop.pattern, pch = 1, col = 'blue')
points(all.results$brier.x*all.results$prop.pattern, all.results$brier.y*all.results$prop.pattern, pch = 1, col = 'darkgreen')


plot(seq(1, length(all.results$Row.names), by=1),all.results$brier.x*all.results$prop.pattern,
     pch = NA, xlab = "Pattern", ylab = "Weighted MSE", xaxt='n', ylim = c(-2,60),
    main = 'sps ~ pafi + meanbp + wblc + alb + resp + temp + hrt + bili + crea + sod', cex.main = 1) #, xlim = c(0,5), ylim = c(0,5))


points(jitter(seq(1, length(all.results$Row.names), by=1)), all.results$metric.MI.noy.brier*all.results$prop.pattern, pch = 20, col = 'red')
points(jitter(seq(1, length(all.results$Row.names), by=1)), all.results$metric.MIMI.noy.brier*all.results$prop.pattern, pch = 20, col = 'blue')
points(jitter(seq(1, length(all.results$Row.names), by=1)), all.results$brier.y*all.results$prop.pattern, pch = 20, col = 'darkgreen')
points(jitter(seq(1, length(all.results$Row.names), by=1)), all.results$brier.x*all.results$prop.pattern, pch = 20, col = 'black')


text(x = seq(1, length(all.results$Row.names), by=1), par("usr")[3] - 2.5, 
     labels = all.results$Row.names
     , srt = 90, pos = 1, xpd = TRUE, cex = 0.5)
abline(v = c(2:length(all.results$Row.names))-0.5, lty = 2, lwd=0.5, col='gray')

legend('topright', legend = c("PMKS", "MI", "MIMI", "CCSM",round(sum(all.results$brier.x*all.results$prop.pattern),1),
                              round(sum(all.results$metric.MI.noy.brier*all.results$prop.pattern),1),
                              round(sum(all.results$metric.MIMI.noy.brier*all.results$prop.pattern),1),
                              round(sum(all.results$brier.y*all.results$prop.pattern),1)), ncol = 2, bty='n',
       title = "Total MSE", col = c('black','red','blue','darkgreen'), pch = c(20,20,20,20,NA,NA,NA,NA))

text(x = seq(1, length(all.results$Row.names), by=1), par("usr")[3] + 2.8, 
     labels = round(all.results$prop.pattern,2)
     ,pos = 1, xpd = TRUE, cex = 0.5)
text(x =-1,  par("usr")[3] + 2.8, labels = "Prop/Pattern",pos = 1, xpd = TRUE, cex = 0.5)


save(all.results, file = '/Users/SFM/Desktop/supportAdd25pafi')







###################################
#Plot the Betas for each pattern
##################################

#################################################
# This is the support dataset example
###############################################

source('/Users/SFM/Dropbox/GraduateSchool/Dissertation/MissingDataWithPredictionModels/Functions/PMKSfunctions.R')
dat <- read.csv(file = "/Users/SFM/Dropbox/GraduateSchool/Dissertation/MissingDataWithPredictionModels/support2/support2.csv")
dd <- datadist(dat)
options(datadist='dd')

#either include dzclass or dzgroup
dat <- subset(dat, select = c(sps,
                              pafi,
                              meanbp,
                              wblc,
                              alb,
                              resp,
                              temp,
                              hrt,
                              bili,
                              crea,
                              sod))


#remove the one person with a missing outcome
dat <- dat[-which(is.na(dat$sps)),]

#remove the few people with only one missing 
dat <- dat[-which(is.na(dat$meanbp)),]
#dat <- dat[-which(is.na(dat$resp)),]
#dat <- dat[-which(is.na(dat$temp)),]
#dat <- dat[-which(is.na(dat$hrt)),]
#dat <- dat[-which(is.na(dat$sod)),]


dat.I <- create.miss.ind(dat) 
dd <- datadist(dat.I)
options(datadist='dd')

imp.alldat <- aregImpute( ~ sps  +
                            pafi + 
                            meanbp + 
                            wblc + 
                            alb +
                            resp + 
                            temp + 
                            hrt + 
                            bili + 
                            crea + 
                            sod , data=dat.I, n.impute=10,pr=FALSE)

mod.alldat <- fit.mult.impute(sps ~  
                                     pafi + 
                                     meanbp + 
                                     wblc + 
                                     alb +
                                     resp + 
                                     temp + 
                                     hrt + 
                                     bili + 
                                     crea + 
                                     sod
                                   , ols,imp.alldat,data=dat.I,pr=FALSE)


mod.MIMI.alldat <- fit.mult.impute(sps ~  
                                     pafi + 
                                     meanbp + 
                                     wblc + 
                                     alb +
                                     resp + 
                                     temp + 
                                     hrt + 
                                     bili + 
                                     crea + 
                                     sod +
                                     m.pafi + 
                                     m.wblc +
                                     m.alb +
                                     m.bili +
                                     m.crea
                                   , ols,imp.alldat,data=dat.I,pr=FALSE)



ex.ccsm <- ccsm(DATA=dat, model=ex.fit, logistic=FALSE)
ex.pmks <- pmks(DATA=dat, model=ex.fit, logistic=FALSE)
  
pred.pmks <- predict.sm(prediction.data=dat, model=ex.fit, pmks.object=ex.pmks, logistic = FALSE) 

mod.betas <- lapply(pred.pmks, function(z) data.frame(t(data.frame(coef(z$mod)))))

all.betas <- do.call('rbind.fill', mod.betas)

par(las = 2)
plot(1:11, rep(NA,11), ylim = c(-5, 172), 
     typ = 'n', ylab = "Beta Coefficient", xlab = "",
     xaxt = 'n', main = "Variation of Beta Coefficients Across PMKS")
for(i in 1:ncol(all.betas)){
points(rep(i,23), all.betas[,i], pch = 20 )
}

points(1:11, coef(mod.alldat), pch = 2, col = 'red')
points(1:11, coef(mod.MIMI.alldat)[1:11], pch = 3, col = 'blue')

text(x = seq(1, length(names(all.betas)), by=1), par("usr")[3] - 3.5, 
     labels = c('Intercept',names(all.betas)[2:11])
     , srt = 90, pos = 1, xpd = TRUE, cex = 0.5)
abline(v = c(2:length(names(all.betas)))-0.5, lty = 2, lwd=0.5, col='gray')
legend('topright', legend = c("PMKS","MIMI","MI"), pch = c(20,2,3), col = c('black','red','blue'), bty = 'n')
