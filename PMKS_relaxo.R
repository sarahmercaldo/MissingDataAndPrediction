
library(glmnet)
library(lars)
library(relaxo)
library(corrplot)
library(cvTools)
library(plyr)
library(caret)
expit <- function(aa) { 
  exp_aa <- exp(aa)
  exp_aa/(1+exp_aa)
}

pmks.relaxlasso <- function(DATA, model, logistic=TRUE){
	    mod.DATA        <- get_all_vars(as.formula(model), data=DATA)
	    SDATA 		 <- mod.DATA[,-1] #remove the outcome 
	    tmp.dat      <- as.data.frame(is.na(SDATA)*1)
		tmp.pattern  <- factor(apply(tmp.dat,1,function(z) paste(z,collapse="")))
		all.patterns <- factor(apply(expand.grid(rep(list(0:1),ncol(SDATA))),1,function(z) paste(z,collapse="")))
		obs.patterns <- unique(tmp.pattern)
		tmp.info     <- split(seq(nrow(SDATA)), tmp.pattern)
		mp.levels    <- levels(tmp.pattern)
		mp.pattern   <- do.call(rbind, lapply(as.list(mp.levels),function(ZZ) strsplit(ZZ,'')[[1]])) 		
		mp.info     <- data.frame(cbind(names(tmp.info), unlist(lapply(tmp.info, length))),
		                          stringsAsFactors= FALSE)
		rownames(mp.info) <- seq(nrow(mp.info))
		colnames(mp.info) <- c('mp','n')
		if(length(setdiff(all.patterns,obs.patterns)) == 0){
		  empty.patterns = NULL
		} else {
		empty.patterns <- data.frame(mp = factor(setdiff(all.patterns,obs.patterns)), n=0)
		}	
		mp.info <- rbind(mp.info,empty.patterns)
		
		mod.rhs <- strsplit(model, ' ~')[[1]][1]
		mod.lhs <- strsplit(model, '~')[[1]][2]
		mod.lhs <- strsplit(mod.lhs,' ')[[1]]
		mod.lhs <- mod.lhs[!(mod.lhs%in%c('','+','~'))]

		cc <- ncol(model.matrix(as.formula(model),mod.DATA))
			
		threshold <- cc*2
		mp.info$use.ptmx <- (as.numeric(mp.info$n)>=threshold)*1
		
		reg.out <- vector('list', length(all.patterns))
		names(reg.out) <- mp.info$mp.info

		
		for(ixx in seq(nrow(mp.info))) {
			col.keep  <- which(strsplit(mp.info$mp[ixx],'')[[1]]=='0')
			if(length(col.keep)==0){
			new.mod <- as.formula(paste(mod.rhs,1,sep='~'))
			reg.out[[ixx]] <- list(pattern = mod.lhs[col.keep],
							mod = glm(new.mod,family = 'binomial',DATA))
			} else {
			new.mod   <- as.formula(paste(mod.rhs,paste(mod.lhs[col.keep],collapse='+'),
				                      sep='~'))
				                      
			if(mp.info$use.ptmx[ixx]==1) {
				x <- get_all_vars(as.formula(new.mod), data=mod.DATA[tmp.info[[ixx]],])[,-1]	
				
				x <- apply(x, 2, function(zzz) {z = NULL; z = cbind(z,zzz)})
				y <- get_all_vars(as.formula(new.mod), data=mod.DATA[tmp.info[[ixx]],])[,1]
				y <- cbind(y)
				
				if(length(col.keep)==1){
					reg.out[[ixx]] <- list(pattern = mod.lhs[col.keep],
							mod = glm(as.formula(paste(mod.rhs,paste(mod.lhs[col.keep],collapse='+'),
				                      sep='~')),family = 'binomial',
				                      data=mod.DATA[tmp.info[[ixx]],]))
				} else {
				l1.cv <- cv.glmnet(x,y,alpha=1,family="binomial")	
				b <- coef(l1.cv,s="lambda.1se")
				colnames(b) <- "lambda.1se"
				subset <- colnames(x)[which(b!=0)-1]
				if(length(subset)==0){
				reg.out[[ixx]] <- list(pattern = mod.lhs[col.keep],
							mod = glm(as.formula(paste(mod.rhs,1,sep='~')),family = 'binomial',
							data=mod.DATA[tmp.info[[ixx]],]))	
				} else {reg.out[[ixx]] <- list(pattern = mod.lhs[col.keep],
							mod = glm(as.formula(paste(mod.rhs,paste(mod.lhs[mod.lhs%in%subset],collapse='+'),
				                      sep='~')),family = 'binomial',
				                      data=mod.DATA[tmp.info[[ixx]],]))}
				}
				
				 } else {
				dat.subset <- as.data.frame(DATA[,c(1,col.keep + 1)])
				colnames(dat.subset) <- c(mod.rhs,mod.lhs[col.keep])
				cc <- dat.subset[complete.cases(dat.subset),]
				
				x <- get_all_vars(as.formula(new.mod), data=cc)[,-1]	
				y <- get_all_vars(as.formula(new.mod), data=cc)[,1]
				y <- cbind(y)

				if(length(col.keep)==1){
					reg.out[[ixx]] <- list(pattern = mod.lhs[col.keep],
							mod = glm(as.formula(paste(mod.rhs,paste(mod.lhs[col.keep],collapse='+'),
				                      sep='~')),family = 'binomial',
				                      data=mod.DATA))
				} else {
				x <- apply(x, 2, function(zzz) {z = NULL; z = cbind(z,zzz)})
				y <- get_all_vars(as.formula(new.mod), data=cc)[,1]
				y <- cbind(y)

				l1.cv <- cv.glmnet(x,y,alpha=1,family="binomial")	
				b <- coef(l1.cv,s="lambda.1se")
				colnames(b) <- "lambda.1se"
				subset <- colnames(x)[which(b!=0)-1]
				if(length(subset)==0){
				reg.out[[ixx]] <- list(pattern = mod.lhs[col.keep],
							mod = glm(as.formula(paste(mod.rhs,1,sep='~')),family = 'binomial',
							data=mod.DATA))	
				} else {reg.out[[ixx]] <- list(pattern = mod.lhs[col.keep],
							mod = glm(as.formula(paste(mod.rhs,paste(mod.lhs[mod.lhs%in%subset],collapse='+'),
				                      sep='~')),family = 'binomial',
				                      data=mod.DATA))
							}
				}			
			} 
		   }
		}
		reg.out
}

#Fit all the complete case submodels
ccsm.relaxlasso <- function(DATA, model, logistic=TRUE){
	     mod.DATA        <- get_all_vars(as.formula(model), data=DATA)
	    SDATA 		 <- mod.DATA[,-1] #remove the outcome 
	    tmp.dat      <- as.data.frame(is.na(SDATA)*1)
		tmp.pattern  <- factor(apply(tmp.dat,1,function(z) paste(z,collapse="")))
		all.patterns <- factor(apply(expand.grid(rep(list(0:1),ncol(SDATA))),1,function(z) paste(z,collapse="")))
		obs.patterns <- unique(tmp.pattern)
		tmp.info     <- split(seq(nrow(SDATA)), tmp.pattern)
		mp.levels    <- levels(tmp.pattern)
		mp.pattern   <- do.call(rbind, lapply(as.list(mp.levels),function(ZZ) strsplit(ZZ,'')[[1]])) 		
		mp.info     <- data.frame(cbind(names(tmp.info), unlist(lapply(tmp.info, length))),
		                          stringsAsFactors= FALSE)
		rownames(mp.info) <- seq(nrow(mp.info))
		colnames(mp.info) <- c('mp','n')
		if(length(setdiff(all.patterns,obs.patterns)) == 0){
		  empty.patterns = NULL
		} else {
		empty.patterns <- data.frame(mp = factor(setdiff(all.patterns,obs.patterns)), n=0)
		}	
		mp.info <- rbind(mp.info,empty.patterns)
		
		mod.rhs <- strsplit(model, ' ~')[[1]][1]
		mod.lhs <- strsplit(model, '~')[[1]][2]
		mod.lhs <- strsplit(mod.lhs,' ')[[1]]
		mod.lhs <- mod.lhs[!(mod.lhs%in%c('','+','~'))]

		cc <- ncol(model.matrix(as.formula(model),mod.DATA))
			
		threshold <- cc*2
		mp.info$use.ptmx <- (as.numeric(mp.info$n)>=threshold)*1
		
		reg.out <- vector('list', length(all.patterns))
		names(reg.out) <- mp.info$mp.info

		
		for(ixx in seq(nrow(mp.info))) {
			col.keep  <- which(strsplit(mp.info$mp[ixx],'')[[1]]=='0')
			if(length(col.keep)==0){
			new.mod <- as.formula(paste(mod.rhs,1,sep='~'))
			reg.out[[ixx]] <- list(pattern = mod.lhs[col.keep],
							mod = glm(new.mod,family = 'binomial',DATA))
			} else {
			new.mod   <- as.formula(paste(mod.rhs,paste(mod.lhs[col.keep],collapse='+'),
				                      sep='~'))
				                      
			
				dat.subset <- as.data.frame(DATA[,c(1,col.keep + 1)])
				colnames(dat.subset) <- c(mod.rhs,mod.lhs[col.keep])
				cc <- dat.subset[complete.cases(dat.subset),]
				
				x <- get_all_vars(as.formula(new.mod), data=cc)[,-1]	
				y <- get_all_vars(as.formula(new.mod), data=cc)[,1]
				y <- cbind(y)

				if(length(col.keep)==1){
					reg.out[[ixx]] <- list(pattern = mod.lhs[col.keep],
							mod = glm(as.formula(paste(mod.rhs,paste(mod.lhs[col.keep],collapse='+'),
				                      sep='~')),family = 'binomial',
				                      data=mod.DATA))
				} else {
				x <- apply(x, 2, function(zzz) {z = NULL; z = cbind(z,zzz)})
				y <- get_all_vars(as.formula(new.mod), data=cc)[,1]
				y <- cbind(y)

				l1.cv <- cv.glmnet(x,y,alpha=1,family="binomial")	
				b <- coef(l1.cv,s="lambda.1se")
				colnames(b) <- "lambda.1se"
				subset <- colnames(x)[which(b!=0)-1]
				if(length(subset)==0){
				reg.out[[ixx]] <- list(pattern = mod.lhs[col.keep],
							mod = glm(as.formula(paste(mod.rhs,1,sep='~')),family = 'binomial',
							data=mod.DATA))	
				} else {reg.out[[ixx]] <- list(pattern = mod.lhs[col.keep],
							mod = glm(as.formula(paste(mod.rhs,paste(mod.lhs[mod.lhs%in%subset],collapse='+'),
				                      sep='~')),family = 'binomial',
				                      data=mod.DATA))
							}
				}			
			} 
		   
		}
		reg.out

}

last <- function(x) { return( x[length(x)] ) }

#Predict for a new individual using a PMKS or CCSM object.  
predict.sm.relaxlasso <- function(prediction.data, model, pmks.object, logistic = TRUE){
		mod.DATA     <- get_all_vars(as.formula(model), data=prediction.data)
	    PDATA 		 <- mod.DATA[,-1] #remove the outcome 
		tmp.dat      <- as.data.frame(is.na(PDATA)*1)
		tmp.pattern  <- factor(apply(tmp.dat,1,function(z) paste(z,collapse="")))
		tmp.info     <- split(seq(nrow(PDATA)), tmp.pattern)
		mp.levels    <- levels(tmp.pattern)
		mp.pattern   <- do.call(rbind, lapply(as.list(mp.levels),function(ZZ) strsplit(ZZ,'')[[1]])) 		
		mp.info     <- data.frame(cbind(names(tmp.info), unlist(lapply(tmp.info, length))),
		                          stringsAsFactors= FALSE)
		rownames(mp.info) <- seq(nrow(mp.info))
		colnames(mp.info) <- c('mp','n')
	    mod.rhs <- strsplit(model, '~')[[1]][1]		
	    mod.rhs <- strsplit(mod.rhs,' ')[[1]]
		mod.lhs <- strsplit(model, '~')[[1]][2]
		mod.lhs <- strsplit(mod.lhs,' ')[[1]]
		mod.lhs <- mod.lhs[!(mod.lhs%in%c('','+','~'))]
		pred.out <- vector('list', nrow(mp.info))
		#For the different patterns 
		for(ixx in seq(length(tmp.info))){
			col.keep  <- which(strsplit(mp.info$mp[ixx],'')[[1]]=='0')
			pattern   <- mod.lhs[col.keep]
			which.mod <- which(lapply(pmks.object, function(z) identical(z$pattern, pattern))==TRUE)
						linear.predictor <- predict(pmks.object[[which.mod]]$mod, 
        										   prediction.data[tmp.info[[ixx]],])
						
      pred.out[[ixx]] <- list(
        		            numeric.pattern = mp.info$mp[ixx],
        		            
        		            pattern = pattern,
        		            
        		            mod = pmks.object[[which.mod]]$mod,
        		            
        								lin.pred = linear.predictor,
        								truth = prediction.data[tmp.info[[ixx]],mod.rhs],
        								strat.auc = auc(expit(linear.predictor),
        										   prediction.data[tmp.info[[ixx]],mod.rhs]),
        								strat.brier = 	brier.score(expit(linear.predictor),
        										   prediction.data[tmp.info[[ixx]],mod.rhs]),
        								strat.logscore = 	logarithmic.scoring.rule(expit(linear.predictor),
        								                           prediction.data[tmp.info[[ixx]],mod.rhs]))		
		
		}
				pred.out
	}
	


source('PMKSfunctions.R')
dat <- read.csv(file = "support2.csv")
dat$sps <- ifelse(dat$sps>=23.90,1,0)

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
  
  ex.ccsm <- ccsm.relaxlasso(DATA=trainingset, model=ex.fit, logistic=TRUE)
  ex.pmks <- pmks.relaxlasso(DATA=trainingset, model=ex.fit, logistic=TRUE)
  
  pred.pmks <- tryCatch(predict.sm.relaxlasso(prediction.data=testset, model=ex.fit, pmks.object=ex.pmks, logistic = TRUE),
                        error = function(e) {NA})
  if(is.na(pred.pmks)) pred.pmks.res <- NA
  if(!is.na(pred.pmks)) pred.pmks.res <- do.call(rbind,lapply(pred.pmks, function(z) data.frame(pattern = z$numeric.pattern, brier = z$strat.brier, logscore = z$strat.logscore,auc=z$strat.auc, prop.pattern = (length(z$truth)/nrow(testset)))))
  
  pred.ccsm <- tryCatch(predict.sm.relaxlasso(prediction.data=testset, model=ex.fit, pmks.object=ex.ccsm, logistic = TRUE),
                        error = function(e) {NA})
  if(is.na(pred.ccsm)) pred.ccsm.res <- NA
  if(!is.na(pred.ccsm)) pred.ccsm.res <- do.call(rbind,lapply(pred.ccsm, function(z) data.frame(pattern = z$numeric.pattern,brier = z$strat.brier, logscore = z$strat.logscore,auc=z$strat.auc,  prop.pattern = (length(z$truth)/nrow(testset)))))
  
  results.pmks[[i]] <- pred.pmks.res
  results.ccsm[[i]] <- pred.ccsm.res
  
}


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

results.pmks.avg$pmks.brier.weight.avg = results.pmks.avg[,1]*results.pmks.avg[,4]
results.pmks.avg$pmks.log.weight.avg = results.pmks.avg[,2]*results.pmks.avg[,4]
results.pmks.avg$pmks.auc.weight.avg = results.pmks.avg[,3]*results.pmks.avg[,4]

results.ccsm.avg$ccsm.brier.weight.avg = results.ccsm.avg[,1]*results.ccsm.avg[,4]
results.ccsm.avg$ccsm.log.weight.avg = results.ccsm.avg[,1]*results.ccsm.avg[,4]
results.ccsm.avg$ccsm.auc.weight.avg = results.ccsm.avg[,1]*results.ccsm.avg[,4]
all.results <- merge(results.pmks.avg,results.ccsm.avg, by = 0)



##################################################################
# Add 25 to pafi before dichotomizing it
##################################################################


dat <- read.csv(file = "support2.csv")

dat[which(is.na(dat[,'pafi'])),'sps'] <- dat[which(is.na(dat[,'pafi'])),'sps'] + rnorm(length(which(is.na(dat[,'pafi']))),mean = 10, sd = 2)

dat$sps <- ifelse(dat$sps >= 23.90,1,0)

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

#####################################################
# Cross Validate within each missing data pattern
####################################################
set.seed(3)
k=10


results.pmks <- vector('list',length = k)
results.ccsm <- vector('list',length = k)



ex.fit <- "sps ~ pafi + meanbp + wblc + alb + resp +  temp + hrt +  bili + crea + sod"

flds <- createFolds(dat$sps, k=k, list=TRUE, returnTrain=FALSE)



ex.fit <- "sps ~ pafi + meanbp + wblc + alb + resp +  temp + hrt +  bili + crea + sod"

flds <- createFolds(dat$sps, k=k, list=TRUE, returnTrain=FALSE)

for( i in 1:k){
  
 
  trainingset <- dat[-flds[[i]],]
  testset <- dat[flds[[i]],]
  
  ex.ccsm <- ccsm.relaxlasso(DATA=trainingset, model=ex.fit, logistic=TRUE)
  ex.pmks <- pmks.relaxlasso(DATA=trainingset, model=ex.fit, logistic=TRUE)
  
  pred.pmks <- tryCatch(predict.sm.relaxlasso(prediction.data=testset, model=ex.fit, pmks.object=ex.pmks, logistic = TRUE),
                        error = function(e) {NA})
  if(is.na(pred.pmks)) pred.pmks.res <- NA
  if(!is.na(pred.pmks)) pred.pmks.res <- do.call(rbind,lapply(pred.pmks, function(z) data.frame(pattern = z$numeric.pattern, brier = z$strat.brier, logscore = z$strat.logscore,auc=z$strat.auc, prop.pattern = (length(z$truth)/nrow(testset)))))
  
  pred.ccsm <- tryCatch(predict.sm.relaxlasso(prediction.data=testset, model=ex.fit, pmks.object=ex.ccsm, logistic = TRUE),
                        error = function(e) {NA})
  if(is.na(pred.ccsm)) pred.ccsm.res <- NA
  if(!is.na(pred.ccsm)) pred.ccsm.res <- do.call(rbind,lapply(pred.ccsm, function(z) data.frame(pattern = z$numeric.pattern,brier = z$strat.brier, logscore = z$strat.logscore,auc=z$strat.auc,  prop.pattern = (length(z$truth)/nrow(testset)))))
  
  results.pmks[[i]] <- pred.pmks.res
  results.ccsm[[i]] <- pred.ccsm.res
}

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

results.pmks.avg$pmks.brier.weight.avg = results.pmks.avg[,1]*results.pmks.avg[,4]
results.pmks.avg$pmks.log.weight.avg = results.pmks.avg[,2]*results.pmks.avg[,4]
results.pmks.avg$pmks.auc.weight.avg = results.pmks.avg[,3]*results.pmks.avg[,4]

results.ccsm.avg$ccsm.brier.weight.avg = results.ccsm.avg[,1]*results.ccsm.avg[,4]
results.ccsm.avg$ccsm.log.weight.avg = results.ccsm.avg[,1]*results.ccsm.avg[,4]
results.ccsm.avg$ccsm.auc.weight.avg = results.ccsm.avg[,1]*results.ccsm.avg[,4]

all.results <- merge(results.pmks.avg,results.ccsm.avg, by = 0)



