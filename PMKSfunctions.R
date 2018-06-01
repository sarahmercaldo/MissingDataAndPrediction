
library(plyr)
library(rms)
library(mice)
library(cvTools)
library(MASS)
library(Hmisc)
library(lme4)
library(mi)
library(doSNOW)

########################
# Functions to use
#######################


#This function creates missing data indicators 
#for your dataset and then appends those indicators
#To the original dataset

create.miss.ind <- function(DATA){
	 tmp.dat <- as.data.frame(is.na(DATA)*1)
	 names(tmp.dat) <- paste('m.',names(tmp.dat),sep="")
    cbind(DATA,tmp.dat)
}

#This function creates missing data indicators 
#for your dataset 

create.miss.ind.only <- function(DATA){
	 tmp.dat <- as.data.frame(is.na(DATA)*1)
	 names(tmp.dat) <- paste('m.',names(tmp.dat),sep="")
    tmp.dat
}




#Function to combine two data frames and match them up by column names
rbind.match.columns <- function(input1, input2) {
    n.input1 <- ncol(input1)
    n.input2 <- ncol(input2)
    if (n.input2 < n.input1) {
        TF.names <- which(names(input2) %in% names(input1))
        column.names <- names(input2[, TF.names])
    } else {
        TF.names <- which(names(input1) %in% names(input2))
        column.names <- names(input1[, TF.names])
    }
    return(rbind(input1[, column.names], input2[, column.names]))
}

#predict an out of sample individual where you have a single model created
#by mice, takes the mice model object and the newdata from a new individual
predict.oos.mice <- function(FIT, NEWDATA){
	mod.call <- as.character(FIT$call1)[3]
	mod.call <- strsplit(mod.call, '\\(')[[1]][2]
	mod.call <- strsplit(mod.call, ',')[[1]][1]
	lp <- NA
	bb <- FIT$qbar
	mm <- model.matrix(as.formula(mod.call), NEWDATA)
	if(length(bb) != ncol(mm)) {
		stop('issue with mm')
	} else {
		lp <- mm %*% matrix(bb, ncol=1)
	}
	lp
}

#AUC function
auc=function(score,status){
######################
## auc version 1.0 
## Compute Area under Rmpirical ROC curve by Trapezoidal Rule
## Author: J. Blume
## Date:   July 2014
######################

	pos=score[status==1]
	neg=score[status==0]

	ct1=sum(outer(pos,neg,">"))
	ct2=sum(outer(pos,neg,"=="))
	den=length(pos)*length(neg)
	auc=(ct1+0.5*ct2)/den
	auc=max(auc,1-auc)
	auc	

}


#Function to calculate the brier score/MSE
brier.score <- function(pred, outcome){
	mean((pred - outcome)^2)
}

#Function to calculate the logarithmic scoring rule
logarithmic.scoring.rule <- function(pred, outcome){
	
mean(outcome*log(pred) + (1-outcome)*log(1-pred))	
	
}

#########################################################################################
#1) Find all the patterns of missingness
#2) For patterns of missingness with membership >= 2p + 2, fit a submodel with all the
#   observed covariates, ussing only the people in that pattern.
#3) For patterns of missingness with mempership < 2p + 2 fit the complete case submodel, which
#   only includes the observed covariates in that patter, but if fitted with ALL the 
#   individuals in the data set
#4) Allow function to handle both linear and logistic models.
#5) Output prediction models in a way that would be usable with a validation 
#   dataset or new person
# For all unobserved patterns do the complete case submodel
#########################################################################################


pmks <- function(DATA, model, logistic=TRUE){
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
		
		mod.rhs <- strsplit(model, '~')[[1]][1]
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
			} else {new.mod   <- as.formula(paste(mod.rhs,paste(mod.lhs[col.keep],collapse='+'),
				                      sep='~'))}
				                      
			if(mp.info$use.ptmx[ixx]==1) {
				if(logistic == TRUE){
				reg.out[[ixx]] <- list(pattern = mod.lhs[col.keep],
							mod = glm(new.mod,data=mod.DATA[tmp.info[[ixx]],],family = 'binomial'))
				} else {
				  reg.out[[ixx]] <- list(pattern = mod.lhs[col.keep],
				                         mod = glm(new.mod,data=mod.DATA[tmp.info[[ixx]],],family = 'gaussian'))  
				}
			} else {
			  if(logistic == TRUE){
			  reg.out[[ixx]] <- list(pattern = mod.lhs[col.keep],
			  						 mod	 = glm(new.mod,data=mod.DATA, family='binomial'))
			  } else {
			  reg.out[[ixx]] <- list(pattern = mod.lhs[col.keep],
			                           mod	 = glm(new.mod,data=mod.DATA, family='gaussian'))
			  }
			}
		}
		reg.out
}


#Fit all the complete case submodels
ccsm <- function(DATA, model, logistic=TRUE){
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
		empty.patterns <- data.frame(mp = factor(setdiff(all.patterns,obs.patterns)), n=0)	
		mp.info <- rbind(mp.info,empty.patterns)
		
		mod.rhs <- strsplit(model, '~')[[1]][1]
		mod.lhs <- strsplit(model, '~')[[1]][2]
		mod.lhs <- strsplit(mod.lhs,' ')[[1]]
		mod.lhs <- mod.lhs[!(mod.lhs%in%c('','+','~'))]

		reg.out <- vector('list', length(all.patterns))
		names(reg.out) <- mp.info$mp.info

		for(ixx in seq(nrow(mp.info))) {
			col.keep  <- which(strsplit(mp.info$mp[ixx],'')[[1]]=='0')
			if(length(col.keep)==0){
			new.mod <- as.formula(paste(mod.rhs,1,sep='~'))
			} else {new.mod   <- as.formula(paste(mod.rhs,paste(mod.lhs[col.keep],collapse='+'),
				                      sep='~'))}			
			if(logistic == TRUE){
			  # Use complete case submodel
			  reg.out[[ixx]] <- list(pattern = mod.lhs[col.keep],
			  						 mod	 = glm(new.mod, data=mod.DATA, family = 'binomial'))
			} else {reg.out[[ixx]] <- list(pattern = mod.lhs[col.keep],
			                               mod	 = glm(new.mod, data=mod.DATA, family = 'gaussian'))}
		}
		reg.out
}


#Predict for a new individual using a PMKS or CCSM object.  
predict.sm <- function(prediction.data, model, pmks.object, logistic = TRUE){
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
			
			
			if(logistic == TRUE){
			
      pred.out[[ixx]] <- list(
        		            numeric.pattern = mp.info$mp[ixx],
        		            
        		            pattern = pattern,
        		            
        		            mod = pmks.object[[which.mod]]$mod,
        		            
        								lin.pred = predict(pmks.object[[which.mod]]$mod, 
        										   prediction.data[tmp.info[[ixx]],]),
        								truth = prediction.data[tmp.info[[ixx]],mod.rhs],
        								strat.auc = auc(expit(predict(pmks.object[[which.mod]]$mod, 
        										   prediction.data[tmp.info[[ixx]],])),
        										   prediction.data[tmp.info[[ixx]],mod.rhs]),
        								strat.brier = 	brier.score(expit(predict(pmks.object[[which.mod]]$mod, 
        										   prediction.data[tmp.info[[ixx]],])),
        										   prediction.data[tmp.info[[ixx]],mod.rhs]),
        								strat.logscore = 	logarithmic.scoring.rule(expit(predict(pmks.object[[which.mod]]$mod, 
        								                                   prediction.data[tmp.info[[ixx]],])),
        								                           prediction.data[tmp.info[[ixx]],mod.rhs]))		
		} else {
		  pred.out[[ixx]] <- list(
		    numeric.pattern = mp.info$mp[ixx],
		    
		    pattern = pattern,
		    
		    mod = pmks.object[[which.mod]]$mod,
		    
		    lin.pred = predict(pmks.object[[which.mod]]$mod, 
		                       prediction.data[tmp.info[[ixx]],]),
		    truth = prediction.data[tmp.info[[ixx]],mod.rhs],
		    
		    strat.brier = 	brier.score(predict(pmks.object[[which.mod]]$mod, 
		                                             prediction.data[tmp.info[[ixx]],]),
		                               prediction.data[tmp.info[[ixx]],mod.rhs]))
			}
		}
				pred.out
	}
	

#Function to find all the observed missing data patterns
which.pattern <- function(DATA, model){
  
  mod.DATA        <- get_all_vars(as.formula(model), data=DATA)
  SDATA 		      <- mod.DATA[,-1] #remove the outcome 
  tmp.dat         <- as.data.frame(is.na(SDATA)*1)
  tmp.pattern     <- factor(apply(tmp.dat,1,function(z) paste(z,collapse="")))
  
   pattern <- tmp.pattern 
   pattern
}


#Get metrics by pattern for a logistic model
get_metric <- function(DATA) {
  colnames(DATA) <- c('pred','true','pattern')
  do.call(rbind, by(DATA, pattern,function(xx) {
    data.frame('log'=logarithmic.scoring.rule(xx$pred,xx$true),
               'auc'=auc(xx$pred,xx$true),
               'brier'=brier.score(xx$pred,xx$true)
               , 'prop.pattern' = (length(xx$true)/nrow(DATA)))}))
}

#Get metrics by pattern for a linear model
get_metric_linear <- function(DATA) {
  colnames(DATA) <- c('pred','true','pattern')
  do.call(rbind, lapply( split(DATA, DATA$pattern), function(xx) data.frame('brier'=brier.score(xx$pred,xx$true), 'prop.pattern' = (length(xx$true)/nrow(DATA))) ))    
}

#MSE from a fitted model
mse <- function(sm) { 
  mse <- mean(sm$residuals^2)
  return(mse)
}

#logit function
logit <- function(p) log(p/(1-p))


#Empirically calculate the intercept for the missing data 
#mechanism 
miss.int <- function(miss.mech, betaM=1,p.miss){
  if(ncol(miss.mech)!=2){
    logit(p.miss) - betaM*mean(miss.mech[,1])
  } else {
    logit(p.miss) - betaM*mean((((miss.mech[,2]/sd(miss.mech[,2]))
                                 + miss.mech[,1])/as.numeric(sqrt(2*(1+cor(miss.mech[,2],miss.mech[,1]))))))	
  }
}

#Function to create a data frame for simulation
dat.frame <- function(n.tr,
                      mu.z,
                      mu.w,
                      mu.q,
                      s.z,
                      s.w,
                      s.q,
                      s.wz,
                      s.zq,
                      s.wq,
                      p.miss,
                      missing.type,
                      b.x,
                      b.z,
                      b.w,
                      b.m,
                      mu.e,
                      s.e,
                      betaM=1){
  x = rep(1,n.tr)
  mu = c(mu.z,mu.w,mu.q)
  Sigma = matrix(c(1,s.wz,s.zq,
                   s.wz,1,s.wq,
                   s.zq,s.wq,1), byrow=TRUE, ncol=3)
  dat =  as.data.frame(mvrnorm(n=n.tr, mu=mu, Sigma=Sigma, empirical=TRUE))
  colnames(dat) = c('z','w','q')
  coef.tru	= rbind( b.x , b.z , b.w)
  dat$y=cbind(x,dat[,'z'], dat[,'w'])%*%coef.tru + rnorm(n.tr,mean=mu.e,sd=s.e)
  
  if(missing.type=='MCAR'){
    m = rbinom(n.tr,1,p.miss)
    y.pmm = b.x*x + b.z*dat[,'z'] +  b.w*dat[,'w'] + b.m*m + rnorm(n.tr,mean=mu.e,sd=s.e)
    data.frame(x=x, z=dat[,'z'], w=dat[,'w'], q=dat[,'q'],m=m, y=dat[,'y'], y.pmm = y.pmm)
  } else if(missing.type=='MAR'){
    m = rbinom(n.tr, 1, expit(miss.int(miss.mech=cbind(dat[,'w']), 
                                       betaM=betaM,p.miss) + betaM*dat[,'w']))
    y.pmm = b.x*x + b.z*dat[,'z'] +  b.w*dat[,'w'] + b.m*m + rnorm(n.tr,mean=mu.e,sd=s.e)
    data.frame(x=x, z=dat[,'z'], w=dat[,'w'], q=dat[,'q'], m=m, y=dat[,'y'], y.pmm = y.pmm)
  } else if(missing.type=='MNAR'){
    m = rbinom(n.tr, 1, expit(miss.int(miss.mech=cbind(dat[,'z']), 
                                       betaM=betaM,p.miss) + betaM*dat[,'z']))
    y.pmm = b.x*x + b.z*dat[,'z'] +  b.w*dat[,'w'] + b.m*m + rnorm(n.tr,mean=mu.e,sd=s.e)
    data.frame(x=x, z=dat[,'z'], w=dat[,'w'], q=dat[,'q'], m=m, y=dat[,'y'], y.pmm = y.pmm)
  } else if(missing.type=='MARY'){
    m = rbinom(n.tr, 1, 
               expit(miss.int(miss.mech=dat[,c('w','y')], betaM=betaM,p.miss)
                     + betaM*(((dat[,'y']/sd(dat[,'y'])) + dat[,'w'])/as.numeric(sqrt(2*(1+cor(dat[,'y'],dat[,'w'])))))))
    data.frame(x=x, z=dat[,'z'], w=dat[,'w'], q=dat[,'q'], m=m, y=dat[,'y'])
  } else if(missing.type=='MNARY'){
    m = rbinom(n.tr, 1, 
               expit(miss.int(miss.mech=dat[,c('z','y')], betaM=betaM,p.miss)
                     + betaM*(((dat[,'y']/sd(dat[,'y'])) + dat[,'z'])/as.numeric(sqrt(2*(1+cor(dat[,'y'],dat[,'z'])))))))
    data.frame(x=x, z=dat[,'z'], w=dat[,'w'],  q=dat[,'q'], m=m, y=dat[,'y'])
  }
  
}



#Create an out of sample data frame, returns a list
# with different kinds of imputations done on the
# out of sample individuals
oos.data <- function(n.tr.new,
                     mu.z.new,
                     mu.w.new,
                     mu.q.new,
                     s.z.new,
                     s.w.new,
                     s.q.new,
                     s.wz.new,
                     s.zq.new,
                     s.wq.new,
                     p.miss.new,
                     missing.type.new,
                     b.x.new,
                     b.z.new,
                     b.w.new,
                     mu.e.new,
                     s.e.new,
                     betaM.new=1,
                     n.imp = 5,
                     original.dat.miss,
                     Y.PM
) {
  
  dat.new <- dat.frame(n.tr=n.tr.new,
                       mu.z = mu.z.new,
                       mu.w = mu.w.new,
                       mu.q = mu.q.new,
                       s.z = s.z.new,
                       s.w = s.w.new,
                       s.q = s.q.new,
                       s.wz = s.wz.new,
                       s.zq = s.zq.new,
                       s.wq = s.wq.new,
                       p.miss = p.miss.new,
                       missing.type = missing.type.new,
                       b.x = b.x.new,
                       b.z = b.z.new,
                       b.w = b.w.new,
                       b.m = betaM.new,
                       mu.e = mu.e.new,
                       s.e = s.e.new,
                       betaM = betaM.new)
  if(Y.PM == TRUE) {dat.new$y <- dat.new$y.pmm}
  
  dat.miss.new <- dat.new
  dat.miss.new$z <- ifelse(dat.new$m==1,NA,dat.new$z)
  
  imputed.dat.cm.mipack <- imputed.dat.cm.Hmisc <- dat.miss.new
  imputed.dat.mi <- imputed.dat.miq <- imputed.dat.miy <- rep(list(dat.miss.new),n.imp)
  
  dat.cc.mean <- dat.miss.new
  dat.cc.mean[which(is.na(dat.cc.mean$z)==TRUE),'z'] <- mean(original.dat.miss$z, na.rm=TRUE) 
  
  dat.cond.mean <- dat.miss.new
  mod.predict.z <- lm(z ~ w, data = original.dat.miss)
  dat.cond.mean[dat.cond.mean$m==1,'z'] <- predict(mod.predict.z, dat.miss.new[dat.miss.new$m==1,])
  
  for(i in which(dat.new$m==1)){
    newperson <- dat.miss.new[i,]
    newperson$y <- NA 
    addNewPatient <- rbind.match.columns(original.dat.miss, newperson)
    
    #Single Conditional Mean Imputation Using MI
    mdf <- suppressMessages(suppressWarnings(missing_data.frame(addNewPatient[,c('x','z','w')])))
    
    #the expectation argument is conditional mean expectation
    mdf <-  suppressMessages(change(mdf, y=c('x','z','w'),what = "imputation_method", to = "expectation"))
    imputations <- mi(mdf,n.iter=1,n.chains=1,verbose=FALSE,parallel=FALSE)
    
    imputed.dat.cm.mipack[i,c('x','z','w')] <- complete(imputations,1)[nrow(addNewPatient),c('x','z','w')]
    
    f <- aregImpute( ~ w + z, data=addNewPatient, n.impute=1,type='regression',pr=FALSE)
    imputed.dat.cm.Hmisc[i,c('z','w')] <- as.data.frame(impute.transcan(f, 
                                                                        imputation=1, data=addNewPatient, list.out=TRUE,
                                                                        pr=FALSE, check=FALSE))[nrow(addNewPatient),c('z','w')]
    
    # imp.pmm  <- aregImpute(~ w + z, n.impute=n.imp, x=TRUE,
    #                        nk=3, tlinear=F, data=addNewPatient, pr=FALSE)
    # imp.pmm.q  <- aregImpute(~ w + z + q, n.impute=n.imp, x=TRUE, 
    #                          nk=3, tlinear=F, data=addNewPatient, pr=FALSE)
    # imp.pmm.y  <- aregImpute(~  w + z + y, n.impute=n.imp, x=TRUE,
    #                          nk=3, tlinear=F, data=addNewPatient, pr=FALSE)
      
    
    imp.pmm  <- aregImpute(~ w + z, n.impute=n.imp, x=TRUE,
                           nk=0, tlinear=TRUE, data=addNewPatient, pr=FALSE)
    imp.pmm.q  <- aregImpute(~ w + z + q, n.impute=n.imp, x=TRUE,
                             nk=0, tlinear=TRUE, data=addNewPatient, pr=FALSE)
    imp.pmm.y  <- aregImpute(~  w + z + y, n.impute=n.imp, x=TRUE,
                             nk=0, tlinear=TRUE, data=addNewPatient, pr=FALSE)
    
    for(j in 1:n.imp) {
      imputed.dat.mi[[j]][i,c('z')] <- as.data.frame(impute.transcan(imp.pmm, imputation=j,
                                                                     data=addNewPatient, list.out=TRUE,
                                                                     pr=FALSE, check=FALSE))[nrow(addNewPatient),c('z')]
      
      imputed.dat.miq[[j]][i,c('z')] <- as.data.frame(impute.transcan(imp.pmm.q, imputation=j,
                                                                      data=addNewPatient, list.out=TRUE,
                                                                      pr=FALSE, check=FALSE))[nrow(addNewPatient),c('z')] 
      
      imputed.dat.miy[[j]][i,c('z')] <- as.data.frame(impute.transcan(imp.pmm.y, imputation=j,
                                                                      data=addNewPatient, list.out=TRUE,
                                                                      pr=FALSE, check=FALSE))[nrow(addNewPatient),c('z')]                               
    }
    
  }
  
  #I'm pretty sure taking the mean of all the imputations and getting the prediction is the same as getting 10 predictions and then taking the mean - This is only true for out of sample individuals!
  
  list(
    dat.miss.new = dat.miss.new,
    imputed.dat.cm.mipack = imputed.dat.cm.mipack,
    imputed.dat.cm.Hmisc = imputed.dat.cm.Hmisc,
    imputed.dat.mi.avg = Reduce("+",imputed.dat.mi)/length(imputed.dat.mi),
    imputed.dat.miq.avg = Reduce("+",imputed.dat.miq)/length(imputed.dat.miq),
    imputed.dat.miy.avg = Reduce("+",imputed.dat.miy)/length(imputed.dat.miy),
    dat.new = dat.new,
    dat.mean = dat.cc.mean,
    dat.cond.mean =  dat.cond.mean
  )
  
  
}

#Simulation
simulation <- function(n.tr,
                       mu.z,
                       mu.w,
                       mu.q,
                       s.z,
                       s.w,
                       s.q,
                       s.wz,
                       s.zq,
                       s.wq,
                       p.miss,
                       missing.type,
                       b.x,
                       b.z,
                       b.w,
                       mu.e,
                       s.e,
                       betaM,
                       n.tr.new,
                       mu.z.new,
                       mu.w.new,
                       mu.q.new,
                       s.z.new,
                       s.w.new,
                       s.q.new,
                       s.wz.new,
                       s.zq.new,
                       s.wq.new,
                       p.miss.new,
                       missing.type.new,
                       b.x.new,
                       b.z.new,
                       b.w.new,
                       mu.e.new,
                       s.e.new,
                       betaM.new,
                       n.imp,
                       n.sim,
                       Y.PM
                       ){
results = list(is.mse = data.frame(mse.truth = rep(NA,n.sim), mse.truth.int = rep(NA,n.sim), 
                                   mse.cc = rep(NA,n.sim), mse.full.cm.mipack = rep(NA,n.sim), 
                                   mse.full.cm.Hmisc = rep(NA, n.sim), mse.marg.cm.mipack = rep(NA, n.sim), 
                                   mse.marg.cm.Hmisc = rep(NA, n.sim), mse.full.mi.Hmisc = rep(NA,n.sim),
                                   mse.marg.mi.Hmisc = rep(NA, n.sim), mse.full.miy.Hmisc = rep(NA, n.sim), 
                                   mse.marg.miy.Hmisc = rep(NA, n.sim), mse.full.miq.Hmisc = rep(NA, n.sim),
                                   mse.marg.miq.Hmisc = rep(NA, n.sim), mse.full.miyq.Hmisc = rep(NA, n.sim), 
                                   mse.marg.miyq.Hmisc = rep(NA, n.sim), 
                                   mse.full.cond.mean = rep(NA, n.sim), mse.marg.cond.mean = rep(NA, n.sim)),
               is.pattern.00 = data.frame(oracle.marg.mse = rep(NA, n.sim), oracle.marg.prop = rep(NA, n.sim),
                                           oracle.full.mse = rep(NA, n.sim), oracle.full.prop = rep(NA, n.sim),
                                           cc.mean.mse = rep(NA, n.sim),cc.mean.prop = rep(NA, n.sim),
                                           full.cm.mipack.mse = rep(NA, n.sim),full.cm.mipack.prop = rep(NA, n.sim),
                                           marg.cm.mipack.mse = rep(NA, n.sim),marg.cm.mipack.prop = rep(NA, n.sim),
                                           full.cm.Hmisc.mse = rep(NA, n.sim), full.cm.Hmisc.prop = rep(NA, n.sim),
                                           marg.cm.Hmisc.mse = rep(NA, n.sim),marg.cm.Hmisc.prop = rep(NA, n.sim),
                                           full.mi.Hmisc.mse = rep(NA, n.sim),full.mi.Hmisc.prop = rep(NA, n.sim),
                                           marg.mi.Hmisc.mse = rep(NA, n.sim),marg.mi.Hmisc.prop = rep(NA, n.sim),
                                           full.miq.Hmisc.mse = rep(NA, n.sim),full.miq.Hmisc.prop = rep(NA, n.sim),
                                           marg.miq.Hmisc.mse = rep(NA, n.sim),marg.miq.Hmisc.prop = rep(NA, n.sim),
                                           full.miy.Hmisc.mse = rep(NA, n.sim),full.miy.Hmisc.prop = rep(NA, n.sim),
                                           marg.miy.Hmisc.mse = rep(NA, n.sim) ,marg.miy.Hmisc.prop = rep(NA, n.sim) ,
                                           metric.pmks.mse = rep(NA, n.sim),metric.pmks.prop = rep(NA, n.sim),
                                           metric.ccsm.mse = rep(NA, n.sim),metric.ccsm.prop = rep(NA, n.sim),
                                          full.cond.mean = rep(NA, n.sim), full.cond.mean.prop = rep(NA, n.sim),
                                          marg.cond.mean = rep(NA,n.sim), marg.cond.mean = rep(NA, n.sim)),
               is.pattern.10 = data.frame(oracle.marg.mse = rep(NA, n.sim), oracle.marg.prop = rep(NA, n.sim),
                                           oracle.full.mse = rep(NA, n.sim), oracle.full.prop = rep(NA, n.sim),
                                           cc.mean.mse = rep(NA, n.sim),cc.mean.prop = rep(NA, n.sim),
                                           full.cm.mipack.mse = rep(NA, n.sim),full.cm.mipack.prop = rep(NA, n.sim),
                                           marg.cm.mipack.mse = rep(NA, n.sim),marg.cm.mipack.prop = rep(NA, n.sim),
                                           full.cm.Hmisc.mse = rep(NA, n.sim), full.cm.Hmisc.prop = rep(NA, n.sim),
                                           marg.cm.Hmisc.mse = rep(NA, n.sim),marg.cm.Hmisc.prop = rep(NA, n.sim),
                                           full.mi.Hmisc.mse = rep(NA, n.sim),full.mi.Hmisc.prop = rep(NA, n.sim),
                                           marg.mi.Hmisc.mse = rep(NA, n.sim),marg.mi.Hmisc.prop = rep(NA, n.sim),
                                           full.miq.Hmisc.mse = rep(NA, n.sim),full.miq.Hmisc.prop = rep(NA, n.sim),
                                           marg.miq.Hmisc.mse = rep(NA, n.sim),marg.miq.Hmisc.prop = rep(NA, n.sim),
                                           full.miy.Hmisc.mse = rep(NA, n.sim),full.miy.Hmisc.prop = rep(NA, n.sim),
                                           marg.miy.Hmisc.mse = rep(NA, n.sim), marg.miy.Hmisc.prop = rep(NA, n.sim),
                                           metric.pmks.mse = rep(NA, n.sim),metric.pmks.prop = rep(NA, n.sim),
                                           metric.ccsm.mse = rep(NA, n.sim),metric.ccsm.prop = rep(NA, n.sim),
                                          full.cond.mean = rep(NA, n.sim), full.cond.mean.prop = rep(NA, n.sim),
                                          marg.cond.mean = rep(NA,n.sim), marg.cond.mean.prop = rep(NA, n.sim)),
               oos.mse = data.frame(mse.cc.mean.new= rep(NA,n.sim),
                                    mse.cc.cm.Hmisc.new= rep(NA,n.sim),
                                    mse.full.cm.mipack.new= rep(NA,n.sim),
                                    mse.marg.cm.mipack.new= rep(NA,n.sim),
                                    mse.full.cm.Hmisc.new= rep(NA,n.sim),
                                    mse.marg.cm.Hmisc.new= rep(NA,n.sim),
                                    mse.full.mi.Hmisc.new= rep(NA,n.sim),
                                    mse.marg.mi.Hmisc.new= rep(NA,n.sim),
                                    mse.full.miq.Hmisc.new= rep(NA,n.sim),
                                    mse.marg.miq.Hmisc.new= rep(NA,n.sim), 
                                    mse.full.miy.Hmisc.new= rep(NA,n.sim),
                                    mse.marg.miy.Hmisc.new= rep(NA,n.sim),
                                    oracle.marg= rep(NA,n.sim),
                                    oracle.full= rep(NA,n.sim),
                                    mse.full.cond.mean = rep(NA,n.sim), 
                                    mse.marg.cond.mean = rep(NA, n.sim)
                                   ),
               oos.pattern.00 = data.frame(oracle.marg.mse = rep(NA, n.sim), oracle.marg.prop = rep(NA, n.sim),
                                           oracle.full.mse = rep(NA, n.sim), oracle.full.prop = rep(NA, n.sim),
                                           cc.mean.mse = rep(NA, n.sim),cc.mean.prop = rep(NA, n.sim),
                                           cc.cm.mipack.mse = rep(NA, n.sim),cc.cm.mipack.prop = rep(NA, n.sim),
                                           cc.cm.Hmisc.mse = rep(NA, n.sim),cc.cm.Hmisc.prop = rep(NA, n.sim),
                                           full.cm.mipack.mse = rep(NA, n.sim),full.cm.mipack.prop = rep(NA, n.sim),
                                           marg.cm.mipack.mse = rep(NA, n.sim),marg.cm.mipack.prop = rep(NA, n.sim),
                                           full.cm.Hmisc.mse = rep(NA, n.sim), full.cm.Hmisc.prop = rep(NA, n.sim),
                                           marg.cm.Hmisc.mse = rep(NA, n.sim),marg.cm.Hmisc.prop = rep(NA, n.sim),
                                           full.mi.Hmisc.mse = rep(NA, n.sim),full.mi.Hmisc.prop = rep(NA, n.sim),
                                           marg.mi.Hmisc.mse = rep(NA, n.sim),marg.mi.Hmisc.prop = rep(NA, n.sim),
                                           full.miq.Hmisc.mse = rep(NA, n.sim),full.miq.Hmisc.prop = rep(NA, n.sim),
                                           marg.miq.Hmisc.mse = rep(NA, n.sim),marg.miq.Hmisc.prop = rep(NA, n.sim),
                                           full.miy.Hmisc.mse = rep(NA, n.sim),full.miy.Hmisc.prop = rep(NA, n.sim),
                                           marg.miy.Hmisc.mse = rep(NA, n.sim) ,marg.miy.Hmisc.prop = rep(NA, n.sim) ,
                                          metric.pmks.mse = rep(NA, n.sim),metric.pmks.prop = rep(NA, n.sim),
                                           metric.ccsm.mse = rep(NA, n.sim),metric.ccsm.prop = rep(NA, n.sim),
                                          full.cond.mean = rep(NA, n.sim), full.cond.mean.prop = rep(NA, n.sim),
                                          marg.cond.mean = rep(NA,n.sim), marg.cond.mean.prop = rep(NA, n.sim)),
               oos.pattern.10 = data.frame(oracle.marg.mse = rep(NA, n.sim), oracle.marg.prop = rep(NA, n.sim),
                                           oracle.full.mse = rep(NA, n.sim), oracle.full.prop = rep(NA, n.sim),
                                           cc.mean.mse = rep(NA, n.sim),cc.mean.prop = rep(NA, n.sim),
                                           cc.cm.mipack.mse = rep(NA, n.sim),cc.cm.mipack.prop = rep(NA, n.sim),
                                           cc.cm.Hmisc.mse = rep(NA, n.sim),cc.cm.Hmisc.prop = rep(NA, n.sim),
                                           full.cm.mipack.mse = rep(NA, n.sim),full.cm.mipack.prop = rep(NA, n.sim),
                                           marg.cm.mipack.mse = rep(NA, n.sim),marg.cm.mipack.prop = rep(NA, n.sim),
                                           full.cm.Hmisc.mse = rep(NA, n.sim), full.cm.Hmisc.prop = rep(NA, n.sim),
                                           marg.cm.Hmisc.mse = rep(NA, n.sim),marg.cm.Hmisc.prop = rep(NA, n.sim),
                                           full.mi.Hmisc.mse = rep(NA, n.sim),full.mi.Hmisc.prop = rep(NA, n.sim),
                                           marg.mi.Hmisc.mse = rep(NA, n.sim),marg.mi.Hmisc.prop = rep(NA, n.sim),
                                           full.miq.Hmisc.mse = rep(NA, n.sim),full.miq.Hmisc.prop = rep(NA, n.sim),
                                           marg.miq.Hmisc.mse = rep(NA, n.sim),marg.miq.Hmisc.prop = rep(NA, n.sim),
                                           full.miy.Hmisc.mse = rep(NA, n.sim),full.miy.Hmisc.prop = rep(NA, n.sim),
                                           marg.miy.Hmisc.mse = rep(NA, n.sim), marg.miy.Hmisc.prop = rep(NA, n.sim),
                                           metric.pmks.mse = rep(NA, n.sim),metric.pmks.prop = rep(NA, n.sim),
                                           metric.ccsm.mse = rep(NA, n.sim),metric.ccsm.prop = rep(NA, n.sim),
                                           full.cond.mean = rep(NA, n.sim), full.cond.mean.prop = rep(NA, n.sim),
                                           marg.cond.mean = rep(NA,n.sim), marg.cond.mean.prop = rep(NA, n.sim)),
               oos.z.mse = data.frame(dat.miss.new= rep(NA, n.sim),
                                      imputed.dat.cm.mipack = rep(NA, n.sim),
                                      imputed.dat.cm.Hmisc = rep(NA, n.sim), 
                                      imputed.dat.mi.avg = rep(NA, n.sim),
                                      imputed.dat.miq.avg = rep(NA, n.sim),
                                      imputed.dat.miy.avg = rep(NA, n.sim), 
                                      dat.new = rep(NA, n.sim),
                                      dat.mean = rep(NA, n.sim),
                                      dat.cond.mean = rep(NA, n.sim)     )
               )
#Replicate Simulation
for(RS in 1:n.sim){
  dat = dat.frame(n.tr=n.tr,
                  mu.z = mu.z,
                  mu.w = mu.w,
                  mu.q = mu.q,
                  s.z = s.z,
                  s.w = s.w,
                  s.q = s.q,
                  s.wz = s.wz,
                  s.zq = s.zq,
                  s.wq = s.wq,
                  p.miss = p.miss,
                  missing.type = missing.type,
                  b.x = b.x,
                  b.z = b.z,
                  b.w = b.w,
                  b.m = betaM,
                  mu.e = mu.e,
                  s.e = s.e,
                  betaM = betaM)
  
  if(Y.PM == TRUE) {dat$y <- dat[,'y.pmm'] }
   
  #Full Model that we could use if we had all the data-we will use this model to get the 'truth'
  mod.truth <- lm(y ~ z + w, data=dat)
  mod.truth.int <- lm(y ~ (z + w)*m, data=dat)

  #From those people make some of their information missing
  dat.miss <- dat
  dat.miss$z <- ifelse(dat$m==1,NA,dat$z)
 
  datdis <<- datadist(dat.miss)
  options(datadist='datdis')
  
  #complete case model    
  mod.cc <- lm(y ~  z + w, data=dat.miss)
  mod.ccsm.w <- lm(y ~ w, data=dat.miss)
  
  #############################################################################
  #Conditional Mean Imputations and Results Using both the MI package and Hmisc
  #############################################################################
  
  #mi package
  mdf <- suppressMessages(suppressWarnings(missing_data.frame(dat.miss[,c('x','z','w')])))
  #the expectation argument is conditional mean expectation but the underlying model is Bayesian
  mdf <-  suppressMessages(change(mdf, y=c('x','z','w'),what = "imputation_method", to = "expectation"))
  imputations <- mi(mdf,n.iter=1,n.chains=1,verbose=FALSE,parallel=FALSE)
  dat.miss.cm.mipack <- dat.miss
  dat.miss.cm.mipack[,c('x','z','w')] <- complete(imputations,1)[,c('x','z','w')]
  
  #regression imputation from Hmisc
  f <- aregImpute( ~ w + z, data=dat.miss, n.impute=1,type='regression',pr=FALSE)
  dat.miss.cm.Hmisc <- dat.miss
  dat.miss.cm.Hmisc[,c('z')] <- as.data.frame(impute.transcan(
    f, imputation=1, data=dat.miss, 
    list.out=TRUE, pr=FALSE, 
    check=FALSE))[,c('z')]
  
  ##My conditional mean imputation with no error
  dat.miss.cond.mean <- dat.miss
  mod.pred.z <- lm(z ~ w, data = dat.miss)
  dat.miss.cond.mean[ dat.miss.cond.mean$m==1,'z'] <- predict(mod.pred.z,  dat.miss.cond.mean[ dat.miss.cond.mean$m==1,])
  mod.full.cond.mean <- lm(y ~ (z + w)*m, data=dat.miss.cond.mean)
  mod.marg.cond.mean <- lm(y ~ z + w, data=dat.miss.cond.mean)
  
  #Fit model with the imputed conditional means indicator for missingness
  mod.full.cm.mipack <- lm(y ~ (z + w)*m, data=dat.miss.cm.mipack)
  mod.full.cm.Hmisc <- lm(y ~ (z + w)*m, data=dat.miss.cm.Hmisc)
  mod.marg.cm.mipack <- lm(y ~  z + w, data=dat.miss.cm.mipack)
  mod.marg.cm.Hmisc <- lm(y ~ z + w, data=dat.miss.cm.Hmisc)
  
  ### PMKS
  mod.zmiss <- lm(y ~ w, data=dat.miss[dat.miss$m==1,])
  pmks.mod <- pmks(DATA = dat.miss, model = "y ~ z + w", logistic = FALSE)
  ccsm.mod <- ccsm(DATA = dat.miss, model = "y ~ z + w", logistic = FALSE)
  
  #Missing data pattern percentages
  mpp <- (table(mdf@patterns)/n.tr)
  
  #####################################################
  #Mulitple Imputation Using a Congenial Model with pmm
  #####################################################
  imp.congenial <- aregImpute( ~  w + z, data=dat.miss, n.impute=10,pr=FALSE)
  mod.full.mi.Hmisc <- fit.mult.impute(y ~ (z + w)*m, ols,imp.congenial,data=dat.miss,
                                       pr=FALSE)
  mod.marg.mi.Hmisc <- fit.mult.impute(y ~ (z + w), ols,imp.congenial,data=dat.miss,pr=FALSE)
  
  ###################################################
  #Mulitple Imputation Using a Congenial Model with Y
  ###################################################
  
  imp.y <- aregImpute( ~ w + z + y, data=dat.miss, n.impute=10,pr=FALSE)
  mod.full.miy.Hmisc <- fit.mult.impute(y ~ (z + w)*m, ols,imp.y,data=dat.miss,pr=FALSE)
  mod.marg.miy.Hmisc <- fit.mult.impute(y ~ (z + w), ols,imp.y,data=dat.miss,pr=FALSE)
  
  #######################################################################
  #Mulitple Imputation Using an Imputation Model with extra information q
  ########################################################################
  
  imp.q <- aregImpute( ~ w + z + q, data=dat.miss, n.impute=10,pr=FALSE)
  mod.full.miq.Hmisc <- fit.mult.impute(y ~ (z + w)*m, ols,imp.q,data=dat.miss,pr=FALSE)
  mod.marg.miq.Hmisc <- fit.mult.impute(y ~ (z + w), ols,imp.q,data=dat.miss,pr=FALSE)
  
  ##############################################################################
  #Mulitple Imputation Using an Imputation Model with extra information y and q
  ##############################################################################
  
  imp.yq <- aregImpute( ~  w + z + q + y, data=dat.miss, n.impute=10,pr=FALSE)
  mod.full.miyq.Hmisc <- fit.mult.impute(y ~ (z + w)*m, ols,imp.yq,data=dat.miss,pr=FALSE)
  mod.marg.miyq.Hmisc <- fit.mult.impute(y ~ (z + w), ols,imp.yq,data=dat.miss,pr=FALSE)
  
  
  
  ###########################
  # In Sample Pattern Metrics
  ###########################
  
  pattern.is <- which.pattern(dat.miss, model = 'y ~ z + w')
  
  oracle.marg.pattern.is <- get_metric_linear(data.frame(predict(mod.truth, dat),
                                                         dat[,'y'],
                                                         pattern.is))
  
  oracle.full.pattern.is <- get_metric_linear(data.frame(predict(mod.truth.int, dat),
                                                         dat[,'y'],
                                                         pattern.is))
  
  dat.cc.pred <- dat.miss
  dat.cc.pred[which(is.na(dat.cc.pred$z)==TRUE),'z'] <- mean(dat.miss$z, na.rm=TRUE) 
  
  cc.mean.pattern.is <- get_metric_linear(data.frame( predict(mod.cc, dat.cc.pred),
                                                         dat.miss[,'y'],
                                                         pattern.is))
  
  full.cond.mean.pattern.is <- get_metric_linear(data.frame(predict(mod.full.cond.mean),
                                                            dat.miss[,'y'],
                                                            pattern.is))
    
  marg.cond.mean.pattern.is <- get_metric_linear(data.frame(predict(mod.marg.cond.mean),
                                                            dat.miss[,'y'],
                                                            pattern.is))
  
  full.cm.mipack.pattern.is <- get_metric_linear(data.frame(predict(mod.full.cm.mipack),
                                                            dat.miss[,'y'],
                                                            pattern.is))
  
  marg.cm.mipack.pattern.is <- get_metric_linear(data.frame(predict(mod.marg.cm.mipack),
                                                            dat.miss[,'y'],
                                                            pattern.is))
  
  #CM with Hmisc				     
  
  full.cm.Hmisc.pattern.is <- get_metric_linear(data.frame(predict(mod.full.cm.Hmisc),
                                                           dat.miss[,'y'],
                                                           pattern.is))
  
  marg.cm.Hmisc.pattern.is <- get_metric_linear(data.frame(predict(mod.marg.cm.Hmisc),
                                                           dat.miss[,'y'],
                                                           pattern.is))
  
  #MI predictions
  ###
  
  full.mi.Hmisc.pattern.is <- get_metric_linear(data.frame(predict(mod.full.mi.Hmisc),
                                                           dat.miss[,'y'],
                                                           pattern.is))
  
  marg.mi.Hmisc.pattern.is <- get_metric_linear(data.frame(predict(mod.marg.mi.Hmisc),
                                                           dat.miss[,'y'],
                                                           pattern.is))
  
  
  #MI q predictions
  
  full.miq.Hmisc.pattern.is <- get_metric_linear(data.frame(predict(mod.full.miq.Hmisc),
                                                            dat.miss[,'y'],
                                                            pattern.is))
  
  marg.miq.Hmisc.pattern.is <- get_metric_linear(data.frame(predict(mod.marg.miq.Hmisc),
                                                            dat.miss[,'y'],
                                                            pattern.is))
  
  
  #MI y predictions
  full.miy.Hmisc.pattern.is <- get_metric_linear(data.frame(predict(mod.full.miy.Hmisc),
                                                            dat.miss[,'y'],
                                                            pattern.is))
  
  marg.miy.Hmisc.pattern.is <- get_metric_linear(data.frame(predict(mod.marg.miy.Hmisc),
                                                            dat.miss[,'y'],
                                                            pattern.is))
  
  #Submodel predictions
  is.pmks <- predict.sm(dat.miss, 'y ~ z + w', pmks.object = pmks.mod, logistic = FALSE)
  is.ccsm <- predict.sm(dat.miss, 'y ~ z + w', pmks.object = ccsm.mod, logistic = FALSE)
  metric.pmks.is <- data.frame(mse = c(is.pmks[[1]]$strat.brier, is.pmks[[2]]$strat.brier), 
                               proportion = c(length(is.pmks[[1]]$truth)/n.tr, length(is.pmks[[2]]$truth)/n.tr))
  rownames(metric.pmks.is) <- c(is.pmks[[1]][1], is.pmks[[2]][1])
  metric.ccsm.is <- data.frame(mse = c(is.ccsm[[1]]$strat.brier, is.ccsm[[2]]$strat.brier), 
                               proportion = c(length(is.ccsm[[1]]$truth)/n.tr, length(is.ccsm[[2]]$truth)/n.tr))
  rownames(metric.ccsm.is) <- c(is.ccsm[[1]][1], is.ccsm[[2]][1])
  
  
  
  ###############################################################################################
  #Get new sample, make people missing, imput their missing data, and then get model predictions
  # - Assume 1 by 1 predictions
  ###############################################################################################
  
  new.sample <- oos.data(n.tr.new,
                         mu.z.new,
                         mu.w.new,
                         mu.q.new,
                         s.z.new,
                         s.w.new,
                         s.q.new,
                         s.wz.new,
                         s.zq.new,
                         s.wq.new,
                         p.miss.new,
                         missing.type.new,
                         b.x.new,
                         b.z.new,
                         b.w.new,
                         mu.e.new,
                         s.e.new,
                         betaM.new=betaM.new,
                         n.imp = n.imp,
                         original.dat.miss = dat.miss,
                         Y.PM = Y.PM
                         )
  
  
  
  pattern <- which.pattern(new.sample[["dat.miss.new"]], model = 'y ~ z + w')
  
  #The oracle model
  oracle.marg <- mean((new.sample[["dat.new"]][,'y'] -
                         predict(mod.truth, new.sample[["dat.new"]]))^2)
  
  oracle.full <- mean((new.sample[['dat.new']][,'y']  - 
                         predict(mod.truth.int, new.sample[["dat.new"]]))^2)
  
  oracle.marg.pattern <- get_metric_linear(data.frame(predict(mod.truth, new.sample[["dat.new"]]),
                               new.sample[["dat.new"]][,'y'],
                               pattern))
  
  oracle.full.pattern <- get_metric_linear(data.frame(predict(mod.truth.int, new.sample[["dat.new"]]),
                                                      new.sample[["dat.new"]][,'y'],
                                                      pattern))
                               
  
  #CM with MI pack on the complete case model
  mse.cc.mean.new <- mean((new.sample[["dat.mean"]][,'y'] - 
                                  predict(mod.cc,new.sample[["dat.mean"]]))^2)
  
  mse.cc.mean.new.pattern <- get_metric_linear(data.frame( predict(mod.cc,new.sample[["dat.mean"]]),
                                                           new.sample[["dat.mean"]][,'y'],
                                                           pattern))
  
  mse.cc.cm.Hmisc.new <- mean((new.sample[["imputed.dat.cm.Hmisc"]][,'y'] - 
                                 predict(mod.cc,new.sample[["imputed.dat.cm.Hmisc"]]))^2)
  
  cc.cm.mipack.pattern <- get_metric_linear(data.frame( predict(mod.cc,new.sample[["imputed.dat.cm.mipack"]]),
                                                       new.sample[["imputed.dat.cm.mipack"]][,'y'],
                                                      pattern))
  
  cc.cm.Hmisc.pattern <- get_metric_linear(data.frame( predict(mod.cc,new.sample[["imputed.dat.cm.Hmisc"]]),
                                                       new.sample[["imputed.dat.cm.Hmisc"]][,'y'],
                                                        pattern))
  #Strict Conditional Mean
  
  mse.full.cond.mean.new <- mean((new.sample[["dat.cond.mean"]][,'y'] -
                                   predict(mod.full.cond.mean,new.sample[["dat.cond.mean"]]))^2)
    
    
  mse.marg.cond.mean.new <- mean((new.sample[["dat.cond.mean"]][,'y'] -
                                    predict(mod.marg.cond.mean,new.sample[["dat.cond.mean"]]))^2)
  
  full.cond.mean.pattern <- get_metric_linear(data.frame(predict(mod.full.cond.mean,new.sample[["dat.cond.mean"]]),
                                                         new.sample[["dat.cond.mean"]][,'y'],
                                                         pattern))
    
  marg.cond.mean.pattern <- get_metric_linear(data.frame(predict(mod.marg.cond.mean,new.sample[["dat.cond.mean"]]),
                                                         new.sample[["dat.cond.mean"]][,'y'],
                                                         pattern))
  #CM with MI pack			
  
  mse.full.cm.mipack.new <-  mean((new.sample[["imputed.dat.cm.mipack"]][,'y'] -
                                     predict(mod.full.cm.mipack,new.sample[["imputed.dat.cm.mipack"]]))^2)
  
  mse.marg.cm.mipack.new <-  mean((new.sample[["imputed.dat.cm.mipack"]][,'y'] - 
                                     predict(mod.marg.cm.mipack,new.sample[["imputed.dat.cm.mipack"]]))^2)
  
  
  full.cm.mipack.pattern <- get_metric_linear(data.frame(predict(mod.full.cm.mipack,new.sample[["imputed.dat.cm.mipack"]]),
                                                        new.sample[["imputed.dat.cm.mipack"]][,'y'],
                                                        pattern))
  
  marg.cm.mipack.pattern <- get_metric_linear(data.frame(predict(mod.marg.cm.mipack,new.sample[["imputed.dat.cm.mipack"]]),
                                                         new.sample[["imputed.dat.cm.mipack"]][,'y'],
                                                         pattern))
  
  #CM with Hmisc				     
  mse.full.cm.Hmisc.new <- mean((new.sample[["imputed.dat.cm.Hmisc"]][,'y'] - 	
                                   predict(mod.full.cm.Hmisc,new.sample[["imputed.dat.cm.Hmisc"]]))^2)
  
  mse.marg.cm.Hmisc.new <- mean((new.sample[["imputed.dat.cm.Hmisc"]][,'y'] - 
                                   predict(mod.marg.cm.Hmisc,new.sample[["imputed.dat.cm.Hmisc"]]))^2)
  
  full.cm.Hmisc.pattern <- get_metric_linear(data.frame(predict(mod.full.cm.Hmisc,new.sample[["imputed.dat.cm.Hmisc"]]),
                                                        new.sample[["imputed.dat.cm.Hmisc"]][,'y'],
                                                         pattern))
  
  marg.cm.Hmisc.pattern <- get_metric_linear(data.frame(predict(mod.marg.cm.Hmisc,new.sample[["imputed.dat.cm.Hmisc"]]),
                                                        new.sample[["imputed.dat.cm.Hmisc"]][,'y'],
                                                        pattern))
  
  #MI predictions
  ###
  
  mse.full.mi.Hmisc.new <-  mean((new.sample[["imputed.dat.mi.avg"]][,'y'] - 
                                    predict(mod.full.mi.Hmisc,new.sample[["imputed.dat.mi.avg"]]))^2)
  
  mse.marg.mi.Hmisc.new <-  mean((new.sample[["imputed.dat.mi.avg"]][,'y'] - 
                                    predict(mod.marg.mi.Hmisc,new.sample[["imputed.dat.mi.avg"]]))^2)
  
  full.mi.Hmisc.pattern <- get_metric_linear(data.frame(predict(mod.full.mi.Hmisc,new.sample[["imputed.dat.mi.avg"]]),
                                                        new.sample[["imputed.dat.mi.avg"]][,'y'],
                                                        pattern))
  
  marg.mi.Hmisc.pattern <- get_metric_linear(data.frame(predict(mod.marg.mi.Hmisc,new.sample[["imputed.dat.mi.avg"]]),
                                                        new.sample[["imputed.dat.mi.avg"]][,'y'],
                                                        pattern))
  
  
  #MI q predictions
  mse.full.miq.Hmisc.new <- mean((new.sample[["imputed.dat.miq.avg"]][,'y'] - 
                                    predict(mod.full.miq.Hmisc,new.sample[["imputed.dat.miq.avg"]]))^2)
  
  mse.marg.miq.Hmisc.new <- mean((new.sample[["imputed.dat.miq.avg"]][,'y'] - 
                                    predict(mod.marg.miq.Hmisc,new.sample[["imputed.dat.miq.avg"]]))^2)
  
  full.miq.Hmisc.pattern <- get_metric_linear(data.frame(predict(mod.full.miq.Hmisc,new.sample[["imputed.dat.miq.avg"]]),
                                                         new.sample[["imputed.dat.miq.avg"]][,'y'],
                                                        pattern))
  
  marg.miq.Hmisc.pattern <- get_metric_linear(data.frame(predict(mod.marg.miq.Hmisc,new.sample[["imputed.dat.miq.avg"]]),
                                                         new.sample[["imputed.dat.miq.avg"]][,'y'],
                                                         pattern))
  
  
  #MI y predictions
  mse.full.miy.Hmisc.new <- mean((new.sample[["imputed.dat.miy.avg"]][,'y'] - 
                                    predict(mod.full.miy.Hmisc,new.sample[["imputed.dat.miy.avg"]]))^2)
  
  mse.marg.miy.Hmisc.new <- mean((new.sample[["imputed.dat.miy.avg"]][,'y'] - 
                                    predict(mod.marg.miy.Hmisc,new.sample[["imputed.dat.miy.avg"]]))^2)
  
  full.miy.Hmisc.pattern <- get_metric_linear(data.frame(predict(mod.full.miy.Hmisc,new.sample[["imputed.dat.miy.avg"]]),
                                                         new.sample[["imputed.dat.miy.avg"]][,'y'],
                                                         pattern))
  
  marg.miy.Hmisc.pattern <- get_metric_linear(data.frame(predict(mod.marg.miy.Hmisc,new.sample[["imputed.dat.miy.avg"]]),
                                                         new.sample[["imputed.dat.miy.avg"]][,'y'],
                                                         pattern))
  
  #Submodel predictions
  oos.pmks <- predict.sm(new.sample[["dat.miss.new"]], 'y ~ z + w', pmks.object = pmks.mod, logistic = FALSE)
  oos.ccsm <- predict.sm(new.sample[["dat.miss.new"]], 'y ~ z + w', pmks.object = ccsm.mod, logistic = FALSE)
  metric.pmks <- data.frame(mse = c(oos.pmks[[1]]$strat.brier, oos.pmks[[2]]$strat.brier), 
                            proportion = c(length(oos.pmks[[1]]$truth)/n.tr.new, length(oos.pmks[[2]]$truth)/n.tr.new))
  rownames(metric.pmks) <- c(oos.pmks[[1]][1], oos.pmks[[2]][1])
  metric.ccsm <- data.frame(mse = c(oos.ccsm[[1]]$strat.brier, oos.ccsm[[2]]$strat.brier), 
                            proportion = c(length(oos.ccsm[[1]]$truth)/n.tr.new, length(oos.ccsm[[2]]$truth)/n.tr.new))
  rownames(metric.ccsm) <- c(oos.ccsm[[1]][1], oos.ccsm[[2]][1])
  
  
  ###########################
  #  Storing all the results
  ###########################
 results[['is.mse']][RS,] <- c(
      #In sample MSE
     mse(mod.truth),
     mse(mod.truth.int),
     mse(mod.cc),
     mse(mod.full.cm.mipack),
     mse(mod.full.cm.Hmisc),
     mse(mod.marg.cm.mipack),
     mse(mod.marg.cm.Hmisc),
     mse(mod.full.mi.Hmisc),
     mse(mod.marg.mi.Hmisc),
     mse(mod.full.miy.Hmisc),
     mse(mod.marg.miy.Hmisc),
     mse(mod.full.miq.Hmisc),
     mse(mod.marg.miq.Hmisc),
     mse(mod.full.miyq.Hmisc),
     mse(mod.marg.miyq.Hmisc),
     mse(mod.full.cond.mean),
     mse(mod.marg.cond.mean)
     )
 
 results[['is.pattern.00']][RS,] <- c(
   oracle.marg.pattern.is[1,] ,
   oracle.full.pattern.is[1,],
   cc.mean.pattern.is[1,],
   full.cm.mipack.pattern.is[1,],
   marg.cm.mipack.pattern.is[1,], 
   full.cm.Hmisc.pattern.is[1,],
   marg.cm.Hmisc.pattern.is[1,], 
   full.mi.Hmisc.pattern.is[1,], 
   marg.mi.Hmisc.pattern.is[1,], 
   full.miq.Hmisc.pattern.is[1,], 
   marg.miq.Hmisc.pattern.is[1,], 
   full.miy.Hmisc.pattern.is[1,], 
   marg.miy.Hmisc.pattern.is[1,],
   metric.pmks.is[1,],
   metric.ccsm.is[1,],
   full.cond.mean.pattern.is[1,],
   marg.cond.mean.pattern.is[1,]
 )
 
 results[['is.pattern.10']][RS,] = c(
   oracle.marg.pattern.is[2,] ,
   oracle.full.pattern.is[2,],
   cc.mean.pattern.is[2,],
   full.cm.mipack.pattern.is[2,],
   marg.cm.mipack.pattern.is[2,], 
   full.cm.Hmisc.pattern.is[2,],
   marg.cm.Hmisc.pattern.is[2,], 
   full.mi.Hmisc.pattern.is[2,], 
   marg.mi.Hmisc.pattern.is[2,], 
   full.miq.Hmisc.pattern.is[2,], 
   marg.miq.Hmisc.pattern.is[2,], 
   full.miy.Hmisc.pattern.is[2,], 
   marg.miy.Hmisc.pattern.is[2,],
   metric.pmks.is[2,],
   metric.ccsm.is[2,],
   full.cond.mean.pattern.is[2,],
   marg.cond.mean.pattern.is[2,]) 
 
 #Out of sample MSE    
 results[['oos.mse']][RS,] <-  c(    
      mse.cc.mean.new,
      mse.cc.cm.Hmisc.new,
      mse.full.cm.mipack.new,
      mse.marg.cm.mipack.new,
      mse.full.cm.Hmisc.new,
      mse.marg.cm.Hmisc.new,
      mse.full.mi.Hmisc.new,
      mse.marg.mi.Hmisc.new,
      mse.full.miq.Hmisc.new,
      mse.marg.miq.Hmisc.new,
      mse.full.miy.Hmisc.new,
      mse.marg.miy.Hmisc.new,
      oracle.marg,
      oracle.full,
      mse.full.cond.mean.new,
      mse.marg.cond.mean.new
    )
    
 results[['oos.pattern.00']][RS,] <- c(
      oracle.marg.pattern[1,] ,
      oracle.full.pattern[1,],
      mse.cc.mean.new.pattern[1,],
      cc.cm.mipack.pattern[1,],
      cc.cm.Hmisc.pattern[1,],
      full.cm.mipack.pattern[1,],
      marg.cm.mipack.pattern[1,], 
      full.cm.Hmisc.pattern[1,],
      marg.cm.Hmisc.pattern[1,], 
      full.mi.Hmisc.pattern[1,], 
      marg.mi.Hmisc.pattern[1,], 
      full.miq.Hmisc.pattern[1,], 
      marg.miq.Hmisc.pattern[1,], 
      full.miy.Hmisc.pattern[1,], 
      marg.miy.Hmisc.pattern[1,],
      metric.pmks[1,],
      metric.ccsm[1,],
      full.cond.mean.pattern[1,],
      marg.cond.mean.pattern[1,]
      )
    
    results[['oos.pattern.10']][RS,] = c(
      oracle.marg.pattern[2,] ,
      oracle.full.pattern[2,],
      mse.cc.mean.new.pattern[2,],
      cc.cm.mipack.pattern[2,],
      cc.cm.Hmisc.pattern[2,],
      full.cm.mipack.pattern[2,],
      marg.cm.mipack.pattern[2,], 
      full.cm.Hmisc.pattern[2,],
      marg.cm.Hmisc.pattern[2,], 
      full.mi.Hmisc.pattern[2,], 
      marg.mi.Hmisc.pattern[2,], 
      full.miq.Hmisc.pattern[2,], 
      marg.miq.Hmisc.pattern[2,], 
      full.miy.Hmisc.pattern[2,], 
      marg.miy.Hmisc.pattern[2,],
      metric.pmks[2,],
      metric.ccsm[2,],
      full.cond.mean.pattern[2,],
      marg.cond.mean.pattern[2,])   
    
    new.sample.z <- data.frame(do.call(cbind,lapply(new.sample, function(z) cbind(z$z))))
    names(new.sample.z) <- names(new.sample)
    results[['oos.z.mse']][RS,] = c(apply(new.sample.z, 2,function(z) mean((z - new.sample.z$dat.new)^2)))
    
}
results
}



###############
#EPE Function
##############
#This function simulates the data for Figure 1 of the supplemetary material, comparing the Expected prediction error for the large and small model

#Install Libraries
library(MASS)
library(DMwR)
library(mice)
library(doSNOW)
library(rms)

#########################################
####  Training data & prediction interval functions
#########################################
expit <- function(x) exp(x)/(1 + exp(x))
logit <- function(p) log(p/(1-p))

miss.int <- function(miss.mech, betaM=1,p.miss){
  if(ncol(miss.mech)!=2){
    logit(p.miss) - betaM*mean(miss.mech[,1])
  } else {
    logit(p.miss) - betaM*mean((((miss.mech[,2]/sd(miss.mech[,2]))
                                 + miss.mech[,1])/as.numeric(sqrt(2*(1+cor(miss.mech[,2],miss.mech[,1]))))))
  }
}


dat.frame <- function(n.tr,mu.z,mu.w,s.z,s.w,s.wz,p.miss,
                      missing.type,
                      b.x,b.z,b.w,
                      mu.e,s.e,betaM=1){
  x = rep(1,n.tr)
  mu = c(mu.z,mu.w)
  Sigma = matrix(c(1,s.wz,s.wz,1), byrow=TRUE, ncol=2)
  dat =  as.data.frame(mvrnorm(n=n.tr, mu=mu, Sigma=Sigma, empirical=TRUE))
  colnames(dat) = c('z','w')
  coef.tru	= rbind( b.x , b.z , b.w)
  dat$y=cbind(x,dat[,'z'], dat[,'w'])%*%coef.tru + rnorm(n.tr,mean=mu.e,sd=s.e)

  if(missing.type=='MCAR'){
    m = rbinom(n.tr,1,p.miss)
    y.pmm = cbind(x,dat[,'z'], dat[,'w'])%*%coef.tru + betaM*m + betaM*dat[,'z']*m + betaM*dat[,'w']*m +  rnorm(n.tr,mean=mu.e,sd=s.e)
    data.frame(x=x, z=dat[,'z'], w=dat[,'w'], m=m, y=dat[,'y'], y.pmm = y.pmm)
  } else if(missing.type=='MAR'){
    m = rbinom(n.tr, 1, expit(miss.int(miss.mech=cbind(dat[,'w']),
                                       betaM=betaM,p.miss) + betaM*dat[,'w']))
    y.pmm = cbind(x,dat[,'z'], dat[,'w'])%*%coef.tru + betaM*m + betaM*dat[,'z']*m + betaM*dat[,'w']*m + rnorm(n.tr,mean=mu.e,sd=(s.e + 2*m))
    data.frame(x=x, z=dat[,1], w=dat[,2], m=m, y=dat[,'y'], y.pmm = y.pmm)
  } else if(missing.type=='MNAR'){
    m = rbinom(n.tr, 1, expit(miss.int(miss.mech=cbind(dat[,'z']),
                                       betaM=betaM,p.miss) + betaM*dat[,'z']))
    y.pmm = cbind(x,dat[,'z'], dat[,'w'])%*%coef.tru + betaM*m + betaM*dat[,'z']*m + betaM*dat[,'w']*m + rnorm(n.tr,mean=mu.e,sd=(s.e + 2*m))
    data.frame(x=x, z=dat[,'z'], w=dat[,'w'], m=m, y=dat[,'y'], y.pmm = y.pmm)
  } else if(missing.type=='MARY'){
    m = rbinom(n.tr, 1,
               expit(miss.int(miss.mech=dat[,c('w','y')], betaM=betaM,p.miss)
                     + betaM*(((dat[,'y']/sd(dat[,'y'])) + dat[,'w'])/as.numeric(sqrt(2*(1+cor(dat[,'y'],dat[,'w'])))))))
    data.frame(x=x, z=dat[,'z'], w=dat[,'w'], m=m, y=dat[,'y'])
  } else if(missing.type=='MNARY'){
    m = rbinom(n.tr, 1,
               expit(miss.int(miss.mech=dat[,c('z','y')], betaM=betaM,p.miss)
                     + betaM*(((dat[,'y']/sd(dat[,'y'])) + dat[,'z'])/as.numeric(sqrt(2*(1+cor(dat[,'y'],dat[,'z'])))))))
    data.frame(x=x, z=dat[,'z'], w=dat[,'w'], m=m, y=dat[,'y'])
  }

}

sims=1000;
n.y.tr = 100;
n.tr=10000;
b.x=0.5; b.z = 3; b.w = 3;
mu.z=3;
mu.w=3;
s.z=1;
s.w=1;
s.wz=0.5;
p.miss=0.5;
missing.type='MCAR';
lo.z=-1; hi.z=7;
acc=0.001;
mu.e=0;
s.e=1;
mu.cond.imp=0;
s.cond.imp=1;
imputed='conditional';
orig.imp='includey';
y.lo=7; y.hi=15;
b.m=1
missing.type.orig='MCAR';
PMMY = FALSE


EPE.T 	= function(sims=10000,
                  n.y.tr = 10,
                  n.tr=50,
                  b.x=bx, b.z = bz, b.w = bw,
                  mu.z=mean.z,
                  mu.w=mean.w,
                  s.z=sd.z,
                  s.w=sd.w,
                  s.wz=sd.wz,
                  p.miss,
                  missing.type.orig,
                  missing.type,
                  lo.z=-1, hi.z=7,
                  acc=0.01,
                  mu.e=mean.e,
                  s.e=sd.e,
                  mu.cond.imp,
                  s.cond.imp,
                  imputed,
                  orig.imp,
                  y.lo=7, y.hi=15,
                  b.m = b.m,
                  PMMY = FALSE) {

  coef.tru	= rbind( b.x , b.z , b.w)

  ## Predictors (Intercept in X)
  xz.dat 		= replicate(n.y.tr,
                       dat.frame(n.tr = n.tr,
                                 mu.z = mu.z,
                                 mu.w = mu.w,
                                 s.z = s.z,
                                 s.wz = s.wz,
                                 p.miss = p.miss,
                                 missing.type = missing.type.orig,
                                 b.x = b.x,
                                 b.z = b.z,
                                 b.w = b.w,
                                 mu.e = mu.e,
                                 s.e = s.e,
                                 betaM=b.m
                       ),
                       simplify=FALSE)

  data.tr.list = lapply(xz.dat, function(q) {data.frame(x=q$x, z=q$z,
                                                        w=q$w, m=q$m,
                                                        y = q$y,
                                                        y.pmm = q$y.pmm)})


  epes =	lapply(data.tr.list, function(q) {

    if(PMMY==FALSE){
      y.tr      	= q$y
    } else {
      y.tr= q$y.pmm}
    x.tr      	= q$x
    z.tr.truth	= q$z
    w.tr		= q$w
    m.z.tr      = q$m
    z.tr        = ifelse(q$m==1,NA,q$z)
    observed.miss    = sum(m.z.tr)/n.tr

    #Complete Case Design Matrix
    Z.cc = cbind(x.tr[m.z.tr==0],z.tr[m.z.tr==0],w.tr[m.z.tr==0],y.tr[m.z.tr==0])
    Z.orig.withMiss = cbind(x.tr,z.tr,w.tr,y.tr)

    #Design Matrix Imputation Model
    Z.s   	= cbind(x.tr[m.z.tr==0],w.tr[m.z.tr==0])
    bhat.z  = chol2inv(chol(t(Z.s)%*%Z.s))%*%t(Z.s)%*%z.tr[m.z.tr==0]

    #Get beta for a imputation model that includes y
    Z.y   	= cbind(x.tr[m.z.tr==0],w.tr[m.z.tr==0], y.tr[m.z.tr==0])
    bhat.zy  = chol2inv(chol(t(Z.y)%*%Z.y))%*%t(Z.y)%*%z.tr[m.z.tr==0]

    #Impute E[z ~ x + w] for those individuals missing
    if(orig.imp=='noy'){
      z.tr[is.na(z.tr)] = cbind(x.tr[m.z.tr==1],w.tr[m.z.tr==1])%*%bhat.z
    } else if(orig.imp=='includey'){
      z.tr[is.na(z.tr)] = cbind(x.tr[m.z.tr==1],w.tr[m.z.tr==1],y.tr[m.z.tr==1])%*%bhat.zy
    }


    ## Design Matrix
    if(PMMY==FALSE){
      X.truth = cbind(x.tr,z.tr.truth, w.tr)
    } else {
      X.truth = cbind(x.tr,z.tr.truth,w.tr,m.z.tr,z.tr.truth*m.z.tr,w.tr*m.z.tr)
    }



    X.l   	= cbind(x.tr,z.tr,w.tr)
    X.s   	= cbind(x.tr,w.tr)
    X.l.cc  = cbind(x.tr[m.z.tr==0],z.tr[m.z.tr==0],w.tr[m.z.tr==0])
    X.s.cc  = cbind(x.tr[m.z.tr==1], w.tr[m.z.tr==1])

    ## Least Squares estimates

    bhat.truth = chol2inv(chol(t(X.truth)%*%X.truth))%*%t(X.truth)%*%y.tr
    bhat.l     = chol2inv(chol(t(X.l)%*%X.l))%*%t(X.l)%*%y.tr
    bhat.s	   = chol2inv(chol(t(X.s)%*%X.s))%*%t(X.s)%*%y.tr
    bhat.l.cc  = chol2inv(chol(t(X.l.cc)%*%X.l.cc))%*%t(X.l.cc)%*%y.tr[m.z.tr==0]
    bhat.s.cc  = chol2inv(chol(t(X.s.cc)%*%X.s.cc))%*%t(X.s.cc)%*%y.tr[m.z.tr==1]

    #NOTE: bhat.s != bhat.s.cc because the lengths are different!

    ##################################################
    # This Section Starts the Part of the Simulation
    # code that includes the OUT OF SAMPLE population
    ##################################################

    ## Create Prediction Set
    z.pred 	  = seq(lo.z, hi.z, acc)
    #should I predict w based on z.pred
    w.pred	  = rnorm(length(z.pred),
                     mean=(mu.w + s.wz*(s.w/s.z)*(z.pred - mu.z)), sd = (1-s.wz^2)*s.w^2)
    x.pred 	  = rep(1,length(z.pred))
    n.pred	  = length(z.pred)

    z.pred.sim 	= rep(z.pred, each = sims)
    x.pred.sim	= rep(x.pred, each = sims)
    w.pred.sim = rep(w.pred, each = sims)
    ## Make new true responses (y) at x.new (sim replicates)


    m.pred    = if(missing.type=='MCAR'){
      rbinom(length(z.pred),1,p.miss)
    } else if(missing.type=='MAR'){
      rbinom(length(z.pred), 1,
             expit(miss.int(
               miss.mech=cbind(w.pred),betaM=b.m,p.miss=p.miss) + b.m*w.pred))
    } else if(missing.type=='MNAR'){
      rbinom(length(z.pred), 1,
             expit(miss.int(miss.mech=cbind(z.pred),
                            betaM=b.m,p.miss=p.miss) + b.m*z.pred))
    } else if(missing.type=='MARY'){
      rbinom(length(z.pred.sim), 1,
             expit(miss.int(miss.mech=cbind(w.pred.sim,y.new),
                            betaM=b.m,p.miss=p.miss)
                   + b.m*(((y.new/sd(y.new))
                           + w.pred.sim)/as.numeric(sqrt(2*(1+cor(y.new,w.pred.sim)))))))
    } else if(missing.type=='MNARY'){
      rbinom(length(z.pred.sim), 1,
             expit(miss.int(miss.mech=cbind(z.pred.sim,y.new),
                            betaM=b.m,p.miss=p.miss)
                   + b.m*(((y.new/sd(y.new))
                           + z.pred.sim)/as.numeric(sqrt(2*(1+cor(y.new,z.pred.sim)))))))
    }

    if(PMMY==FALSE){
      y.new = cbind(x.pred.sim,z.pred.sim,w.pred.sim)%*%coef.tru + rnorm(n.pred*sims,mean=mu.e,sd=s.e)
    } else {
      y.new = cbind(x.pred.sim,z.pred.sim,w.pred.sim)%*%coef.tru + b.m*rep(m.pred,each=sims) +  b.m*z.pred.sim*rep(m.pred,each=sims) + b.m*rep(m.pred,each=sims)*w.pred.sim + rnorm(n.pred*sims,mean=mu.e,sd=(s.e + 2*rep(m.pred,each=sims)))
    }


    Z.pred.miss = cbind(x.pred[m.pred==1],NA,w.pred[m.pred==1],NA)

    prop.pred.miss		 = 	sum(m.pred)/length(z.pred)
    n.pred.miss	         = 	length(z.pred) - sum(m.pred)
    n.pred.notmiss       = 	length(z.pred) - sum(abs(m.pred-1))
    n.pred.sim		  	 = 	length(z.pred.sim)
    n.pred				 =	length(z.pred)

    z.pred.sim.miss 		= rep(z.pred[m.pred==1], each = sims)
    z.pred.sim.notmiss 	= rep(z.pred[m.pred==0], each = sims)

    x.pred.sim.miss = rep(x.pred[m.pred==1], each = sims)
    x.pred.sim.notmiss = rep(x.pred[m.pred==0], each = sims)

    w.pred.sim.miss = rep(w.pred[m.pred==1], each = sims)
    w.pred.sim.notmiss = rep(w.pred[m.pred==0], each = sims)

    y.new.miss  	   = y.new[rep(m.pred,each=sims)==1]
    y.new.notmiss  = y.new[rep(m.pred,each=sims)==0]

    #Impute for the missing z's
    z.pred.imputed = z.pred
    z.pred.imputed.sim = z.pred.sim

    z.pred.imputed[m.pred==1] = cbind(x.pred[m.pred==1],w.pred[m.pred==1])%*%bhat.z +
      rnorm(sum(m.pred==1),mean=mu.cond.imp,sd=s.cond.imp)
    z.pred.sim.imputed = rep(z.pred.imputed, each = sims)

    ## Model predictions fhat(x.pred,z.pred) ; out-of-sample
    yhat.new.truth 		 = cbind(x.pred.sim,z.pred.sim,w.pred.sim)%*%bhat.truth
    yhat.new.l 	    	 = cbind(x.pred.sim,z.pred.sim,w.pred.sim)%*%bhat.l
    yhat.new.s    		 = cbind(x.pred.sim,w.pred.sim)%*%bhat.s

    #These take into account missing data
    yhat.new.l.pmks   	 = cbind(x.pred.sim.notmiss,z.pred.sim.notmiss,w.pred.sim.notmiss)%*%bhat.l.cc
    yhat.new.s.pmks  		 = cbind(x.pred.sim.miss,w.pred.sim.miss)%*%bhat.s.cc
    yhat.new.l.imputed   = cbind(x.pred.sim,z.pred.sim.imputed,w.pred.sim)%*%bhat.l
    yhat.new.l.cc  	 = cbind(x.pred.sim.notmiss,z.pred.sim.notmiss,w.pred.sim.notmiss)%*%bhat.l.cc
    yhat.new.s.cc 		 = cbind(x.pred.sim.miss,w.pred.sim.miss)%*%bhat.s


    ## Prediction Errors ( (y.new-fhat)^2 ) ; out-of-sample
    err.truth 	  = matrix(((y.new-yhat.new.truth)^2),ncol=n.pred,nrow=sims,byrow=FALSE)
    err.l 		  = matrix(((y.new-yhat.new.l)^2),ncol=n.pred,nrow=sims,byrow=FALSE)
    err.s  		  = matrix(((y.new-yhat.new.s)^2),ncol=n.pred,nrow=sims,byrow=FALSE)
    err.l.imputed = matrix(((y.new-yhat.new.l.imputed)^2),ncol=n.pred,nrow=sims,byrow=FALSE)


    err.l.pmks 	 = matrix(((y.new.notmiss-yhat.new.l.pmks)^2),
                          ncol=n.pred-sum(m.pred),nrow=sims,byrow=FALSE)
    err.s.pmks = matrix(((y.new.miss-yhat.new.s.pmks)^2),
                        ncol=sum(m.pred),nrow=sims,byrow=FALSE)
    err.l.cc 	 = matrix(((y.new.notmiss-yhat.new.l.cc)^2),
                        ncol=n.pred-sum(m.pred),nrow=sims,byrow=FALSE)
    err.s.cc	 = matrix(((y.new.miss - yhat.new.s.cc)^2),
                       ncol=sum(m.pred),nrow=sims,byrow=FALSE)


    # Mean prediction error (over sims) out-of-sample & overall
    avg.pe.truth		= colMeans(err.truth)		# square matrix	only
    avg.pe.l	    = colMeans(err.l)		# square matrix	only
    avg.pe.s    		= colMeans(err.s)	# square matrix	only
    avg.pe.l.pmks	 = colMeans(err.l.pmks)		# square matrix	only
    avg.pe.s.pmks  	= colMeans(err.s.pmks) 	# square matrix	only
    avg.pe.l.imputed = colMeans(err.l.imputed)	# square matrix	only
    avg.pe.l.cc	 	= colMeans(err.l.cc)		# square matrix	only
    avg.pe.s.cc  	= colMeans(err.s.cc) 	# square matrix	only


    merged.ls.sm  = cbind(c(z.pred[m.pred==0],z.pred[m.pred==1]),c(avg.pe.l.pmks,avg.pe.s.pmks))
    merged.ls = cbind(c(z.pred[m.pred==0],z.pred[m.pred==1]),c(avg.pe.l.cc,avg.pe.s.cc))

    avg.pe.pmks  = merged.ls.sm[sort(merged.ls.sm[,1],index.return=TRUE)$ix,2]
    avg.pe.cc  	 = merged.ls[sort(merged.ls[,1],index.return=TRUE)$ix,2]

    ###########
    # Exact Equations
    ##########
    X.p		    = cbind(x.pred,z.pred, w.pred)
    X.p.imputed = cbind(x.pred, z.pred.imputed, w.pred)
    X.p.small		= cbind(x.pred, w.pred)

    ## Hat Matricies for prediction
    hat.truth	= diag(X.p%*%chol2inv(chol(t(X.truth)%*%X.truth))%*%t(X.p))
    hat.l   	= diag(X.p%*%chol2inv(chol(t(X.l)%*%X.l))%*%t(X.p))
    hat.s   	= diag(X.p.small%*%chol2inv(chol(t(X.s)%*%X.s))%*%t(X.p.small))
    hat.l.cc	= diag(X.p[m.pred==0,]%*%chol2inv(chol(t(X.l.cc)%*%X.l.cc))%*%t(X.p[m.pred==0,]))
    hat.s.cc 	= diag(X.p.small[m.pred==1,]%*%chol2inv(chol(t(X.s.cc)%*%X.s.cc))%*%t(X.p.small[m.pred==1,]))
    hat.l.imputed = diag(X.p.imputed%*%chol2inv(chol(t(X.l)%*%X.l))%*%t(X.p.imputed))


    var.truth.tr	= s.e^2
    var.l.tr    	= s.e^2
    var.l.x	    	= s.e^2

    bias.truth.tr	= cbind(x.pred,z.pred,w.pred)%*%(bhat.truth-coef.tru)
    bias.l.tr   	= cbind(x.pred,z.pred,w.pred)%*%(bhat.l-coef.tru)

    bias.truth.x	= 0 	# Function assumes large model is true
    bias.l.x    	= 0 	# ONLY TRUE WITH MCAR/MAR

    epe.truth.tr	= var.truth.tr * (1) + bias.truth.tr^2
    epe.truth.x		= var.l.x * (1+hat.truth) + bias.l.x^2

    epe.l.tr	= var.l.tr * (1) + bias.l.tr^2
    epe.l.x		= var.l.x * (1+hat.l) + bias.l.x^2

    # Verify var and bias of small model for generalized matrix z
    var.s.tr	= s.e^2
    bias.s.tr	= cbind(x.pred,z.pred,w.pred)%*%(rbind(bhat.s,0)-coef.tru)

    var.s.x		= s.e^2*diag(cbind(x.pred,w.pred)%*%chol2inv(chol(t(X.s)%*%X.s))%*%t(cbind(x.pred,w.pred)))
    bias.s.x	= (cbind(x.pred,w.pred)%*%chol2inv(chol(t(X.s)%*%X.s))%*%t(X.s)%*%y.tr)- cbind(x.pred,z.pred,w.pred)%*%coef.tru

    epe.s.tr	= as.numeric(var.s.tr) * (1) + bias.s.tr^2
    epe.s.x		=  bias.s.x^2 + var.s.x + s.e^2



    ############
    # End Exact Equations
    ############

    data.frame(avg.pe.truth,
               avg.pe.l,
               avg.pe.s,
               avg.pe.pmks,
               avg.pe.cc,
               avg.pe.l.imputed,
               weighted = avg.pe.l*(1-prop.pred.miss) +
                 avg.pe.s*(prop.pred.miss),
               z.pred,
               epe.truth.tr,
               epe.truth.x,
               epe.l.tr,
               epe.l.x,
               epe.s.tr,
               epe.s.x,
               weighted.epe = epe.truth.x*(1-prop.pred.miss) +
                 epe.s.x*(prop.pred.miss)

    )
  })

  mean.epes = Reduce("+", epes) / length(epes)
  save(mean.epes, '~/SimulatedMCAR.Rda')
}




