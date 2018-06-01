sourcedir = "Sims"

stats <- function(sim.result){
  means <- apply(sim.result,2,mean, na.rm=TRUE)
  ci <- apply(sim.result,2,function(x) quantile(x, c(0.025, 0.975), na.rm=TRUE))
  out <- c()
  for(i in 1:ncol(sim.result)){
    out <- c(out,paste(format(round(means[i],2),nsmall = 2), " (",format(round(ci[1,i],2), nsmall = 2),",",format(round(ci[2,i],2), nsmall = 2),")", sep=""))
  }
  out <- t(as.data.frame(out))
  colnames(out) <- colnames(sim.result)
  out
}

appendList <- function (x, val) 
{
  stopifnot(is.list(x), is.list(val))
  xnames <- names(x)
  for (v in names(val)) {
    x[[v]] <- if (v %in% xnames && is.list(x[[v]]) && is.list(val[[v]])) 
      appendList(x[[v]], val[[v]])
    else rbind(x[[v]], val[[v]])
  }
  x
}

for(missmech in c('MAR','MNAR')){
        for(seed in 1:150){
          rdaname=paste(missmech,seed,'.Rda',sep="")
          LoadMe <-load(paste(sourcedir,'/Miss',missmech,'_seed',seed,"/", rdaname, sep=""))
          assign(paste(missmech,seed,'nopm',sep=""), res)
        }
}

for(missmech in c('MAR','MNAR')){
  for(seed in 1:150){
    rdaname=paste(missmech,seed,'_PMY.Rda',sep="")
    LoadMe <-load(paste(sourcedir,'/Miss',missmech,'_seed',seed,"/", rdaname, sep=""))
    assign(paste(missmech,seed,'pmy',sep=""), res.PMY)
  }
}

MAROutput <- MAR1nopm
MNAROutput <- MNAR1nopm
MARPMYOutput <-  MAR1pmy
MNARPMYOutput <-  MNAR1pmy

for(seed in 2:150){
  assign('tmpMAR',get(paste('MAR',seed,'nopm',sep="")))
  assign('tmpMNAR',get(paste('MNAR',seed,'nopm',sep="")))
  assign('tmpMARPMY',get(paste('MAR',seed,'pmy',sep="")))
  assign('tmpMNARPMY',get(paste('MNAR',seed,'pmy',sep="")))
  MAROutput <- mapply(rbind, MAROutput,tmpMAR)
  MNAROutput <- mapply(rbind, MNAROutput,tmpMNAR)
  MARPMYOutput <- mapply(rbind, MARPMYOutput,tmpMARPMY)
  MNARPMYOutput <- mapply(rbind, MNARPMYOutput,tmpMNARPMY)
}


summarize.results.unweighted <- function(results){
  means <- lapply(results,colMeans, na.rm=TRUE)

  is.patt1 <- c()
  is.patt2 <- c()
 
  for(i in 1:(length(means[[2]])/2)){
    is.patt1   <- c(is.patt1,means[[2]][c(2*i - 1)])  
    is.patt2   <- c(is.patt2,means[[3]][c(2*i - 1)])  
   }
  
  is.avg.mse <- (is.patt1 + is.patt2)/2
 
  oos.patt1 <- c()
  oos.patt2 <- c()

  for(i in 1:(length(means[[5]])/2)){
    oos.patt1   <- c(oos.patt1,means[[5]][c(2*i - 1)])   
    oos.patt2   <- c(oos.patt2,means[[6]][c(2*i - 1)])  
    
    }
  
  oos.avg.mse <- (oos.patt1 + oos.patt2)/2
 
    ###Unweighted MSE
  oos.patt1 <- c(oos.patt1['oracle.marg.mse'],
                 oos.patt1['cc.mean.mse'],
                 oos.patt1['marg.cm.mipack.mse'],
                 oos.patt1['metric.ccsm.mse'],
                 oos.patt1['marg.mi.Hmisc.mse'],
                 oos.patt1['marg.miy.Hmisc.mse'],
                 oos.patt1['marg.cond.mean'],
                 oos.patt1['full.cond.mean'],
                 oos.patt1['full.cm.mipack.mse'],
                 oos.patt1['full.mi.Hmisc.mse'],
                 oos.patt1['full.miy.Hmisc.mse'],
                 oos.patt1['metric.PMKS.mse']
  )
  oos.patt2 <- c(oos.patt2['oracle.marg.mse'],
                 oos.patt2['cc.mean.mse'],
                 oos.patt2['marg.cm.mipack.mse'],
                 oos.patt2['metric.ccsm.mse'],
                 oos.patt2['marg.mi.Hmisc.mse'],
                 oos.patt2['marg.miy.Hmisc.mse'],
                 oos.patt2['marg.cond.mean'],
                 oos.patt2['full.cond.mean'],
                 oos.patt2['full.cm.mipack.mse'],
                 oos.patt2['full.mi.Hmisc.mse'],
                 oos.patt2['full.miy.Hmisc.mse'],
                 oos.patt2['metric.PMKS.mse']
  )
  
  oos.avg.mse <- c( oos.avg.mse['oracle.marg.mse'],
                    oos.avg.mse['cc.mean.mse'],
                    oos.avg.mse['marg.cm.mipack.mse'],
                    oos.avg.mse['metric.ccsm.mse'],
                    oos.avg.mse['marg.mi.Hmisc.mse'],
                    oos.avg.mse['marg.miy.Hmisc.mse'],
                    oos.avg.mse['marg.cond.mean'],
                    oos.avg.mse['full.cond.mean'],
                    oos.avg.mse['full.cm.mipack.mse'],
                    oos.avg.mse['full.mi.Hmisc.mse'],
                    oos.avg.mse['full.miy.Hmisc.mse'],
                    oos.avg.mse['metric.PMKS.mse']
  )
    
    
rbind( 
                   oos.patt1,
                  oos.patt2,
                  oos.avg.mse
                  )  
  
}
  
  

summarize.results.weighted <- function(results, weight){
  means <- lapply(results,colMeans, na.rm=TRUE)
  
  is.patt1 <- c()
  is.patt2 <- c()
  
  for(i in 1:(length(means[[2]])/2)){
    is.patt1   <- c(is.patt1,means[[2]][c(2*i - 1)]*weight) 
    is.patt2   <- c(is.patt2,means[[3]][c(2*i - 1)]*weight ) 
  }
  
  is.avg.mse <- (is.patt1 + is.patt2)
  
  oos.patt1 <- c()
  oos.patt2 <- c()
  
  for(i in 1:(length(means[[5]])/2)){
    oos.patt1   <- c(oos.patt1,means[[5]][c(2*i - 1)]*weight) 
    oos.patt2   <- c(oos.patt2,means[[6]][c(2*i - 1)]*weight ) 
    
  }
  
  oos.avg.mse <- (oos.patt1 + oos.patt2)
    ###Weighted MSE
  oos.patt1 <-  c(oos.patt1['oracle.marg.mse'],
                  oos.patt1['cc.mean.mse'],
                  oos.patt1['marg.cm.mipack.mse'],
                  oos.patt1['metric.ccsm.mse'],
                  oos.patt1['marg.mi.Hmisc.mse'],
                  oos.patt1['marg.miy.Hmisc.mse'],
                  oos.patt1['marg.cond.mean'],
                  oos.patt1['full.cond.mean'],
                  oos.patt1['full.cm.mipack.mse'],
                  oos.patt1['full.mi.Hmisc.mse'],
                  oos.patt1['full.miy.Hmisc.mse'],
                  oos.patt1['metric.PMKS.mse']
  )
  oos.patt2 <- c(oos.patt2['oracle.marg.mse'],
                 oos.patt2['cc.mean.mse'],
                 oos.patt2['marg.cm.mipack.mse'],
                 oos.patt2['metric.ccsm.mse'],
                 oos.patt2['marg.mi.Hmisc.mse'],
                 oos.patt2['marg.miy.Hmisc.mse'],
                 oos.patt2['marg.cond.mean'],
                 oos.patt2['full.cond.mean'],
                 oos.patt2['full.cm.mipack.mse'],
                 oos.patt2['full.mi.Hmisc.mse'],
                 oos.patt2['full.miy.Hmisc.mse'],
                 oos.patt2['metric.PMKS.mse']
  )
  oos.avg.mse <- c( oos.avg.mse['oracle.marg.mse'],
                    oos.avg.mse['cc.mean.mse'],
                    oos.avg.mse['marg.cm.mipack.mse'],
                    oos.avg.mse['metric.ccsm.mse'],
                    oos.avg.mse['marg.mi.Hmisc.mse'],
                    oos.avg.mse['marg.miy.Hmisc.mse'],
                    oos.avg.mse['marg.cond.mean'],
                    oos.avg.mse['full.cond.mean'],
                    oos.avg.mse['full.cm.mipack.mse'],
                    oos.avg.mse['full.mi.Hmisc.mse'],
                    oos.avg.mse['full.miy.Hmisc.mse'],
                    oos.avg.mse['metric.PMKS.mse']
  )
  
 rbind(              #oos.patt1.diff,
#                     oos.patt2.diff,
#                     oos.avg.mse.diff,
#                     oos.patt1.ratio,
#                     oos.patt2.ratio,
#                     oos.avg.mse.ratio,
                    oos.patt1,
                    oos.patt2,
                    oos.avg.mse
  )  
  
}



MAR.unweighted <- summarize.results.unweighted(MAROutput)
MNAR.unweighted <- summarize.results.unweighted(MNAROutput)
MARPMY.unweighted <- summarize.results.unweighted(MARPMYOutput)
MNARPMY.unweighted <- summarize.results.unweighted(MNARPMYOutput)

MAR.weighted <- summarize.results.weighted(MAROutput, weight = 0.5)
MNAR.weighted <- summarize.results.weighted(MNAROutput, weight = 0.5)
MARPMY.weighted <- summarize.results.weighted(MARPMYOutput, weight = 0.5)
MNARPMY.weighted <- summarize.results.weighted(MNARPMYOutput, weight = 0.5)


##############
#Weighted Raw
#############
#Save as and 8x10 pdf
  XLAB = 'Pattern Specific Contribution to Prediction Error'
  TITLE = 'Comparison of Pattern Prediction Error Among Missingness Mechanisms'
  
  plot(MAR.weighted['oos.patt1' ,c(  "cc.mean.mse",
                                     "marg.cond.mean",
                                     "marg.cm.mipack.mse",
                                     "marg.mi.Hmisc.mse",
                                     "marg.miy.Hmisc.mse",
                                     "metric.ccsm.mse",
                                     "full.cond.mean",
                                     "full.cm.mipack.mse",
                                     "full.mi.Hmisc.mse",
                                     "full.miy.Hmisc.mse",
                                     "metric.PMKS.mse")],
       2:12, 
       pch = 19, 
       xlab = '',
       ylab = "", 
       yaxt = 'n',
       xlim = c(-3, 12), main = "", col = 'blue', ylim = c(0, 48), cex = 1,
       xaxt = 'n', bty = 'n') 
  axis(1,at = c(0,1,2,3,4,5,6,7,8,9,10),line =-1.95)
  segments(x0 = rep(0,6),x1=rep(12,6), y0= c(2:13)-0.5,y1= c(2:13)-0.5, lty = 2, col='gray')
  segments(x0 = rep(0,6),x1=rep(12,6), y0= c(14:25)-0.5,y1= c(14:25)-0.5, lty = 2, col='gray')
  segments(x0 = rep(0,6),x1=rep(12,6), 
           y0= c(26:36)-0.5,y1= c(26:36)-0.5, lty = 2, col='gray')
  
  segments(x0 = rep(0,6),x1=rep(12,6), 
           y0= c(37:49)-0.5,y1= c(37:49)-0.5, lty = 2, col='gray')
  
  points(MAR.weighted['oos.patt1' ,c(  "cc.mean.mse",
                                       "marg.cond.mean",
                                       "marg.cm.mipack.mse",
                                       "marg.mi.Hmisc.mse",
                                       "marg.miy.Hmisc.mse",
                                       "metric.ccsm.mse",
                                       "full.cond.mean",
                                       "full.cm.mipack.mse",
                                       "full.mi.Hmisc.mse",
                                       "full.miy.Hmisc.mse",
                                       "metric.PMKS.mse")],
         2:12, 
         pch = 19, 
         col = 'blue',cex = 1)
  
  points(MAR.weighted['oos.patt2',c("cc.mean.mse",
                                    "marg.cond.mean",
                                    "marg.cm.mipack.mse",
                                    "marg.mi.Hmisc.mse",
                                    "marg.miy.Hmisc.mse",
                                    "metric.ccsm.mse",
                                    "full.cond.mean",
                                    "full.cm.mipack.mse",
                                    "full.mi.Hmisc.mse",
                                    "full.miy.Hmisc.mse",
                                    "metric.PMKS.mse")],
         2:12, pch = 19, col = 'darkgreen', cex = 1) 
  points(MAR.weighted['oos.avg.mse',c("cc.mean.mse",
                                      "marg.cond.mean",
                                      "marg.cm.mipack.mse",
                                      "marg.mi.Hmisc.mse",
                                      "marg.miy.Hmisc.mse",
                                      "metric.ccsm.mse",
                                      "full.cond.mean",
                                      "full.cm.mipack.mse",
                                      "full.mi.Hmisc.mse",
                                      "full.miy.Hmisc.mse",
                                      "metric.PMKS.mse")],
         2:12, pch = 17, col = 'firebrick3', cex = 1)
  text(x= 0.1,y = seq(1, 12, by=1) , 
       labels = c("",
                  "Complete Case Model",
                  "Cond. Mean Imp. (Freq.)",
                  "Cond. Mean Imp. (Bayes.)",
                  "Multiple Imp. (no y)",
                  "Multiple Imp. (include y)",
                  "CCS",
                  "Cond. Mean Imp. (Freq: MIMI)",
                  "Cond. Mean Imp. (Bayes: MIMI)",
                  "Multiple Imp. (no y: MIMI)",
                  "Multiple Imp. (include y: MIMI)",
                  "PS"
       ), font = 2,  pos = 2, xpd = TRUE, cex = 0.65)
  
  #MAR PMY
  points(MARPMY.weighted['oos.patt1',c("cc.mean.mse",
                                       "marg.cond.mean",
                                       "marg.cm.mipack.mse",
                                       "marg.mi.Hmisc.mse",
                                       "marg.miy.Hmisc.mse",
                                       "metric.ccsm.mse",
                                       "full.cond.mean",
                                       "full.cm.mipack.mse",
                                       "full.mi.Hmisc.mse",
                                       "full.miy.Hmisc.mse",
                                       "metric.PMKS.mse")],
         14:24, pch = 19, col = 'blue', cex = 1) 
  
  points(MARPMY.weighted['oos.patt2',c("cc.mean.mse",
                                       "marg.cond.mean",
                                       "marg.cm.mipack.mse",
                                       "marg.mi.Hmisc.mse",
                                       "marg.miy.Hmisc.mse",
                                       "metric.ccsm.mse",
                                       "full.cond.mean",
                                       "full.cm.mipack.mse",
                                       "full.mi.Hmisc.mse",
                                       "full.miy.Hmisc.mse",
                                       "metric.PMKS.mse")],
         14:24, pch = 19, col = 'darkgreen', cex = 1) 
  points(MARPMY.weighted['oos.avg.mse',c("cc.mean.mse",
                                         "marg.cond.mean",
                                         "marg.cm.mipack.mse",
                                         "marg.mi.Hmisc.mse",
                                         "marg.miy.Hmisc.mse",
                                         "metric.ccsm.mse",
                                         "full.cond.mean",
                                         "full.cm.mipack.mse",
                                         "full.mi.Hmisc.mse",
                                         "full.miy.Hmisc.mse",
                                         "metric.PMKS.mse")], 14:24, pch = 17, 
         col = 'firebrick3', cex = 1)
  text(x= 0.1,y = seq(13, 24, by=1) , 
       labels = c("",
                  "Complete Case Model",
                  "Cond. Mean Imp. (Freq.)",
                  "Cond. Mean Imp. (Bayes.)",
                  "Multiple Imp. (no y)",
                  "Multiple Imp. (include y)",
                  "CCS",
                  "Cond. Mean Imp. (Freq: MIMI)",
                  "Cond. Mean Imp. (Bayes: MIMI)",
                  "Multiple Imp. (no y: MIMI)",
                  "Multiple Imp. (include y: MIMI)",
                  "PS"
       ), font = 2,  pos = 2, xpd = TRUE, cex = 0.65)
  
  
  #MNAR
  points(MNAR.weighted['oos.patt1',c("cc.mean.mse",
                       "marg.cond.mean",
                       "marg.cm.mipack.mse",
                       "marg.mi.Hmisc.mse",
                       "marg.miy.Hmisc.mse",
                       "metric.ccsm.mse",
                       "full.cond.mean",
                       "full.cm.mipack.mse",
                       "full.mi.Hmisc.mse",
                       "full.miy.Hmisc.mse",
                       "metric.PMKS.mse")],
         26:36, pch = 19, col = 'blue', cex = 1) 
  
  points(MNAR.weighted['oos.patt2',c("cc.mean.mse",
                                     "marg.cond.mean",
                                     "marg.cm.mipack.mse",
                                     "marg.mi.Hmisc.mse",
                                     "marg.miy.Hmisc.mse",
                                     "metric.ccsm.mse",
                                     "full.cond.mean",
                                     "full.cm.mipack.mse",
                                     "full.mi.Hmisc.mse",
                                     "full.miy.Hmisc.mse",
                                     "metric.PMKS.mse")],
         26:36, pch = 19, col = 'darkgreen', cex = 1) 
  points(MNAR.weighted['oos.avg.mse',c("cc.mean.mse",
                                       "marg.cond.mean",
                                       "marg.cm.mipack.mse",
                                       "marg.mi.Hmisc.mse",
                                       "marg.miy.Hmisc.mse",
                                       "metric.ccsm.mse",
                                       "full.cond.mean",
                                       "full.cm.mipack.mse",
                                       "full.mi.Hmisc.mse",
                                       "full.miy.Hmisc.mse",
                                       "metric.PMKS.mse")],
         26:36, pch = 17, col = 'firebrick3', cex = 1)
    text(x= 0.1,y = seq(25, 36, by=1) , 
       labels = c("",
                  "Complete Case Model",
                  "Cond. Mean Imp. (Freq.)",
                  "Cond. Mean Imp. (Bayes.)",
                  "Multiple Imp. (no y)",
                  "Multiple Imp. (include y)",
                  "CCS",
                  "Cond. Mean Imp. (Freq: MIMI)",
                  "Cond. Mean Imp. (Bayes: MIMI)",
                  "Multiple Imp. (no y: MIMI)",
                  "Multiple Imp. (include y: MIMI)",
                  "PS"
       ), font = 2,  pos = 2, xpd = TRUE, cex = 0.65)
  
  #MNAR PMY
  points(MNARPMY.weighted['oos.patt1',c("cc.mean.mse",
                                        "marg.cond.mean",
                                        "marg.cm.mipack.mse",
                                        "marg.mi.Hmisc.mse",
                                        "marg.miy.Hmisc.mse",
                                        "metric.ccsm.mse",
                                        "full.cond.mean",
                                        "full.cm.mipack.mse",
                                        "full.mi.Hmisc.mse",
                                        "full.miy.Hmisc.mse",
                                        "metric.PMKS.mse")],
         38:48, pch = 19, col = 'blue', cex = 1) 
  
  points(MNARPMY.weighted['oos.patt2',c("cc.mean.mse",
                                        "marg.cond.mean",
                                        "marg.cm.mipack.mse",
                                        "marg.mi.Hmisc.mse",
                                        "marg.miy.Hmisc.mse",
                                        "metric.ccsm.mse",
                                        "full.cond.mean",
                                        "full.cm.mipack.mse",
                                        "full.mi.Hmisc.mse",
                                        "full.miy.Hmisc.mse",
                                        "metric.PMKS.mse")],
         38:48, pch = 19, col = 'darkgreen', cex = 1) 
  points(MNARPMY.weighted['oos.avg.mse',c("cc.mean.mse",
                                          "marg.cond.mean",
                                          "marg.cm.mipack.mse",
                                          "marg.mi.Hmisc.mse",
                                          "marg.miy.Hmisc.mse",
                                          "metric.ccsm.mse",
                                          "full.cond.mean",
                                          "full.cm.mipack.mse",
                                          "full.mi.Hmisc.mse",
                                          "full.miy.Hmisc.mse",
                                          "metric.PMKS.mse")],
         38:48, pch = 17, col = 'firebrick3', cex = 1)
   text(x= 0.1,y = seq(37, 48, by=1) , 
       labels = c("",
                  "Complete Case Model",
                  "Cond. Mean Imp. (Freq.)",
                  "Cond. Mean Imp. (Bayes.)",
                  "Multiple Imp. (no y)",
                  "Multiple Imp. (include y)",
                  "CCS",
                  "Cond. Mean Imp. (Freq: MIMI)",
                  "Cond. Mean Imp. (Bayes: MIMI)",
                  "Multiple Imp. (no y: MIMI)",
                  "Multiple Imp. (include y: MIMI)",
                  "PS"
       ), font = 2,  pos = 2, xpd = TRUE, cex = 0.65)
  
  abline(h = 13, lwd = 2)
  abline(h = 25, lwd = 2)
  abline(h = 37, lwd = 2)
  abline(h = 49, lwd = 2)
  abline(h = 0.51, lwd = 2)
  abline(v = -5, lwd = 2)
  #abline(h = 0.1, lwd = 2)
  
  
 segments(x0 = 0, x1 = 0, y0 = 0, y1 = 49, lwd= 2)
  

  mtext(bquote(.("Red: Total Prediction Error, Blue: Pattern 1 - No Missing Data, Green: Pattern 2 - Missing ") ~ X[1]),
        side =3, line = 0.1, cex = 0.85, col = 'black')
 
  mtext(TITLE,
        side =3, line = 2.3, cex = 1.25, font = 2)
  
  text(x=-4.0,y = c(8,21.5,33,46.75), labels = c('MAR','MAR PMY', 'MNAR', 'MNAR PMY'),
       srt = 90, pos =2, xpd = TRUE, font = 2)
  
  text(x = 5, y = -3.9, labels = XLAB, xpd = TRUE)
  
  
#################
# EPE Final Plot
################
load('SimulatedMCAR.Rda')
lo.z=-1; hi.z=7; 
z.pred 	  = seq(lo.z, hi.z, acc)
s.e = 1
sorted.index = sort(z.pred,index.return=TRUE)$ix
with(mean.epes,	{
  par(las=1)
  plot(z.pred, avg.pe.s,pch=20,type="n",
       ylab="Avg. Prediction Error",
       xlab=expression(paste("Out of Sample ", X[1])),
       ylim=c(0,100),xlim=c(-1,7),
       main='Simulated MCAR',xaxt = 'n')
  axis(1, at = c(-1,0,1,2,3,4,5,6,7),
       labels = c(expression(paste("-4",sigma,sep = "")),
                  expression(paste("-3",sigma,sep = "")),
                  expression(paste("-2",sigma,sep = "")),
                  expression(paste("-",sigma,sep = "")),
                  expression(paste("E[", X[1],"]", sep = "")),
                  expression(paste("",sigma,sep = "")),
                  expression(paste("2",sigma,sep = "")),
                  expression(paste("3",sigma,sep = "")),
                  expression(paste("-4",sigma,sep = ""))))
  abline(h=(s.e^2),lty=2,lwd=0.5)
  points(z.pred,avg.pe.l,pch=20,col="chartreuse4",cex=0.8)
  points(z.pred,avg.pe.s,pch=20,col="blueviolet",cex=0.8)
  points(z.pred,avg.pe.PS, pch=20, col="goldenrod1",cex=0.8)
  points(z.pred, weighted, pch=20, col='red')
  legend(1.75,95,c("Large" ,"Small","PS","Weighted Avg."),lty=1,lwd=3,
         col=c("chartreuse4","blueviolet","goldenrod1","red"),bty="n")
})


##################################
#SUPPORT plot with original data
#Unweighted MSE
#################################
#Save as pdf 8 x 6
par(mfrow = c(2,2))
load("supportSM.Rda")
par(mar= c(1.95,1.5,3,0.5), xpd=TRUE)

#layout(matrix(c(1,1,2,2,1,1,2,2), 2, 4, byrow = TRUE), width=c(1,1), respect=FALSE)

#Start with the grid do show which patterns are present
plot(1:10, 1:10, ylim = c(0,23), xlim = c(0,800), pch = NA, bty= 'n', yaxt = 'n', xaxt = 'n', xlab = "", ylab = "", main = 'SUPPORT Example')
sorted.n <- sort(round(all.results$prop.pattern*9103,2),index.return=TRUE)$ix

text(x = seq(0,180,20)+15 , par("usr")[3] + 0.95 , 
     labels = c('pafi', 'meanbp', 'wblc', 'alb', 'resp', 'temp', 'hrt', 'bili', 'crea', 'sod')
     , srt = 90, pos = 2, xpd = TRUE, cex = 0.75)

#"0000000000" #1
points(seq(0,180,20),rep(23,10), 
       pch=c('','','','','','','','','',''), cex = 0.75)
#"0000000100" #2
points(seq(0,180,20),rep(18,10),
       pch=c('','','','','','','','X','',''), cex = 0.75)
#"0001000000" #3
points(seq(0,180,20),rep(20,10), 
       pch=c('','','','X','','','','','',''), cex = 0.75)
#"0001000010" #4
points(seq(0,180,20),rep(1,10), 
       pch=c('','','','X','','','','','X',''), cex = 0.75)
#"0001000100" #5
points(seq(0,180,20),rep(22,10), 
       pch=c('','','','X','','','','X','',''), cex = 0.75)
#"0001000110" #6
points(seq(0,180,20),rep(5,10), 
       pch=c('','','','X','','','','X','X',''), cex = 0.75)
#"0010000000" #7
points(seq(0,180,20),rep(13,10), 
       pch=c('','','X','','','','','','',''), cex = 0.75)
#"0010000100" #8
points(seq(0,180,20),rep(7,10), 
       pch=c('','','X','','','','','X','',''), cex = 0.75)
#"0011000000" #9
points(seq(0,180,20),rep(6,10), 
       pch=c('','','X','X','','','','','',''), cex = 0.75)

#"0011000100" #10
points(seq(0,180,20),rep(12,10), 
       pch=c('','','X','X','','','','X','',''), cex = 0.75)

#"0011000110" #11
points(seq(0,180,20),rep(9,10), 
       pch=c('','','X','X','','','','X','X',''), cex = 0.75)

#"1000000000" #12
points(seq(0,180,20),rep(21,10), 
       pch=c('X','','','','','','','','',''), cex = 0.75)

#"1000000100" #13
points(seq(0,180,20),rep(16,10), 
       pch=c('X','','','','','','','X','',''), cex = 0.75)

#"1001000000" #14
points(seq(0,180,20),rep(17,10), 
       pch=c('X','','','X','','','','','',''), cex = 0.75)

#"1001000010" #15
points(seq(0,180,20),rep(2,10), 
       pch=c('X','','','X','','','','','X',''), cex = 0.75)

#"1001000100" #16
points(seq(0,180,20),rep(19,10), 
       pch=c('X','','','X','','','','X','',''), cex = 0.75)

#"1001000110" #17
points(seq(0,180,20),rep(10,10), 
       pch=c('X','','','X','','','','X','X',''), cex = 0.75)

#"1010000000" #18
points(seq(0,180,20),rep(11,10), 
       pch=c('X','','X','','','','','','',''), cex = 0.75)

#"1010000100" #19
points(seq(0,180,20),rep(3,10), 
       pch=c('X','','X','','','','','X','',''), cex = 0.75)

#"1010000110" #20
points(seq(0,180,20),rep(4,10), 
       pch=c('X','','X','','','','','X','X',''), cex = 0.75)

#"1011000000" #21
points(seq(0,180,20),rep(8,10), 
       pch=c('X','','X','X','','','','','',''), cex = 0.75)

#"1011000100" #22
points(seq(0,180,20),rep(15,10), 
       pch=c('X','','X','X','','','','X','',''), cex = 0.75)

#"1011000110" #23
points(seq(0,180,20),rep(14,10), 
       pch=c('X','','X','X','','','','X','X',''), cex = 0.75)

segments(x0 = 200, x1 = 200, y0 =0, y1 = 23.5, xpd = TRUE)

points(((all.results$metric.MI.withy.brier) + 200)[sorted.n],
       (1:23), pch = 19, col = 'red')
points(((all.results$metric.MIMI.withy.brier) + 200)[sorted.n],
       (1:23), pch = 19, col = 'blue')
points(((all.results$brier.y) + 200)[sorted.n],
       (1:23), pch = 19, col = 'darkgreen')
points(((all.results$brier.x) + 200)[sorted.n],
       (1:23),pch = 19, col = 'black')
abline(h = c(2:length(all.results$Row.names))-0.5, lty = 2, lwd=0.5, col='gray')


legend(590,22, legend = c("PS", "MIMI", "CCS", "MI"),
       col = c('black','blue','darkgreen','red'),
       pch = c(19,19,19,19), bty = 'n', cex = 0.75,
       x.intersp = 0.2
)

legend(660,23.8, legend = c(round(sum(all.results$brier.x*all.results$prop.pattern),1),
                            round(sum(all.results$metric.MIMI.withy.brier*all.results$prop.pattern),1),
                            round(sum(all.results$brier.y*all.results$prop.pattern),1),
                            round(sum(all.results$metric.MI.withy.brier*all.results$prop.pattern),1)),
       bty='n',
       col = c(NA,NA,NA,NA),
       pch = c(NA,NA,NA,NA), cex = 0.75, title = "Model PE")

text(x= par("usr")[3] - 5.5,y = 1:23 - 0.1 , 
     labels = round(all.results$prop.pattern*9103,2)[sorted.n]
     ,pos = 2, xpd = TRUE, cex = 0.75)

text(par("usr")[3]-25, y = 25.5, labels = "N",pos = 1, xpd = TRUE, cex = 0.75)
text(par("usr")[3] + 105, y = 25.5, labels = "Pattern",pos = 1, xpd = TRUE, cex = 0.75)
axis(1,at = seq(200,800,20), line = -0.60, tick = seq(200,800,20), labels = seq(0,600,20))
segments(x0 =-12, y0 = 0, x1 = -12, y1 = 24.5 , xpd = TRUE)
segments(x0 = -12, y0 = 0.32, x1 = 200, y1 = 0.32)
segments(x0 = seq(0,180,20) + 10, x1 = seq(0,180,20) + 10, y0 = rep(0.5,10), y1 = rep(23.5,10), lty = 2, lwd=0.5, col='gray')
segments(x0 = -35, y0=23.5,y1=23.5,x1=200, xpd = TRUE)

mtext('Prediction Error', side = 1, line = 1.55, at = 450, cex =  0.75)


####################################
#SUPPORT plot with manipulated data
#Unweighted
###################################

load("supportAdd25pafi.Rda")
par(mar= c(1.95,1.3,3,0.5), xpd=TRUE)

#par(mar= c(2,0.5,4,0.5), xpd=TRUE)
#Start with the grid do show which patterns are present
plot(1:10, 1:10, ylim = c(0,23), xlim = c(0,800), pch = NA,
     bty= 'n', yaxt = 'n', xaxt = 'n', xlab = "", ylab = "",
     main = 'SUPPORT Example (Induced MNAR)')

sorted.n <- sort(round(all.results$prop.pattern*9103,2),index.return=TRUE)$ix

text(x = seq(0,180,20)+15 , par("usr")[3] + 0.95 , 
     labels = c('pafi', 'meanbp', 'wblc', 'alb', 'resp', 'temp', 'hrt', 'bili', 'crea', 'sod')
     , srt = 90, pos = 2, xpd = TRUE, cex = 0.75)

#"0000000000" #1
points(seq(0,180,20),rep(23,10), 
       pch=c('','','','','','','','','',''), cex = 0.75)
#"0000000100" #2
points(seq(0,180,20),rep(18,10),
       pch=c('','','','','','','','X','',''), cex = 0.75)
#"0001000000" #3
points(seq(0,180,20),rep(20,10), 
       pch=c('','','','X','','','','','',''), cex = 0.75)
#"0001000010" #4
points(seq(0,180,20),rep(1,10), 
       pch=c('','','','X','','','','','X',''), cex = 0.75)
#"0001000100" #5
points(seq(0,180,20),rep(22,10), 
       pch=c('','','','X','','','','X','',''), cex = 0.75)
#"0001000110" #6
points(seq(0,180,20),rep(5,10), 
       pch=c('','','','X','','','','X','X',''), cex = 0.75)
#"0010000000" #7
points(seq(0,180,20),rep(13,10), 
       pch=c('','','X','','','','','','',''), cex = 0.75)
#"0010000100" #8
points(seq(0,180,20),rep(7,10), 
       pch=c('','','X','','','','','X','',''), cex = 0.75)
#"0011000000" #9
points(seq(0,180,20),rep(6,10), 
       pch=c('','','X','X','','','','','',''), cex = 0.75)

#"0011000100" #10
points(seq(0,180,20),rep(12,10), 
       pch=c('','','X','X','','','','X','',''), cex = 0.75)

#"0011000110" #11
points(seq(0,180,20),rep(9,10), 
       pch=c('','','X','X','','','','X','X',''), cex = 0.75)

#"1000000000" #12
points(seq(0,180,20),rep(21,10), 
       pch=c('X','','','','','','','','',''), cex = 0.75)

#"1000000100" #13
points(seq(0,180,20),rep(16,10), 
       pch=c('X','','','','','','','X','',''), cex = 0.75)

#"1001000000" #14
points(seq(0,180,20),rep(17,10), 
       pch=c('X','','','X','','','','','',''), cex = 0.75)

#"1001000010" #15
points(seq(0,180,20),rep(2,10), 
       pch=c('X','','','X','','','','','X',''), cex = 0.75)

#"1001000100" #16
points(seq(0,180,20),rep(19,10), 
       pch=c('X','','','X','','','','X','',''), cex = 0.75)

#"1001000110" #17
points(seq(0,180,20),rep(10,10), 
       pch=c('X','','','X','','','','X','X',''), cex = 0.75)

#"1010000000" #18
points(seq(0,180,20),rep(11,10), 
       pch=c('X','','X','','','','','','',''), cex = 0.75)

#"1010000100" #19
points(seq(0,180,20),rep(3,10), 
       pch=c('X','','X','','','','','X','',''), cex = 0.75)

#"1010000110" #20
points(seq(0,180,20),rep(4,10), 
       pch=c('X','','X','','','','','X','X',''), cex = 0.75)

#"1011000000" #21
points(seq(0,180,20),rep(8,10), 
       pch=c('X','','X','X','','','','','',''), cex = 0.75)

#"1011000100" #22
points(seq(0,180,20),rep(15,10), 
       pch=c('X','','X','X','','','','X','',''), cex = 0.75)

#"1011000110" #23
points(seq(0,180,20),rep(14,10), 
       pch=c('X','','X','X','','','','X','X',''), cex = 0.75)

segments(x0 = 200, x1 = 200, y0 =0, y1 = 23.5, xpd = TRUE)

points(((all.results$metric.MI.withy.brier) + 200)[sorted.n],
       (1:23), pch = 19, col = 'red')
points(((all.results$metric.MIMI.withy.brier) + 200)[sorted.n],
       (1:23), pch = 19, col = 'blue')
points(((all.results$brier.y) + 200)[sorted.n],
       (1:23), pch = 19, col = 'darkgreen')
points(((all.results$brier.x) + 200)[sorted.n],
       (1:23),pch = 19, col = 'black')
points(800,4,pch = 8)
arrows(x0=800,x1=820, y0=4,y1=4, length = 0.05)
#which(all.results$brier.x == all.results$brier.y)

abline(h = c(2:length(all.results$Row.names))-0.5, lty = 2, lwd=0.5, col='gray')

legend(590,22, legend = c("PS", "MIMI", "CCS", "MI"),
       col = c('black','blue','darkgreen','red'),
       pch = c(19,19,19,19), bty = 'n', cex = 0.75,
       x.intersp = 0.2
)

legend(660,23.8, legend = c(round(sum(all.results$brier.x*all.results$prop.pattern),1),
                            round(sum(all.results$metric.MIMI.withy.brier*all.results$prop.pattern),1),
                            round(sum(all.results$brier.y*all.results$prop.pattern),1),
                            round(sum(all.results$metric.MI.withy.brier*all.results$prop.pattern),1)),
       bty='n',
       col = c(NA,NA,NA,NA),
       pch = c(NA,NA,NA,NA), cex = 0.75, title = 'Model PE')


text(x= par("usr")[3] - 5.5,y = (1:23 - 0.1) , 
     labels = round(all.results$prop.pattern*9103,2)[sorted.n]
     ,pos = 2, xpd = TRUE, cex = 0.75)

text(par("usr")[3]-25, y = 25.5, labels = "N",pos = 1, xpd = TRUE, cex = 0.75)
text(par("usr")[3] + 105, y = 25.5, labels = "Pattern",pos = 1, xpd = TRUE, cex = 0.75)
axis(1,at = seq(200,800,20), line = -0.60, tick = seq(200,800,20), labels = seq(0,600,20))
segments(x0 =-12, y0 = 0, x1 = -12, y1 = 24.5 , xpd = TRUE)
segments(x0 = -12, y0 = 0.32, x1 = 200, y1 = 0.32)
segments(x0 = seq(0,180,20) + 10, x1 = seq(0,180,20) + 10, y0 = rep(0.5,10), y1 = rep(23.5,10), lty = 2, lwd=0.5, col='gray')
segments(x0 = -35, y0=23.5,y1=23.5,x1=200, xpd = TRUE)

mtext('Prediction Error', side = 1, line = 1.55, at = 450, cex =  0.75)

###################
#SUPPORT
#Weighted by prop
###################
#layout(matrix(c(1,1,2,2,1,1,2,2), 2, 4, byrow = TRUE), width=c(1,1), respect=FALSE)

load("supportSM.Rda")
par(mar= c(2.6,1.3,3,0.5), xpd=TRUE)

#Start with the grid do show which patterns are present
plot(1:10, 1:10, ylim = c(0,23), xlim = c(0,80), pch = NA, bty= 'n', yaxt = 'n', xaxt = 'n', xlab = "", ylab = "", main = 'SUPPORT Example')

text(x = seq(0,18,2)+1.5 , par("usr")[3] + 0.95 , 
     labels = c('pafi', 'meanbp', 'wblc', 'alb', 'resp', 'temp', 'hrt', 'bili', 'crea', 'sod')
     , srt = 90, pos = 2, xpd = TRUE, cex = 0.75)

#"0000000000" #1
points(seq(0,18,2),rep(23,10), 
       pch=c('','','','','','','','','',''), cex = 0.75)
#"0000000100" #2
points(seq(0,18,2),rep(18,10),
       pch=c('','','','','','','','X','',''), cex = 0.75)
#"0001000000" #3
points(seq(0,18,2),rep(20,10), 
       pch=c('','','','X','','','','','',''), cex = 0.75)
#"0001000010" #4
points(seq(0,18,2),rep(1,10), 
       pch=c('','','','X','','','','','X',''), cex = 0.75)
#"0001000100" #5
points(seq(0,18,2),rep(22,10), 
       pch=c('','','','X','','','','X','',''), cex = 0.75)
#"0001000110" #6
points(seq(0,18,2),rep(5,10), 
       pch=c('','','','X','','','','X','X',''), cex = 0.75)
#"0010000000" #7
points(seq(0,18,2),rep(13,10), 
       pch=c('','','X','','','','','','',''), cex = 0.75)
#"0010000100" #8
points(seq(0,18,2),rep(7,10), 
       pch=c('','','X','','','','','X','',''), cex = 0.75)
#"0011000000" #9
points(seq(0,18,2),rep(6,10), 
       pch=c('','','X','X','','','','','',''), cex = 0.75)

#"0011000100" #10
points(seq(0,18,2),rep(12,10), 
       pch=c('','','X','X','','','','X','',''), cex = 0.75)

#"0011000110" #11
points(seq(0,18,2),rep(9,10), 
       pch=c('','','X','X','','','','X','X',''), cex = 0.75)

#"1000000000" #12
points(seq(0,18,2),rep(21,10), 
       pch=c('X','','','','','','','','',''), cex = 0.75)

#"1000000100" #13
points(seq(0,18,2),rep(16,10), 
       pch=c('X','','','','','','','X','',''), cex = 0.75)

#"1001000000" #14
points(seq(0,18,2),rep(17,10), 
       pch=c('X','','','X','','','','','',''), cex = 0.75)

#"1001000010" #15
points(seq(0,18,2),rep(2,10), 
       pch=c('X','','','X','','','','','X',''), cex = 0.75)

#"1001000100" #16
points(seq(0,18,2),rep(19,10), 
       pch=c('X','','','X','','','','X','',''), cex = 0.75)

#"1001000110" #17
points(seq(0,18,2),rep(10,10), 
       pch=c('X','','','X','','','','X','X',''), cex = 0.75)

#"1010000000" #18
points(seq(0,18,2),rep(11,10), 
       pch=c('X','','X','','','','','','',''), cex = 0.75)

#"1010000100" #19
points(seq(0,18,2),rep(3,10), 
       pch=c('X','','X','','','','','X','',''), cex = 0.75)

#"1010000110" #20
points(seq(0,18,2),rep(4,10), 
       pch=c('X','','X','','','','','X','X',''), cex = 0.75)

#"1011000000" #21
points(seq(0,18,2),rep(8,10), 
       pch=c('X','','X','X','','','','','',''), cex = 0.75)

#"1011000100" #22
points(seq(0,18,2),rep(15,10), 
       pch=c('X','','X','X','','','','X','',''), cex = 0.75)

#"1011000110" #23
points(seq(0,18,2),rep(14,10), 
       pch=c('X','','X','X','','','','X','X',''), cex = 0.75)

segments(x0 = 20, x1 = 20, y0 =0, y1 = 23.5, xpd = TRUE)

points(((all.results$metric.MI.withy.brier*all.results$prop.pattern) + 20)[sorted.n],
       (1:23), pch = 19, col = 'red')
points(((all.results$metric.MIMI.withy.brier*all.results$prop.pattern) + 20)[sorted.n],
       (1:23), pch = 19, col = 'blue')
points(((all.results$brier.y*all.results$prop.pattern) + 20)[sorted.n],
       (1:23), pch = 19, col = 'darkgreen')
points(((all.results$brier.x*all.results$prop.pattern) + 20)[sorted.n],
       (1:23),pch = 19, col = 'black')
abline(h = c(2:length(all.results$Row.names))-0.5, lty = 2, lwd=0.5, col='gray')

legend(59,12, legend = c("PS", "MIMI", "CCS", "MI"),
       col = c('black','blue','darkgreen','red'),
       pch = c(19,19,19,19), bty = 'n', cex = 0.75, x.intersp = 0.2
)


legend(66,13.8, legend = c(round(sum(all.results$brier.x*all.results$prop.pattern),1),
                           round(sum(all.results$metric.MIMI.withy.brier*all.results$prop.pattern),1),
                           round(sum(all.results$brier.y*all.results$prop.pattern),1),
                           round(sum(all.results$metric.MI.withy.brier*all.results$prop.pattern),1)),
       bty='n',
       col = c(NA,NA,NA,NA),
       pch = c(NA,NA,NA,NA), cex = 0.75,
       title = "Model PE")


text(x= par("usr")[3] + .55,y = 1:23 - 0.1 , 
     labels = round(all.results$prop.pattern*9103,2)[sorted.n]
     ,pos = 2, xpd = TRUE, cex = 0.75)

text(par("usr")[3]-2.3, y = 25.7, labels = "N",pos = 1, xpd = TRUE, cex = 0.75)
text(par("usr")[3] + 10.5, y = 25.7, labels = "Pattern",pos = 1, xpd = TRUE, cex = 0.75)
axis(1,at = seq(20,80,2), line =-0.60, tick = seq(20,80,2), labels = seq(0,60,2))
segments(x0 = -1.2, y0 = 0, x1 = -1.2, y1 = 24.5 , xpd = TRUE)
segments(x0 = -1.2, y0 = 0.32, x1 = 20, y1 = 0.32)
segments(x0 = seq(0,18,2) + 1, x1 = seq(0,18,2) + 1, y0 = rep(0.5,10), y1 = rep(23.5,10), lty = 2, lwd=0.5, col='gray')
segments(x0 = -3.5, y0=23.5,y1=23.5,x1=20, xpd = TRUE)
mtext('Pattern Specific Contribution (PE*Pattern Weight)', side = 1, line = 1.50, at = 50, cex =  0.75)




#Manipulated Data
load("supportAdd25pafi.Rda")
par(mar= c(2.6,1.3,3,0.5), xpd=TRUE)

#Start with the grid do show which patterns are present
plot(1:10, 1:10, ylim = c(0,23), xlim = c(0,80),
     pch = NA, bty= 'n', yaxt = 'n', xaxt = 'n', 
     xlab = "", ylab = "", main = 'SUPPORT Example (Induced MNAR)')

text(x = seq(0,18,2)+1.5 , par("usr")[3] + 0.95 , 
     labels = c('pafi', 'meanbp', 'wblc', 'alb', 'resp', 'temp', 'hrt', 'bili', 'crea', 'sod')
     , srt = 90, pos = 2, xpd = TRUE, cex = 0.75)

#"0000000000" #1
points(seq(0,18,2),rep(23,10), 
       pch=c('','','','','','','','','',''), cex = 0.75)
#"0000000100" #2
points(seq(0,18,2),rep(18,10),
       pch=c('','','','','','','','X','',''), cex = 0.75)
#"0001000000" #3
points(seq(0,18,2),rep(20,10), 
       pch=c('','','','X','','','','','',''), cex = 0.75)
#"0001000010" #4
points(seq(0,18,2),rep(1,10), 
       pch=c('','','','X','','','','','X',''), cex = 0.75)
#"0001000100" #5
points(seq(0,18,2),rep(22,10), 
       pch=c('','','','X','','','','X','',''), cex = 0.75)
#"0001000110" #6
points(seq(0,18,2),rep(5,10), 
       pch=c('','','','X','','','','X','X',''), cex = 0.75)
#"0010000000" #7
points(seq(0,18,2),rep(13,10), 
       pch=c('','','X','','','','','','',''), cex = 0.75)
#"0010000100" #8
points(seq(0,18,2),rep(7,10), 
       pch=c('','','X','','','','','X','',''), cex = 0.75)
#"0011000000" #9
points(seq(0,18,2),rep(6,10), 
       pch=c('','','X','X','','','','','',''), cex = 0.75)

#"0011000100" #10
points(seq(0,18,2),rep(12,10), 
       pch=c('','','X','X','','','','X','',''), cex = 0.75)

#"0011000110" #11
points(seq(0,18,2),rep(9,10), 
       pch=c('','','X','X','','','','X','X',''), cex = 0.75)

#"1000000000" #12
points(seq(0,18,2),rep(21,10), 
       pch=c('X','','','','','','','','',''), cex = 0.75)

#"1000000100" #13
points(seq(0,18,2),rep(16,10), 
       pch=c('X','','','','','','','X','',''), cex = 0.75)

#"1001000000" #14
points(seq(0,18,2),rep(17,10), 
       pch=c('X','','','X','','','','','',''), cex = 0.75)

#"1001000010" #15
points(seq(0,18,2),rep(2,10), 
       pch=c('X','','','X','','','','','X',''), cex = 0.75)

#"1001000100" #16
points(seq(0,18,2),rep(19,10), 
       pch=c('X','','','X','','','','X','',''), cex = 0.75)

#"1001000110" #17
points(seq(0,18,2),rep(10,10), 
       pch=c('X','','','X','','','','X','X',''), cex = 0.75)

#"1010000000" #18
points(seq(0,18,2),rep(11,10), 
       pch=c('X','','X','','','','','','',''), cex = 0.75)

#"1010000100" #19
points(seq(0,18,2),rep(3,10), 
       pch=c('X','','X','','','','','X','',''), cex = 0.75)

#"1010000110" #20
points(seq(0,18,2),rep(4,10), 
       pch=c('X','','X','','','','','X','X',''), cex = 0.75)

#"1011000000" #21
points(seq(0,18,2),rep(8,10), 
       pch=c('X','','X','X','','','','','',''), cex = 0.75)

#"1011000100" #22
points(seq(0,18,2),rep(15,10), 
       pch=c('X','','X','X','','','','X','',''), cex = 0.75)

#"1011000110" #23
points(seq(0,18,2),rep(14,10), 
       pch=c('X','','X','X','','','','X','X',''), cex = 0.75)

segments(x0 = 20, x1 = 20, y0 =0, y1 = 23.5, xpd = TRUE)

points(((all.results$metric.MI.withy.brier*all.results$prop.pattern) + 20)[sorted.n],
       (1:23), pch = 19, col = 'red')
points(((all.results$metric.MIMI.withy.brier*all.results$prop.pattern) + 20)[sorted.n],
       (1:23), pch = 19, col = 'blue')
points(((all.results$brier.y*all.results$prop.pattern) + 20)[sorted.n],
       (1:23), pch = 19, col = 'darkgreen')
points(((all.results$brier.x*all.results$prop.pattern) + 20)[sorted.n],
       (1:23),pch = 19, col = 'black')

abline(h = c(2:length(all.results$Row.names))-0.5, lty = 2, lwd=0.5, col='gray')

legend(59,12, legend = c("PS", "MIMI", "CCS", "MI"),
       col = c('black','blue','darkgreen','red'),
       pch = c(19,19,19,19), bty = 'n', cex = 0.75,
       x.intersp = 0.2
)

legend(66,13.8, legend = c(round(sum(all.results$brier.x*all.results$prop.pattern),1),
                           round(sum(all.results$metric.MIMI.withy.brier*all.results$prop.pattern),1),
                           round(sum(all.results$brier.y*all.results$prop.pattern),1),
                           round(sum(all.results$metric.MI.withy.brier*all.results$prop.pattern),1)),
       bty='n',
       col = c(NA,NA,NA,NA),
       pch = c(NA,NA,NA,NA), cex = 0.75,
       title = "Model PE")

text(x= par("usr")[3] + .55,y = 1:23 - 0.1 , 
     labels = round(all.results$prop.pattern*9103,2)[sorted.n]
     ,pos = 2, xpd = TRUE, cex = 0.75)

text(par("usr")[3]-2.3, y = 25.7, labels = "N",pos = 1, xpd = TRUE, cex = 0.75)
text(par("usr")[3] + 10.5, y = 25.7, labels = "Pattern",pos = 1, xpd = TRUE, cex = 0.75)
axis(1,at = seq(20,80,2), line = -0.60, tick = seq(20,80,2), labels = seq(0,60,2))
segments(x0 = -1.2, y0 = 0, x1 = -1.2, y1 = 24.5 , xpd = TRUE)
segments(x0 = -1.2, y0 = 0.32, x1 = 20, y1 = 0.32)
segments(x0 = seq(0,18,2) + 1, x1 = seq(0,18,2) + 1, y0 = rep(0.5,10), y1 = rep(23.5,10), lty = 2, lwd=0.5, col='gray')
segments(x0 = -3.5, y0=23.5,y1=23.5,x1=20, xpd = TRUE)
mtext('Pattern Specific Contribution (PE*Pattern Weight)', side = 1, line = 1.55, at = 50, cex =  0.75)


load("supportSM-LogisticAdd10.Rda")

all.results1 <- all.results
load("supportSM-LogisticAdd10Relaxo.Rda")

all.results.relaxo <- all.results
all.results <- all.results1
par(mfrow = c(2,2))
par(mar= c(3,1.5,2.5,0.5), xpd=TRUE)

#layout(matrix(c(1,1,2,2,1,1,2,2), 2, 4, byrow = TRUE), width=c(1,1), respect=FALSE)

#Start with the grid do show which patterns are present
plot(1:10, 1:10, ylim = c(0,23), xlim = c(0,3), pch = NA, bty= 'n', yaxt = 'n', xaxt = 'n', xlab = "", ylab = "", main = 'Support Example (Induced MNAR)   ')
#sorted.n <- sort(round(all.results$metric.MI.withy.prop.pattern*9103,2),index.return=TRUE)$ix
sorted.n <- sort(round(prop.pattern*9103,2),index.return=TRUE)$ix


text(x = seq(0,0.9,.1)+.065 , par("usr")[3] + 0.95 , 
     labels = c('pafi', 'meanbp', 'wblc', 'alb', 'resp', 'temp', 'hrt', 'bili', 'crea', 'sod')
     , srt = 90, pos = 2, xpd = TRUE, cex = 0.75)

#"0000000000" #1
points(seq(0,0.9,.1),rep(23,10), 
       pch=c('','','','','','','','','',''), cex = 0.75)
#"0000000100" #2
points(seq(0,0.9,.1),rep(18,10),
       pch=c('','','','','','','','X','',''), cex = 0.75)
#"0001000000" #3
points(seq(0,0.9,.1),rep(20,10), 
       pch=c('','','','X','','','','','',''), cex = 0.75)
#"0001000010" #4
points(seq(0,0.9,.1),rep(1,10), 
       pch=c('','','','X','','','','','X',''), cex = 0.75)
#"0001000100" #5
points(seq(0,0.9,.1),rep(22,10), 
       pch=c('','','','X','','','','X','',''), cex = 0.75)
#"0001000110" #6
points(seq(0,0.9,.1),rep(5,10), 
       pch=c('','','X','','','','','X','',''), cex = 0.75)
#"0010000000" #7
points(seq(0,0.9,.1),rep(13,10), 
       pch=c('X','','X','X','','','','X','X',''), cex = 0.75)
#"0010000100" #8
points(seq(0,0.9,.1),rep(7,10), 
       pch=c('','','','X','','','','X','X',''), cex = 0.75)
#"0011000000" #9
points(seq(0,0.9,.1),rep(6,10), 
       pch=c('','','X','X','','','','','',''), cex = 0.75)

#"0011000100" #10
points(seq(0,0.9,.1),rep(12,10), 
       pch=c('','','X','X','','','','X','',''), cex = 0.75)

#"0011000110" #11
points(seq(0,0.9,.1),rep(9,10), 
       pch=c('','','X','X','','','','X','X',''), cex = 0.75)

#"1000000000" #12
points(seq(0,0.9,.1),rep(21,10), 
       pch=c('X','','','','','','','','',''), cex = 0.75)

#"1000000100" #13
points(seq(0,0.9,.1),rep(16,10), 
       pch=c('X','','','','','','','X','',''), cex = 0.75)

#"1001000000" #14
points(seq(0,0.9,.1),rep(17,10), 
       pch=c('X','','','X','','','','','',''), cex = 0.75)

#"1001000010" #15
points(seq(0,0.9,.1),rep(2,10), 
       pch=c('X','','','X','','','','','X',''), cex = 0.75)

#"1001000100" #16
points(seq(0,0.9,.1),rep(19,10), 
       pch=c('X','','','X','','','','X','',''), cex = 0.75)

#"1001000110" #17
points(seq(0,0.9,.1),rep(10,10), 
       pch=c('X','','','X','','','','X','X',''), cex = 0.75)

#"1010000000" #18
points(seq(0,0.9,.1),rep(11,10), 
       pch=c('X','','X','','','','','','',''), cex = 0.75)

#"1010000100" #19
points(seq(0,0.9,.1),rep(3,10), 
       pch=c('X','','X','','','','','X','',''), cex = 0.75)

#"1010000110" #20
points(seq(0,0.9,.1),rep(4,10), 
       pch=c('X','','X','','','','','X','X',''), cex = 0.75)

#"1011000000" #21
points(seq(0,0.9,.1),rep(8,10), 
       pch=c('X','','X','X','','','','','',''), cex = 0.75)

#"1011000100" #22
points(seq(0,0.9,.1),rep(15,10), 
       pch=c('X','','X','X','','','','X','',''), cex = 0.75)

#"1011000110" #23
points(seq(0,0.9,.1),rep(14,10), 
       pch=c('','','X','','','','','','',''), cex = 0.75)

segments(x0 = 1, x1 = 1, y0 =0, y1 = 23.5, xpd = TRUE)

points(((all.results$metric.MI.withy.brier) + 1)[sorted.n],
       (1:23), pch = 19, col = 'red',cex = 0.75)
points(((all.results$metric.MIMI.withy.brier) + 1)[sorted.n],
       (1:23), pch = 19, col = 'blue', cex = 0.75)
points(((all.results$brier.y) + 1)[sorted.n],
       (1:23), pch = 19, col = 'darkgreen', cex = 0.75)
points(((all.results.relaxo$brier.y) + 1)[sorted.n],
       (1:23),pch = 19, col = 'darkolivegreen3', cex = 0.75)
points(((all.results.relaxo$brier.x) + 1)[sorted.n],
       (1:23),pch = 19, col = 'dodgerblue', cex = 0.75)

points(((all.results$brier.x) + 1)[sorted.n],
       (1:23),pch = 19, col = 'black', cex = 0.75)



abline(h = c(2:length(all.results$Row.names))-0.5, lty = 2, lwd=0.5, col='gray')

legend(2.05,22.2, legend = c("PS","PS Relaxed Lasso" ,"MIMI", "CCS", "CCS Relaxed Lasso", "MI"),
       col = c('black','dodgerblue','blue','darkgreen','darkolivegreen3','red'),
       pch = c(19,19,19,19,19), bty = 'n', cex = 0.75,
       x.intersp = 0.4
)

legend(2.7,23.7, legend = c(round(sum(all.results$brier.x*all.results$prop.pattern.x),3),
                            round(sum(all.results.relaxo$brier.x*all.results.relaxo$prop.pattern.x),3),
                            round(sum(all.results$metric.MIMI.withy.brier*all.results$metric.MIMI.withy.prop.pattern),3),
                            round(sum(all.results$brier.y*all.results$prop.pattern.y),3),
                            round(sum(all.results.relaxo$brier.y*all.results.relaxo$prop.pattern.y),3),
                            round(sum(all.results$metric.MI.withy.brier*all.results$metric.MI.withy.prop.pattern),3)
),
bty='n',
col = c(NA,NA,NA,NA),
pch = c(NA,NA,NA,NA), cex = 0.75, title = "Model BS")


text(x= -0.05 ,y = 1:23 - 0.1 , 
     labels = round(prop.pattern*9103,0)[sorted.n]
     ,pos = 2, xpd = TRUE, cex = 0.75)

text(-.25, y = 25.5, labels = "N",pos = 1, xpd = TRUE, cex = 0.75)
text(0.5, y = 25.5, labels = "Pattern",pos = 1, xpd = TRUE, cex = 0.75)
axis(1,at = c(-0.5,1), line = -0.80, tick = TRUE, labels = c("",""), cex =0.75)
axis(1,at = seq(1,3,0.2), line = -0.80, tick = seq(1,3,0.2), labels = seq(0,2,0.2), cex =0.75)
segments(x0 =-0.5, y0 = 23.5, x1 = 1, y1 = 23.5 , xpd = TRUE)
#segments(x0 = 0.09, y0 = 0.32, x1 = 1, y1 = 0.32)
segments(x0 = seq(-0.1,0.9,.1) + .05, x1 = seq(-0.1,0.9,.1) + .05, y0 = rep(0.5,10), y1 = rep(23.5,10), lty = 2, lwd=0.5, col='gray')
#segments(x0 = -1, y0=23.5,y1=23.5,x1=1, xpd = TRUE)

mtext('Brier Score', side = 1, line = 1.55, at = 2, cex =  0.75)
####################################################################

#Start with the grid do show which patterns are present
plot(1:10, 1:10, ylim = c(0,23), xlim = c(0,3), pch = NA, bty= 'n', yaxt = 'n', xaxt = 'n', xlab = "", ylab = "", main = 'Support Example (Induced MNAR)   ')
sorted.n <- sort(round(prop.pattern*9103,2),index.return=TRUE)$ix

text(x = seq(0,0.9,.1)+.065 , par("usr")[3] + 0.95 , 
     labels = c('pafi', 'meanbp', 'wblc', 'alb', 'resp', 'temp', 'hrt', 'bili', 'crea', 'sod')
     , srt = 90, pos = 2, xpd = TRUE, cex = 0.75)

#"0000000000" #1
points(seq(0,0.9,.1),rep(23,10), 
       pch=c('','','','','','','','','',''), cex = 0.75)
#"0000000100" #2
points(seq(0,0.9,.1),rep(18,10),
       pch=c('','','','','','','','X','',''), cex = 0.75)
#"0001000000" #3
points(seq(0,0.9,.1),rep(20,10), 
       pch=c('','','','X','','','','','',''), cex = 0.75)
#"0001000010" #4
points(seq(0,0.9,.1),rep(1,10), 
       pch=c('','','','X','','','','','X',''), cex = 0.75)
#"0001000100" #5
points(seq(0,0.9,.1),rep(22,10), 
       pch=c('','','','X','','','','X','',''), cex = 0.75)
#"0001000110" #6
points(seq(0,0.9,.1),rep(5,10), 
       pch=c('','','X','','','','','X','',''), cex = 0.75)
#"0010000000" #7
points(seq(0,0.9,.1),rep(13,10), 
       pch=c('X','','X','X','','','','X','X',''), cex = 0.75)
#"0010000100" #8
points(seq(0,0.9,.1),rep(7,10), 
       pch=c('','','','X','','','','X','X',''), cex = 0.75)
#"0011000000" #9
points(seq(0,0.9,.1),rep(6,10), 
       pch=c('','','X','X','','','','','',''), cex = 0.75)

#"0011000100" #10
points(seq(0,0.9,.1),rep(12,10), 
       pch=c('','','X','X','','','','X','',''), cex = 0.75)

#"0011000110" #11
points(seq(0,0.9,.1),rep(9,10), 
       pch=c('','','X','X','','','','X','X',''), cex = 0.75)

#"1000000000" #12
points(seq(0,0.9,.1),rep(21,10), 
       pch=c('X','','','','','','','','',''), cex = 0.75)

#"1000000100" #13
points(seq(0,0.9,.1),rep(16,10), 
       pch=c('X','','','','','','','X','',''), cex = 0.75)

#"1001000000" #14
points(seq(0,0.9,.1),rep(17,10), 
       pch=c('X','','','X','','','','','',''), cex = 0.75)

#"1001000010" #15
points(seq(0,0.9,.1),rep(2,10), 
       pch=c('X','','','X','','','','','X',''), cex = 0.75)

#"1001000100" #16
points(seq(0,0.9,.1),rep(19,10), 
       pch=c('X','','','X','','','','X','',''), cex = 0.75)

#"1001000110" #17
points(seq(0,0.9,.1),rep(10,10), 
       pch=c('X','','','X','','','','X','X',''), cex = 0.75)

#"1010000000" #18
points(seq(0,0.9,.1),rep(11,10), 
       pch=c('X','','X','','','','','','',''), cex = 0.75)

#"1010000100" #19
points(seq(0,0.9,.1),rep(3,10), 
       pch=c('X','','X','','','','','X','',''), cex = 0.75)

#"1010000110" #20
points(seq(0,0.9,.1),rep(4,10), 
       pch=c('X','','X','','','','','X','X',''), cex = 0.75)

#"1011000000" #21
points(seq(0,0.9,.1),rep(8,10), 
       pch=c('X','','X','X','','','','','',''), cex = 0.75)

#"1011000100" #22
points(seq(0,0.9,.1),rep(15,10), 
       pch=c('X','','X','X','','','','X','',''), cex = 0.75)

#"1011000110" #23
points(seq(0,0.9,.1),rep(14,10), 
       pch=c('','','X','','','','','','',''), cex = 0.75)

segments(x0 = 1, x1 = 1, y0 =0, y1 = 23.5, xpd = TRUE)

points(((all.results$metric.MI.withy.brier*all.results$metric.MI.withy.prop.pattern) + 1)[sorted.n],
       (1:23), pch = 19, col = 'red',cex = 0.75)
points(((all.results$metric.MIMI.withy.brier*all.results$metric.MIMI.withy.prop.pattern) + 1)[sorted.n],
       (1:23), pch = 19, col = 'blue', cex = 0.75)
points(((all.results$brier.y*all.results$prop.pattern.y) + 1)[sorted.n],
       (1:23), pch = 19, col = 'darkgreen', cex = 0.75)
points(((all.results.relaxo$brier.y*all.results$prop.pattern.y) + 1)[sorted.n],
       (1:23),pch = 19, col = 'darkolivegreen3', cex = 0.75)
points(((all.results.relaxo$brier.x*all.results$prop.pattern.x) + 1)[sorted.n],
       (1:23),pch = 19, col = 'dodgerblue', cex = 0.75)
points(((all.results$brier.x*all.results$prop.pattern.x) + 1)[sorted.n],
       (1:23),pch = 19, col = 'black', cex = 0.75)

abline(h = c(2:length(all.results$Row.names))-0.5, lty = 2, lwd=0.5, col='gray')


legend(2.05,22.2, legend = c("PS","PS Relaxed Lasso" ,"MIMI", "CCS", "CCS Relaxed Lasso", "MI"),
       col = c('black','dodgerblue','blue','darkgreen','darkolivegreen3','red'),
       pch = c(19,19,19,19,19), bty = 'n', cex = 0.75,
       x.intersp = 0.4
)

legend(2.7,23.7, legend = c(round(sum(all.results$brier.x*all.results$prop.pattern.x),3),
                            round(sum(all.results.relaxo$brier.x*all.results.relaxo$prop.pattern.x),3),
                            round(sum(all.results$metric.MIMI.withy.brier*all.results$metric.MIMI.withy.prop.pattern),3),
                            round(sum(all.results$brier.y*all.results$prop.pattern.y),3),
                            round(sum(all.results.relaxo$brier.y*all.results.relaxo$prop.pattern.y),3),
                            round(sum(all.results$metric.MI.withy.brier*all.results$metric.MI.withy.prop.pattern),3)
),
bty='n',
col = c(NA,NA,NA,NA),
pch = c(NA,NA,NA,NA), cex = 0.75, title = "Model BS")



text(x= -0.05 ,y = 1:23 - 0.1 , 
     labels = round(prop.pattern*9103,0)[sorted.n]
     ,pos = 2, xpd = TRUE, cex = 0.75)

text(-.25, y = 25.5, labels = "N",pos = 1, xpd = TRUE, cex = 0.75)
text(0.5, y = 25.5, labels = "Pattern",pos = 1, xpd = TRUE, cex = 0.75)
axis(1,at = c(-0.5,1), line = -0.80, tick = TRUE, labels = c("",""), cex =0.75)
axis(1,at = seq(1,3,0.2), line = -0.80, tick = seq(1,3,0.2), labels = seq(0,2,0.2), cex =0.75)
segments(x0 =-0.5, y0 = 23.5, x1 = 1, y1 = 23.5 , xpd = TRUE)
#segments(x0 = 0.09, y0 = 0.32, x1 = 1, y1 = 0.32)
segments(x0 = seq(-0.1,0.9,.1) + .05, x1 = seq(-0.1,0.9,.1) + .05, y0 = rep(0.5,10), y1 = rep(23.5,10), lty = 2, lwd=0.5, col='gray')
#segments(x0 = -1, y0=23.5,y1=23.5,x1=1, xpd = TRUE)

mtext('Pattern Specific Contribution (BS*Pattern Weight)', side = 1, line = 1.55, at = 2, cex =  0.75)

########################################################

#Start with the grid do show which patterns are present
plot(1:10, 1:10, ylim = c(0,23), xlim = c(0,3), pch = NA, bty= 'n', yaxt = 'n', xaxt = 'n', xlab = "", ylab = "", main = 'Support Example (Induced MNAR)')
sorted.n <- sort(round(prop.pattern*9103,2),index.return=TRUE)$ix

text(x = seq(0,0.9,.1)+.065 , par("usr")[3] + 0.95 , 
     labels = c('pafi', 'meanbp', 'wblc', 'alb', 'resp', 'temp', 'hrt', 'bili', 'crea', 'sod')
     , srt = 90, pos = 2, xpd = TRUE, cex = 0.75)

#"0000000000" #1
points(seq(0,0.9,.1),rep(23,10), 
       pch=c('','','','','','','','','',''), cex = 0.75)
#"0000000100" #2
points(seq(0,0.9,.1),rep(18,10),
       pch=c('','','','','','','','X','',''), cex = 0.75)
#"0001000000" #3
points(seq(0,0.9,.1),rep(20,10), 
       pch=c('','','','X','','','','','',''), cex = 0.75)
#"0001000010" #4
points(seq(0,0.9,.1),rep(1,10), 
       pch=c('','','','X','','','','','X',''), cex = 0.75)
#"0001000100" #5
points(seq(0,0.9,.1),rep(22,10), 
       pch=c('','','','X','','','','X','',''), cex = 0.75)
#"0001000110" #6
points(seq(0,0.9,.1),rep(5,10), 
       pch=c('','','X','','','','','X','',''), cex = 0.75)
#"0010000000" #7
points(seq(0,0.9,.1),rep(13,10), 
       pch=c('X','','X','X','','','','X','X',''), cex = 0.75)
#"0010000100" #8
points(seq(0,0.9,.1),rep(7,10), 
       pch=c('','','','X','','','','X','X',''), cex = 0.75)
#"0011000000" #9
points(seq(0,0.9,.1),rep(6,10), 
       pch=c('','','X','X','','','','','',''), cex = 0.75)

#"0011000100" #10
points(seq(0,0.9,.1),rep(12,10), 
       pch=c('','','X','X','','','','X','',''), cex = 0.75)

#"0011000110" #11
points(seq(0,0.9,.1),rep(9,10), 
       pch=c('','','X','X','','','','X','X',''), cex = 0.75)

#"1000000000" #12
points(seq(0,0.9,.1),rep(21,10), 
       pch=c('X','','','','','','','','',''), cex = 0.75)

#"1000000100" #13
points(seq(0,0.9,.1),rep(16,10), 
       pch=c('X','','','','','','','X','',''), cex = 0.75)

#"1001000000" #14
points(seq(0,0.9,.1),rep(17,10), 
       pch=c('X','','','X','','','','','',''), cex = 0.75)

#"1001000010" #15
points(seq(0,0.9,.1),rep(2,10), 
       pch=c('X','','','X','','','','','X',''), cex = 0.75)

#"1001000100" #16
points(seq(0,0.9,.1),rep(19,10), 
       pch=c('X','','','X','','','','X','',''), cex = 0.75)

#"1001000110" #17
points(seq(0,0.9,.1),rep(10,10), 
       pch=c('X','','','X','','','','X','X',''), cex = 0.75)

#"1010000000" #18
points(seq(0,0.9,.1),rep(11,10), 
       pch=c('X','','X','','','','','','',''), cex = 0.75)

#"1010000100" #19
points(seq(0,0.9,.1),rep(3,10), 
       pch=c('X','','X','','','','','X','',''), cex = 0.75)

#"1010000110" #20
points(seq(0,0.9,.1),rep(4,10), 
       pch=c('X','','X','','','','','X','X',''), cex = 0.75)

#"1011000000" #21
points(seq(0,0.9,.1),rep(8,10), 
       pch=c('X','','X','X','','','','','',''), cex = 0.75)

#"1011000100" #22
points(seq(0,0.9,.1),rep(15,10), 
       pch=c('X','','X','X','','','','X','',''), cex = 0.75)

#"1011000110" #23
points(seq(0,0.9,.1),rep(14,10), 
       pch=c('','','X','','','','','','',''), cex = 0.75)

segments(x0 = 1, x1 = 1, y0 =0, y1 = 23.5, xpd = TRUE)

points(((-all.results$metric.MI.withy.log) + 1)[sorted.n],
       (1:23), pch = 19, col = 'red',cex = 0.75)
points(((-all.results$logscore.y) + 1)[sorted.n],
       (1:23), pch = 19, col = 'darkgreen', cex = 0.75)

points(((-all.results.relaxo$logscore.y) + 1)[sorted.n],
       (1:23), pch = 19, col = 'darkolivegreen3', cex = 0.75)
points(((-all.results$metric.MIMI.withy.log) + 1)[sorted.n],
       (1:23), pch = 19, col = 'blue', cex = 0.75)

points(((-all.results.relaxo$logscore.x) + 1)[sorted.n],
       (1:23),pch = 19, col = 'dodgerblue', cex = 0.75)

points(((-all.results$logscore.x) + 1)[sorted.n],
       (1:23),pch = 19, col = 'black', cex = 0.75)


abline(h = c(2:length(all.results$Row.names))-0.5, lty = 2, lwd=0.5, col='gray')


legend(2.05,22.2, legend = c("PS","PS Relaxed Lasso", "MIMI", "CCS", "CC Relax Lasso","MI"),
       col = c('black','dodgerblue','blue','darkgreen','darkolivegreen3','red'),
       pch = c(19,19,19,19,19), bty = 'n', cex = 0.75,
       x.intersp = 0.4
)

legend(2.7,23.7, legend = c(round(sum(-all.results$logscore.x*all.results$prop.pattern.x),3),
                            round(sum(-all.results.relaxo$logscore.x*all.results.relaxo$prop.pattern.x),3),
                            
                            round(sum(-all.results$metric.MIMI.withy.log*all.results$metric.MIMI.withy.prop.pattern),3),
                            round(sum(-all.results$logscore.y*all.results$prop.pattern.y),3),
                            round(sum(-all.results.relaxo$logscore.y*all.results.relaxo$prop.pattern.y),3),
                            round(sum(-all.results$metric.MI.withy.log*all.results$metric.MI.withy.prop.pattern),3)
),
bty='n',
col = c(NA,NA,NA,NA),
pch = c(NA,NA,NA,NA), cex = 0.75, title = "Model LS")


text(x= -0.05 ,y = 1:23 - 0.1 , 
     labels = round(prop.pattern*9103,0)[sorted.n]
     ,pos = 2, xpd = TRUE, cex = 0.75)

text(-.25, y = 25.5, labels = "N",pos = 1, xpd = TRUE, cex = 0.75)
text(0.5, y = 25.5, labels = "Pattern",pos = 1, xpd = TRUE, cex = 0.75)
axis(1,at = c(-0.5,1), line = -0.80, tick = TRUE, labels = c("",""), cex =0.75)
axis(1,at = seq(1,3,0.2), line = -0.80, tick = seq(1,3,0.2), labels = seq(0,2,0.2), cex =0.75)
segments(x0 =-0.5, y0 = 23.5, x1 = 1, y1 = 23.5 , xpd = TRUE)
#segments(x0 = 0.09, y0 = 0.32, x1 = 1, y1 = 0.32)
segments(x0 = seq(-0.1,0.9,.1) + .05, x1 = seq(-0.1,0.9,.1) + .05, y0 = rep(0.5,10), y1 = rep(23.5,10), lty = 2, lwd=0.5, col='gray')
#segments(x0 = -1, y0=23.5,y1=23.5,x1=1, xpd = TRUE)

mtext('Log Score', side = 1, line = 1.55, at = 2, cex =  0.75)




##############################################################################
#Start with the grid do show which patterns are present
plot(1:10, 1:10, ylim = c(0,23), xlim = c(0,3), pch = NA, bty= 'n', yaxt = 'n', xaxt = 'n', xlab = "", ylab = "", main = 'Support Example (Induced MNAR)   ')
sorted.n <- sort(round(prop.pattern*9103,2),index.return=TRUE)$ix

text(x = seq(0,0.9,.1)+.065 , par("usr")[3] + 0.95 , 
     labels = c('pafi', 'meanbp', 'wblc', 'alb', 'resp', 'temp', 'hrt', 'bili', 'crea', 'sod')
     , srt = 90, pos = 2, xpd = TRUE, cex = 0.75)

#"0000000000" #1
points(seq(0,0.9,.1),rep(23,10), 
       pch=c('','','','','','','','','',''), cex = 0.75)
#"0000000100" #2
points(seq(0,0.9,.1),rep(18,10),
       pch=c('','','','','','','','X','',''), cex = 0.75)
#"0001000000" #3
points(seq(0,0.9,.1),rep(20,10), 
       pch=c('','','','X','','','','','',''), cex = 0.75)
#"0001000010" #4
points(seq(0,0.9,.1),rep(1,10), 
       pch=c('','','','X','','','','','X',''), cex = 0.75)
#"0001000100" #5
points(seq(0,0.9,.1),rep(22,10), 
       pch=c('','','','X','','','','X','',''), cex = 0.75)
#"0001000110" #6
points(seq(0,0.9,.1),rep(5,10), 
       pch=c('','','X','','','','','X','',''), cex = 0.75)
#"0010000000" #7
points(seq(0,0.9,.1),rep(13,10), 
       pch=c('X','','X','X','','','','X','X',''), cex = 0.75)
#"0010000100" #8
points(seq(0,0.9,.1),rep(7,10), 
       pch=c('','','','X','','','','X','X',''), cex = 0.75)
#"0011000000" #9
points(seq(0,0.9,.1),rep(6,10), 
       pch=c('','','X','X','','','','','',''), cex = 0.75)

#"0011000100" #10
points(seq(0,0.9,.1),rep(12,10), 
       pch=c('','','X','X','','','','X','',''), cex = 0.75)

#"0011000110" #11
points(seq(0,0.9,.1),rep(9,10), 
       pch=c('','','X','X','','','','X','X',''), cex = 0.75)

#"1000000000" #12
points(seq(0,0.9,.1),rep(21,10), 
       pch=c('X','','','','','','','','',''), cex = 0.75)

#"1000000100" #13
points(seq(0,0.9,.1),rep(16,10), 
       pch=c('X','','','','','','','X','',''), cex = 0.75)

#"1001000000" #14
points(seq(0,0.9,.1),rep(17,10), 
       pch=c('X','','','X','','','','','',''), cex = 0.75)

#"1001000010" #15
points(seq(0,0.9,.1),rep(2,10), 
       pch=c('X','','','X','','','','','X',''), cex = 0.75)

#"1001000100" #16
points(seq(0,0.9,.1),rep(19,10), 
       pch=c('X','','','X','','','','X','',''), cex = 0.75)

#"1001000110" #17
points(seq(0,0.9,.1),rep(10,10), 
       pch=c('X','','','X','','','','X','X',''), cex = 0.75)

#"1010000000" #18
points(seq(0,0.9,.1),rep(11,10), 
       pch=c('X','','X','','','','','','',''), cex = 0.75)

#"1010000100" #19
points(seq(0,0.9,.1),rep(3,10), 
       pch=c('X','','X','','','','','X','',''), cex = 0.75)

#"1010000110" #20
points(seq(0,0.9,.1),rep(4,10), 
       pch=c('X','','X','','','','','X','X',''), cex = 0.75)

#"1011000000" #21
points(seq(0,0.9,.1),rep(8,10), 
       pch=c('X','','X','X','','','','','',''), cex = 0.75)

#"1011000100" #22
points(seq(0,0.9,.1),rep(15,10), 
       pch=c('X','','X','X','','','','X','',''), cex = 0.75)

#"1011000110" #23
points(seq(0,0.9,.1),rep(14,10), 
       pch=c('','','X','','','','','','',''), cex = 0.75)

segments(x0 = 1, x1 = 1, y0 =0, y1 = 23.5, xpd = TRUE)

points(((-all.results$metric.MI.withy.log*all.results$metric.MI.withy.prop.pattern) + 1)[sorted.n],
       (1:23), pch = 19, col = 'red',cex = 0.75)
points(((-all.results$logscore.y*all.results$prop.pattern.y) + 1)[sorted.n],
       (1:23), pch = 19, col = 'darkgreen', cex = 0.75)

points(((-all.results.relaxo$logscore.y*all.results.relaxo$prop.pattern.y) + 1)[sorted.n],
       (1:23), pch = 19, col = 'darkolivegreen3', cex = 0.75)
points(((-all.results$metric.MIMI.withy.log*all.results$metric.MIMI.withy.prop.pattern) + 1)[sorted.n],
       (1:23), pch = 19, col = 'blue', cex = 0.75)

points(((-all.results.relaxo$logscore.x*all.results.relaxo$prop.pattern.x) + 1)[sorted.n],
       (1:23),pch = 19, col = 'dodgerblue', cex = 0.75)
points(((-all.results$logscore.x*all.results$prop.pattern.x) + 1)[sorted.n],
       (1:23),pch = 19, col = 'black', cex = 0.75)


abline(h = c(2:length(all.results$Row.names))-0.5, lty = 2, lwd=0.5, col='gray')


legend(2.05,22.2, legend = c("PS","PS Relaxed Lasso", "MIMI", "CCS", "CC Relax Lasso","MI"),
       col = c('black','dodgerblue','blue','darkgreen','darkolivegreen3','red'),
       pch = c(19,19,19,19,19), bty = 'n', cex = 0.75,
       x.intersp = 0.4
)

legend(2.7,23.7, legend = c(round(sum(-all.results$logscore.x*all.results$prop.pattern.x),3),
                            round(sum(-all.results.relaxo$logscore.x*all.results.relaxo$prop.pattern.x),3),
                            
                            round(sum(-all.results$metric.MIMI.withy.log*all.results$metric.MIMI.withy.prop.pattern),3),
                            round(sum(-all.results$logscore.y*all.results$prop.pattern.y),3),
                            round(sum(-all.results.relaxo$logscore.y*all.results.relaxo$prop.pattern.y),3),
                            round(sum(-all.results$metric.MI.withy.log*all.results$metric.MI.withy.prop.pattern),3)
),
bty='n',
col = c(NA,NA,NA,NA),
pch = c(NA,NA,NA,NA), cex = 0.75, title = "Model LS")

text(x= -0.05 ,y = 1:23 - 0.1 , 
     labels = round(prop.pattern*9103,0)[sorted.n]
     ,pos = 2, xpd = TRUE, cex = 0.75)

text(-.25, y = 25.5, labels = "N",pos = 1, xpd = TRUE, cex = 0.75)
text(0.5, y = 25.5, labels = "Pattern",pos = 1, xpd = TRUE, cex = 0.75)
axis(1,at = c(-0.5,1), line = -0.80, tick = TRUE, labels = c("",""), cex =0.75)
axis(1,at = seq(1,3,0.2), line = -0.80, tick = seq(1,3,0.2), labels = seq(0,2,0.2), cex =0.75)
segments(x0 =-0.5, y0 = 23.5, x1 = 1, y1 = 23.5 , xpd = TRUE)
#segments(x0 = 0.09, y0 = 0.32, x1 = 1, y1 = 0.32)
segments(x0 = seq(-0.1,0.9,.1) + .05, x1 = seq(-0.1,0.9,.1) + .05, y0 = rep(0.5,10), y1 = rep(23.5,10), lty = 2, lwd=0.5, col='gray')

mtext('Pattern Specific Contribution (LS*Pattern Weight)', side = 1, line = 1.55, at = 2, cex =  0.75)



