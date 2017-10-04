#this code aiming at improving Mugnier dynamics data


pdf('Mungnier data.pdf')
#PART2 mugnier dynamics 

#pdf('Mugnier dynamics.pdf')

#get the mugnier dataset 

# source("https://bioconductor.org/biocLite.R")
# biocLite('GenomicRanges')
# biocLite('rtracklayer')

# library(GenomicRanges)
# library(rtracklayer)

BedData<-read.table(file = 'tb144398.bed',header = F)
head(BedData)
par(mfrow=c(1,1))
VSGLib<-BedData$V8-BedData$V7

hist(VSGLib,breaks=20)
abline(v = median(VSGLib),col='red',lty=2)
boxplot(VSGLib)
summary(VSGLib)
# hist(VSGLib, xlim = c(900, 1900), breaks = seq(from=900, to=1900, by=50),
#      col=adjustcolor("blue", alpha.f = 0.5), freq = FALSE, ylim = c(0, 0.0065))
# hist(CrossLib-70, xlim = c(900, 1900), breaks = seq(from=900, to=1900, by=50),
#      col=adjustcolor("red", alpha.f = 0.5), freq = FALSE, add=TRUE)
# 
# t.test(VSGLib,CrossLib)
# wilcox.test(VSGLib,CrossLib)
# qqplot(VSGLib,CrossLib)
# qqline()
#boxplot(VSGLib)
#plot(density(VSGLib))
load('./CompleteVSGlengthCross.R')
CrossLib<-as.numeric(as.character(Output.Complete$VSGLength))

##function to get the length 
GetVSGLength<-function(Target){#function to get Target length from the name and add additional colum to the datasets
        Target<-cbind(Target,as.numeric(unlist(lapply(strsplit(x=as.character(Target[,1]),split ='_' ),FUN =function(x) {unlist(x)[5]}))))
        names(Target)<-c(names(Target)[1:4],'VSGlength')
        return(Target)
}

##get the VSG data from file 
M1D6to30<-read.csv('Mugnier_Database_S1.csv')
M2D6to30<-read.csv('Mugnier_Database_S2.csv')
M3D7to30<-read.csv('Mugnier_Database_S3.csv')
M4D7to30<-read.csv('Mugnier_Database_S4.csv')
M3D96to105<-read.csv('Mugnier_Database_S5.csv')



M1D6to30<-GetVSGLength(M1D6to30)
M2D6to30<-GetVSGLength(M2D6to30)
M3D7to30<-GetVSGLength(M3D7to30)
M4D7to30<-GetVSGLength(M4D7to30)
M3D96to105<-GetVSGLength(M3D96to105)

AllMice<-list(M1D6to30,M2D6to30,M3D7to30,M4D7to30,M3D96to105)


#pdf('Mugnier data.pdf')
hist(M3D7to30$VSGlength[M3D7to30$pct>0.001], freq = FALSE, col = "red")
hist(M3D96to105$VSGlength[M3D96to105$pct>0.001], add=TRUE, freq = FALSE)                            

Beg <- M3D7to30$VSGlength[M3D7to30$pct>0.001]             
End <- M3D96to105$VSGlength[M3D96to105$pct>0.001]
Beg <- Beg[!is.na(Beg)]
End <- End[!is.na(End)]

T1 <- table(cut(Beg, breaks = c(0, 1200, 1400, 1600, Inf)))
T2 <- table(cut(End, breaks = c(0, 1200, 1400, 1600, Inf)))

chisq.test(rbind(T1, T2))

sum(Beg < mean(Beg))/length(Beg)
sum(End < mean(End))/length(End)

#dominating clone 

for (i in 1:length(AllMice)){
        Target<-AllMice[[i]]
        head(Target)
        Target<-Target[Target$VSG!='parasitemia',]  
        TargetVec<-split(Target,f = Target$day)
        VecPlot1<-lapply(TargetVec,function(x){x[which.max(x$pct),'VSGlength']})
        CorVec<-cor(unlist(VecPlot1[-1]),as.numeric(names(unlist(VecPlot1[-1]))))
        
        par(mfrow=c(2,1))
        
        plot(as.numeric(names(unlist(VecPlot1))),unlist(VecPlot1),xlab='day',ylab ='Length of dominating VSG' ,pch=4,ylim=c(1200,1800),
             main=paste('Sample',i,'Dominating VSG clones, ','Cor',signif(x=CorVec,digits = 3)))
        
        
        LM1<-lm(unlist(VecPlot1[-1])~as.numeric(names(unlist(VecPlot1[-1]))))
        abline(LM1,lty=2,lwd=2,col='red')
        
        LM2<-lm(unlist(VecPlot1[c(1:3)])~as.numeric(names(unlist(VecPlot1[c(1:3)]))))
        abline(LM2,lty=2,lwd=2,col='blue')
        
        abline(h = mean(VSGLib),col='brown',lty=3)
        
        TargetParasitemia<-AllMice[[i]][AllMice[[i]]$VSG=='parasitemia',]
        plot(TargetParasitemia$day,TargetParasitemia$parasites,xlab='day',ylab ='parasitemia' ,
             main=paste('mouse',i,'parasitemia'),type = 'l',col='red',lwd=2)
        
        par(mfrow=c(1,1))
        
}
par(mfrow=c(1,1))

#vsg length of dominating clones
#pdf('Mugnier Domiating clone.pdf')
par(mfrow=c(2,2),mar=c(2,2.5,2,2))    
for (i in 1:4){
        Target<-AllMice[[i]]
        head(Target)
        Target<-Target[Target$VSG!='parasitemia'&Target$pct>0.001,]  
        TargetVec<-split(Target,f = Target$day)
        VecPlot1<-lapply(TargetVec,function(x){x[which.max(x$pct),'VSGlength']})
        CorVec<-cor(unlist(VecPlot1[-1]),as.numeric(names(unlist(VecPlot1[-1]))))
        VecDays<-as.numeric(names(unlist(VecPlot1)))
        CorVec<-cor(unlist(VecPlot1[c(1:3)]),VecDays[c(1:3)])
        CorVec2<-cor(unlist(VecPlot1[-c(1:3)]),VecDays[-c(1:3)])
        
        
        
        plot(as.numeric(names(unlist(VecPlot1))),unlist(VecPlot1),type='p',cex.main=2,cex.axis=1.5 ,pch=19,ylim=c(1200,1800),
             main=paste('Mouse',i),xlab='',ylab='')
        
        abline(h = median(VSGLib),col='brown',lty=3,lwd=3) 
        
        if (i!=3){# mice 3 is different from others 
                LM1<-lm(unlist(VecPlot1)[c(1:3)]~VecDays[c(1:3)])
                
                clip(x1 = 2,y1 = 0,x2 = 16,y2 = 2000) # resitrict the length of line 
                
                abline(LM1,lty=2,lwd=2,col='red')
                
                
                clip(x1 = 16,y1 = 0,x2 = 32,y2 = 2000) 
                LM2<-lm(unlist(VecPlot1)[-c(1:3)]~VecDays[-c(1:3)])
                abline(LM2,lty=2,lwd=2,col='blue')
                
                abline(h = median(VSGLib),col='brown',lty=3) }
        
        
        if (i==3){# mice 3 is different from others 
                LM1<-lm(unlist(VecPlot1)[c(1:3)]~VecDays[c(1:3)])
                
                clip(x1 = 2,y1 = 0,x2 = 16,y2 = 2000) # resitrict the length of line 
                
                abline(LM1,lty=2,lwd=2,col='red')
                
                
                clip(x1 = 16,y1 = 0,x2 = 24,y2 = 2000) 
                LM2<-lm(unlist(VecPlot1)[c(4:6)]~VecDays[c(4:6)])
                abline(LM2,lty=2,lwd=2,col='blue')
                
                abline(h = median(VSGLib),col='brown',lty=3) 
                
                clip(x1 = 24,y1 = 0,x2 = 32,y2 = 2000) 
                LM3<-lm(unlist(VecPlot1)[c(6:8)]~VecDays[c(6:8)])
                abline(LM3,lty=2,lwd=2,col='red')
                
                
                
        }
        
        
        
}
par(mfrow=c(1,1))
#dev.off()


#distribution of VSG length on each mouse

DesLib<-density(VSGLib)
DesLib$y<-80*DesLib$y

for (i in 1:length(AllMice)){
        Target<-AllMice[[i]]
        Target<-GetVSGLength(Target)
        Target<-Target[Target$VSG!='parasitemia'&Target$pct>0.001,] 
        
        AllVSG<-sapply(as.character(unique(Target$VSG)),FUN = function(x){unique(Target[Target$VSG==x[1],]$VSGlength)})
        par(mfrow=c(3,3))
        # plot(density(AllVSG))
        boxplot(VSGLib,horizontal = T)
        
        sapply(split(Target,Target$day),FUN = function(x){DesVecT<-density(x$VSGlength,weights=x$pct);
        plot(DesVecT,xlim=range(VSGLib),col='blue',ylim=c(0,1),
             main=paste('Sample',i,'Day',unique(x$day)),xlab='VSGlength',ylab='pct');
        abline(v = median(VSGLib),col='red',lty=2);lines(DesLib,col='red')} )
        
        par(mfrow=c(1,1))
}


##histogram version


for (i in 1:4){
        Target<-AllMice[[i]]
        Target<-GetVSGLength(Target)
        Target<-Target[Target$VSG!='parasitemia'&Target$pct>0.001,] 
        
        AllVSG<-sapply(as.character(unique(Target$VSG)),FUN = function(x){unique(Target[Target$VSG==x[1],]$VSGlength)})
        par(mfrow=c(3,3))
        # plot(density(AllVSG))
        boxplot(VSGLib,horizontal = T)
        
        sapply(split(Target,Target$day),FUN = function(x){
                plot(x = x$VSGlength ,y = x$pct,xlim=range(VSGLib),col='blue',type='h',ylim=c(0,100),
                     main=paste('Mouse',i,'Day',unique(x$day)),xlab='VSGlength',ylab='pct');
                abline(v = median(CrossLib),col='red',lty=2)} )
        
        par(mfrow=c(1,1))
}


#density version
library(sfsmisc)
for (i in 1:4){
        Target<-AllMice[[i]]
        Target<-GetVSGLength(Target)
        Target<-Target[Target$VSG!='parasitemia'&Target$pct>0.001,] 
        
        AllVSG<-sapply(as.character(unique(Target$VSG)),FUN = function(x){unique(Target[Target$VSG==x[1],]$VSGlength)})
        par(mfrow=c(3,3),mar=c(4.2,5,4,4))
        # plot(density(AllVSG))
        boxplot(VSGLib,horizontal = T,ylim = range(VSGLib),xlab='VSG length (bp)',cex.lab=1.5,cex.axis=1.5)
        head(Target)
        VecPlot<-lapply(split(Target,f = Target$day),FUN = function(x){rep(x$VSGlength,x$pct/min(x$pct))})
        
        #pdf('Mouse1density.pdf')
       # par(mfrow=c(3,3))
        VecDens<-lapply(VecPlot,density,adjust=5)
        VecDens<-lapply(VecDens,function(x){x$y=x$y*100;return(x)})
        for (j in 1:length(VecPlot)){
         
        plot(VecDens[[j]],main=paste0("Mouse ",i,' Day ',names(VecDens)[j]),log='',
             xlab='VSG length (bp)',ylab='Kernel density /1e-2',xlim = range(VSGLib),cex.lab=1.5,cex.axis=1.5,cex.main=2)
        
          polygon(VecDens[[j]],col='red',border = 'red')
        abline(v = median(VSGLib),col='blue',lty=2,lwd=2)
        abline(v = median(VecPlot[[j]]),col='green',lty=2,lwd=2)}
        
      
        #dev.off()
      
        par(mfrow=c(1,1))
}




#look at weighted mean Aand improve the plot 

#setwd("~/Box Sync/research/Third year 3/multiple gene copies/Record/20160801/Mugnier data")

#pdf('WeightedMean fited line.pdf')
##this is the function to calculate std of weighted value o
WeightedSD<-function(x,w){
        WtSd<-sqrt(sum(var(x)*(w^2)))
        return(WtSd)        
}
Wt.mean<-function(x,weight){
        if (!is.na(weight)&sum(weight)<=1){return(sum(x*weight))} else {return(NaN);warning('error with wtmean')}
        
}



par(mfrow=c(1,1))    
for (i in 1:4){
        Target<-AllMice[[i]]
        head(Target)
        Target<-Target[Target$VSG!='parasitemia'&Target$pct>0.001,]  
        TargetVec<-split(Target,f = Target$day)
        VecPlot1<-lapply(TargetVec,function(x){weighted.mean(x =x$VSGlength,w =x$pct/100)})#this is where it is different 
        VecDays<-as.numeric(names(unlist(VecPlot1)))
        
        CorVec<-cor(unlist(VecPlot1[c(1:3)]),VecDays[c(1:3)])
        CorVec2<-cor(unlist(VecPlot1[-c(1:3)]),VecDays[-c(1:3)])
        
        
        
        plot(as.numeric(names(unlist(VecPlot1))),unlist(VecPlot1),type='p',xlab='day',ylab ='Weighted mean of VSG' ,pch=19,ylim=c(1200,1800),
             main=paste('Mouse',i))
        WtSd<-sapply(TargetVec,FUN = function(x){WeightedSD(x =x$VSGlength,w =x$pct/100)})#weighted sd
        
        arrows(as.numeric(names(unlist(VecPlot1))), unlist(VecPlot1)-WtSd,
               as.numeric(names(unlist(VecPlot1))), unlist(VecPlot1)+WtSd, length=0.05, angle=90, code=3)
        abline(h = median(VSGLib),col='brown',lty=3) 
        
        if (i!=3){# mice 3 is different from others 
                LM1<-lm(unlist(VecPlot1)[c(1:3)]~VecDays[c(1:3)])
                
                clip(x1 = 2,y1 = 0,x2 = 16,y2 = 2000) # resitrict the length of line 
                
                abline(LM1,lty=2,lwd=2,col='red')
                
                
                clip(x1 = 16,y1 = 0,x2 = 32,y2 = 2000) 
                LM2<-lm(unlist(VecPlot1)[-c(1:3)]~VecDays[-c(1:3)])
                abline(LM2,lty=2,lwd=2,col='blue')
                
                abline(h = median(VSGLib),col='brown',lty=3) }
        
        
        if (i==3){# mice 3 is different from others 
                LM1<-lm(unlist(VecPlot1)[c(1:3)]~VecDays[c(1:3)])
                
                clip(x1 = 2,y1 = 0,x2 = 16,y2 = 2000) # resitrict the length of line 
                
                abline(LM1,lty=2,lwd=2,col='red')
                
                
                clip(x1 = 16,y1 = 0,x2 = 24,y2 = 2000) 
                LM2<-lm(unlist(VecPlot1)[c(4:6)]~VecDays[c(4:6)])
                abline(LM2,lty=2,lwd=2,col='blue')
                
                abline(h = median(VSGLib),col='brown',lty=3) 
                
                clip(x1 = 24,y1 = 0,x2 = 32,y2 = 2000) 
                LM3<-lm(unlist(VecPlot1)[c(6:8)]~VecDays[c(6:8)])
                abline(LM3,lty=2,lwd=2,col='red')
                
                
                
        }
        
        
        
}
par(mfrow=c(1,1))

#dev.off()
        #this part generate the growth model simulation 
        #VSGLib<-VSGLib[VSGLib>1000]# get rid of those potential incomplete VSG     
#         
#         VecTiSimul<-Target[Target$day==min(Target$day),c('pct','VSGlength')]
#         VecTiSimul[,'pct']<-ceiling(VecTiSimul[,'pct']*1)
#         
#         setwd("~/Box Sync/research/Third year 3/multiple gene copies/Record/20160427")
#         source('./20160328 growth model function2.R')
#         SimuGrowthResult<-SimuDifGrowth(VecTiSimul,VSGLib)#simulate grwoth model 
#         
#         WtMean<-mapply(1:nrow(SimuGrowthResult),
#                        FUN = function(x){Wt.mean(VSGLib,weight = SimuGrowthResult[x,-c(1,2)]/SimuGrowthResult[x,2])})
#         names(WtMean)<-SimuGrowthResult[,1]
#         
#         plot(SimuGrowthResult[,1]/4,WtMean,type='l',xlab = 'days',ylab='Weigthed mean of VSG length',
#              main='Differential growth model',ylim=c(1200,1800))
#         abline(h = median(VSGLib),col='blue',lty=2)
#         
#         #this part generate the null model simulation 
#         
#         VecTiSimul<-Target[Target$day==min(Target$day),c('pct','VSGlength')]
#         VecTiSimul[,'pct']<-ceiling(VecTiSimul[,'pct']*10)
#         
#         setwd("~/Box Sync/research/Third year 3/multiple gene copies/Record/20160427")
#         source('./20160328 null model function.R')
#         SimuNUllResult<-SimuGrowthNULL(VecTiSimul,VSGLib)#simulate grwoth model 
#         
#         WtMean<-mapply(1:nrow(SimuNUllResult),
#                        FUN = function(x){Wt.mean(VSGLib,weight = SimuNUllResult[x,-c(1,2)]/SimuNUllResult[x,2])})
#         names(WtMean)<-SimuNUllResult[,1]
#         
#         plot(SimuNUllResult[,1]/4,WtMean,type='l',xlab = 'days',ylab='Weigthed mean of VSG length',
#              main='Non-differential growth model',ylim=c(1200,1800))
#         abline(h = median(VSGLib),col='brown',lty=2)
        
        
        
        #TargetParasitemia<-AllMice[[i]][AllMice[[i]]$VSG=='parasitemia',]
        #plot(TargetParasitemia$day,TargetParasitemia$parasites,xlab='day',ylab ='parasitemia' ,
        #main=paste('mouse',i,'parasitemia'),type = 'l',col='red',lwd=2)
        


#dev.off()

#explore detailed balance 

# Target<-AllMice[[1]]
# Target<-Target[Target$VSG!='parasitemia',]
# Target<-Target[Target$pct>0.001,]
# unique(Target$day)
# Target<-Target[Target$day==30,]
# unique(Target$day)
# 
# par(mfrow=c(1,1))
# 
# Vec1<-hist(Target$VSGlength,col = adjustcolor('red',alpha.f = 0.5),freq = F,xlim=range(VSGLib))
# Vec2<-hist(VSGLib,col = adjustcolor('blue',alpha.f = 0.5),add=T,freq = F)
# wilcox.test(Target$VSGlength,VSGLib)




#this is to see whether other central tendency trend will do better than weighter mean 











#explore whether there is a weak heirachy , I could not get the answer yet. 

Target<-AllMice[[1]]
Target<-Target[Target$VSG!='parasitemia',]
Target<-Target[Target$day<=16,]
Target<-Target[Target$pct>0.001,]


head(Target)
Target<-split(Target,Target$VSG)
Target<-Target[-length(Target)]
#Target<-Target[which(!sapply(Target,function(x){30%in%x$day}))]



#sapply(Target,function(x){x$day[which.max(x$pct)]}
#sapply(Target,function(x){max(x$pct)})
#summary(sapply(Target,function(x){max(x$parasites)}))
# sapply(Target,function(x){min(x$day[which(x$pct>0.5)])})
# 
# sapply(Target,function(x){min(x$day[(x$pct>0.01)])})
# sapply(Target,function(x){x$day[which.max(x$parasites)]})
# sapply(Target,function(x){min(x$day[(x$pct>0.01)])})
# sapply(Target,function(x){min(x$day[which(x$parasites>1000000)])})
# 
# 
# VecPlot<-cbind(sapply(Target,function(x){unique(x$VSGlength)}),sapply(Target,function(x){x$day[which.max(x$parasites)]}))
# plot(VecPlot)
# 
# cor(VecPlot[,1],VecPlot[,2],method = 'spearman')
# cor.test(VecPlot[,1],VecPlot[,2],method = 'spearman')

# # look at number of VSG per day 
# #pdf('Number of detected VSG.pdf')
# par(mfrow=c(2,2))
# for (i in 1:4){
#         Target<-AllMice[[i]]
#         head(Target)
#         Target<-Target[Target$VSG!='parasitemia',]
#         Target<-Target[Target$pct>0.001,]
#         Target<-split(Target,Target$day)
#         names(Target)
#         VecPlot<-data.frame(Time=as.numeric(names(Target)),NVSG=sapply(Target,nrow))
#         plot(VecPlot$Time,VecPlot$NVSG,pch=3,cex=3,xlab='Days after infection',ylab='Number of different VSG detected ',
#              main=paste('Mouse',i))
#         
#         LM1<-lm(VecPlot$NVSG[VecPlot$Time<=21]~VecPlot$Time[VecPlot$Time<=21])
#         clip(x1 = 0,y1 = 10,x2 = 21,y2 = 100)
#         abline(LM1,col='brown',lty=2,lwd=3)
#         
#         LM2<-lm(VecPlot$NVSG[VecPlot$Time==21|VecPlot$Time==24]~c(21,24))
#         clip(x1 = 21,y1 = 100,x2 = 24,y2 = 10)
#         abline(LM2,col='blue',lty=2,lwd=3)}
# par(mfrow=c(1,1))
# segments(x0 = 21,y=VecPlot$NVSG[VecPlot$Time==21],x1 = 24,y1=VecPlot$NVSG[VecPlot$Time==24])





#loess version of Wtmean 
#look at weighted mean Aand improve the plot 

#setwd("~/Box Sync/research/Third year 3/multiple gene copies/Record/20160427")
##this is the function to calculate std of weighted value o
WeightedSD<-function(x,w){
        WtSd<-sqrt(sum(var(x)*(w^2)))
        return(WtSd)        
}
Wt.mean<-function(x,weight){
        if (all(!is.na(weight))&sum(weight)<=1){return(sum(x*weight)/sum(weight))} else {return(NaN);warning('error with wtmean')}
        
}



par(mfrow=c(2,2))    
for (i in 1:4){
        Target<-AllMice[[i]]
        head(Target)
        Target<-Target[Target$VSG!='parasitemia'&Target$pct>0.001,]  
        TargetVec<-split(Target,f = Target$day)
        VecPlot1<-unlist(lapply(TargetVec,function(x){weighted.mean(x =x$VSGlength,w =x$pct/100)}))#this is where it is different 
        VecDays<-as.numeric(names(unlist(VecPlot1)))
        LVec<-loess(VecPlot1~VecDays)
       
        CorVec<-cor(unlist(VecPlot1[c(1:3)]),VecDays[c(1:3)])
        CorVec2<-cor(unlist(VecPlot1[-c(1:3)]),VecDays[-c(1:3)])
        
        
        
        plot(as.numeric(names(unlist(VecPlot1))),unlist(VecPlot1),type='p',xlab='day',ylab ='Weighted mean of VSG' ,pch=19,ylim=c(1200,1800),
             main=paste('Mouse',i),cex.lab=1,cex.main=4)
        lines( predict(LVec,1:30),lty=2, lwd=2,col='brown')
        WtSd<-sapply(TargetVec,FUN = function(x){WeightedSD(x =x$VSGlength,w =x$pct/100)})#weighted sd
        
        arrows(as.numeric(names(unlist(VecPlot1))), unlist(VecPlot1)-WtSd,
               as.numeric(names(unlist(VecPlot1))), unlist(VecPlot1)+WtSd, length=0.05, angle=90, code=3)
        abline(h = median(VSGLib),col='brown',lty=3) 
        
        
        
        
        
}
par(mfrow=c(1,1))


#loess version dominating clone 
#vsg length of dominating clones
#setwd("~/Box Sync/research/Third year 3/multiple gene copies/Record/20160323")
#pdf('Mugnier Domiating clone.pdf')
par(mfrow=c(2,2))    
for (i in 1:4){
        Target<-AllMice[[i]]
        head(Target)
        Target<-Target[Target$VSG!='parasitemia'&Target$pct>0.001,]  
        TargetVec<-split(Target,f = Target$day)
        VecPlot1<-unlist(lapply(TargetVec,function(x){x[which.max(x$pct),'VSGlength']}))
        CorVec<-cor(unlist(VecPlot1[-1]),as.numeric(names(unlist(VecPlot1[-1]))))
        VecDays<-as.numeric(names(unlist(VecPlot1)))
        LVec<-loess(VecPlot1~VecDays)
        
        CorVec<-cor(unlist(VecPlot1[c(1:3)]),VecDays[c(1:3)])
        CorVec2<-cor(unlist(VecPlot1[-c(1:3)]),VecDays[-c(1:3)])
        
        
        
        plot(as.numeric(names(unlist(VecPlot1))),unlist(VecPlot1),type='p',xlab='day',ylab ='VSG length of dominating clone' ,pch=19,ylim=c(1200,1800),
             main=paste('Mouse',i),cex.lab=1,cex.main=4)
        lines( predict(LVec,1:30),lty=2, lwd=2,col='brown')
        
        abline(h = median(VSGLib),col='brown',lty=3,lwd=3) 
        
       
        
        
        
}
par(mfrow=c(1,1))

#dev.off()


#this is to see whether other central tendency trend will do better than weighter mean 
#Firstly lt's try weighted mean of the top 10% types of the patasited 

#("~/Box Sync/research/Third year 3/multiple gene copies/Record/20160801/Mugnier data")
#pdf('wt mean of top VSG.pdf')
par(mfrow=c(2,2))    
for (i in 1:4){
        Target<-AllMice[[i]]
        head(Target)
        Target<-Target[Target$VSG!='parasitemia'&Target$pct>0.001,]  
        TargetVec<-split(Target,f = Target$day)
        round(length(unique(Target$VSG))/10)
       Target$pct
        VecPlot1<-unlist(lapply(TargetVec,function(x){VecID<-order(x$pct)[ 1:round(length(unique(x$VSG))/10)]
                                                      Vecresult<-Wt.mean(x = x$VSGlength[VecID],weight =x$pct[VecID]/100)}))
        #VecPlot1<-unlist(lapply(TargetVec,function(x){x[which.max(x$pct),'VSGlength']}))
        CorVec<-cor(unlist(VecPlot1[-1]),as.numeric(names(unlist(VecPlot1[-1]))))
        VecDays<-as.numeric(names(unlist(VecPlot1)))
        LVec<-loess(VecPlot1~VecDays)
        
        CorVec<-cor(unlist(VecPlot1[c(1:3)]),VecDays[c(1:3)])
        CorVec2<-cor(unlist(VecPlot1[-c(1:3)]),VecDays[-c(1:3)])
        
        
        
        plot(as.numeric(names(unlist(VecPlot1))),unlist(VecPlot1),type='p',xlab='day',ylab ='WT mean of top VSGs' ,pch=19,ylim=c(1200,1800),
             main=paste('Mouse',i),cex.lab=1,cex.main=4)
        lines( predict(LVec,1:30),lty=2, lwd=2,col='brown')
        
        abline(h = median(VSGLib),col='brown',lty=3,lwd=3) }
        
        
#    dev.off()   
        
        
        



#now try median length 

#fuction to ge the median 

TargetMedian<-function(Target){
        Target<-Target[order(Target$VSGlength,decreasing = F),]
        MLength<-Target$VSGlength[which(cumsum(Target$pct/100)>0.5)[1]]
        return(MLength)
}




#setwd("~/Box Sync/research/Third year 3/multiple gene copies/Record/20160801/Mugnier data")
#pdf('median VSG.pdf')
par(mfrow=c(2,2))    
for (i in 1:4){
        Target<-AllMice[[i]]
        head(Target)
        Target<-Target[Target$VSG!='parasitemia'&Target$pct>0.001,]  
        TargetVec<-split(Target,f = Target$day)
        round(length(unique(Target$VSG))/10)
        Target$pct
        VecPlot1<-unlist(lapply(TargetVec,function(x){TargetMedian(x)}))
        #VecPlot1<-unlist(lapply(TargetVec,function(x){x[which.max(x$pct),'VSGlength']}))
        CorVec<-cor(unlist(VecPlot1[-1]),as.numeric(names(unlist(VecPlot1[-1]))))
        VecDays<-as.numeric(names(unlist(VecPlot1)))
        LVec<-loess(VecPlot1~VecDays)
        
        CorVec<-cor(unlist(VecPlot1[c(1:3)]),VecDays[c(1:3)])
        CorVec2<-cor(unlist(VecPlot1[-c(1:3)]),VecDays[-c(1:3)])
        
        
        
        plot(as.numeric(names(unlist(VecPlot1))),unlist(VecPlot1),type='p',xlab='day',ylab ='median VSG length' ,pch=19,ylim=c(1200,1800),
             main=paste('Mouse',i),cex.lab=1,cex.main=4)
        lines( predict(LVec,1:30),lty=2, lwd=2,col='brown')
        
        abline(h = median(VSGLib),col='brown',lty=3,lwd=3) }

#dev.off()



#plot scatter plot but use circle size to suggest weight 


dfx = data.frame(ev1=1:10, ev2=sample(10:99, 10), ev3=10:1)

with(dfx, symbols(x=ev1, y=ev2, circles=ev3, inches=1/3,
                  ann=F, bg="steelblue2", fg=NULL))


##this is the function to calculate std of weighted value o
WeightedSD<-function(x,w){
        WtSd<-sqrt(sum(var(x)*(w^2)))
        return(WtSd)        
}
Wt.mean<-function(x,weight){
        if (!is.na(weight)&sum(weight)<=1){return(sum(x*weight))} else {return(NaN);warning('error with wtmean')}
        
}





#pdf('symbol.pdf')

par(mfrow=c(2,2),mar=c(2,2.5,2,2))    
for (i in 1:4){
        Target<-AllMice[[i]]
        head(Target)
        Target<-Target[Target$VSG!='parasitemia'&Target$pct>0.001,]  
        TargetVec<-split(Target,f = Target$day)
        #VecPlot1<-lapply(TargetVec,function(x){weighted.mean(x =x$VSGlength,w =x$pct/100)})#this is where it is different 
        VecDays<-as.numeric(names(unlist(VecPlot1)))
        
        CorVec<-cor(unlist(VecPlot1[c(1:3)]),VecDays[c(1:3)])
        CorVec2<-cor(unlist(VecPlot1[-c(1:3)]),VecDays[-c(1:3)])
        
        #plot(as.numeric(unlist(Target$VSGlength)),unlist(Target$day),type='p',xlab='day',ylab ='Weighted mean of VSG' ,pch=19,
          #   main=paste('Mouse',i))
        CircleSize<-Target$pct#ensure people can see the small populations 
        CircleSize[CircleSize<5]<-5
        with(Target, symbols(x=day, y=VSGlength, circles=CircleSize, inches=1/3,cex.axis=1.5,cex.main=2,
                          bg=adjustcolor(col = "steelblue2",alpha.f = 0.5), fg=NULL,xlab='',ylab ='',main=paste('Mouse',i),ylim=c(1200,1800)))
       
        # plot(as.numeric(names(unlist(VecPlot1))),unlist(VecPlot1),type='p',xlab='day',ylab ='Weighted mean of VSG' ,pch=19,ylim=c(1200,1800),
             #main=paste('Mouse',i))
        WtSd<-sapply(TargetVec,FUN = function(x){WeightedSD(x =x$VSGlength,w =x$pct/100)})#weighted sd
        
        #arrows(as.numeric(names(unlist(VecPlot1))), unlist(VecPlot1)-WtSd,
               #as.numeric(names(unlist(VecPlot1))), unlist(VecPlot1)+WtSd, length=0.05, angle=90, code=3)
        abline(h = median(VSGLib),col='brown',lty=3) 
        
        # if (i!=3){# mice 3 is different from others 
                LM1<-lm(unlist(VecPlot1)[c(1:3)]~VecDays[c(1:3)])
                
                clip(x1 = 2,y1 = 0,x2 = 16,y2 = 2000) # resitrict the length of line 
                
                abline(LM1,lty=2,lwd=2,col='red')
                
                
                clip(x1 = 16,y1 = 0,x2 = 32,y2 = 2000) 
                LM2<-lm(unlist(VecPlot1)[-c(1:3)]~VecDays[-c(1:3)])
                abline(LM2,lty=2,lwd=2,col='blue')
                
                abline(h = median(VSGLib),col='brown',lty=3) #}
        
        
#         if (i==3){# mice 3 is different from others 
#                 LM1<-lm(unlist(VecPlot1)[c(1:3)]~VecDays[c(1:3)])
#                 
#                 clip(x1 = 2,y1 = 0,x2 = 16,y2 = 2000) # resitrict the length of line 
#                 
#                 abline(LM1,lty=2,lwd=2,col='red')
#                 
#                 
#                 clip(x1 = 16,y1 = 0,x2 = 24,y2 = 2000) 
#                 LM2<-lm(unlist(VecPlot1)[c(4:6)]~VecDays[c(4:6)])
#                 abline(LM2,lty=2,lwd=2,col='blue')
#                 
#                 abline(h = median(VSGLib),col='brown',lty=3) 
#                 
#                 clip(x1 = 24,y1 = 0,x2 = 32,y2 = 2000) 
#                 LM3<-lm(unlist(VecPlot1)[c(6:8)]~VecDays[c(6:8)])
#                 abline(LM3,lty=2,lwd=2,col='red')
                
                
                
       # }
        
        
        
}
par(mfrow=c(1,1))

dev.off()

