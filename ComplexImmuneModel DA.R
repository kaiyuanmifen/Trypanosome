#this version include a demo of growth rate and immune detection time of all clones 
#this version tries 3 different maximum load and put back in the maximumload 
#this version include the immune data
#Look at VSG dependent Differential growth 
#this version explore parameters 

#Use DeSolve , chekc whether it works 
#install.packages('deSolve')
library(deSolve)
#install.packages('deSolve')

#Trial<-function(t,x,parms){return(list(1))}
#Trial<-function(t,x,parms){return(list(exp(x = -x)))}

# Output<-ode(times = seq(0,1,0.1),func = Trial,y=c(0,0),parm=1,method = 'rk4')
# plot(Output[,2]~Output[,1],type='l')
#pdf('ProMaxDelay.pdf')
#load('Parameters.Rdata')
load('MugnierVSGLib.Rdata')
LAll=LAll+rnorm(n = length(LAll))#introduction of a little bit random noise to avoid messing up of
hist(LAll,col='brown',xlab='VSG length(bp)',ylab='Counts' ,main='',
     cex.lab=2,cex.main=2,cex.axis=1.5)

#class(ParaSpace)
#ParaSpace<-data.frame(t(ParaSpace))
#ParaSpace[1100,]
pdf('DA CM.pdf')


DS<-0
DG<-0
plot(2-DG*(LAll-median(LAll))~LAll)
hist(2-DG*(LAll-median(LAll)))


DST<-10
1/length(LAll)+DST*(LAll-median(LAll))
Threshold<-1e3
KillingRate<-8
ControlV<-1
MaximumLoad=1e6
#the required input include LAll(Library), DS,DG, DST, Threshold,KillingRate


ImmundetectionRecord<-NULL #vector to dectect length of VSG that has been dectected, this inreoduce an saasumption that each VSG length is unique 
DetectionTime<-NULL
ImmuneKillStarted=NULL
#ImmuneRecordMatrix<-NULL#matrix to record which VSG got detected each time point 


PopuDynamic<-function(t,Nl,l){
        
        DS<-DS#diffrentiation factor for switch, range  range c(0,2e-5)
        Swi<-function(l){return(4e-4 + DS*(l-median(LAll)))}
        #Swi<-function(l){1e-5}
        #Swi<-function(l){(1e-7)*l}#dfferential switch rate 
        #Swi<-function(l){return(0)}
        #
        
        
        #MaxLoad<-10^99
        DG<-DG#factor to describe effect from length of VSG in growth range c(0,0.008)
        R<-function(l){2-DG*(l-median(LAll))}#new version of differentiation growth 
        #R<-function(l){(1-(sum(Nl)/MaxLoad))*4*(median(LAll)/l)}
        #R<-function(l){return(4+Diffac*(median(LAll)-l))}#I abondoned the max load 
        #R<-function(l){4}
        
        DST<-DST #diffrentiation activation for switch, range  range c(0,1/500*N)
        SwTo<-function(l){return((DST+max(LAll)-l)/sum((DST+max(LAll)-LAll)))}
        #SwTo<-function(l){return(1/length(LAll))}
        #SwTo<-function(l){return(0)}
       #SwTo<-function(l){return((1+Diffac*(median(LAll)-l))/length(LAll))}#switch preference
        
        #Die<-function(l){2.5}
        #Die<-function(l){2.5-Diffac*(median(LAll)-l)}#differential deatch
        
        Threshold<-Threshold
        KillingRate<-KillingRate
        
        #deterministic version with delay
        # ImKill<-function(Nl){#immune dection and kill,threshold version
        #         #VecPro=rexp(length(l),rate = Nl)<1e-3&(!l%in%ImmundetectionRecord)
        #          VecPro=(Nl>Threshold)&(!l%in%ImmundetectionRecord)
        #          if(any(Nl>Threshold&(!(l%in%ImmundetectionRecord)))){
        #                  DetectionTime<<-c(DetectionTime,rep(t,sum(VecPro)))#record detection time
        #                 ImmundetectionRecord<<-c(ImmundetectionRecord,l[VecPro])
        #                }
        #      if(any(l%in%ImmundetectionRecord)){ OutputVec<-rep(0,length(l))
        #                  ImmuneKillStarted<<-ImmundetectionRecord[t-DetectionTime>=5]
        #                 OutputVec[l%in%ImmuneKillStarted]<-KillingRate
        #                  #OutputVec[l%in%ImmundetectionRecord]<- -Diffac*(l[l%in%ImmundetectionRecord]-median(LAll))+KillingRate #differential immune kill
        #                    return(OutputVec)} else {return(rep(0,length(l)))} }
         
        # 
        
        
         #Problistic version 
        ImKill<-function(Nl){#immune dection and kill,threshold version
                VecPro=rnorm(n = length(l),mean = Nl)>Threshold&(!l%in%ImmundetectionRecord)
                if(any((VecPro)&(!(l%in%ImmundetectionRecord)))){
                        DetectionTime<<-c(DetectionTime,rep(t,sum(VecPro)))#record detection time
                        ImmundetectionRecord<<-c(ImmundetectionRecord,l[VecPro])
                }
                if(any(l%in%ImmundetectionRecord)){ OutputVec<-rep(0,length(l))
                ImmuneKillStarted<<-ImmundetectionRecord[t-DetectionTime>=5]
                OutputVec[l%in%ImmuneKillStarted]<-KillingRate
                #OutputVec[l%in%ImmundetectionRecord]<- -Diffac*(l[l%in%ImmundetectionRecord]-median(LAll))+KillingRate #differential immune kill
                return(OutputVec)} else {return(rep(0,length(l)))} }
        

         

       #  ImmuneRecordVec<-rep(0,length(LAll))
       # ImmuneRecordVec[LAll%in%ImmundetectionRecord]<-1
       #  #ImmuneRecordMatrix<<-rbind(ImmuneRecordMatrix,ImmuneRecordVec)
       # 
       #   ImKill<-function(Nl){#immune dection and kill,probabolity version
       #           if(any(!l%in%ImmundetectionRecord)){
       #                   ImmundetectionRecord<<-c(ImmundetectionRecord,l[rexp(length(l),rate = Nl)<1e-8])
       #                   ImmundetectionRecord<<-unique(ImmundetectionRecord)}
       #           if(any(l%in%ImmundetectionRecord)){ OutputVec<-rep(0,length(l));OutputVec[l%in%ImmundetectionRecord]<-KillingRate
       #           return(OutputVec)} else {return(rep(0,length(l)))} }
       #  ImKill<-function(Nl){return(0)}
       # 
        
        #dNl<- (1-(sum(Nl)/MaximumLoad))*R(l)*Nl-Swi(l)*Nl+SwTo(l)*sum(Nl*Swi(l))-ImKill(Nl)*Nl
        dNl<- (1-sum(Nl)/MaximumLoad)*R(l)*Nl-Swi(l)*Nl+SwTo(l)*sum(Nl*Swi(l))-ImKill(Nl)*Nl
          #dNl<- R(l)*Nl-Swi(l)*Nl+SwTo(l)*sum(Nl*Swi(l))-ImKill(Nl)*Nl
        
        return(list(c(dNl)))        
}


YINI<-rep(0,length(LAll))
YINI[218]<-1
#MuGa<-ode(times = seq(from=0,to=30,by=0.01),y = 0,func = PopuDynamic,parms = 1500,method='rk4')  
Output<-ode(times = seq(from=0,to=30,by=0.01),y= YINI,func =PopuDynamic,parms = LAll)

#plot(ImmuneKillStarted~DetectionTime)
#Output[,2:ncol(Output)][Output[,2:ncol(Output)]<1]<-0


#pdf('DG8.pdf')
Parasitemia<-apply(Output[,-1],MARGIN = 1,sum)

plot(Parasitemia~Output[,1],log='y')
#functin to calculate weighted mean 

WTmean<-apply(Output[,-1],MARGIN = 1,FUN = function(x){if(sum(x)>1) return(sum(LAll[x>1]*(x[x>1]/sum(x[x>1])))) else {return(0)}})

LengthDom<-apply(Output[,-1],MARGIN = 1,FUN = function(x){LAll[which.max(x)]})                                                                                                      

#Get the 5 responsive variable needed 

WtMeanNTime<-data.frame(WTmean,Time=Output[,1],Parasitemia,LengthDom)
head(WtMeanNTime)
min(WtMeanNTime$WTmean)

WtMeanNTime<-WtMeanNTime[apply(Output[,-1],MARGIN = 1,FUN = function(x){any(x>=1)}),]
range(WtMeanNTime$Time)
range(WtMeanNTime$WTmean)

T1<-WtMeanNTime$Time[which.min(WtMeanNTime$WTmean)]#Time when lowest WTmean 
T2<-WtMeanNTime$Time[which.max(WtMeanNTime$Time)]#end time 
EndLength<-WtMeanNTime$WTmean[which.max(WtMeanNTime$Time)]#end point WTmean length

LowestLength<-WtMeanNTime$WTmean[which.min(WtMeanNTime$WTmean)]#lowest Wtmean VSG length

CorDecrease<-cor(WtMeanNTime[WtMeanNTime$Time<=T1,1],WtMeanNTime[WtMeanNTime$Time<=T1,2])#whether first stage decrease 
CorIncrease<-cor(WtMeanNTime[WtMeanNTime$Time>T1,1],WtMeanNTime[WtMeanNTime$Time>T1,2])


#Immune infor
ImmuneDist<-matrix(data = 0,ncol = length(LAll),nrow = nrow(Output))#matrix to store which VSG detected each time point 
colnames(ImmuneDist)<-LAll

ImmuneVec<-data.frame(VSGlength=ImmundetectionRecord,DetectionTime)
head(ImmuneVec)

for (i in 1:nrow(ImmuneDist)){
        ImmuneDist[i,which(LAll%in%ImmuneVec$VSGlength[ImmuneVec$DetectionTime<=Output[i,1]-5])]<-1
}
#apply(ImmuneDist, MARGIN = 1, sum)


#VecToSave<-unlist(c(ParaSpace[SimuNumber,],T1=T1,T2=T2,EndLength=EndLength,
#                    LowestLength=LowestLength,CorDecrease=CorDecrease,CorIncrease=CorIncrease))

#look at the results 
#look at the results 


par(mfrow=c(2,2))
for (i in c(1,10,15,20)){
        Target<-which.min(abs(Output[,1]-i))        
        plot(Output[Target,-1]/Parasitemia[Target]~LAll,type='h',xlab='VSGlength',ylab='% in population',main=paste0('Day',i))
}
par(mfrow=c(1,1))

# 
# par(mfrow=c(2,2))
# for (i in c(1,10,15,20)){
#         Target<-which.min(abs(Output[,1]-i))
#         VecPlot<-data.frame(VSGL=LAll,Freq=Output[Target,-1])
#         MyHist<-list(breaks=seq(from=1000,to=2000,by=100),counts=VecPlot$Freq, density=VecPlot$Freq/diff(c(1000,VecPlot$VSGL),col='red'),
#                      xname="VSGLdistr")
#         class(MyHist)<-'histogram'
#         
#         plot(MyHist,col='red',main=paste0('Day',i))
#         
# }
# par(mfrow=c(1,1))
# par(mfrow=c(2,2))
# for (i in c(1,10,15,20)){
#         Target<-which.min(abs(Output[,1]-i))
#         VecPlot<-data.frame(VSGL=LAll,Freq=Output[Target,-1])
#         VecPlot<-rep(LAll,round(Output[Target,-1]/min(Output[Target,-1]))) 
#         hist(VecPlot,col='red',main=paste0('Day',i),probability = T)
#         
# }
# par(mfrow=c(1,1))

par(mfrow=c(2,3),mar=c(4,4,4,4))
for (i in c(1,5,10,15,20,30)){
        Target<-which.min(abs(Output[,1]-i))
        VecPlot<-data.frame(VSGL=LAll,Freq=Output[Target,-1])
        VecParasitemia<-sum(Output[Target,-1])
        VecPlot$Freq[VecPlot$Freq<1]<-0
        VecPlot<-rep(LAll,round(Output[Target,-1]/min(Output[Target,-1])))
        DenVec<-density(VecPlot,adjust = 5)#density plot 
        
        if(any(Output[Target,-1]>=1)){
                plot(DenVec,col='red',main=paste0('Day',i),lwd=2,xlab='',ylab='',xlim=c(min(LAll),max(LAll)),
                     cex.lab=2,cex.main=2,cex.axis=2)
                
                if ((range(DenVec$y)[2]-range(DenVec$y)[1])<=0.2){polygon(DenVec,col='red',border = 'red')}}#density plot for a single line looks odd 
        
        
        if(all(Output[Target,-1]<1)){
                plot(1,col='red',main=paste0('Day',i,' (died out)'),lwd=2,xlab='',ylab='',xlim=c(min(LAll),max(LAll)),
                     ylim = range(DenVec$y),
                     cex.lab=2,cex.main=2,cex.axis=2)}
        
        abline(v = median(LAll),col='blue',lwd=2,lty=2)
        
        #hist(VecPlot,col='red',main=paste0('Day',i),probability = T)
        
        if (sum(ImmuneDist[Target,]==1)>1){
                ImmuneDens<-density(LAll[ImmuneDist[Target,]==1])
                ImmuneDens$y<-ImmuneDens$y*( max(DenVec$y)/max(ImmuneDens$y))#need an additional matrix 
                lines(ImmuneDens)
                polygon(ImmuneDens,col=adjustcolor( "blue", alpha.f = 0.5),border = adjustcolor( "blue", alpha.f = 0.5))}
}
par(mfrow=c(1,1))


apply(Output[,-1],MARGIN = 1,FUN = function(x){any(x>=1)})

TimeBeforeDieout<-max(Output[apply(Output[,-1],MARGIN = 1,FUN = function(x){any(x>=1)}),1])

#plot(WtMeanNTime$WTmean~WtMeanNTime$Time,xlim=c(0,30),type='s',xlab='Days after infection',ylab='WtMean of VSG length')
plot(WtMeanNTime$WTmean~WtMeanNTime$Time,xlim=c(0,30),type='s',xlab='',ylab='',main='DA',cex.main=3,cex.axis=2,ylim=c(min(LAll),max(LAll)))
polygon(x = c(TimeBeforeDieout,TimeBeforeDieout,32,32),y = c(800,2000,2000,800),col='grey')

#plot(WtMeanNTime$LengthDom~WtMeanNTime$Time,xlim=c(0,30),type='s',xlab='Days after infection',ylab='VSG length of dominating clone')
plot(WtMeanNTime$LengthDom~WtMeanNTime$Time,xlim=c(0,30),type='s',xlab='',ylab='',main='DA',cex.main=3,cex.axis=2,ylim=c(min(LAll),max(LAll)))

polygon(x = c(TimeBeforeDieout,TimeBeforeDieout,32,32),y = c(800,2000,2000,800),col='grey')



#Growth rate vs. VSG length 
plot((2-DG*(LAll-median(LAll)))~LAll,ylab='Replications/day',xlab='VSG lengths(bp)',type='l',lwd=2)

#immune detection time vs length 
plot(ImmuneVec$DetectionTime~ImmuneVec$VSGlength,xlab='VSG lengths(bp)',ylab='Time detected by adaptive immune system(days)')

dev.off()
