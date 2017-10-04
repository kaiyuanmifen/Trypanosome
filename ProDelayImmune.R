#improved grwoth model , unified code with other model 


par(mfcol=c(1,1))

#set.seed(131)
#setwd("/Users/dliu/Box Sync/research/Third year 3/multiple gene copies/Record/20160323/")
load(file = 'MugnierVSGLib.Rdata')
#hist(VSGData,breaks = 20)
#VSGData<-sort(sample(BedData$V8-BedData$V7,size = 100,replace = F))


# VSGData<-sort(sample(x = VSGData, size =300 , replace = TRUE))#try this example 
VSGData<-LAll

range(VSGData)
#length(VSGData)

#set.seed(113)
# ActiveVSGVec<-rep(F,length(VSGData))#  A Vec to store which VSGs are actively expressed  
Pro.Vec<-rep(x = 0,length(VSGData))
Pro.Vec[which.min(abs(VSGData-1600))]<-1#start simulation from median

names(Pro.Vec)<-as.character(VSGData)
names(Pro.Vec)<-1:length(VSGData)

# ActiveVSGVec[which.min(abs(VSGData-1700))]<-T
#names(ActiveVSGVec)<-as.character(VSGData)
#names(ActiveVSGVec)<-1:length(VSGData)


##all VSG are asigned identical switch rates 
Qv<-10^-6#stall rate per base. increase it to speed up simulation, I realise having slightly lower switch rate will give more similar results 
Nv<-1e6#average distance from replicaiton origin
#Swrate<-((1-Qv)^Nv)*(1-(1-Qv)^VSGData)#siwtch rate per replication
Swrate<-rep(1e-4,length(Pro.Vec))

##generate the transformation matrix 
SwitchFactor<-1#this decides the bia od swtich , no switch preference 

Trans.Matrix<-matrix(0,nrow = length(Pro.Vec),ncol =length(Pro.Vec) )
rownames(Trans.Matrix)<-as.character(VSGData)
colnames(Trans.Matrix)<-as.character(VSGData)

for (i in 1:length(Pro.Vec)){
        Trans.Matrix[i,i]<-1-Swrate[i]
        Trans.Matrix[i,]<-Trans.Matrix[i,]+(Swrate[i]/length(VSGData))*(SwitchFactor/(sum(SwitchFactor)))
      
}


###check it , OI think I got some problem here Ans: no , this is fine 
# Trans.Matrix
# t(Trans.Matrix)%*%Pro.Vec
#now run the simulation

output<-c(0,sum(Pro.Vec),Pro.Vec)
Count<-0
TotalRounds<-4*30
ImmuneThreshold<-10^5
DetectionFactor<-10^(8)
KillDelay<-4*5
MaxLoad<-10^9 #max number of parasite the system can obtain 
KillRate<-4
ActiFactor<-10^(-5)

DectectedVSG<-data.frame(VSGID=names(Pro.Vec),Dectection=rep(F,length(Pro.Vec)),
                         DectTime=rep(-1,length(Pro.Vec)),ActiveKill=rep(F,length(Pro.Vec)))#a data frame to record detected state and killing stateof the VSG


ProgVec<-txtProgressBar(min = Count<-0,max =TotalRounds )
repeat{
        Count<-Count+1  
        if (Count>TotalRounds){break}
        
        # Pro.Vec[!ActiveVSGVec]<-0
        
        #ActiveVSGVec[sample(x = c(1:length(ActiveVSGVec)),
        #                    replace = T,size = floor(ActiFactor*sum(Swrate*Pro.Vec)))]<-TRUE#randomly sample VSG to activate
        
        
        
        Pro.Vec[DectectedVSG$ActiveKill==T]<-Pro.Vec[DectectedVSG$ActiveKill==T]/KillRate#immune kill
        Pro.Vec <- (Pro.Vec%*%Trans.Matrix)#switch
        
        ##this where differential growth rate take places 
        #GrowthFactor<-1#no differential growth 
        GrowthFactor<-(median(VSGData)/VSGData)^0.2
        GrowthVec<-GrowthFactor*Pro.Vec*(1 - sum(Pro.Vec)/MaxLoad)#growth vector
        Pro.Vec<-Pro.Vec+GrowthVec #growth 
        
        #Pro.Vec[!ActiveVSGVec]<-0
        
    
        #   if ( sum(DectectedVSG$ActiveKill==T)>0){#This is the point, kill all the detected cells 
        #     Pro.Vec[DectectedVSG$ActiveKill==T]<-0
        #   }  
        
        Parasitemia<-sum(Pro.Vec)
        UpdateVec<-c(TimeStep=Count,Parasitemia,Pro.Vec)
        
        output<-rbind(output,UpdateVec)
        
        
          ##random detection of VSG by immune system
          VSGID <- runif(length(Pro.Vec)) < 1-exp(-Pro.Vec/(DetectionFactor))
          DectectedVSG$DectTime[DectectedVSG$DectTime == -1 & VSGID] <- Count
          DectectedVSG$Dectection[VSGID]<-T
           
#           #detection be threshold 
#           VSGID<-DectectedVSG$Dectection==F&Pro.Vec>ImmuneThreshold
#           DectectedVSG$DectTime[DectectedVSG$DectTime == -1 & VSGID] <- Count
#           DectectedVSG$Dectection[VSGID]<-T
#         
#         #just kill the dominating clone 
#         
#         if (Parasitemia>ImmuneThreshold){
#                 #VSGID<- which(order(Pro.Vec,decreasing = T)%in%c(1:1))
#                 VSGID<- which(Pro.Vec==max(Pro.Vec[DectectedVSG$Dectection==F]))
#                 if(any(DectectedVSG$Dectection[VSGID]==F)) {DectectedVSG$Dectection[VSGID]<-T}
#                 if(any(DectectedVSG$DectTime[VSGID] ==-1)){DectectedVSG$DectTime[VSGID][DectectedVSG$DectTime[VSGID] ==-1] <- Count}
#         }
#         
        if(sum(DectectedVSG$Dectection==T&DectectedVSG$ActiveKill==F&DectectedVSG$DectTime<=Count-KillDelay)>0){
                DectectedVSG$ActiveKill[DectectedVSG$Dectection==T&DectectedVSG$ActiveKill==F&DectectedVSG$DectTime<=Count-KillDelay]<-T
                #activate active kill for those dectedfor more than killding delay 
        }  
        
        
        
        
        
        setTxtProgressBar(pb = ProgVec,value = Count)
        
}
colnames(output)<-c('TimeStep','Parasitemia',as.character(VSGData))

#use a way , consistent with other figures to plot
Output=output[,-2]
close(ProgVec)



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
        ImmuneDist[i,which(LAll%in%ImmuneVec$VSGlength[ImmuneVec$DetectionTime<=Output[i,1]])]<-1
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
plot(WtMeanNTime$WTmean~WtMeanNTime$Time,xlim=c(0,30),type='s',xlab='',ylab='',main='DG',cex.main=3,cex.axis=2,ylim=c(min(LAll),max(LAll)))
polygon(x = c(TimeBeforeDieout,TimeBeforeDieout,32,32),y = c(800,2000,2000,800),col='grey')

#plot(WtMeanNTime$LengthDom~WtMeanNTime$Time,xlim=c(0,30),type='s',xlab='Days after infection',ylab='VSG length of dominating clone')
plot(WtMeanNTime$LengthDom~WtMeanNTime$Time,xlim=c(0,30),type='s',xlab='',ylab='',main='DG',cex.main=3,cex.axis=2,ylim=c(min(LAll),max(LAll)))

polygon(x = c(TimeBeforeDieout,TimeBeforeDieout,32,32),y = c(800,2000,2000,800),col='grey')



#Growth rate vs. VSG length 
plot((2-DG*(LAll-median(LAll)))~LAll,ylab='Replications/day',xlab='VSG lengths(bp)',type='l',lwd=2)

#immune detection time vs length 
plot(ImmuneVec$DetectionTime~ImmuneVec$VSGlength,xlab='VSG lengths(bp)',ylab='Time detected by adaptive immune system(days)')
