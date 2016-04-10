# R-script for running simulations used for estimating defaunation effects on aboveground carbon storage,
#stand wood density and stand volume in Osuri et al. 'Contrasting effects of defaunation on aboveground carbon storage
#across the global tropics'. Nature Communications (2016)
#Authors: Anand M Osuri & Varun Varma

library(scales)
library(sampling)
library(matrixStats)
library(Hmisc)

dat<-read.csv("C:/Users/Uttara & Anand/Dropbox/analysis/test.csv") #input plot dataframe 'dat'
dat<-dat[complete.cases(dat),]
dat$basal<-pi*(dat$diameter^2)/4
#'dat' needs to be composed of the following columns:
#'site'= spatial unit for simulation
#'accpt.nam'= species name
#'diameter'= tree basal diameter (DBH)
#'basal'= tree basal area
#'volume'= tree estimated volume
#'carbon'= tree estimated carbon stock
#'seed'= species seed length
#'dmode'= species seed dispersal mode
#'wden'= species wood density
# In this study, carbon and volume are calculated using the equation for moist forests in 
# Chave et al. (2005) Oecologia 145: 87-99 (see Methods and below).
# dat$carbon<-0.5*dat$wden*exp(-1.499+2.148*log(dat$diameter)+0.207*log(dat$diameter)^2-0.0281*log(dat$diameter)^3)
# dat$volume<-exp(-1.499+2.148*log(dat$diameter)+0.207*log(dat$diameter)^2-0.0281*log(dat$diameter)^3)



site<-"A"
X1<-NA
X2<-NA
X3<-NA
X4<-NA
X5<-NA
baseframe<-data.frame(site,X1,X2,X3,X4,X5)

# X1, X2, X3, X4 and X5 correspond to results at 0, 25, 50, 75 and 100% extirpation
# of large-seeded animal-dispersed tree species
S.SEED.N<-baseframe # Community-weighted seed size (defaunation scenario)
S.BASAL<-baseframe # Total basal area (defaunation scenario)
S.VOLUME<-baseframe # Stand volume (defaunation scenario)
S.WDEN.BA<-baseframe # Community-weighted wood density (defaunation scenario)
S.LARGE70<-baseframe # Relative abundance of large trees (defaunation scenario)
S.CARBON<-baseframe # Carbon stock (defaunation scenario)

RI.SEED.N<-baseframe # Community-weighted seed size (control 1 scenario)
RI.BASAL<-baseframe # Total basal area (control 1 scenario)
RI.VOLUME<-baseframe  # Stand volume (control 1 scenario)
RI.WDEN.BA<-baseframe # Community-weighted wood density (control 1 scenario)
RI.LARGE70<-baseframe # Relative abundance of large trees (control 1 scenario)
RI.CARBON<-baseframe # Carbon stock (control 1 scenario)

RS.SEED.N<-baseframe # Community-weighted seed size (control 2 scenario)
RS.BASAL<-baseframe # Total basal area (control 2 scenario)
RS.VOLUME<-baseframe  # Stand volume (control 2 scenario)
RS.WDEN.BA<-baseframe # Community-weighted wood density (control 2 scenario)
RS.LARGE70<-baseframe # Relative abundance of large trees (control 2 scenario)
RS.CARBON<-baseframe # Carbon stock (control 2 scenario)

# Each data frame to which 'baseframe' is assigned will grow to hold 1000 simulation outputs for each site

###

sites<-sort(unique(dat$site))
Seedquan<-matrix(0,nrow=length(sites),ncol=5)
rownames(Seedquan)<-sites
colnames(Seedquan)<-c("Q0","Q25","Q50","Q75","Q100")
for(k in sites){
  
  data<-subset(dat,dat$site==as.character(k))
  
  species<-sort(unique(data$accpt.nam))
  seed<-data$seed[match(species,data$accpt.nam)]
  dmode<-data$dmode[match(species,data$accpt.nam)]
  wden<-data$wden[match(species,data$accpt.nam)]
  basal<-tapply(data$basal,data$accpt.nam,sum)
  basal<-basal[complete.cases(basal)]
  
  Seedquan[k,]<-quantile(seed[dmode=="A"|dmode=="MULTI"],na.rm=T)
  upquan<-Seedquan[k,4]
  #upquan corresponds to the 75th percentile of seed size distribution of site 'k'
  large<-species[seed>=upquan]
  
  INDEX<-match(data$accpt.nam,large)
  data$large<-ifelse(is.na(INDEX)==T,0,1)
  Nlarge<-sum(data$large)
  
  ITER<-2 # No. of iterations to be used
  use.iter<-1:ITER
  
  use.perc<-c(0,0.25,0.5,0.75,1)
  use.seq<-round(use.perc*Nlarge,0)
  # use.seq contains the numbers of individuals to 0-100% of 
  #large-seeded trees in site 'k'
  
  S.Seed.n<-matrix(0,nrow=length(use.iter),ncol=length(use.seq))
  S.Basal<-matrix(0,nrow=length(use.iter),ncol=length(use.seq))
  S.Volume<-matrix(0,nrow=length(use.iter),ncol=length(use.seq))
  S.Wden.ba<-matrix(0,nrow=length(use.iter),ncol=length(use.seq))
  S.Large70<-matrix(0,nrow=length(use.iter),ncol=length(use.seq))
  S.Carbon<-matrix(0,nrow=length(use.iter),ncol=length(use.seq))
  
  # Initial (zero removal) values set based on original community of site 'k'
  S.Seed.n[,1]<-mean(data$seed[data$dmode=="A"|data$dmode=="MULTI"])
  S.Basal[,1]<-sum(data$basal)
  S.Volume[,1]<-sum(data$volume)
  S.Wden.ba[,1]<-weighted.mean(wden,basal)
  S.Large70[,1]<-nrow(data[data$diameter>=70,])/nrow(data)
  S.Carbon[,1]<-sum(data$carbon)/1000
  
  RI.Seed.n<-S.Seed.n
  RI.Basal<-S.Basal
  RI.Volume<-S.Volume
  RI.Wden.ba<-S.Wden.ba
  RI.Large70<-S.Large70
  RI.Carbon<-S.Carbon
  
  RS.Seed.n<-S.Seed.n
  RS.Basal<-S.Basal
  RS.Volume<-S.Volume
  RS.Wden.ba<-S.Wden.ba
  RS.Large70<-S.Large70
  RS.Carbon<-S.Carbon
  

  test.basal<-sum(data$basal,na.rm=T) # original basal area
  
  for(i in use.iter){
    
    
    for(j in 2:length(use.seq)){
      ##Defaunation scenario
      INDEX<-1:nrow(data)
      #drawing and removing a random sample of large-seeded animal-dispersed trees
      REM<-sample(INDEX[data$large==1],use.seq[j],replace=F)
      MATCH<-match(INDEX,REM)
      seeddat<-subset(data,is.na(MATCH)==T)
      BA.loss<-test.basal-sum(seeddat$basal)
      test.buffer<-0.01*BA.loss
      # in the loop below, individuals from the remaining pool are sampled at random
      # to make up lost basal area. This step is repeated until basal area is recovered
      # to within 1% of the original
      count<-0
      repeat{
        count<-count+1
        newframe<-seeddat[sample(nrow(seeddat),25*nrow(data),replace=T),]
        cumba<-cumsum(newframe$basal)
        test<-abs(BA.loss-cumba)
        min.diff<-min(test)
        if(min.diff<=test.buffer){break}
        if(count==100){break}
      }
      
      newframe<-newframe[1:match(min.diff,test),] # cutting the recovery data frame at the target basal area
      seeddat<-rbind(seeddat,newframe)
      
      S.Seed.n[i,j]<-ifelse(count<100,mean(seeddat$seed[seeddat$seed>0]),NA)
      S.Basal[i,j]<-ifelse(count<100,sum(seeddat$basal),NA)
      S.Volume[i,j]<-ifelse(count<100,sum(seeddat$volume),NA)
      SP<-sort(unique(seeddat$accpt.nam))
      WD<-seeddat$wden[match(SP,seeddat$accpt.nam)]
      BA<-tapply(seeddat$basal,seeddat$accpt.nam,sum)
      BA<-BA[complete.cases(BA)]
      S.Wden.ba[i,j]<-ifelse(count<100,weighted.mean(WD,BA,na.rm=T),NA)
      S.Large70[i,j]<-ifelse(count<100,nrow(seeddat[seeddat$diameter>=70,])/nrow(seeddat),NA)
      S.Carbon[i,j]<-ifelse(count<100,sum(seeddat$carbon,na.rm=T)/1000,NA)
      print(paste("Iter",i,"Seq",j,"Count",count,k))
      sploss<-length(unique(data$accpt.nam))-length(unique(seeddat$accpt.nam))
      # sploss is the number of species extirpated in iteration i,j of the defaunation
      # scenario. This number is used subsequently in control 2 scenario
      
      #Control scenario 1: no. of individuals removed = no. of individuals removed in defaunation scenario
      INDEX<-1:nrow(data)
      # drawing and removing a random sample of large-seeded, small-seeded and
      # abiotically-dispersed trees, such that the number of individuals removed is same
      # as in the defaunation scenario
      REM<-sample(INDEX,use.seq[j],replace=F)
      MATCH<-match(INDEX,REM)
      randdat<-subset(data,is.na(MATCH)==T)
      BA.loss<-test.basal-sum(randdat$basal)
      test.buffer<-0.01*BA.loss    
      count<-0
      repeat{
        count<-count+1
        newframe<-randdat[sample(nrow(randdat),25*nrow(data),replace=T),]
        cumba<-cumsum(newframe$basal)
        test<-abs(BA.loss-cumba)
        min.diff<-min(test)
        if(min.diff<=test.buffer){break}
        if(count==100){break}
      }
      newframe<-newframe[1:match(min.diff,test),] # cutting the recovery data frame at the target basal area
      randdat<-rbind(randdat,newframe)
      randdat<-randdat[sample(nrow(randdat),nrow(randdat),replace=F),]
      RI.Seed.n[i,j]<-ifelse(count<100,mean(randdat$seed[randdat$seed>0]),NA)
      RI.Basal[i,j]<-ifelse(count<100,sum(randdat$basal),NA)
      RI.Volume[i,j]<-ifelse(count<100,sum(randdat$volume),NA)
      SP<-sort(unique(randdat$accpt.nam))
      WD<-randdat$wden[match(SP,randdat$accpt.nam)]
      BA<-tapply(randdat$basal,randdat$accpt.nam,sum)
      BA<-BA[complete.cases(BA)]
      RI.Wden.ba[i,j]<-ifelse(count<100,weighted.mean(WD,BA,na.rm=T),NA)
      RI.Large70[i,j]<-ifelse(count<100,nrow(randdat[randdat$diameter>=70,])/nrow(randdat),NA)
      RI.Carbon[i,j]<-ifelse(count<100,sum(randdat$carbon,na.rm=T)/1000,NA)
      
      
      ##Control scenario 2: no. of species removed = no. of species lost in defaunation scenario
      # drawing and removing a random sample of large-seeded, small-seeded and
      # abiotically-dispersed trees, such that the number of species lost is the same as
      # in the defaunation scenario
      
      REM<-sample(species,sploss,replace=F)
      INDEX<-match(data$accpt.nam,REM)
      INDEX<-ifelse(is.na(INDEX)==T,1,0)
      randdat<-data[INDEX==1,]
      BA.loss<-test.basal-sum(randdat$basal)
      test.buffer<-0.01*BA.loss    
      count<-0
      repeat{
        count<-count+1
        newframe<-randdat[sample(nrow(randdat),25*nrow(data),replace=T),]
        cumba<-cumsum(newframe$basal)
        test<-abs(BA.loss-cumba)
        min.diff<-min(test)
        if(min.diff<=test.buffer){break}
        if(count==100){break}
      }
      newframe<-newframe[1:match(min.diff,test),] # cutting the recovery data frame at the target basal area
      randdat<-rbind(randdat,newframe)
      randdat<-randdat[sample(nrow(randdat),nrow(randdat),replace=F),]
      RS.Seed.n[i,j]<-ifelse(count<100,mean(randdat$seed[randdat$seed>0]),NA)
      RS.Basal[i,j]<-ifelse(count<100,sum(randdat$basal),NA)
      RS.Volume[i,j]<-ifelse(count<100,sum(randdat$volume),NA)
      SP<-sort(unique(randdat$accpt.nam))
      WD<-randdat$wden[match(SP,randdat$accpt.nam)]
      BA<-tapply(randdat$basal,randdat$accpt.nam,sum)
      BA<-BA[complete.cases(BA)]
      RS.Wden.ba[i,j]<-ifelse(count<100,weighted.mean(WD,BA,na.rm=T),NA)
      RS.Large70[i,j]<-ifelse(count<100,nrow(randdat[randdat$diameter>=70,])/nrow(randdat),NA)
      RS.Carbon[i,j]<-ifelse(count<100,sum(randdat$carbon,na.rm=T)/1000,NA)
      
      
    }
  }
  s.Seed.n<-data.frame("site"=rep(k,length(use.iter)),S.Seed.n)
  s.Basal<-data.frame("site"=rep(k,length(use.iter)),S.Basal)
  s.Volume<-data.frame("site"=rep(k,length(use.iter)),S.Volume)
  s.Wden.ba<-data.frame("site"=rep(k,length(use.iter)),S.Wden.ba)
  s.Large70<-data.frame("site"=rep(k,length(use.iter)),S.Large70)
  s.Carbon<-data.frame("site"=rep(k,length(use.iter)),S.Carbon)
  
  
  rs.Seed.n<-data.frame("site"=rep(k,length(use.iter)),RS.Seed.n)
  rs.Basal<-data.frame("site"=rep(k,length(use.iter)),RS.Basal)
  rs.Volume<-data.frame("site"=rep(k,length(use.iter)),RS.Volume)
  rs.Wden.ba<-data.frame("site"=rep(k,length(use.iter)),RS.Wden.ba)
  rs.Large70<-data.frame("site"=rep(k,length(use.iter)),RS.Large70)
  rs.Carbon<-data.frame("site"=rep(k,length(use.iter)),RS.Carbon)
  
  ri.Seed.n<-data.frame("site"=rep(k,length(use.iter)),RI.Seed.n)
  ri.Basal<-data.frame("site"=rep(k,length(use.iter)),RI.Basal)
  ri.Volume<-data.frame("site"=rep(k,length(use.iter)),RI.Volume)
  ri.Wden.ba<-data.frame("site"=rep(k,length(use.iter)),RI.Wden.ba)
  ri.Large70<-data.frame("site"=rep(k,length(use.iter)),RI.Large70)
  ri.Carbon<-data.frame("site"=rep(k,length(use.iter)),RI.Carbon)  
  
  #output dataframes follow
  #each output dataframe comprises 5 columns, corresponding to
  #'site' = site name
  #'X1'= Value in original community
  #'X2'= Value at 25% extirpation of large-seeded animal-dispersed species (or, corresponding value for control simulations)
  #'X3'= Value at 50% extirpation
  #'X4'= Value at 75% extirpation
  #'X5'= Value at 100% extirpation
  
# Note that each dataframe contains a dummy first row, which can be ignored/deleted.
  # Defaunation scenario
  S.SEED.N<-rbind(S.SEED.N,s.Seed.n) # Community-weighted seed size
  S.BASAL<-rbind(S.BASAL,s.Basal) # Total basal area
  S.VOLUME<-rbind(S.VOLUME,s.Volume) # Total volume
  S.WDEN.BA<-rbind(S.WDEN.BA,s.Wden.ba) # Basal area-weighted wood density
  S.LARGE70<-rbind(S.LARGE70,s.Large70) # Relative abundance of large trees
  S.CARBON<-rbind(S.CARBON,s.Carbon) # Total carbon stocks
  
  # Control scenario 1
  RI.SEED.N<-rbind(RI.SEED.N,ri.Seed.n) # Community-weighted seed size
  RI.BASAL<-rbind(RI.BASAL,ri.Basal) # Total basal area
  RI.VOLUME<-rbind(RI.VOLUME,ri.Volume) # Total volume
  RI.WDEN.BA<-rbind(RI.WDEN.BA,ri.Wden.ba) # Basal area-weighted wood density
  RI.LARGE70<-rbind(RI.LARGE70,ri.Large70) # Relative abundance of large trees
  RI.CARBON<-rbind(RI.CARBON,ri.Carbon) # Total carbon stocks
  
  # Control scenario 2
  RS.SEED.N<-rbind(RS.SEED.N,rs.Seed.n) # Community-weighted seed size
  RS.BASAL<-rbind(RS.BASAL,rs.Basal) # Total basal area
  RS.VOLUME<-rbind(RS.VOLUME,rs.Volume) # Total volume
  RS.WDEN.BA<-rbind(RS.WDEN.BA,rs.Wden.ba) # Basal area-weighted wood density
  RS.LARGE70<-rbind(RS.LARGE70,rs.Large70) # Relative abundance of large trees
  RS.CARBON<-rbind(RS.CARBON,rs.Carbon) # Total carbon stocks

}

