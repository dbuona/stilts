rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())
graphics.off()


setwd("~/Documents/git/stilts/")

d<-read.csv("data/phenology survey 2025.xlsx - full_datasheet.csv")

demos<-d[,1:16]

cols_to_keep <- grep("114", names(d), value = TRUE)

# Subset the data frame to keep only those columns
timechange <- d[ , cols_to_keep]
colnames(timechange)<-c("Retnoutria","Artemesia","Berberis","Microstegium","Pueria","Celastrus","Pinus","Vinetoxicum","Rhamnus",
                        "Elaeagnus","impatiens","Lonicera","Miscanthus","Rosa","Acer","Phragmotes","Cytisus","Centaurea","Ailanthis","Pastinaca","Ligustrum")


#timechange<-timechange[,4]

cols_to_keep <- grep("117", names(d), value = TRUE)
freqchange <- d[ , cols_to_keep]

colnames(freqchange)<-c("Retnoutria","Artemesia","Berberis","Microstegium","Pueria","Celastrus","Pinus","Vinetoxicum","Rhamnus",
  "Elaeagnus","impatiens","Lonicera","Miscanthus","Rosa","Acer","Phragmotes","Cytisus","Centaurea","Ailanthis","Pastinaca","Ligustrum")


#freqchange<-freqchange[,4]

cols_to_keep <- grep("118", names(d), value = TRUE)
dur<-d[ , cols_to_keep]
colnames(dur)<-c("Retnoutria","Artemesia","Berberis","Microstegium","Pueria","Celastrus","Pinus","Vinetoxicum","Rhamnus",
                   "Elaeagnus","impatiens","Lonicera","Miscanthus","Rosa","Acer","Phragmotes","Cytisus","Centaurea","Ailanthis","Pastinaca","Ligustrum")

d$Q35
d$state
timechange<-as.data.frame(cbind(timechange,d$Q35,d$state))
freqchange<-as.data.frame(cbind(freqchange))
dur<-as.data.frame(cbind(dur))


timechange<-tidyr::gather(timechange,"species","timechange",1:21)
freqchange<-tidyr::gather(freqchange,"species","freqchange",1:21)
dur<-tidyr::gather(dur,"species","dur",1:21)

freqchange<-dplyr::select(freqchange,-species)
dur<-dplyr::select(dur,-species)

dat<-cbind(freqchange,timechange,dur)
dat[dat == ""] <- NA


dat$phenshift<-NA
dat$phenshift<-ifelse(dat$timechange %in% c("At the same time due to no change in phenology","Other","Earlier or later for other reasons"),"no","yes")
dat$phenshift<-ifelse(dat$timechange=="I don't know","unknown",dat$phenshift)

dat$phenfreq<-NA
dat$phenfreq<-ifelse(dat$freqchange %in% c("At the same frequency due to no change in phenology","More or less frequently for other reasons","Other" ),"no","yes")
dat$phenfreq<-ifelse(dat$freqchange=="I don’t know" ,"unknown",dat$phenfreq)

dat$region<-ifelse(dat$`d$state` %in% c("Maine","Massachusetts", "Rhode Island"),"New England","Mid-Atlantic")
dat$region<-ifelse(dat$`d$state` %in% c("Michigan","Indiana","Ohio"),"Mid-West",dat$region)
dat$region<-ifelse(dat$`d$state` %in% c("Missouri","North Carolina","Tennessee"),"Southeast",dat$region)


dat1<-dplyr::filter(dat,phenshift!="unknown")  
dat2<-dplyr::filter(dat,phenfreq!="unknown")  

stilt1<-filter(dat1,species=="Microstegium")
stilt2<-filter(dat2,species=="Microstegium")


ggplot(stilt1,aes(phenshift))+geom_bar(aes(fill=timechange))+xlab("Is phenology shifting?")+ggthemes::theme_few()+scale_fill_viridis_d()

ggplot(stilt2,aes(phenfreq))+geom_bar(aes(fill=freqchange))+xlab("Is phenology shifting?")+ggthemes::theme_few()+scale_fill_viridis_d()






ggplot(dat1,aes(phenshift))+geom_bar(aes(fill=timechange),position = "dodge")+facet_grid(region~species)

ggplot(dat2,aes(phenfreq))+geom_bar(aes(fill=phenfreq),position = "dodge")+facet_grid(~species)


ggplot(stilt1,aes(phenshift))+geom_bar()


dat1$resp<-ifelse(dat1$phenshift=="yes",1,0)
dat2$resp<-ifelse(dat2$phenfreq=="yes",1,0)

dat1<-dplyr::filter(dat1,!species %in% c("Pinus","Cytisus"))
dat2<-dplyr::filter(dat2,!species %in% c("Pinus","Cytisus"))



modstilt<-brm(resp~1,data=dat1,family = "bernoulli",warmup = 3000,iter = 4000,control=list(adapt_delta=.99))

library(brms)
mod<-brm(resp~dur+species+(dur+species|region),data=dat1,family = "bernoulli",warmup = 3000,iter = 4000,control=list(adapt_delta=.99))
mod1<-brm(resp~dur+species+(dur+species|region),data=dat2,family = "bernoulli",warmup = 3000,iter = 4000,control=list(adapt_delta=.99))


modz<-brm(resp~dur+species+(1|region),data=dat1,family = "bernoulli",warmup = 3000,iter = 4000,control=list(adapt_delta=.99))
modz1<-brm(resp~dur+species+(1|region),data=dat2,family = "bernoulli",warmup = 3000,iter = 4000,control=list(adapt_delta=.99))


new.dat<-data.frame(species=rep(c("Retnoutria","Artemesia","Berberis","Microstegium","Pueria","Celastrus","Vinetoxicum","Rhamnus",
                     "Elaeagnus","impatiens","Lonicera","Miscanthus","Rosa","Acer","Phragmotes"
                     ,"Centaurea","Ailanthis","Pastinaca","Ligustrum"),each=4),region=rep(c("Mid-Atlantic","New England", "Mid-West","Southeast"),19),dur=median(dat1$dur,na.rm=TRUE))


library(tidybayes)
goober<-epred_draws(modz1,newdata = new.dat)
goober2<-epred_draws(modz,newdata = new.dat)
pd<-position_dodge(width = 0.8)

ggplot(goober,aes(reorder(species,.epred),.epred))+stat_pointinterval(.width = c(.5,.9),aes(color=region),alpha=0.4,shape=0,position=pd)+
  stat_pointinterval(.width = c(.5,.9))+ylim(0,1)+ggthemes::theme_few()

ggplot(goober2,aes(reorder(species,.epred),.epred))+stat_pointinterval(.width = c(.5,.9),aes(color=region),alpha=0.4,shape=0,position=pd)+
  stat_pointinterval(.width = c(.5,.9))+ylim(0,1)+ggthemes::theme_few()


dat11<-dplyr::filter(dat1,resp==1)                               
unique(dat11$timechange)
dat11<-dplyr::filter(dat11,timechange %in% c("Later due to changing phenology","Earlier due to changing phenology", "At the same time despite changing phenology" ))

dat22<-dplyr::filter(dat2,resp==1) 
unique(dat22$freqchange)
dat22<-dplyr::filter(dat22,freqchange %in% c("At the same frequency despite a change in phenology","More frequently due to changing phenology", "Less frequently due to changing phenology"   ))

dat22$nimble<-ifelse(dat22$freqchange=="At the same frequency despite a change in phenology",0,1)
dat11$nimble<-ifelse(dat11$timechange=="At the same time despite changing phenology" ,0,1)

table(dat11$nimble)
dat11<-dplyr::filter(dat11,!species %in% c("impatiens","Miscanthus"))

mod.nim<-brm(nimble~dur+species+(1|region),data=dat11,family = "bernoulli",warmup = 3000,iter = 4000,control=list(adapt_delta=.99))


dat22<-dplyr::filter(dat22,!species %in% c("Pueria"))
mod.nim2<-brm(nimble~dur+species+(1|region),data=dat22,family = "bernoulli",warmup = 3000,iter = 4000,control=list(adapt_delta=.99))



new.dat2<-data.frame(species=rep(c("Retnoutria","Artemesia","Berberis","Microstegium","Pueria","Celastrus","Vinetoxicum","Rhamnus",
                                  "Elaeagnus","Lonicera","Rosa","Acer","Phragmotes"
                                  ,"Centaurea","Ailanthis","Pastinaca","Ligustrum"),each=4),region=rep(c("Mid-Atlantic","New England", "Mid-West","Southeast"),17),dur=median(dat1$dur,na.rm=TRUE))


new.dat3<-data.frame(species=rep(c("Retnoutria","Artemesia","Berberis","Microstegium","Celastrus","Vinetoxicum","Rhamnus",
                                  "Elaeagnus","impatiens","Lonicera","Miscanthus","Rosa","Acer","Phragmotes"
                                  ,"Centaurea","Ailanthis","Pastinaca","Ligustrum"),each=4),region=rep(c("Mid-Atlantic","New England", "Mid-West","Southeast"),18),dur=median(dat1$dur,na.rm=TRUE))



goober3<-epred_draws(mod.nim,newdata = new.dat2)
goober4<-epred_draws(mod.nim2,newdata = new.dat3)
pd<-position_dodge(width = 0.8)

ggplot(goober3,aes(reorder(species,.epred),.epred))+stat_pointinterval(.width = c(.5,.9),aes(color=region),alpha=0.4,shape=0,position=pd)+
  stat_pointinterval(.width = c(.5,.9))+ylim(0,1)+ggthemes::theme_few()

ggplot(goober4,aes(reorder(species,.epred),.epred))+stat_pointinterval(.width = c(.5,.9),aes(color=region),alpha=0.4,shape=0,position=pd)+
  stat_pointinterval(.width = c(.5,.9))+ylim(0,1)+ggthemes::theme_few()

dat111<-dplyr::filter(dat11,nimble==1)
dat222<-dplyr::filter(dat22,nimble==1)

ggplot(dat111,aes(species))+geom_bar(aes(fill=timechange),position="dodge",color="black")+scale_fill_manual(values=c("yellow","black"))+
  ggthemes::theme_few()

ggplot(dat222,aes(species))+geom_bar(aes(fill=freqchange),position="dodge",color="black")+scale_fill_manual(values=c("yellow","black"))+
  ggthemes::theme_few()


dat111$ealier<-ifelse(dat111$timechange=="Earlier due to changing phenology" ,1,0)

mod.Early<-brm(ealier~dur+species+(1|region),data=dat111,family = "bernoulli",warmup = 3000,iter = 4000,control=list(adapt_delta=.99))


dat$frequnknown<-ifelse(dat$freqchange=="I don’t know",1,0)
dat$phenunknown<-ifelse(dat$phenshift=="unknown",1,0)

umod1<-brm(phenunknown~dur+(1|species),data=dat,family="bernoulli")
umod2<-brm(frequnknown~dur+(1|species),data=dat,family="bernoulli")
conditional_effects(umod2)


Mv1<-dplyr::filter(dat1,species=="Microstegium")
Mv2<-dplyr::filter(dat2,species=="Microstegium")

Mv.mod<-brm(resp~dur+(1|region),data=Mv1,family="bernoulli")
Mv.mod2<-brm(resp~dur+(1|region),data=Mv2,family="bernoulli")
