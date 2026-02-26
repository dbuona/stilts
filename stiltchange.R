####started by Dan Feb 2026
###goal is to analyze the the climate change part of the RISCC phenology survey
##collaborators include Fitz and Bethany

###house keeping
rm(list=ls()) #removes everything in R environment
options(stringsAsFactors = FALSE) #stops R from automatically treating string as factors
options(mc.cores = parallel::detectCores()) # for Bayesian models, paralleled chain sampling
graphics.off() # removes anything in graphic environemnt



if(length(grep("dbuona", getwd()) > 0)) { # a fancy way to set your working direcyotu
  setwd("~/Documents/git/stilts/")
} else if(length(grep("fitz", getwd()) > 0)) {
  setwd("")}

d<-read.csv("data/phenology survey 2025.xlsx - full_datasheet.csv")

demos<-d[,1:16] ##this seperates out the 16 columns at are demographics

cols_to_keep <- grep("114", names(d), value = TRUE) ##this makes a list of all columns names related to question 114


timechange <- d[ , cols_to_keep] # Subset the data frame to keep only those columns

colnames(timechange)<-c("Retnoutria","Artemesia","Berberis","Microstegium","Pueria","Celastrus",  # here I rename the columns based on their order in the Qualtrics table
                        "Pinus","Vinetoxicum","Rhamnus",
                        "Elaeagnus","impatiens","Lonicera","Miscanthus",
                        "Rosa","Acer","Phragmotes","Cytisus","Centaurea",
                        "Ailanthis","Pastinaca","Ligustrum")



#same as above but for questions related to frequency of management
cols_to_keep <- grep("117", names(d), value = TRUE) 
freqchange <- d[ , cols_to_keep]

colnames(freqchange)<-c("Retnoutria","Artemesia","Berberis","Microstegium","Pueria","Celastrus","Pinus","Vinetoxicum","Rhamnus",
  "Elaeagnus","impatiens","Lonicera","Miscanthus","Rosa","Acer","Phragmotes","Cytisus","Centaurea","Ailanthis","Pastinaca","Ligustrum")


##same as above, but related to how long you've been managing the specific species
cols_to_keep <- grep("118", names(d), value = TRUE)
dur<-d[ , cols_to_keep]
colnames(dur)<-c("Retnoutria","Artemesia","Berberis","Microstegium","Pueria","Celastrus","Pinus","Vinetoxicum","Rhamnus",
                   "Elaeagnus","impatiens","Lonicera","Miscanthus","Rosa","Acer","Phragmotes","Cytisus","Centaurea","Ailanthis","Pastinaca","Ligustrum")

d$Q35 ##yhis is the question
d$state ## this is the state they arwe part of
timechange<-as.data.frame(cbind(timechange,d$Q35,d$state))
freqchange<-as.data.frame(cbind(freqchange,d$Q35,d$state))
dur<-as.data.frame(cbind(dur,d$Q35,d$state))


timechange<-tidyr::gather(timechange,"species","timechange",1:21) ##switches data to long formate for regression analysis 
freqchange<-tidyr::gather(freqchange,"species","freqchange",1:21)
dur<-tidyr::gather(dur,"species","dur",1:21)

freqchange<-dplyr::select(freqchange,-species,-`d$Q35`,-`d$state`) ### get ride of species colummn
dur<-dplyr::select(dur,-species,-`d$Q35`,-`d$state`)

dat<-cbind(freqchange,timechange,dur) ## bind the three datasets together
dat[dat == ""] <- NA # convert blanks to NAs

#The data is now formated

### new column to assess whether on not climate change timing shifts are reported
dat$phenshift<-NA
dat$phenshift<-ifelse(dat$timechange %in% c("At the same time due to no change in phenology","Other","Earlier or later for other reasons"),"no","yes")
dat$phenshift<-ifelse(dat$timechange=="I don't know","unknown",dat$phenshift)

### new column to assess whether on not climate change frequency shifts are reported
dat$phenfreq<-NA
dat$phenfreq<-ifelse(dat$freqchange %in% c("At the same frequency due to no change in phenology","More or less frequently for other reasons","Other" ),"no","yes")
dat$phenfreq<-ifelse(dat$freqchange=="I don’t know" ,"unknown",dat$phenfreq)

###since there are no lat long trends, I combined states into regions
dat$region<-ifelse(dat$`d$state` %in% c("Maine","Massachusetts", "Rhode Island"),"New England","Mid-Atlantic")
dat$region<-ifelse(dat$`d$state` %in% c("Michigan","Indiana","Ohio"),"Mid-West",dat$region)
dat$region<-ifelse(dat$`d$state` %in% c("Missouri","North Carolina","Tennessee"),"Southeast",dat$region)

## make two subsheet with the unknowns removed in phenshift and freqshift respectively
dat1<-dplyr::filter(dat,phenshift!="unknown")  
dat2<-dplyr::filter(dat,phenfreq!="unknown")  

#this was for stiltgrass only anayses
#stilt1<-filter(dat1,species=="Microstegium")
#stilt2<-filter(dat2,species=="Microstegium")




ggplot(dat1,aes(phenshift))+geom_bar(aes(fill=timechange),position = "dodge")+facet_grid(region~species)

ggplot(dat2,aes(phenfreq))+geom_bar(aes(fill=phenfreq),position = "dodge")+facet_grid(~species)



dat1$resp<-ifelse(dat1$phenshift=="yes",1,0) # convert to numeric
dat2$resp<-ifelse(dat2$phenfreq=="yes",1,0)

dat1<-dplyr::filter(dat1,!species %in% c("Pinus","Cytisus")) # remove species with little data
dat2<-dplyr::filter(dat2,!species %in% c("Pinus","Cytisus"))


####regressions are found below, stopping here for now because I am not sure if fitz
### will want to use bayesian or lmer
##Note, we should add respondent as a random effect, so need to go back and merge


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
