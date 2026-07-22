####started by Dan Feb 2026 updated June 25
###goal is to analyze the the climate change part of the RISCC phenology survey
##collaborators include Fitz and Bethany

###house keeping
rm(list=ls()) #
removes everything in R environment
options(stringsAsFactors = FALSE) #stops R from automatically treating string as factors
options(mc.cores = parallel::detectCores()) # for Bayesian models, paralleled chain sampling
graphics.off() # removes anything in graphic environemnt
library(dplyr)
library(ggplot2)
library(brms)
library(tidybayes)

if(length(grep("dbuona", getwd()) > 0)) { # a fancy way to set your working direcyotu
  setwd("~/Documents/git/stilts/")
} else if(length(grep("fitz", getwd()) > 0)) {
  setwd("")}

d<-read.csv("data/phenology survey 2025.xlsx - full_datasheet.csv") ##read in data

demos<-d[,1:16] ##this seperates out the 16 columns at are demographics
colnames(demos) #which do we want
demos<-dplyr::select(demos,ResponseId,Q1,Q2,Q3)
colnames(demos)<-c("ResponseId","field","career","agency")

###now seperate phenology shift questions
cols_to_keep <- grep("114", names(d), value = TRUE) ##this makes a list of all columns names related to question 114
  timechange <- d[ , cols_to_keep] # Subset the data frame to keep only those columns
  colnames(timechange)<-c("Retnoutria","Artemesia","Berberis","Microstegium","Pueria","Celastrus",  # here I rename the columns based on their order in the Qualtrics table
                        "Pinus","Vinetoxicum","Rhamnus",
                        "Elaeagnus","impatiens","Lonicera","Miscanthus",
                        "Rosa","Acer","Phragmites","Cytisus","Centaurea",
                        "Ailanthis","Pastinaca","Ligustrum")


##bindthem with the demographics we care about
#timechange<-cbind(demos,timechange)  
if(FALSE){
#same as above but for questions related to frequency of management
cols_to_keep <- grep("117", names(d), value = TRUE) 
  freqchange <- d[ , cols_to_keep]
  colnames(freqchange)<-c("Retnoutria","Artemesia","Berberis",
                          "Microstegium","Pueria","Celastrus","Pinus",
                          "Vinetoxicum","Rhamnus",
                          "Elaeagnus","impatiens","Lonicera","Miscanthus","Rosa","Acer",
                          "Phragmotes","Cytisus","Centaurea","Ailanthis","Pastinaca","Ligustrum")
#freqchange<-cbind(demos,freqchange)  

  

}
##same as above, but related to how long you've been managing the specific species
cols_to_keep <- grep("118", names(d), value = TRUE)
  dur<-d[ , cols_to_keep]
  colnames(dur)<-c("Retnoutria","Artemesia","Berberis","Microstegium","Pueria","Celastrus","Pinus","Vinetoxicum","Rhamnus",
                   "Elaeagnus","impatiens","Lonicera","Miscanthus","Rosa","Acer","Phragmites","Cytisus","Centaurea","Ailanthis","Pastinaca","Ligustrum")

dur$ResponseId  
dur$ResponseId<-d$ResponseId #append dur with a response id for indexing


#now try and pull in 2026 data



#timechange<-tidyr::gather(timechange,"species","timechange",1:21) ##switches data to long formate for regression analysis 
#freqchange<-tidyr::gather(freqchange,"species","freqchange",1:21)

#dur<-tidyr::gather(dur,"species","dur",1:21)

#timechange<-dplyr::left_join(timechange,dur)
#freqchange<-dplyr::left_join(freqchange,dur)


timechange[timechange == ""] <- NA # convert blanks to NAs
#freqchange[freqchange == ""] <- NA # convert blanks to NAs


d2<-read.csv("data/phenology survey 2026.xlsx - full_datasheet.csv")

cols_to_keep2 <- grep("114", names(d2), value = TRUE) ##this makes a list of all columns names related to question 114
timechange2 <- d2[ , cols_to_keep2] # Subset the data frame to keep only those columns
colnames(timechange2)<-c("Celastrus","Lonicera","Retnoutria","Pueria","Rosa","Microstegium", "Elaeagnus","Berberis","Pinus", 
                         "Vinetoxicum","Rhamnus","Miscanthus","impatiens","Lonicera japonica","Artemisia","Acer",
                         "Phragmites","Ligustrum","Cytisus", "Centaurea","Ailanthis","Pastinaca")


timechange2[timechange2 == ""] <- NA # convert blanks to NAs



demos2<-d2[,1:16] ##this seperates out the 16 columns at are demographics
colnames(demos2) #which do we want
demos2<-dplyr::select(demos2,ResponseId,Q1,Q2,Q3)
colnames(demos2)<-c("ResponseId","field","career","agency")


cols_to_keep <- grep("118", names(d2), value = TRUE)
dur2<-d2[ , cols_to_keep]
colnames(dur2)<-c("Celastrus","Lonicera","Retnoutria","Pueria","Rosa","Microstegium", "Elaeagnus","Berberis","Pinus", 
                         "Vinetoxicum","Rhamnus","Miscanthus","impatiens","Lonicera japonica","Artemisia","Acer",
                         "Phragmites","Ligustrum","Cytisus", "Centaurea","Ailanthis","Pastinaca")


dur2$ResponseId<-d2$ResponseId #append dur with a response id for indexing
dur$state<-d$state
dur2$state<-d2$state

dur$juris<-d$Q35
dur2$juris<-d2$Q35


timechange$ResponseId<-d$ResponseId
timechange2$ResponseId<-d2$ResponseId

timechange<-tidyr::gather(timechange,"species","timechange",1:21)
timechange2<-tidyr::gather(timechange2,"species","timechange",1:22)
dur<-tidyr::gather(dur,"species","dur",1:21)
dur2<-tidyr::gather(dur2,"species","dur",1:22)

twen5<-left_join(dur,timechange)
twen6<-left_join(dur2,timechange2)

twen5<-left_join(twen5,demos)
twen6<-left_join(twen6,demos2)



twen5$surveyYear<-2025
twen6$surveyYear<-2026



dat<-rbind(twen5,twen6)
unique(dat$timechange)

#idk<-filter(dat,timechange %in% c("I don't know","Other,I don't know", "Other"))

###remove NA's frome response
dat<-filter(dat,!is.na(timechange))
dat$minCareer<-ifelse(dat$career=="10+ years",10,0.1)
dat$minCareer<-ifelse(dat$career=="2 to 5 years",2,dat$minCareer)
dat$minCareer<-ifelse(dat$career=="6 to 10 years",6,dat$minCareer)

unique(dat$career)

dat$idk<-ifelse(dat$timechange%in% c("I don't know","Other,I don't know", "Other"),0,1)

cor(dat$minCareer,dat$dur,use = "pairwise.complete.obs")
summary(lm(idk~dur+minCareer,data=dat))

dater<-filter(dat,idk==1)
unique(dater$timechange)

noCC<-c("At the same time due to no change in phenology","Earlier or later for other reasons",
        "Earlier or later for other reasons,At the same time due to no change in phenology",
        "Starting earlier or ending later for other reasons",
        "Starting earlier or ending later for other reasons,At the same time due to no change in phenology" )

dater$CC<-ifelse(dater$timechange %in% noCC,0,1)

dater$duration.cent<-scale(dater$dur,center = TRUE)
dater$career.cent<-scale(dater$minCareer,center = TRUE)

dater<-filter(dater,!state %in% c("","Wyoming","Alaska","Nova Scotia"))
table(dater$state)

quantile(dater$dur,na.rm=TRUE)
model.cc1<-brm(CC~dur+minCareer+(1|state)+(1|species)+(1|agency)+(1|ResponseId), iter=5000,warmup=4000, control = list(adapt_delta=0.99), data=dater,family = "bernoulli")

draws <- as.data.frame(brms::as_draws_df(model.cc1, variable = "b_Intercept"))

# 2. Transform to Probability Scale using the inverse-logit function
draws$prob_Intercept <- plogis(draws$b_Intercept)

# 3. Calculate mean and 50% Credible Intervals (25th and 75th percentiles)
mean(draws$prob_Intercept)
quantile(draws$prob_Intercept, probs = c(0.25, 0.75))



summary(model.cc1)
conditional_effects(model.cc1)


spdata<-data.frame(species=unique(dater$species),dur=mean(dater$dur,na.rm=TRUE),minCareer=mean(dater$minCareer))
statedata<-data.frame(state=unique(dater$state),dur=mean(dater$dur,na.rm=TRUE),minCareer=mean(dater$minCareer))
agedata<-data.frame(agency=unique(dater$agency),dur=mean(dater$dur,na.rm=TRUE),minCareer=mean(dater$minCareer))

range(dater$dur,na.rm=TRUE)

newdaters<-data.frame(dur=mean(dater$dur,na.rm=TRUE),minCareer=c(0.1,10))
newdaters2<-data.frame(minCareer=mean(dater$minCareer,na.rm=TRUE),dur=c(0,40))

predo<-epred_draws(model.cc1,newdata=newdaters,re_formula = NA)
predo2<-epred_draws(model.cc1,newdata=newdaters2,re_formula = NA)

oneA<-ggpubr::ggarrange(ggplot(predo,aes(minCareer,.epred))+
                          #geom_line(aes(group=.draw),size=0.01,color = "grey25")+
                          stat_lineribbon(.width = c(.5), color = "black", alpha=0.5)+scale_fill_brewer()+
                          
                          
  #geom_smooth(method="lm")+
                            xlab("Career length")+ylab( "Likelihood of observing \nphenological shifts")+
  ggthemes::theme_few()+coord_cartesian(ylim=c(.1,.8)),#+#ylim(0,1),


ggplot(predo2,aes(dur,.epred))+ #geom_line(aes(group=.draw),size=0.01,color = "grey25")+
  #geom_smooth(method="lm")+
  stat_lineribbon(.width = c(.5),alpha=0.5, color = "black")+scale_fill_brewer()+
  xlab("Years of management")+ylab( "Likelihood of observing \nphenological shifts")+
  ggthemes::theme_few()+
  coord_cartesian(ylim=c(.1,.8)),labels=c("a)","b)"),common.legend = TRUE)


library(tidybayes)
preds<-epred_draws(model.cc1,newdata=spdata,re_formula = ~(1|species))
predstate<-epred_draws(model.cc1,newdata=statedata,re_formula = ~(1|state))
predage<-epred_draws(model.cc1,newdata=agedata,re_formula = ~(1|agency))

unique(predstate$state)
states<-c("NY","ME","VT","MD","ON","MA","MN",
  "IN","MI","WI","VA","PE","NC","RI","MO","OH","TN","WV","NH","NJ","CT","KY","QC")

statedata$states<-states

oneB<-ggpubr::ggarrange(ggplot(preds,aes(reorder(species,.epred),.epred))+
  stat_pointinterval(.width = .5)+#ylim(0,1)+
    ggthemes::theme_few()+coord_cartesian(ylim=c(.1,.8))+
  ylab( "Likelihood of observing \nphenological shifts")+xlab("species")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)),

ggplot(predstate,aes(reorder(state,.epred),.epred))+
  stat_pointinterval(.width = .5)+#ylim(0,1)+
  ggthemes::theme_few()+coord_cartesian(ylim=c(.1,.8))+
  ylab( "Likelihood of observing \nphenological shifts")+xlab("state/province")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)),

ggplot(predage,aes(reorder(agency,.epred),.epred))+
  stat_pointinterval(.width = .5)+coord_cartesian(ylim=c(.1,.8))+
  #ylim(0,1)
ggthemes::theme_few()+
  ylab( "Likelihood of observing \nphenological shifts")+xlab("agency")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)),nrow=1,ncol=3,labels=c("c)","d)","e)"))

ggpubr::ggarrange(oneA,oneB,ncol=1,nrow=2) ###FIgure 1

#question2 if you see cliamte change what determines whether or not you shift?

daterCC<-filter(dater,CC==1)

unique(daterCC$timechange)

NoChangers<-c("At the same time despite changing phenology","At the same time due to no change in phenology,At the same time despite changing phenology",
              "Starting earlier or ending later for other reasons,At the same time despite changing phenology")

daterCC$management_shift<-ifelse(daterCC$timechange %in% NoChangers,0,1)

model.mc1<-brm(management_shift~dur+minCareer+(1|state)+(1|species)+(1|agency)+(1|ResponseId), iter=5000,warmup=4000, control = list(adapt_delta=0.99), data=daterCC,family = "bernoulli")
summary(model.mc1,prob = .5)

# 1. Extract the posterior draws for the intercept
draws <- as.data.frame(brms::as_draws_df(model.mc1, variable = "b_Intercept"))

# 2. Transform to Probability Scale using the inverse-logit function
draws$prob_Intercept <- plogis(draws$b_Intercept)

# 3. Calculate mean and 50% Credible Intervals (25th and 75th percentiles)
mean(draws$prob_Intercept)
quantile(draws$prob_Intercept, probs = c(0.25, 0.75))



spdata2<-data.frame(species=unique(daterCC$species),dur=mean(daterCC$dur,na.rm=TRUE),minCareer=mean(daterCC$minCareer))
statedata2<-data.frame(state=unique(daterCC$state),dur=mean(daterCC$dur,na.rm=TRUE),minCareer=mean(daterCC$minCareer))
agedata2<-data.frame(agency=unique(daterCC$agency),dur=mean(daterCC$dur,na.rm=TRUE),minCareer=mean(daterCC$minCareer))
#jurdata2<-data.frame(juris=unique(daterCC$juris),dur=mean(daterCC$dur,na.rm=TRUE),minCareer=mean(daterCC$minCareer))




predskis<-epred_draws(model.mc1,newdata=spdata2,re_formula = ~(1|species))
predstate2<-epred_draws(model.mc1,newdata=statedata2,re_formula = ~(1|state))
predage2<-epred_draws(model.mc1,newdata=agedata2,re_formula = ~(1|agency))
#jurdage2<-epred_draws(model.mc1,newdata=jurdata2,re_formula = ~(1|juris))

unique(daterCC$species)
annualBiennial<-c("Microstegium","Centaurea", "Pastinaca","impatiens")
Vine<-c("Pueria","Celastrus","Lonicera japonica")
perennial<-c("Retnoutria","Artemesia","Miscanthus","Phragmites","Vinetoxicum")
shrub<-c("Berberis","Rhamnus","Elaeagnus","Lonicera","Rosa","Ligustrum","Cytisus")
tree<-c("Acer","Ailanthis")
predsp2<-predskis
predsp2$LHC<-NA
predsp2$LHC<-ifelse(predsp2$species %in% c(annualBiennial,perennial),"herbacious","tree")
predsp2$LHC<-ifelse(predsp2$species %in% c(Vine),"vine",predsp2$LHC)
#predsp2$LHC<-ifelse(predsp2$species %in% c(perennial),"perennial",predsp2$LHC)  
predsp2$LHC<-ifelse(predsp2$species %in% c(shrub),"shrub",predsp2$LHC) 
predsp2$woody<-NA
predsp2$woody<-ifelse(predsp2$species %in% c(shrub,tree,Vine),"woody","non-woody")

p2a<-ggplot(predsp2,aes(reorder(species,.epred),.epred))+
  stat_eye(.width=0.5)+coord_cartesian(ylim = c(.5,1))+ggthemes::theme_few()+
  ylab( "Likelihood of management shifts")+xlab("species")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(predsp2,aes(reorder(LHC,.epred),.epred))+
  stat_eye(.width=0.5)+ylim(0,1)+ggthemes::theme_few()+
  ylab( "Likelihood of management shifts")+xlab("life-history class")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p2b<-ggplot(predsp2,aes(reorder(LHC,.epred),.epred))+
  stat_eye(.width = c(.5))+coord_cartesian(ylim = c(.5,1))+ggthemes::theme_few()+
  ylab( "Likelihood of management shifts")+xlab("life-history class")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p2c<-ggplot(predstate2,aes(reorder(state,.epred),.epred))+
  stat_eye(.width=0.5)+coord_cartesian(ylim = c(.7,1))+ggthemes::theme_few()+
  ylab( "Likelihood of management shifts")+xlab("state/province")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p2d<-ggplot(predage2,aes(reorder(agency,.epred),.epred))+
  stat_eye(.width=0.5)+coord_cartesian(ylim = c(.7,1))+ggthemes::theme_few()+
  ylab( "Likelihood of management shifts")+xlab("agency")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggpubr::ggarrange(p2a,p2b,p2c,p2d,widths = c(.6,.4))

#ggplot(jurdage2,aes(reorder(juris,.epred),.epred))+
#  stat_pointinterval()+ylim(0,1)+ggthemes::theme_few()+
#  ylab( "Likelihood of management shifts")+xlab("jurisdiction")+
#  theme(axis.text.x = element_text(angle = 45, hjust = 1))


daterDir<-filter(daterCC,management_shift==1)
unique(daterDir$timechange)

earlier<-c( "Earlier due to changing phenology","Treatment starting earlier due to changing phenology",
            "Treatment starting earlier due to changing phenology,Starting earlier or ending later for other reasons",
            "Treatment starting earlier due to changing phenology,At the same time due to no change in phenology")

later<-c("Later due to changing phenology","Treatment ending later due to changing phenology",
         "Treatment ending later due to changing phenology,Starting earlier or ending later for other reasons")

earlLate<-c("Earlier due to changing phenology,Later due to changing phenology","Earlier due to changing phenology,Later due to changing phenology,At the same time despite changing phenology",
            "Treatment starting earlier due to changing phenology,Treatment ending later due to changing phenology",
            "Treatment starting earlier due to changing phenology,Treatment ending later due to changing phenology,Starting earlier or ending later for other reasons"
            )
daterDir$howshift<-NA
daterDir$howshift<-ifelse(daterDir$timechange %in% earlier,"earlier","later")
daterDir$howshift<-ifelse(daterDir$timechange %in% earlLate,"earlier & later",daterDir$howshift)
daterDir$howshift <- relevel(as.factor(daterDir$howshift), ref = "earlier & later")

modhow1 <- brm(
  formula = howshift ~ (1|species)+(1|ResponseId), iter=5000,warmup=4000, control = list(adapt_delta=0.99),
  data = daterDir,
  family = "categorical")

ranef(modhow1,probs =c(0.25,0.75))

predskis<-epred_draws(modhow1,newdata = spdata2,re_formula = ~(1|species))
pd<-position_dodge(width=0.3)
p3a<-ggplot(predskis,aes(reorder(species,.epred),.epred))+
  stat_eye(aes(fill=.category,),position=pd,.width=0.5,alpha=0.8)+facet_wrap(~species,nrow=2,scales="free_x")+
  ggthemes::theme_few()+ylab("mangement shift")+xlab("species")+scale_fill_viridis_d(option="C")+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

predskis$LHC<-NA
predskis$LHC<-ifelse(predskis$species %in% c(annualBiennial,perennial),"herbacious","tree")
predskis$LHC<-ifelse(predskis$species %in% c(Vine),"vine",predskis$LHC)
#predskis$LHC<-ifelse(predskis$species %in% c(perennial),"perennial",predskis$LHC)  
predskis$LHC<-ifelse(predskis$species %in% c(shrub),"shrub",predskis$LHC) 
predskis$woody<-NA
predskis$woody<-ifelse(predskis$species %in% c(shrub,tree,Vine),"woody","non-woody")

p3ba<-ggplot(predskis,aes(reorder(LHC,.epred),.epred))+stat_eye(aes(fill=.category),.width=0.5,position=pd,alpha=0.8)+
  ggthemes::theme_few()+
  ylab("mangement shift")+xlab("functional type")+scale_fill_viridis_d(option="C")
p3b<-ggplot(predskis,aes(reorder(woody,.epred),.epred))+stat_eye(aes(fill=.category),.width=0.5,position=pd,alpha=0.8)+ggthemes::theme_few()+
  ylab("mangement shift")+xlab("functional type")+scale_fill_viridis_d(option="C")


ggpubr::ggarrange(p3a,p3ba,common.legend = TRUE,widths = c(.7,.3))
stop()

if(FALSE){
conditional_effects(model.cc1,prob = .89)

timechange2$phenshift<-NA
table(timechange2$timechange)


timechange2<-filter(timechange2,!is.na(timechange))
timechange2$phenshift<-ifelse(timechange2$timechange %in% c("Treatment starting earlier due to changing phenology", "Treatment ending later due to changing phenology",
                                                            "Treatment starting earlier due to changing phenology,Treatment ending later due to changing phenology",
                                                            "Treatment starting earlier due to changing phenology,Treatment ending later due to changing phenology,Starting earlier or ending later for other reasons",
                                                            "Treatment ending later due to changing phenology,Starting earlier or ending later for other reasons",
                                                            "Starting earlier or ending later for other reasons,At the same time due to no change in phenology"),"yes","no")  

timechange2$phenshift<-ifelse(timechange2$timechange %in% c("I don't know","Other","Other,I don't know"),"unknown",timechange2$phenshift)
ggplot(timechange2,aes(species))+geom_bar(aes(fill=phenshift),position="dodge")
timechange2<-dplyr::filter(timechange2,phenshift!="unknown")  
ggplot(timechange2,aes(phenshift))+geom_bar(aes(fill=timechange))+facet_grid(~species)





###useful columns
#timechange<-as.data.frame(cbind(timechange,d$Q35,d$state))
#freqchange<-as.data.frame(cbind(freqchange,d$Q35,d$state))
#colnames(timechange)[26:27]<-c("managementPeriod","state") ## give some columns better names
#colnames(freqchange)[26:27]<-c("managementPeriod","state")

#timechange[timechange == ""] <- NA # convert blanks to NAs
timechange<-filter(timechange,!is.na(timechange))
timechange$phenshift<-NA
table(timechange$timechange)
timechange$phenshift<-ifelse(timechange$timechange %in% c("At the same time due to no change in phenology","Earlier or later for other reasons","At the same time despite changing phenology",
                                                          "Earlier or later for other reasons,At the same time due to no change in phenology"),"no","yes")
timechange$phenshift<-ifelse(timechange$timechange %in% c("I don't know","Other"),"unknown",timechange$phenshift)

timechange<-dplyr::filter(timechange,phenshift!="unknown")  
oldy<-dplyr::select(timechange,species,timechange,phenshift)
oldy$survey<-"2025"
timechange2$survey<-"2026"
com1<-rbind(oldy,timechange2)
ggplot(com1,aes(survey))+geom_bar(aes(color=phenshift,fill=phenshift),position="fill")+facet_wrap(~species,nrow=3)+geom_hline(yintercept=.5)+
  scale_fill_viridis_d()
com1$species<-ifelse(com1$species=="Artemesia","Artemisia",com1$species)
ggplot(com1,aes(survey))+geom_bar(aes(color=phenshift,fill=phenshift),position="dodge")+facet_wrap(~species,nrow=3)+scale_fill_viridis_d()



cols_to_keep22 <- grep("117", names(d2), value = TRUE) 
freqchange2 <- d2[ , cols_to_keep22]
colnames(freqchange2)<-c("Celastrus","Lonicera","Retnoutria","Pueria","Rosa","Microstegium", "Elaeagnus","Berberis","Pinus", 
                         "Vinetoxicum","Rhamnus","Miscanthus","impatiens","Lonicera japonica","Artemisia","Acer",
                         "Phragmotes","Ligustrum","Cytisus", "Centaurea","Ailanthis","Pastinaca")

freqchange2<-tidyr::gather(freqchange2,"species","freqchange",1:22)
freqchange2[freqchange2 == ""] <- NA # convert blanks to NAs
freqchange2<-filter(freqchange2,!is.na(freqchange))
freqchange2$freqshift<-NA
table(freqchange2$freqchange)
freqchange2$freqshift<-ifelse(freqchange2$freqchange %in% c("Less frequently due to changing phenology","Less frequently due to changing phenology,Other",
                                                            "More frequently due to changing phenology","More frequently due to changing phenology,Less frequently due to changing phenology"),"yes","no")
freqchange2$freqshift<-ifelse(freqchange2$freqchange %in% c("I don't know","Other","Other,I don't know"),"unknown",freqchange2$freqshift)

freqchange<-filter(freqchange,!is.na(freqchange))
freqchange$freqshift<-NA
table(freqchange$freqchange)
freqchange$freqshift<-ifelse(freqchange$freqchange %in% c("Less frequently due to changing phenology","More frequently due to changing phenology",
                                                          "More frequently due to changing phenology,At the same frequency despite a change in phenology"),"yes","no")                                                          
freqchange$freqshift<-ifelse(freqchange$freqchange %in% c("I don't know","Other","Other,I don't know"),"unknown",freqchange$freqshift)



freqchange$survey<-"2025"
freqchange2$survey<-"2026"
com2<-rbind(freqchange,freqchange2)
com2<-dplyr::filter(com2,freqshift!="unknown")
com2<-filter(com2,freqchange!="I don’t know" )
ggplot(com2,aes(survey))+geom_bar(aes(color=freqshift,fill=freqshift),position="fill")+facet_wrap(~species,nrow=3)+geom_hline(yintercept=.5)+
  scale_fill_viridis_d()
com2$species<-ifelse(com2$species=="Artemesia","Artemisia",com2$species)
ggplot(com2,aes(survey))+geom_bar(aes(color=freqshift,fill=freqshift),position="dodge")+facet_wrap(~species,nrow=3)+scale_fill_viridis_d()

com1n<-filter(com1,timechange %in% c("At the same time due to no change in phenology","At the same time despite changing phenology"))
unique(com1n$timechange)
ggplot(com1n,aes(survey))+geom_bar(aes(color=timechange,fill=timechange))+facet_wrap(~species,nrow=3)+scale_fill_viridis_d()
ggplot(com1n,aes(survey))+geom_bar(aes(color=timechange,fill=timechange),position="fill")+facet_wrap(~species,nrow=3)+scale_fill_viridis_d()


  
 unique(com2$freqchange)
 com2n<-filter(com2, freqchange %in%  c("At the same frequency due to no change in phenology","At the same frequency despite a change in phenology"))

 
ggplot(com1n,aes(survey))+geom_bar(aes(color=timechange,fill=timechange))+facet_wrap(~species,nrow=3)+scale_fill_viridis_d()
ggplot(com2n,aes(survey))+geom_bar(aes(color=freqchange,fill=freqchange),position="fill")+facet_wrap(~species,nrow=3)+scale_fill_viridis_d()
unique(com1$timechange)
com1$cc<-ifelse(com1$timechange %in% c("At the same time due to no change in phenology","Earlier or later for other reasons", "Earlier or later for other reasons,At the same time due to no change in phenology"),"no","yes") 

ggplot(com1,aes(cc))+geom_bar(aes(fill=phenshift),position="fill")+facet_wrap(~species,nrow=3)+scale_fill_viridis_d(option = "C")
ggplot(com1,aes(phenshift))+geom_bar(aes(fill=cc))+facet_wrap(~species,nrow=3)+scale_fill_viridis_d(option = "C")


ggplot(com1,aes(species))+geom_bar(aes(fill=cc,color=phenshift),position="fill")+facet_wrap(~species,nrow=3,scale="free")+
  scale_fill_viridis_d()+scale_colour_manual(values=c("grey","black"))
jpeg("longtermphenshift.jpeg")
ggplot(com1,aes(species))+geom_bar(aes(fill=cc,color=phenshift))+facet_wrap(~species,nrow=3,scale="free")+
  scale_fill_viridis_d()+scale_colour_manual(values=c("grey","black"))
dev.off()

unique(com2$freqchange)
com2$cc<-ifelse(com2$freqchange %in%c("At the same frequency due to no change in phenology","More or less frequently for other reasons"),"no","yes")
ggplot(com2,aes(species))+geom_bar(aes(fill=cc,color=freqshift),position="fill")+facet_wrap(~species,nrow=3,scale="free")+
  scale_fill_viridis_d(option="C")+scale_colour_manual(values=c("grey","black"))
jpeg("longtermfreqshift.jpeg")
ggplot(com2,aes(species))+geom_bar(aes(fill=cc,color=freqshift))+facet_wrap(~species,nrow=3,scale="free")+
  scale_fill_viridis_d(option="C")+scale_colour_manual(values=c("grey","black"))
dev.off()

commytime<-filter(com1,phenshift=="yes")
table(commytime$timechange)

commytime<-filter(commytime,!timechange %in% c("At the same time due to no change in phenology,At the same time despite changing phenology"," Earlier due to changing phenology,Later due to changing phenology,At the same time despite changing phenology",
                                               "Starting earlier or ending later for other reasons,At the same time due to no change in phenology", "Treatment ending later due to changing phenology,Starting earlier or ending later for other reasons"))
commytime$timechange<-ifelse(commytime$timechange=="Treatment starting earlier due to changing phenology","Earlier due to changing phenology",commytime$timechange)
commytime$timechange<-ifelse(commytime$timechange=="Treatment ending later due to changing phenology","Later due to changing phenology",commytime$timechange)
commytime$timechange<-ifelse(commytime$timechange=="Treatment starting earlier due to changing phenology,Treatment ending later due to changing phenology","Earlier due to changing phenology,Later due to changing phenology",commytime$timechange)

commytime$timechange<-ifelse(commytime$timechange=="Earlier due to changing phenology,Later due to changing phenology,At the same time despite changing phenology","Earlier due to changing phenology,Later due to changing phenology",commytime$timechange)
commytime$timechange<-ifelse(commytime$timechange=="Treatment starting earlier due to changing phenology,Treatment ending later due to changing phenology,Starting earlier or ending later for other reasons","Earlier due to changing phenology,Later due to changing phenology",commytime$timechange)

ggplot(commytime,aes(species))+geom_bar(aes(fill=timechange))+facet_wrap(~species,nrow=3,scale="free")+
  scale_fill_viridis_d(option="B")

commyfreq<-filter(com2,freqshift=="yes")

table(commyfreq$freqchange)
commyfreq$freqchange<-ifelse(commyfreq$freqchange=="Less frequently due to changing phenology,Other","Less frequently due to changing phenology",commyfreq$freqchange)
commyfreq$freqchange<-ifelse(commyfreq$freqchange=="More frequently due to changing phenology,At the same frequency despite a change in phenology","More frequently due to changing phenology",commyfreq$freqchange)

commyfreq<-filter(commyfreq,freqchange!="More frequently due to changing phenology,Less frequently due to changing phenology")

ggplot(commyfreq,aes(species))+geom_bar(aes(fill=freqchange))+facet_wrap(~species,nrow=3,scale="free")+
  scale_fill_viridis_d(option="B")

stilt1<-filter(com1,species=="Microstegium")
stilt2<-filter(com2,species=="Microstegium")

ggpubr::ggarrange(
ggplot(stilt1,aes(species))+geom_bar(aes(fill=cc,color=phenshift))+facet_wrap(~species,nrow=3,scale="free")+
  scale_fill_viridis_d()+scale_colour_manual(values=c("grey","black")),
ggplot(stilt2,aes(species))+geom_bar(aes(fill=cc,color=freqshift))+facet_wrap(~species,nrow=3,scale="free")+
  scale_fill_viridis_d(option="C")+scale_colour_manual(values=c("grey","black")))

########Pause the data are now formated (i.e., clean!)##################Pause the data are now formated scale_colour_grey()(i.e., clean!)##########
####you could write out csv's at this stage to use for analyses and visualizations##### 
############I might consider renaming this file about cleaning###################
stop("not an error, I think this a a good place to pause in terms of understanding the cleaing work flow and preparing to integrate new survey data")




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
conditional_effects(modz)



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
conditional_effects(mod.nim)


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
}
