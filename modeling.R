####started by Dan Feb 2026 updated June 25
###goal is to analyze the the climate change part of the RISCC phenology survey
##collaborators include Fitz and Bethany

###house keeping
rm(list=ls()) #removes everything in R environment
options(stringsAsFactors = FALSE) #stops R from automatically treating string as factors
options(mc.cores = parallel::detectCores()) # for Bayesian models, paralleled chain sampling
graphics.off() # removes anything in graphic environemnt
options(device = "RStudioGD")
library(dplyr)
library(ggplot2)
library(rstan)
library(brms)
library(tidybayes)
graphics.off()

if(length(grep("danielbuonaiuto", getwd()) > 0)) { # a fancy way to set your working direcyotu
  setwd("~/Documents/git/stilts/")
} else if(length(grep("fitz", getwd()) > 0)) {
  setwd("")}

timer<-read.csv("Input/timing_clean.csv")
freqer<-read.csv("Input/frequency_clean.csv")
###do managers see climate change?
class(timer$CC)

merge1<-select(timer,ResponseId,CC,dur,minCareer,state,species,agency)
merge2<-select(freqer,ResponseId,CC,dur,minCareer,state,species,agency)

merge1$dataset<-"time"
merge2$dataset<-"frequency"

d<-rbind(merge1,merge2)
rm(merge1)
rm(merge2)

d$minCareer_s<-scale(d$minCareer)
d$dur_s<-scale(d$dur)


model.1<-brm(CC~dur_s+minCareer_s+(1|state)+(1|species)+(1|agency)+(1|ResponseId),
             iter=5000,warmup=4000, control = list(adapt_delta=0.99), data=d,family = "bernoulli")

fixef(model.1)
fixef(model.1)

spdata<-data.frame(species=unique(d$species),dur_s=mean(d$dur_s,na.rm=TRUE),minCareer_s=mean(d$minCareer_s))
statedata<-data.frame(state=unique(d$state),dur_s=mean(d$dur_s,na.rm=TRUE),minCareer_s=mean(d$minCareer_s))
agedata<-data.frame(agency=unique(d$agency),dur_s=mean(d$dur_s,na.rm=TRUE),minCareer_s=mean(d$minCareer_s))


newdatery<-data.frame(dur_s=mean(d$dur_s,na.rm=TRUE),minCareer_s=c(-.5,.7))
newdatery2<-data.frame(minCareer_s=mean(d$minCareer_s,na.rm=TRUE),dur_s=c(-.8,.62))



predCC.car<-epred_draws(model.1,newdata=newdatery,re_formula = NA,ndraws = 100)
predCC.dur<-epred_draws(model.1,newdata=newdatery2,re_formula = NA,ndraws = 100)



oneA<-ggpubr::ggarrange(ggplot(predCC.car,aes(minCareer_s,.epred))+
  stat_lineribbon(.width = c(.5),alpha=0.5, color = "black")+scale_fill_brewer()+
  xlab("Career length")+ylab( "Likelihood of observing \nphenological shifts")+
  ggthemes::theme_few()+coord_cartesian(ylim=c(0.1,.8)),


ggplot(predCC.dur,aes(dur_s,.epred))+ #geom_line(aes(group=.draw),size=0.01,color = "grey25")+
  #geom_smooth(method="lm")+
  stat_lineribbon(.width = c(.5),alpha=0.5, color = "black")+scale_fill_brewer()+
  xlab("Years of management")+ylab( "Likelihood of observing \nphenological shifts")+
  ggthemes::theme_few()+
  coord_cartesian(ylim=c(.1,.8)),labels=c("a)","b)"),common.legend = TRUE)



preds<-epred_draws(model.1,newdata=spdata,re_formula = ~(1|species))
predstate<-epred_draws(model.1,newdata=statedata,re_formula = ~(1|state))
predage<-epred_draws(model.1,newdata=agedata,re_formula = ~(1|agency))



unique(predstate$state)
states<-c("NY","ME","VT","MD","ON","MA","MN",
          "IN","MI","WI","VA","PE","NC","RI","MO","OH","TN","WV","NH","NJ","CT","KY","QC","KS")

statedata$states<-states

oneB<-ggpubr::ggarrange(ggplot(preds,aes(reorder(species,.epred),.epred))+
                          stat_pointinterval(.width = .5)+#ylim(0,1)+
                          ggthemes::theme_few()+coord_cartesian(ylim=c(.1,.9))+
                          ylab( "Likelihood of observing \nphenological shifts")+xlab("species")+
                          theme(axis.text.x = element_text(angle = 45, hjust = 1)),
                        
                        ggplot(predstate,aes(reorder(state,.epred),.epred))+
                          stat_pointinterval(.width = .5)+#ylim(0,1)+
                          ggthemes::theme_few()+coord_cartesian(ylim=c(.1,.9))+
                          ylab( "Likelihood of observing \nphenological shifts")+xlab("state/province")+
                          theme(axis.text.x = element_text(angle = 45, hjust = 1)),
                        
                        ggplot(predage,aes(reorder(agency,.epred),.epred))+
                          stat_pointinterval(.width = .5)+coord_cartesian(ylim=c(.1,.9))+
                          #ylim(0,1)
                          ggthemes::theme_few()+
                          ylab( "Likelihood of observing \nphenological shifts")+xlab("agency")+
                          theme(axis.text.x = element_text(angle = 45, hjust = 1)),nrow=1,ncol=3,labels=c("c)","d)","e)"))

ggpubr::ggarrange(oneA,oneB,ncol=1,nrow=2) ###FIgure 1


####### next
timer2<-read.csv("Input/timing_cleanCC.csv")
freqer2<-read.csv("Input/frequency_cleanCC.csv")
timer2$minCareer_s<-scale(timer2$minCareer)
timer2$dur_s<-scale(timer2$dur)

freqer2$minCareer_s<-scale(freqer2$minCareer)
freqer2$dur_s<-scale(freqer2$dur)


model.mc1<-brm(management_shift~dur_s+minCareer_s+(1|state)+(1|species)+(1|agency)+(1|ResponseId),
               iter=5000,warmup=4000, control = list(adapt_delta=0.99),
               data=timer2,family = "bernoulli")
summary(model.mc1,prob = .5)

model.mc2<-brm(management_shift~dur_s+minCareer_s+(1|state)+(1|species)+(1|agency)+(1|ResponseId),
               iter=5000,warmup=4000, control = list(adapt_delta=0.99),
               data=freqer2,family = "bernoulli")

summary(model.mc2,prob = .5)

merge11<-select(timer2,ResponseId,management_shift,dur,minCareer,state,species,agency)
merge22<-select(freqer2,ResponseId,management_shift,dur,minCareer,state,species,agency)

merge11$dataset<-"time"
merge22$dataset<-"frequency"
d2<-rbind(merge11,merge22)
d2$minCareer_s<-scale(d2$minCareer)
d2$dur_s<-scale(d2$dur)


model.mc1<-brm(management_shift~dur_s+minCareer_s+dataset+(dataset|state)+(dataset|species)+(dataset|agency)+(1|ResponseId),
               iter=5000,warmup=4000, control = list(adapt_delta=0.99),
               data=d2,family = "bernoulli")

conditional_effects(model.mc1,prob = .5)

spdata2<-data.frame(species=rep(unique(d2$species),2),dataset=rep(c("time","frequency"),each=22),dur_s=mean(d2$dur_s,na.rm=TRUE),minCareer_s=mean(d2$minCareer_s))
statedata2<-data.frame(state=rep(unique(d2$state),2),dataset=rep(c("time","frequency"),each=24),dur_s=mean(d2$dur_s,na.rm=TRUE),minCareer_s=mean(d2$minCareer_s))
agedata2<-data.frame(agency=rep(unique(d2$agency),2),dataset=rep(c("time","frequency"),each=9),dur_s=mean(d2$dur_s,na.rm=TRUE),minCareer_s=mean(d2$minCareer_s))


preds2<-epred_draws(model.mc1,newdata=spdata2,re_formula = ~(1|species))
predstate2<-epred_draws(model.mc1,newdata=statedata2,re_formula = ~(1|state))


predage2<-epred_draws(model.mc1,newdata=agedata2,re_formula = ~(1|agency))
options(device = "quartz")

ggplot(preds2,aes(reorder(species,.epred),.epred))+
  stat_pointinterval(.width = .5,aes(color=dataset))+#ylim(0,1)+
  ggthemes::theme_few()+coord_cartesian(ylim=c(0,1))+
  ylab( "Likelihood adjusting management")+xlab("species")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(preds2,
       aes(x = reorder(species, .epred),
           y = .epred,
           fill = dataset)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_errorbar(
    aes(ymin = .lower, ymax = .upper),
    position = position_dodge(width = 0.9),
    width = 0.2
  ) +
  ggthemes::theme_few() +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    y = "Likelihood adjusting management",
    x = "species"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggplot(predstate2,aes(reorder(state,.epred),.epred))+
  stat_pointinterval(.width = .5,aes(color=dataset))+#ylim(0,1)+
  ggthemes::theme_few()+coord_cartesian(ylim=c(0,1))+
  ylab( "Likelihood adjusting management")+xlab("species")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(predage2,aes(reorder(agency,.epred),.epred))+
  stat_pointinterval(.width = .5,aes(color=dataset))+#ylim(0,1)+
  ggthemes::theme_few()+coord_cartesian(ylim=c(0,1))+
  ylab( "Likelihood adjusting management")+xlab("species")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

