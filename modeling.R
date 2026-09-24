####started by Dan Feb 2026 updated June 25
###goal is to analyze the the climate change part of the RISCC phenology survey
##collaborators include Fitz and Bethany

###house keeping
rm(list=ls()) #removes everything in R environment
options(stringsAsFactors = FALSE) #stops R from automatically treating string as factors
options(mc.cores = parallel::detectCores()) # for Bayesian models, paralleled chain sampling
graphics.off() # removes anything in graphic environemnt
#options(device = "RStudioGD")
library(dplyr)
library(ggplot2)
library(rstan)
library(brms)
library(tidybayes)
library(geosphere)
library(rnaturalearth)
library(sf)

#graphics.off()

if(length(grep("dbuona", getwd()) > 0)) { # a fancy way to set your working direcyotu
  setwd("~/Documents/git/stilts/")
} else if(length(grep("fitz", getwd()) > 0)) {
  setwd("")}

timer<-read.csv("Input/timing_clean.csv")
freqer<-read.csv("Input/frequency_clean.csv")
###do managers see climate change?
class(timer$CC)
table(timer$juris)


merge1<-select(timer,ResponseId,CC,dur,minCareer,state,species,agency)
merge2<-select(freqer,ResponseId,CC,dur,minCareer,state,species,agency)

merge1$dataset<-"time"
merge2$dataset<-"frequency"

d<-rbind(merge1,merge2)
rm(merge1)
rm(merge2)

table(d$CC)
1103/(1648+1103)

if(FALSE){
temp<-read.csv("Input/climate_dif.csv")
temp<-temp[,c(2,4)]
colnames(temp)<-c("state","anomoly")
d<-left_join(d,temp)  
}
gren<-read.csv("Input/Spring_Shift_2001_2024.csv")
colnames(gren)[1]<-"state"
d<-left_join(d,gren)  
unique(d$state)
NECASC<-c("Maine", "New Hampshire", "Vermont","Massachusetts", "Rhode Island",
          "Connecticut","New York", "New Jersey","Pennsylvania","Delaware",
          "Maryland","Virginia", "Kentucky", "West Virginia")
d<-filter(d,state %in% NECASC)
d$minCareer_s<-scale(d$minCareer)
d$dur_s<-scale(d$dur)
d$greenup_s<-scale(d$slope)
d$green_var<-scale(d$se_doy)

d$species<-ifelse(d$species=="Artemesia","Artemisia",d$species) ##fix species name
table(d$CC)
759/(759+1125)

if(FALSE){
# Pick your list of states
my_states <- unique(d$state)

# Get map data
s_data <- map_data("state", region = my_states)

# Loop or aggregate to find the centroid of each state
results <- do.call(rbind, lapply(unique(s_data$region), function(st) {
  subset_data <- subset(s_data, region == st)
  # centroid() returns c(lon, lat)
  cent <- geosphere::centroid(subset_data[, c("long", "lat")]) 
  data.frame(state = st, lon = cent[1], lat = cent[2])
}))

print(results)
results$state <- tools::toTitleCase(results$state)


# 1. Fetch Canadian map data (SF Spatial Object)
canada_provs <- rnaturalearth::ne_states(country = "canada", returnclass = "sf")

# 2. Compute true centroids and pull out lat/lon columns
ca_centroids <- canada_provs %>%
  # Extract the true spatial center point
  st_centroid() %>% 
  # Bind coordinates as numeric columns
  mutate(
    lon = st_coordinates(.)[,1],
    lat = st_coordinates(.)[,2]
  ) %>%
  # Keep only the name and coordinate columns, drop spatial geometries
  as.data.frame() %>%
  select(state = name, lon, lat)

# 3. Clean names to match your Title Case requirement
ca_centroids$state <- tools::toTitleCase(tolower(ca_centroids$state))

print(ca_centroids)

results<-rbind(results,ca_centroids)



d<-left_join(d,results)
d$lat_sc<-scale(d$lat)

}
d$career<-as.factor(d$minCareer)
cor(d$greenup_s,d$green_var)
model.1<-brm(CC~dur+slope+(1|species)+(1|ResponseId),
             iter=8000,warmup=7000, control = list(adapt_delta=0.99), data=d,family = "bernoulli")

model.1a<-brm(CC~anom_s+dur_s+minCareer_s+anom_s:dur_s+anom_s:minCareer_s+minCareer_s:dur_s+(1|species)+(1|agency)+(1|ResponseId),
             iter=8000,warmup=7000, control = list(adapt_delta=0.99), data=d,family = "bernoulli")

range(d$anomoly,na.rm=TRUE)
model.1b<-brm(CC~dur_s*minCareer_sanom_s+(anom_s|species)+(lat_sc|agency)+(1|ResponseId),
              iter=8000,warmup=7000, control = list(adapt_delta=0.99), data=d,family = "bernoulli")

fixef(model.1,probs = c(.25,.75))
fixef(model.1a,probs = c(.25,.75))

conditional_effects(model.1,prob = .5)

conditional_effects(model.1b,prob = .5)
fixef(model.1b,probs = c(.25,.75))

spdata<-data.frame(species=unique(d$species),dur=mean(d$dur,na.rm=TRUE),
                   slope=mean(d$slope,na.rm=TRUE),green_var=mean(d$green_var,nm.rm=TRUE))


statedata<-data.frame(state=unique(d$state),dur_s=mean(d$dur_s,na.rm=TRUE),minCareer_s=mean(d$minCareer_s))
agedata<-data.frame(agency=unique(d$agency),dur_s=mean(d$dur_s,na.rm=TRUE),minCareer_s=mean(d$minCareer_s),sd_doy=mean(d$sd_doy,na.rm=TRUE))
agedata<-filter(agedata,agency!='Tribal/First Nation')###not in model


range(d$anom_s,na.rm=TRUE)



preds<-epred_draws(model.1,newdata=spdata,re_formula = ~(1|species))
predage<-epred_draws(model.1,newdata=agedata,re_formula = ~(1|agency))


ggplot(preds,aes(reorder(species,.epred),.epred))+
  stat_pointinterval(.width = .5)+#ylim(0,1)+
  ggthemes::theme_few()+coord_cartesian(ylim=c(.1,.9))+
  ylab( "Likelihood of observing \nphenological shifts")+xlab("species")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(predage,aes(reorder(agency,.epred),.epred))+
  stat_pointinterval(.width = .5)+#ylim(0,1)+
  ggthemes::theme_few()+coord_cartesian(ylim=c(.1,.9))+
  ylab( "Likelihood of observing \nphenological shifts")+xlab("agency")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


predstate<-epred_draws(model.1,newdata=statedata,re_formula = ~(1|state))
predstate<-left_join(predstate,temp)

ggplot(predstate,aes(anomoly,.epred))+
  geom_smooth(method="lm")+

  ggthemes::theme_few()+coord_cartesian(ylim=c(.1,.9))+
  ylab( "Likelihood of observing \nphenological shifts")+xlab("lat")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

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

if(FALSE){
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
}
merge11<-select(timer2,ResponseId,management_shift,dur,minCareer,state,species,agency)
merge22<-select(freqer2,ResponseId,management_shift,dur,minCareer,state,species,agency)

merge11$dataset<-"time"
merge22$dataset<-"frequency"
d2<-rbind(merge11,merge22)
d2$minCareer_s<-scale(d2$minCareer)
d2$dur_s<-scale(d2$dur)
d2$species<-ifelse(d2$species=="Artemesia","Artemisia",d2$species) ##fix species name

d2<-filter(d2, state %in% NECASC)
d2<-left_join(d2,gren)

d2$minCareer_s<-scale(d2$minCareer)
d2$dur_s<-scale(d2$dur)

d2$greenup_s<-scale(d2$slope)
d2$green_var<-scale(d2$se_doy)
cor(d2$green_var,d2$greenup_s)

model.mc1<-brm(management_shift~dur+slope+dataset+(dataset|species)+(1|ResponseId),
               iter=5000,warmup=4000, control = list(adapt_delta=0.99),
               data=d2,family = "bernoulli")
fixef(model.mc1,prob=c(.25,.75))
conditional_effects(model.mc1,prob = .5)

spdata2<-data.frame(species=rep(unique(d2$species),2),dataset=rep(c("time","frequency"),each=21),
                    dur=mean(d2$dur,na.rm=TRUE),slope=mean(d2$slope))
statedata2<-data.frame(state=rep(unique(d2$state),2),dataset=rep(c("time","frequency"),each=24),dur_s=mean(d2$dur_s,na.rm=TRUE),minCareer_s=mean(d2$minCareer_s))
agedata2<-data.frame(agency=rep(unique(d2$agency),2),dataset=rep(c("time","frequency"),each=9),dur_s=mean(d2$dur_s,na.rm=TRUE),minCareer_s=mean(d2$minCareer_s))


preds2<-epred_draws(model.mc1,newdata=spdata2,re_formula = ~(dataset|species))
predstate2<-epred_draws(model.mc1,newdata=statedata2,re_formula = ~(1|state))


predage2<-epred_draws(model.mc1,newdata=agedata2,re_formula = ~(1|agency))
#options(device = "quartz")

ggpubr::ggarrange(ggplot(preds2,aes(reorder(species,.epred),.epred))+
  stat_pointinterval(.width = c(.01,.5),aes(color=dataset,shape=dataset))+#ylim(0,1)+
  ggthemes::theme_few()+coord_cartesian(ylim=c(0,1))+
  ylab( "Likelihood adjusting management")+xlab("species")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+scale_color_manual(values=c("firebrick4","seagreen3")),


ggplot(predstate2,aes(reorder(state,.epred),.epred))+
  stat_pointinterval(.width = c(.01,.5),aes(color=dataset,shape=dataset))+#ylim(0,1)+
  ggthemes::theme_few()+coord_cartesian(ylim=c(0,1))+
  ylab( "Likelihood adjusting management")+xlab("state")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+scale_color_manual(values=c("firebrick4","seagreen3")),

ggplot(predage2,aes(reorder(agency,.epred),.epred))+
  stat_pointinterval(.width = c(.01,.5),aes(color=dataset,shape=dataset))+#ylim(0,1)+
  ggthemes::theme_few()+coord_cartesian(ylim=c(0,1))+
  ylab( "Likelihood adjusting management")+xlab("agency")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+scale_color_manual(values=c("firebrick4","seagreen3")),nrow=1,common.legend = TRUE)

######how is management shifting
timer3<-read.csv("Input/timing_cleanDir.csv")
freqer3<-read.csv("Input/frequency_cleanDir.csv")

earlier<-c( "Earlier due to changing phenology","Treatment starting earlier due to changing phenology",
            "Treatment starting earlier due to changing phenology,Starting earlier or ending later for other reasons",
            "Treatment starting earlier due to changing phenology,At the same time due to no change in phenology")

later<-c("Later due to changing phenology","Treatment ending later due to changing phenology",
         "Treatment ending later due to changing phenology,Starting earlier or ending later for other reasons")

earlLate<-c("Earlier due to changing phenology,Later due to changing phenology","Earlier due to changing phenology,Later due to changing phenology,At the same time despite changing phenology",
            "Treatment starting earlier due to changing phenology,Treatment ending later due to changing phenology",
            "Treatment starting earlier due to changing phenology,Treatment ending later due to changing phenology,Starting earlier or ending later for other reasons"
)

timer3$howshift<-NA
timer3$howshift<-ifelse(timer3$timechange %in% earlier,"earlier","later")
timer3$howshift<-ifelse(timer3$timechange %in% earlLate,"earlier & later",timer3$howshift)
timer3$howshift <- relevel(as.factor(timer3$howshift), ref = "earlier & later")


table(freqer3$freqchange)
freqer3$howshift<-ifelse(freqer3$freqchange %in% c("Less frequently due to changing phenology","Less frequently due to changing phenology,Other"),"less","more")

merge111<-select(timer3,ResponseId,dur,minCareer,state,species,agency,howshift)
merge222<-select(freqer3,ResponseId,dur,minCareer,state,species,agency,howshift)


timer3$species<-ifelse(timer3$species=="Artemesia","Artemisia",timer3$species) ##fix species name

timer3<-filter(timer3,state %in% NECASC)

timer3<-left_join(timer3,gren)
timer3$minCareer_s<-scale(timer3$minCareer)
timer3$dur_s<-scale(timer3$dur)

timer3$greenup_s<-scale(timer3$slope)
timer3$green_var<-scale(timer3$se_doy)

timer3<-filter(timer3,howshift!="earlier & later")
timer3$how<-ifelse(timer3$howshift=="earlier",0,1)
modhow1 <- brm(
  formula = how~(1|species)+(1|ResponseId), iter=5000,warmup=4000, control = list(adapt_delta=0.99),
  data = timer3,
  family = "bernoulli")


conditional_effects(modhow1,prob = .5)

spdata3<-filter(spdata2,species!="impatiens")
pred3<-epred_draws(modhow1,newdata = spdata3,re_formula = ~(1|species))


statedata3<-filter(statedata2,state!="Kansas")
pred3a<-epred_draws(modhow1,newdata = statedata3,re_formula = ~(1|state))

pd<-position_dodge(width=0.3)
ggplot(pred3,aes(reorder(species,.epred),.epred))+
  stat_pointinterval(.width = .5)+#ylim(0,1)+
  ggthemes::theme_few()+coord_cartesian(ylim=c(.1,.9))+
  ylab( "<-eariler-later->")+xlab("species")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



ggplot(pred3a,aes(reorder(state,.epred),.epred))+facet_wrap(~state,scale="free_x",nrow=3)+
  stat_pointinterval(aes(color=.category,shape=.category),position=pd,.width=c(0.01,0.5),alpha=0.8)+
  ggthemes::theme_few()+ylab("mangement shift")+xlab("state")+scale_color_viridis_d(option="H")+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())


#####playing with categories

annualBiennial<-c("Microstegium","Centaurea", "Pastinaca","impatiens")
Vine<-c("Pueria","Celastrus","Lonicera japonica")
perennial<-c("Retnoutria","Artemesia","Miscanthus","Phragmites","Vinetoxicum")
shrub<-c("Berberis","Rhamnus","Elaeagnus","Lonicera","Rosa","Ligustrum","Cytisus")
tree<-c("Acer","Ailanthis")

preds2$woody<-NA
preds2$woody<-ifelse(preds2$species %in% c(shrub,tree,Vine),"woody","non-woody")

preds2$LHC<-NA
preds2$LHC<-ifelse(preds2$species %in% c(annualBiennial,perennial),"herbacious","tree")
preds2$LHC<-ifelse(preds2$species %in% c(Vine),"vine",preds2$LHC)

preds2$LHC<-ifelse(preds2$species %in% c(shrub),"shrub",preds2$LHC) 

ggplot(preds2,aes(reorder(woody,.epred),.epred))+
  stat_pointinterval(.width = c(.01,.5),aes(color=dataset,shape=dataset))+#ylim(0,1)+
  ggthemes::theme_few()+coord_cartesian(ylim=c(0,1))+
  ylab( "Likelihood adjusting management")+xlab("functional type")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+scale_color_manual(values=c("firebrick4","seagreen3"))

preds$woody<-NA
preds$woody<-ifelse(preds$species %in% c(shrub,tree,Vine),"woody","non-woody")

preds$LHC<-NA
preds$LHC<-ifelse(preds$species %in% c(annualBiennial,perennial),"herbacious","tree")
preds$LHC<-ifelse(preds$species %in% c(Vine),"vine",preds$LHC)
preds$LHC<-ifelse(preds$species %in% c(shrub),"shrub",preds$LHC) 

ggplot(preds,aes(reorder(LHC,.epred),.epred))+
  stat_pointinterval(.width = .5)+#ylim(0,1)+
  ggthemes::theme_few()+coord_cartesian(ylim=c(.1,.9))+
  ylab( "Likelihood of observing \nphenological shifts")+xlab("species")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

