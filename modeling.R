####started by Dan Feb 2026 updated June 25
###goal is to analyze the the climate change part of the RISCC phenology survey
##collaborators include Fitz and Bethany

###house keeping
rm(list=ls()) #removes everything in R environment
options(stringsAsFactors = FALSE) #stops R from automatically treating string as factors
options(mc.cores = parallel::detectCores()) # for Bayesian models, paralleled chain sampling
graphics.off() # removes anything in graphic environemnt
dev.off()
options(device = "RStudioGD")
library(dplyr)
library(ggplot2)
library(rstan)
library(brms)
library(tidybayes)

graphics.off()

if(length(grep("danielbuona", getwd()) > 0)) { # a fancy way to set your working direcyotu
  setwd("~/Documents/git/stilts/")
} else if(length(grep("fitz", getwd()) > 0)) {
  setwd("")}

timeidk<-read.csv("Input/idktime.csv")
freqidk<-read.csv("Input/idkfreq.csv")

mergeidk1<-select(timeidk,ResponseId,idk,dur,career,minCareer,state,species,agency,juris)
mergeidk2<-select(freqidk,ResponseId,idk,dur,career,minCareer,state,species,agency,juris)

mergeidk1$dataset<-"time"
mergeidk2$dataset<-"frequency"

IDK<-rbind(mergeidk1,mergeidk2)

model.idk<-brm(idk~dur+minCareer+(1|ResponseId),
              iter=4000,warmup=3000, control = list(adapt_delta=0.99), data=IDK,family = "bernoulli")

summary(model.idk)
conditional_effects(model.idk)
timer<-read.csv("Input/timing_clean.csv")
freqer<-read.csv("Input/frequency_clean.csv")
###do managers see climate change?
class(timer$CC)
table(timer$juris)


merge1<-select(timer,ResponseId,CC,dur,minCareer,state,species,agency,juris)
merge2<-select(freqer,ResponseId,CC,dur,minCareer,state,species,agency,juris)

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

#d<-filter(d,state %in% NECASC)
d$minCareer_s<-scale(d$minCareer)
d$dur_s<-scale(d$dur)
#d$greenup_s<-scale(d$slope)
#d$green_var<-scale(d$se_doy)

d$species<-ifelse(d$species=="Artemesia","Artemisia",d$species) ##fix species name
table(d$CC)


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
#model.1<-brm(CC~dur+slope+(1|species)+(1|ResponseId),
#             iter=8000,warmup=7000, control = list(adapt_delta=0.99), data=d,family = "bernoulli")

unique(d$juris)
d.prime<-filter(d, juris!=" Other (please specify)")
model.1a<-brm(CC~dur+(1|species)+(1|juris)+(1|agency)+(1|state)+(1|ResponseId),
             iter=8000,warmup=7000, control = list(adapt_delta=0.99), data=d,family = "bernoulli")

model.1ab<-brm(CC~dur+(1|species)+(1|juris)+(1|agency)+(1|ResponseId),
              iter=8000,warmup=7000, control = list(adapt_delta=0.99), data=d,family = "bernoulli")

model.1ac<-brm(CC~dur+(1|species)+(1|juris)+(1|ResponseId),
               iter=8000,warmup=7000, control = list(adapt_delta=0.99), data=d,family = "bernoulli")

model.1ad<-brm(CC~dur+(1|species)+(1|juris)+(1|state)+(1|ResponseId),
               iter=8000,warmup=7000, control = list(adapt_delta=0.99), data=d,family = "bernoulli")

model.1ae<-brm(CC~dur+(1|species)+(1|state)+(1|ResponseId),
               iter=8000,warmup=7000, control = list(adapt_delta=0.99), data=d,family = "bernoulli")

model.1af<-brm(CC~dur+(1|species)+(1|agency)+(1|ResponseId),
               iter=8000,warmup=7000, control = list(adapt_delta=0.99), data=d,family = "bernoulli")

model.1ag<-brm(CC~dur+(1|species)+(1|ResponseId),
               iter=8000,warmup=7000, control = list(adapt_delta=0.99), data=d,family = "bernoulli")

summary(model.1ag)

model.1a<-add_criterion(model.1a, "loo",moment_match = TRUE)
model.1ab<-add_criterion(model.1ab, "loo",moment_match = TRUE)
model.1ac<-add_criterion(model.1ac, "loo",moment_match = TRUE)
model.1ad<-add_criterion(model.1ad, "loo",moment_match = TRUE)
model.1ae<-add_criterion(model.1ae, "loo",moment_match = TRUE)
model.1af<-add_criterion(model.1af, "loo",moment_match = TRUE)
model.1ag<-add_criterion(model.1ag, "loo",moment_match = TRUE)
loo_compare(model.1a,model.1ab,model.1ac,model.1ad,model.1ae,model.1af,model.1ag,criterion = "loo")

### now change
timer2<-read.csv("Input/timing_cleanCC.csv")
freqer2<-read.csv("Input/frequency_cleanCC.csv")

merge11<-select(timer2,ResponseId,management_shift,dur,minCareer,state,species,agency,juris)
merge22<-select(freqer2,ResponseId,management_shift,dur,minCareer,state,species,agency,juris)

merge11$dataset<-"time"
merge22$dataset<-"frequency"
d2<-rbind(merge11,merge22)
#d2$minCareer_s<-scale(d2$minCareer)
#d2$dur_s<-scale(d2$dur)
d2$species<-ifelse(d2$species=="Artemesia","Artemisia",d2$species) ##fix species name

#d2<-filter(d2, state %in% NECASC)
#d2<-left_join(d2,gren)

#d2$minCareer_s<-scale(d2$minCareer)
#d2$dur_s<-scale(d2$dur)


model.mc1<-brm(management_shift~dur+dataset+(dataset|species)+(dataset|juris)+(1|ResponseId),
               iter=5000,warmup=4000, control = list(adapt_delta=0.99),
               data=d2,family = "bernoulli")

summary(model.mc1)


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

merge111<-select(timer3,ResponseId,dur,minCareer,state,species,agency,juris,howshift)
merge222<-select(freqer3,ResponseId,dur,minCareer,state,species,agency,juris,howshift)


timer3$species<-ifelse(timer3$species=="Artemesia","Artemisia",timer3$species) ##fix species name

#timer3<-filter(timer3,state %in% NECASC)


#timer3<-filter(timer3,howshift!="earlier & later")
#timer3$how<-ifelse(timer3$howshift=="earlier",0,1)
modhow1 <- brm(
  formula = howshift~(1|species)+(1|juris)+(1|ResponseId), iter=5000,warmup=4000, control = list(adapt_delta=0.99),
  data = timer3,
  family = "categorical")

summary(modhow1)

save.image("phenmods.Rda")
newdaterz<-data.frame(dur=10,dataset=unique(d$dataset))
newdaterz2<-data.frame(dur=seq(1,40,by=1) ,dataset=unique(d$dataset))

predscc<-epred_draws(model.1ac,newdata=newdaterz,re_formula = NA)
predscc2<-epred_draws(model.1ac,newdata=newdaterz2,re_formula = NA)

predmc<-epred_draws(model.mc1,newdata=newdaterz,re_formula = NA)
predmc2<-epred_draws(model.mc1,newdata=newdaterz2,re_formula = NA)

predhow<-epred_draws(modhow1,newdata=newdaterz,re_formula = NA)

predcc1<- predscc %>% mean_qi(.epred, .width = c(0.5))
predmc1<- predmc %>% mean_qi(.epred, .width = c(0.5))
predhow1<- predhow %>% group_by(.category) %>% mean_qi(.epred, .width = c(0.5))




#figure 1
quartz()
pdf("cclike.pdf",width = 3,height = 7)
ggplot(predcc1,aes(0,.epred*100))+
  geom_bar(stat="identity",fill="tan",color="black")+
  geom_errorbar(aes(ymin = .lower*100, 
                    ymax = .upper*100, 
                    width = 0.2))+
  ylab("probability of observing phenological shifts")+
  scale_x_discrete(name="")+ylim(0,70)+ggthemes::theme_few()
dev.off()
head(d)
durdat<-dplyr::select(d,ResponseId,dur)
durdat<-unique(durdat)
quartz()

gumby<-ggpubr::ggarrange(
ggplot(predscc2,aes(dur,.epred*100))+stat_lineribbon(.width = c(.5),fill="tan")+
  ylab("probability of observing phenological shifts")+
  ggthemes::theme_few()+xlab("years of management experience")+ylim(0,100),
ggplot(predmc2,aes(dur,.epred*100))+stat_lineribbon(aes(fill=dataset),.width = c(.5))+
  ylab("probability of adjusting management timing")+theme(legend.position = "top")+
  ggthemes::theme_few()+xlab("years of management experience")+ylim(0,100)+
  scale_fill_manual(values=c("maroon4","skyblue4")),widths=c(1,1.5),labels=c("b)","c)"),common.legend = TRUE)

jpeg("manage_duration.jpeg",width =8,height=7, unit='in',res=100)
ggpubr::ggarrange(ggplot(durdat,aes(dur))+geom_histogram(fill="lightgray",color="black",bins=20)+
                    ggthemes::theme_few()+xlab("years of management experience"),gumby,ncol=1,heights=c(.75,1.5))
dev.off()

quartz()
predsccSP<-preds %>% group_by(species) %>% mean_qi(.epred, .width = c(0.5))


ggplot(preds,aes(reorder(species,.epred),.epred*100))+
  stat_pointinterval(.width = c(0.01,.5),color="tan",shape=15)+#ylim(0,1)+
  ggthemes::theme_few()+coord_cartesian(ylim=c(0,100))+
  ylab( "probability of of observing \nphenological shifts")+xlab("species")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


spdata2<-data.frame(species=rep(unique(d2$species),2),dataset=rep(c("time","frequency"),each=21),
                    dur=10)
pd=position_dodge(width = .8)
preds2<-epred_draws(model.mc1,newdata=spdata2,re_formula = ~(1|species))

predsmcSP<-preds2 %>% group_by(species,dataset) %>% mean_qi(.epred, .width = c(0.5))

spdata3<-filter(spdata2,species!="impatiens")
pred3<-epred_draws(modhow1,newdata = spdata3,re_formula = ~(1|species))
predhow3<- pred3 %>% group_by(.category,species) %>% mean_qi(.epred, .width = c(0.5))

annualBiennial<-c("Microstegium","Centaurea", "Pastinaca","impatiens")
Vine<-c("Pueria","Celastrus","Lonicera japonica")
perennial<-c("Retnoutria","Artemesia","Miscanthus","Phragmites","Vinetoxicum")
shrub<-c("Berberis","Rhamnus","Elaeagnus","Lonicera","Rosa","Ligustrum","Cytisus")
tree<-c("Acer","Ailanthis")

pred3$woody<-NA
pred3$woody<-ifelse(pred3$species %in% c(shrub,tree,Vine),"woody","non-woody")

preds2$woody<-NA
preds2$woody<-ifelse(preds2$species %in% c(shrub,tree,Vine),"woody","non-woody")
predsmcW<-preds2 %>% group_by(woody,dataset) %>% mean_qi(.epred, .width = c(0.5))



predhowW<- pred3 %>% group_by(.category,woody) %>% mean_qi(.epred, .width = c(0.5))
# Get the ordering based on the "later" category
library(forcats)
species_order <- predhow3 %>%
  filter(.category == "later") %>%
  arrange(-.epred) %>%
  pull(species)


####Figure 3a



gumbers<-ggpubr::ggarrange(ggplot(predsmcSP,aes(reorder(species,.epred),.epred*100))+
  geom_bar(stat="identity",position="dodge",aes(fill=dataset),color="black")+
  geom_errorbar(aes(ymin = .lower*100, 
                    ymax = .upper*100, 
                    width = 0.2, 
                    group=dataset),position=pd)+
  ylab("probability of adjusting management timing")+
  scale_x_discrete(name="")+ylim(0,100)+ggthemes::theme_few()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(name="",values=c("maroon4","skyblue4"))+theme(legend.position = "top"),

predhow3 %>%
  mutate(species = factor(species, levels = species_order)) %>%
  ggplot(aes(x = species, y = .epred * 100)) +
  geom_bar(stat = "identity", aes(fill = .category),color="black") +
  ylab("direction of management shift") +
  scale_x_discrete(name = "") +
  ylim(0, 100) +
  ggthemes::theme_few() +
  scale_fill_manual(name="",
    values = c("lightskyblue3", "lightskyblue1", "royalblue3")
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "top"),ncol=2,labels=c("b)","c)"))

jpeg("species_var.jpeg",width =13,height=11, unit='in',res=100)
ggpubr::ggarrange(ggplot(predsccSP,aes(reorder(species,.epred),.epred*100))+
  geom_bar(stat="identity",fill="tan",color="black")+
  geom_errorbar(aes(ymin = .lower*100, 
                    ymax = .upper*100, 
                    width = 0.2))+
  ylab("probability of observing phenological shifts")+
  scale_x_discrete(name="")+ylim(0,100)+ggthemes::theme_few()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)),gumbers,ncol=1,labels=c("a)",""),
  heights=c(.85,1))
dev.off()


jpeg("woody.jpeg",width =5,height=9, unit='in',res=100)
ggpubr::ggarrange(
  ggplot(predsmcW,aes(woody,.epred*100))+
    geom_bar(stat="identity",position="dodge",aes(fill=dataset),color="black")+
    geom_errorbar(aes(ymin = .lower*100, 
                      ymax = .upper*100, 
                      width = 0.2, 
                      group=dataset),position=pd)+
    ylab("probability of adjusting management timing")+
    scale_x_discrete(name="")+ylim(0,100)+ggthemes::theme_few()+
    scale_fill_manual(name="",values=c("maroon4","skyblue4")),
  ggplot(predhowW,aes(x = woody, y = .epred * 100)) +
  geom_bar(stat = "identity", aes(fill = .category),color="black") +
  ylab("direction of management shift") +
  scale_x_discrete(name = "") +
  ylim(0, 100) +
  ggthemes::theme_few() +
  scale_fill_manual(name="",
                    values = c("lightskyblue3", "lightskyblue1", "royalblue3")
  ),ncol=1,heights=c(1.5,1),labels=c("a)","b)"))

dev.off()
  

preds$LHC<-NA
preds$LHC<-ifelse(preds$species %in% c(annualBiennial,perennial),"herbacious","tree")
preds$LHC<-ifelse(preds$species %in% c(Vine),"vine",preds$LHC)
preds$LHC<-ifelse(preds$species %in% c(shrub),"shrub",preds$LHC) 




predscc2$dataset<-"observe"
pred2<-rbind(predscc2,predmc2)


pdf("cclike3.pdf",width = 5,height=5)
ggplot(predhow1,aes(reorder(.category,-.epred),.epred*100))+
  geom_bar(stat="identity",aes())+
  geom_errorbar(aes(ymin = .lower*100, 
                    ymax = .upper*100, 
                    width = 0.2))+
  ylab("direction of management shift")+
  scale_x_discrete(name="")+ylim(0,100)+ggthemes::theme_few()
dev.off()

pdf("cclike2.pdf",width = 4,height=5)
ggplot(predmc1,aes(dataset,.epred*100))+
  geom_bar(stat="identity",aes(fill=dataset))+
  geom_errorbar(aes(ymin = .lower*100, 
                    ymax = .upper*100, 
                width = 0.2))+
  ylab("probability of adjusting management timing")+
  ylim(0,100)+ggthemes::theme_few()+scale_fill_manual(values=c("maroon4","skyblue4"))+
  theme(legend.position ="none")
dev.off()

quartz()
ggpubr::ggarrange(pa,pb,pc,ncol=3,widths=c(1,2,2))
ggplot(predz,aes(response,.epred))+stat_pointinterval(aes(color=dataset))

quartz()
conditional_effects(model.1a,prob = .5)
summary(model.1a,probs = 0.5)
#d.site<-filter(d, juris=="Single site")

#model.1a.site<-brm(CC~dur_s+minCareer_s+(1|species)+(1|state)+(1|agency)+(1|ResponseId),
#              iter=8000,warmup=7000, control = list(adapt_delta=0.99), data=d.site,family = "bernoulli")



#range(d$anomoly,na.rm=TRUE)
#model.1b<-brm(CC~dur_s*minCareer_sanom_s+(anom_s|species)+(lat_sc|agency)+(1|ResponseId),
 #             iter=8000,warmup=7000, control = list(adapt_delta=0.99), data=d,family = "bernoulli")

#fixef(model.1,probs = c(.25,.75))
#fixef(model.1a,probs = c(.25,.75))
quartz()
conditional_effects(model.1,prob = .5)

conditional_effects(model.1a,prob = .5)
fixef(model.1b,probs = c(.25,.75))

spdata<-data.frame(species=unique(d$species),dur=10)

jurisdata<-data.frame(juris=unique(d$juris),dur=10)
agedata<-data.frame(agency=unique(d$agency),dur=15)
agedata<-filter(agedata,agency!='Tribal/First Nation')###not in model

quantile(d$dur,na.rm=TRUE)
range(d$anom_s,na.rm=TRUE)



preds<-epred_draws(model.1ac,newdata=spdata,re_formula = ~(1|species))
#predage<-epred_draws(model.1a,newdata=agedata,re_formula = ~(1|agency))
predjuris<-epred_draws(model.1ac,newdata=jurisdata,re_formula = ~(1|juris))


quartz()

                        
                
                        ggplot(predjuris,aes(reorder(juris,.epred),.epred))+
                          stat_pointinterval(.width = c(.01,.5))+coord_cartesian(ylim=c(.1,.9))+
                          #ylim(0,1)
                          ggthemes::theme_few()+
                          ylab( "Likelihood of observing \nphenological shifts")+xlab("agency")+
                          theme(axis.text.x = element_text(angle = 45, hjust = 1))
                        


ggpubr::ggarrange(oneA,oneB,ncol=1,nrow=2) ###FIgure 1






statedata2<-data.frame(state=rep(unique(d2$state),2),dataset=rep(c("time","frequency"),each=24),dur_s=mean(d2$dur_s,na.rm=TRUE),minCareer_s=mean(d2$minCareer_s))

agedata2<-data.frame(juris=rep(unique(d2$juris),2),dataset=rep(c("time","frequency"),each=4),dur=mean(d2$dur,na.rm=TRUE))



predstate2<-epred_draws(model.mc1,newdata=statedata2,re_formula = ~(1|state))


predage2<-epred_draws(model.mc1,newdata=agedata2,re_formula = ~(1|juris))
#options(device = "quartz")
quartz()
ggplot(preds2,aes(reorder(species,.epred),.epred*100))+
  stat_pointinterval(.width = c(.01,.5),aes(color=dataset,shape=dataset))+#ylim(0,1)+
  ggthemes::theme_few()+coord_cartesian(ylim=c(0,100))+
  ylab( "probability of adjusting \n management timing")+xlab("species")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_color_manual(values=c("maroon4","skyblue4"))+scale_shape_manual(values=c(15,17))


#ggplot(predstate2,aes(reorder(state,.epred),.epred))+
 # stat_pointinterval(.width = c(.01,.5),aes(color=dataset,shape=dataset))+#ylim(0,1)+
#  ggthemes::theme_few()+coord_cartesian(ylim=c(0,1))+
#  ylab( "Likelihood adjusting management")+xlab("state")+
 # theme(axis.text.x = element_text(angle = 45, hjust = 1))+scale_color_manual(values=c("firebrick4","seagreen3")),

ggplot(predage2,aes(reorder(juris,.epred),.epred))+
  stat_pointinterval(.width = c(.01,.5),aes())+#ylim(0,1)+
  ggthemes::theme_few()+coord_cartesian(ylim=c(0,1))+
  ylab( "Likelihood adjusting management")+xlab("agency")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_color_manual(values=c("firebrick4","seagreen3"))

######how is management shifting


conditional_effects(modhow1,prob = .5)




pd<-position_dodge(width=0.3)
quartz()
ggplot(pred3,aes(reorder(species,.epred),.epred))+
  stat_pointinterval(.width = .5,aes(color=.category))+#ylim(0,1)+
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

