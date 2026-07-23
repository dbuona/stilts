rm(list=ls()) #removes everything in R environment
options(stringsAsFactors = FALSE) #stops R from automatically treating string as factors
options(mc.cores = parallel::detectCores()) # for Bayesian models, paralleled chain sampling
graphics.off() # removes anything in graphic environemnt
library(dplyr)
library(ggplot2)
library(brms)
library(tidybayes)

if(length(grep("danielbuonaiuto", getwd()) > 0)) { # a fancy way to set your working direcyotu
  setwd("~/Documents/git/stilts/")
} else if(length(grep("fitz", getwd()) > 0)) {
  setwd("")}

time<-read.csv("Input/timechange.csv")
freq<-read.csv("Input/freqchange.csv")
colnames(time)
colnames(freq)


###reclassify
unique(time$timechange)
###
unique(time$timechange)

###remove NA's frome response
time<-filter(time,!is.na(timechange))
time$minCareer<-ifelse(time$career=="10+ years",10,0.1)
time$minCareer<-ifelse(time$career=="2 to 5 years",2,time$minCareer)
time$minCareer<-ifelse(time$career=="6 to 10 years",6,time$minCareer)

###remove NA's frome response
freq<-filter(freq,!is.na(freqchange))
freq$minCareer<-ifelse(freq$career=="10+ years",10,0.1)
freq$minCareer<-ifelse(freq$career=="2 to 5 years",2,freq$minCareer)
freq$minCareer<-ifelse(freq$career=="6 to 10 years",6,freq$minCareer)

time$idk<-ifelse(time$timechange%in% c("I don't know","Other,I don't know", "Other", "NA"),0,1)
freq$idk<-ifelse(freq$freqchange%in% c("I don't know","Other,I don't know", "Other", "NA", "I don’t know",""),0,1)

### get rid of idk for now
timer<-filter(time,idk==1)
freqer<-filter(freq,idk==1)

unique(timer$timechange)
noCCt<-c("At the same time due to no change in phenology","Earlier or later for other reasons",
        "Earlier or later for other reasons,At the same time due to no change in phenology",
        "Starting earlier or ending later for other reasons",
        "Starting earlier or ending later for other reasons,At the same time due to no change in phenology" )

unique(freqer$freqchange)
noCCf<-c("At the same frequency due to no change in phenology",
         "More or less frequently for other reasons")

timer$CC<-ifelse(timer$timechange %in% noCCt,0,1)
freqer$CC<-ifelse(freqer$freqchange %in% noCCf,0,1)

timer$duration.cent<-scale(timer$dur,center = TRUE)
timer$career.cent<-scale(timer$minCareer,center = TRUE)

timer<-filter(timer,!state %in% c("","Wyoming","Alaska","Nova Scotia"))

freqer$duration.cent<-scale(freqer$dur,center = TRUE)
freqer$career.cent<-scale(freqer$minCareer,center = TRUE)

freqer<-filter(freqer,!state %in% c("","Wyoming","Alaska","Nova Scotia"))

###
timerCC<-filter(timer,CC==1)
freqerCC<-filter(freqer,CC==1)

unique(timerCC$timechange)
NoChangersT<-c("At the same time despite changing phenology",
               "At the same time due to no change in phenology,At the same time despite changing phenology",
              "Starting earlier or ending later for other reasons,At the same time despite changing phenology")

unique(freqerCC$freqchange)
NoChangersF<-c("At the same frequency despite a change in phenology",
               "At the same frequency due to no change in phenology,At the same frequency despite a change in phenology",
               "More or less frequently for other reasons,At the same frequency despite a change in phenology")

timerCC$management_shift<-ifelse(timerCC$timechange %in% NoChangersT,0,1)
freqerCC$management_shift<-ifelse(freqerCC$freqchange %in% NoChangersF,0,1)



timerDir<-filter(timerCC,management_shift==1)
freqerDir<-filter(freqerCC,management_shift==1)


unique(timerDir$timechange)

earlier<-c( "Earlier due to changing phenology","Treatment starting earlier due to changing phenology",
            "Treatment starting earlier due to changing phenology,Starting earlier or ending later for other reasons",
            "Treatment starting earlier due to changing phenology,At the same time due to no change in phenology")

later<-c("Later due to changing phenology","Treatment ending later due to changing phenology",
         "Treatment ending later due to changing phenology,Starting earlier or ending later for other reasons")

earlLate<-c("Earlier due to changing phenology,Later due to changing phenology","Earlier due to changing phenology,Later due to changing phenology,At the same time despite changing phenology",
            "Treatment starting earlier due to changing phenology,Treatment ending later due to changing phenology",
            "Treatment starting earlier due to changing phenology,Treatment ending later due to changing phenology,Starting earlier or ending later for other reasons"
)
timerDir$howshift<-NA
timerDir$howshift<-ifelse(timerDir$timechange %in% earlier,"earlier","later")
timerDir$howshift<-ifelse(timerDir$timechange %in% earlLate,"earlier & later",timerDir$howshift)
timerDir$howshift <- relevel(as.factor(timerDir$howshift), ref = "earlier & later")

unique(freqerDir$freqchange)

freqerDir<-filter(freqerDir,freqchange!="More frequently due to changing phenology,Less frequently due to changing phenology")

more<-c("More frequently due to changing phenology","More frequently due to changing phenology,At the same frequency despite a change in phenology" )

freqerDir$howshift<-ifelse(freqerDir$freqchange %in% more,"more","less")

###write out cleandatasheet for modeling

write.csv(freqer,"Input/frequency_clean.csv",row.names = FALSE)
write.csv(freqerCC,"Input/frequency_cleanCC.csv",row.names = FALSE)
write.csv(freqerDir,"Input/frequency_cleanDir.csv",row.names = FALSE)


write.csv(timer,"Input/timing_clean.csv",row.names = FALSE)
write.csv(timerCC,"Input/timing_cleanCC.csv",row.names = FALSE)
write.csv(timerDir,"Input/timing_cleanDir.csv",row.names = FALSE)
