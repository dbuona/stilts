rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())
graphics.off()
library(ggplot2)


setwd("~/Documents/git/stilts/")

d<-read.csv("data/phenology survey 2025.xlsx - full_datasheet.csv")


p1<-d %>% dplyr::select(Q26,Q27_1)
p2<-d %>% dplyr::select(Q137,Q138_1)
p3<-d %>% dplyr::select(Q156,Q157_1)
p4<-d %>% dplyr::select(Q175,Q176_1)
p5<-d %>% dplyr::select(Q194,Q195_1)

colnames(p1)<-c("goal","efficacy")
colnames(p2)<-c("goal","efficacy")
colnames(p3)<-c("goal","efficacy")
colnames(p4)<-c("goal","efficacy")
colnames(p5)<-c("goal","efficacy")

p1<- p1 %>%filter(efficacy!="")
p2<- p2 %>%filter(efficacy!="")
p3<- p3 %>%filter(efficacy!="")
p4<- p4 %>%filter(efficacy!="")
p5<- p5 %>%filter(efficacy!="")

p1$sp<-"s1"
p2$sp<-"s2"
p3$sp<-"s3"
p4$sp<-"s4"
p5$sp<-"s5"



allps<-rbind(p1,p2,p3,p4,p5)
allps<-filter(allps,goal!="Other")
library(ggplot2)
order<-c("Not effective at all","Slightly effective",
"Moderately effective","Very effective","Extremely effective")

allps$efficacy <- factor(allps$efficacy, levels = order)



ggplot(allps,aes(goal))+geom_bar()



ggplot(allps,aes(goal))+geom_bar(aes(fill=efficacy),position="fill")+
  scale_y_continuous(labels = scales::percent_format()) +scale_fill_brewer(palette = "Greens")+
  ggthemes::theme_few(base_size = 18)+ylab("")

library(dplyr)
library(forcats)

# Calculate proportion of Extremely + Very effective
goal_order <- allps %>%
  mutate(high_eff = efficacy %in% c("Extremely effective", "Very effective")) %>%
  group_by(goal) %>%
  summarise(prop_high = mean(high_eff, na.rm = TRUE)) %>%
  arrange(desc(prop_high))

# Apply ordering
allps$goal <- factor(allps$goal, levels = goal_order$goal)

# Plot
jpeg("updated_goals.jpeg",width = 12,height=6, units = "in",res=200)
ggplot(allps, aes(goal)) +
  geom_bar(aes(fill = efficacy), position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_brewer(palette = "Greens") +xlab("Management goal")+
  ggthemes::theme_few(base_size = 16) +
  ylab("")+ theme(legend.position = "top",legend.title = element_blank())
dev.off()
