# Script used to generate different model free plots of the discrimination and rating datasets
# Author: Arthur S. Courtin
# License: MIT (see LICENSE file)
# Edited with the assistance of Claude Code (Anthropic).

#### Set-up environment ####
library(tidyverse)

#### Extract discrimination data ####
data <-
  read_csv("data/d_at_2ifc.csv")%>%
  mutate(
    relative_adapting_temperature =
      round(adapting - baseline),
    adapting_temperature_idx = 3 + relative_adapting_temperature
  )
participant<-unique(data$participant)
P<-length(participant)
for(pdx in 1:P){
  data$participant[data$participant==participant[pdx]]<-pdx
}

#Group level
data %>% 
  group_by(temperature,adapting_temperature_idx,task) %>%
  summarise(m=mean(accuracy),n=sum(!is.na(accuracy))) %>% 
  ggplot(aes(x=temperature,y=m,alpha=n/250))+
    geom_point()+
    theme_classic()+
    facet_grid(cols=vars(adapting_temperature_idx),rows=vars(task))

#Participant level
data %>% 
  filter(participant<10) %>% 
  group_by(temperature,adapting_temperature_idx,task,participant) %>%
  summarise(m=mean(accuracy),n=sum(!is.na(accuracy))) %>% 
  ggplot(aes(x=temperature,y=m,alpha=n/30,color=as.factor(task)))+
  geom_point()+
  theme_classic()+
  facet_grid(cols=vars(adapting_temperature_idx),rows=vars(participant))
data %>%
  filter(participant>11) %>%
  group_by(temperature,adapting_temperature_idx,task,participant) %>%
  summarise(m=mean(accuracy),n=sum(!is.na(accuracy))) %>%
  ggplot(aes(x=temperature,y=m,alpha=n/30,color=as.factor(task)))+
  geom_point()+
  theme_classic()+
  facet_grid(cols=vars(adapting_temperature_idx),rows=vars(participant))

#### Extract rating data ####
rating_data <-
  read_csv("data/d_at_ratings.csv")%>%
  filter(baseline_flag == 0, confirmed == 1) %>%
  mutate(
    relative_adapting_temperature =
      round(adapting - baseline),
    adapting_temperature_idx = 3 + relative_adapting_temperature,
    rating = rating/100
  )
rating_participant<-unique(rating_data$participant)
rating_P<-length(rating_participant)
for(pdx in 1:rating_P){
  rating_data$participant[rating_data$participant==rating_participant[pdx]]<-pdx
}

#Group level
rating_data %>%
  group_by(temperature,adapting_temperature_idx,task) %>%
  ggplot(aes(x=temperature,y=rating,alpha=.10))+
    geom_point()+
    theme_classic()+
    facet_grid(cols=vars(adapting_temperature_idx),rows=vars(task))

#Participant level
rating_data %>%
  filter(participant<10) %>%
  group_by(temperature,adapting_temperature_idx,task,participant) %>%
  ggplot(aes(x=temperature,y=rating,alpha=.1,color=as.factor(task)))+
  geom_point()+
  theme_classic()+
  facet_grid(cols=vars(adapting_temperature_idx),rows=vars(participant))
rating_data %>%
  filter(participant>11) %>%
  group_by(temperature,adapting_temperature_idx,task,participant) %>%
  ggplot(aes(x=temperature,y=rating,alpha=.1,color=as.factor(task)))+
  geom_point()+
  theme_classic()+
  facet_grid(cols=vars(adapting_temperature_idx),rows=vars(participant))

