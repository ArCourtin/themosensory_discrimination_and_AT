model_data <-
  read_csv("recovery_analysis/mixed_model_data.csv") %>%
  mutate(
    relative_adapting_temperature =
      absolute_adapting_temperature - recorded_baseline_temperature,
    adapting_temperature_idx = 3 + relative_adapting_temperature,
    dataset=dataset
  )

model_data %>% 
  filter(dataset==3) %>% 
  group_by(participant,absolute_adapting_temperature,absolute_target_temperature) %>% 
  summarize(y=sum(choice_accuracy),n=sum(!is.na(choice_accuracy))) %>% 
  mutate(cy=cumsum(y),cn=cumsum(n),ratio=cy/cn) %>% 
  ggplot(aes(x=absolute_target_temperature-absolute_adapting_temperature,y=ratio,color=as.factor(absolute_adapting_temperature)))+
  geom_line()+
  facet_wrap(participant~.)+
  theme_classic()+
  guides(color='none')

model_data %>% 
  filter(dataset==3) %>% 
  ggplot()+
  geom_point(aes(x=trial,y=absolute_target_temperature-absolute_adapting_temperature,color=as.factor(choice_accuracy)))+
  facet_grid(rows=vars(absolute_adapting_temperature),cols=vars(participant))+
  guides(color='none')+
  theme_classic()

