env_numeric <- env[sapply(env,is.numeric)]
env_numeric <- env_numeric [,2:85]

env_numeric_long <- env_numeric %>% 
	pivot_longer(everything(), names_to = "variable", values_to = "value")

env_numeric_long %>% 
	ggplot(aes(x = value))+
	geom_histogram()+
	facet_wrap(~variable, scales = "free")
