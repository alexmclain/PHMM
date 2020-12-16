
library(tidyverse)

date <- "Dec 2020 multi"  # To label files
plot_data <- Estimation$pred_data

### Back transforming the estimates and intervals
gamma <- exp(Estimation$model$modelStruct$varStruct[1])
sigma2 <- Estimation$model$sigma^2
plot_data$SE_var <- sqrt(sigma2*(1 + gamma*plot_data$SE_var^2))
plot_data$Point.Estimate <- exp(plot_data$Y)/(1+exp(plot_data$Y))
plot_data$Std.Err <- plot_data$Point.Estimate*(1-plot_data$Point.Estimate)*(plot_data$SE_var)
plot_data$Prediction <- exp(plot_data$pred)/(1+exp(plot_data$pred))
plot_data$SE_mean_pred <- plot_data$Prediction*(1-plot_data$Prediction)*(plot_data$sigma_T_est)
plot_data$SE_pred <- plot_data$Prediction*(1-plot_data$Prediction)*(plot_data$sigma_Y_est)
plot_data$lower_CI <- exp(plot_data$lower_CI)/(1+exp(plot_data$lower_CI))
plot_data$upper_CI <- exp(plot_data$upper_CI)/(1+exp(plot_data$upper_CI))
plot_data$lower_PI <- exp(plot_data$lower_PI)/(1+exp(plot_data$lower_PI))
plot_data$upper_PI <- exp(plot_data$upper_PI)/(1+exp(plot_data$upper_PI))

plot_data <- plot_data %>% select(-c("Y","pred","sigma_T_est",
                                     "sigma_Y_est","SE_var","pred_fixed","pred_fixpen")) %>% 
  select(c(country, year, Point.Estimate,    Std.Err, Prediction, SE_mean_pred,    SE_pred,  
           lower_CI,  upper_CI),everything())

##### plot_data can be exported, it contains all the predictions, confidence intervals
##### and prediction intervals
write.csv(plot_data,paste("Table of results ",date,".csv",sep=""),row.names = FALSE)






################################# Plotting the results #################################

# Reformatting the data
P_plot_data <- plot_data
P_plot_data <- rename(P_plot_data,"Y"="Point.Estimate","SE_var"="Std.Err","pred"="Prediction","lower_CI2"="lower_PI","upper_CI2"="upper_PI")

# Getting the region values
reg_vals <- sort(as.character(unique(all_data$Region)))

# Getting the maximum upper limit
Lim_U <- max(P_plot_data$upper_CI2)



### The following will create a pdf file of the predictions with confidence intervals (in 
#   blue) and prediction intervals (in green). Each page will contain a different region.

pdf(paste("Prevalence estimates ",date," .pdf",sep = ""),width = 15, height = 8)
for(k in reg_vals){
  reg_nice <- reg_vals[k==reg_vals]
  t_country <- unique(all_data$country[all_data$Region==k])
  plot_data <- P_plot_data[P_plot_data$country %in% t_country,]
  
  ymn <- plot_data$Y-2*(plot_data$SE_var)
  ymn[ymn<0] <- 0
  limits <- aes(y = Y, x = year, ymax = Y+2*(SE_var), ymin=ymn) 
  
  p <- ggplot(data=plot_data, aes(x=year,y=(pred))) + geom_line() + facet_wrap(~country,scales="fixed")
  
  p <- p +  labs(x="Year", y=paste("Prevalence"),title = paste("Prevalence estimates for",reg_nice),size=60) + theme_bw() + theme(axis.text=element_text(family = "Helvetica", color="#666666",size=10), axis.title = element_text(family = "Helvetica", color="#666666", face="bold", size=22)) 
  p <- p+geom_ribbon(data = plot_data, aes(ymin=lower_CI, ymax=upper_CI),fill="blue", linetype=2, alpha=0.2) + geom_line(data=plot_data,aes(x=year,y=(pred)),color="blue")  + geom_line(data=plot_data,aes(x=year,y=lower_CI),color="blue",linetype = 2)+ geom_line(data=plot_data,aes(x=year,y=upper_CI),color="blue",linetype = 2) + geom_line(data=plot_data,aes(x=year,y=lower_CI2),color="green",linetype = 2)+ geom_line(data=plot_data,aes(x=year,y=upper_CI2),color="green",linetype = 2) +  geom_pointrange(limits,fatten=0.5,size = 0.7)
  
  print(p) # This will give a Warning message "Removed XXX rows contiaining...." ignore this.
}
dev.off()





