library(cosinor)
#data should not have any factor

fit_cosinor_cosinor <- function(data, variable='time', observation='Y',replicated_obs=F, group= NULL, period= 24,
                                color_use= dittoSeq::dittoColors(),jitter_width=0,stat_summary=T,two_axis=T){
  variable=sym(variable) #The idiomatic way now would be to convert to a symbol the string that the variable contains, using sym()
  observation_sym=sym(observation)
  
  data_avg= data
  data_avg$mean= data_avg[,observation]
  if(is.null(group)){
    if(all(data[,observation]==0)){ #if all==0, returns an empty plot
      p1= ggplot() + theme_void()
      scale_factor=1
      fit=NULL
    }else{
      if(replicated_obs){
        data_avg= data %>% group_by(!!variable) %>% summarise_at(vars(!!observation_sym), list(mean = mean)) %>% as.data.frame() # nolint
      }
      if(two_axis){ #the second axis is for the cosine curve: magnify the cosine data and devide the axis
        data_avg_max= data %>% group_by(!!variable) %>% summarise_at(vars(!!observation_sym), list(mean = mean)) %>% 
          as.data.frame() %>% pull(mean) %>% max()
        data_max= data %>% summarise_at(vars(!!observation_sym), list(max = max)) %>% pull(max)
        scale_factor= data_max %/% data_avg_max /2
        scale_factor= ifelse(scale_factor<1,1,scale_factor)  #if scale_factor<1, there is no point rescale the y axis
        data_avg$mean= data_avg$mean * scale_factor
      }
      cosinor_formula= paste0('mean ~ time(',variable,')') %>% as.formula()
      fit <- cosinor::cosinor.lm(cosinor_formula, data = data_avg, period = period)
      p1= cosinor::ggplot_cosinor.lm(fit)+geom_jitter(data=data, aes(!!variable,!!observation_sym, colour= factor(!!variable)),width= jitter_width)+
        theme(legend.position = 'none')+scale_color_manual(values= color_use) #!!: unquote  
      
      tmp= summary(fit)
      p1$data$y.upper= tmp$transformed.table$upper.CI[2]/tmp$transformed.table$estimate[2] *p1$data$Y.hat
      p1$data$y.lower= tmp$transformed.table$lower.CI[2]/tmp$transformed.table$estimate[2] *p1$data$Y.hat
      p1= p1+geom_ribbon(aes(ymin = y.lower, ymax = y.upper), 
                         alpha=0.1, 
                         linetype="dashed",
                         color="grey")
    }       
  }else{
    if(all(data[,observation]==0)){ #if all==0, returns an empty plot
      p1= ggplot() + theme_void()
      scale_factor=1
      fit=NULL
    }else{
      if(replicated_obs){
        data_avg= data %>% group_by(group,!!variable) %>% summarise_at(vars(!!observation_sym), list(mean = mean))
      }
      if(two_axis){ #the second axis is for the cosine curve: magnify the cosine data and devide the axis
        data_avg_max= data %>% group_by(!!variable) %>% summarise_at(vars(!!observation_sym), list(mean = mean)) %>% 
          as.data.frame() %>% pull(mean) %>% max()
        data_max= data %>% summarise_at(vars(!!observation_sym), list(max = max))%>% pull(max)
        scale_factor= data_max %/% data_avg_max
        scale_factor= ifelse(scale_factor<1,1,scale_factor)  #if scale_factor<1, there is no point rescale the y axis
        data_avg$mean= data_avg$mean * scale_factor
      }
      cosinor_formula= paste0('mean ~ time(',variable,') + ',group,' + amp.acro(',group,')') %>% as.formula()
      fit <- cosinor::cosinor.lm(cosinor_formula, data = data_avg, period = period)
      p1= ggplot.cosinor.lm(fit, x_str = group)+geom_jitter(data=data, aes(!!variable,!!observation_sym, colour= factor(!!variable),shape=factor(!!sym(group))),width= jitter_width)+
        theme(legend.position = 'none')+scale_color_manual(values= color_use)
      
      tmp= summary(fit)
      p1$data$y.upper= tmp$transformed.table$upper.CI[2]/tmp$transformed.table$estimate[2] *p1$data$Y.hat
      p1$data$y.lower= tmp$transformed.table$lower.CI[2]/tmp$transformed.table$estimate[2] *p1$data$Y.hat
      p1= p1+geom_ribbon(aes(ymin = y.lower, ymax = y.upper), 
                         alpha=0.1, 
                         linetype="dashed",
                         color="grey")
    }       
  }
  
  if(two_axis){ #transform the second y-axis
    data[,observation]= data[,observation] * scale_factor #this has the last step involving data, exept for stat_summary. 
    p1= p1+scale_y_continuous("Norm_expr",sec.axis = sec_axis(~ . /scale_factor, name = "Avg_expr"))
  }else{p1= p1+ylab('Average expression')}
  
  if(stat_summary){
    p1= p1+ stat_summary(mapping=aes(!!variable,!!observation_sym),data=data, fun = mean,geom='point', size = 8, colour = "black", shape = 95,inherit.aes=F)
  }
  
  p1= p1+ theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5,size = 12)) + xlab("Zeitgeber time") +
    theme(text = element_text(size=12))+ scale_x_discrete(limits=c(1,7,13,19))+ 
    theme(legend.position= "none")+guides(colour=FALSE) +theme(aspect.ratio=1)
  
  return(list(p= p1, fit= fit))
}



# seems for psych::cosinor.plot, the fitting is not as good as cosinor::cosinor.lm() fit. Therefore, I will use cosinor for the fitting.


fit_cosinor_cosinor_no_jitter <- function(data, variable='time', observation='Y',replicated_obs=F, group= NULL, period= 24,bar=T,errorbar=T, CI=T,
                                color_use= dittoSeq::dittoColors(),stat_summary=T,two_axis=T,size_ratio=1){ #size ratio specify the height/width ratio of figure
  variable=sym(variable) #The idiomatic way now would be to convert to a symbol the string that the variable contains, using sym()
  observation_sym=sym(observation)
  
  data_avg= data
  data_avg$mean= data_avg[,observation]
  if(is.null(group)){
    if(all(data[,observation]==0)){ #if all==0, returns an empty plot
      p1= ggplot() + theme_void()
      scale_factor=1
      fit=NULL
    }else{
      if(replicated_obs){
        data_avg= data %>% group_by(!!variable) %>% summarise_at(vars(!!observation_sym), list(mean = mean)) %>% as.data.frame() # nolint
      }
      if(two_axis){ #the second axis is for the cosine curve: magnify the cosine data and devide the axis
        data_avg_max= data %>% group_by(!!variable) %>% summarise_at(vars(!!observation_sym), list(mean = mean)) %>% 
          as.data.frame() %>% pull(mean) %>% max()
        data_max= data %>% summarise_at(vars(!!observation_sym), list(max = max)) %>% pull(max)
        scale_factor= data_max %/% data_avg_max /2
        scale_factor= ifelse(scale_factor<1,1,scale_factor)  #if scale_factor<1, there is no point rescale the y axis
        data_avg$mean= data_avg$mean * scale_factor
      }
      cosinor_formula= paste0('mean ~ time(',variable,')') %>% as.formula()
      fit <- cosinor::cosinor.lm(cosinor_formula, data = data_avg, period = period)
      p1= cosinor::ggplot_cosinor.lm(fit)+theme(aspect.ratio = size_ratio)+
        theme(legend.position = 'none')#+scale_color_manual(values= color_use) #!!: unquote  
      
      tmp= summary(fit)
      p1$data$y.upper= tmp$transformed.table$upper.CI[2]/tmp$transformed.table$estimate[2] *p1$data$Y.hat
      p1$data$y.lower= tmp$transformed.table$lower.CI[2]/tmp$transformed.table$estimate[2] *p1$data$Y.hat
    }       
  }else{
    if(all(data[,observation]==0)){ #if all==0, returns an empty plot
      p1= ggplot() + theme_void()
      scale_factor=1
      fit=NULL
    }else{
      if(replicated_obs){
        data_avg= data %>% group_by(group,!!variable) %>% summarise_at(vars(!!observation_sym), list(mean = mean))
      }
      if(two_axis){ #the second axis is for the cosine curve: magnify the cosine data and devide the axis
        data_avg_max= data %>% group_by(!!variable) %>% summarise_at(vars(!!observation_sym), list(mean = mean)) %>% 
          as.data.frame() %>% pull(mean) %>% max()
        data_max= data %>% summarise_at(vars(!!observation_sym), list(max = max))%>% pull(max)
        scale_factor= data_max %/% data_avg_max
        scale_factor= ifelse(scale_factor<1,1,scale_factor)  #if scale_factor<1, there is no point rescale the y axis
        data_avg$mean= data_avg$mean * scale_factor
      }
      cosinor_formula= paste0('mean ~ time(',variable,') + ',group,' + amp.acro(',group,')') %>% as.formula()
      fit <- cosinor::cosinor.lm(cosinor_formula, data = data_avg, period = period)
      p1= ggplot.cosinor.lm(fit, x_str = group)+theme(aspect.ratio = size_ratio)+
        theme(legend.position = 'none')+scale_color_manual(values= color_use)
      
      tmp= summary(fit)
      p1$data$y.upper= tmp$transformed.table$upper.CI[2]/tmp$transformed.table$estimate[2] *p1$data$Y.hat
      p1$data$y.lower= tmp$transformed.table$lower.CI[2]/tmp$transformed.table$estimate[2] *p1$data$Y.hat
    }       
  }
  if(bar){
    data_avg= data %>% group_by(!!variable) %>% summarise_at(vars(!!observation_sym), list(mean = mean))
    p1=p1+geom_bar(data=data_avg, aes(!!variable,mean, fill=factor(!!variable)),stat = 'identity',width = 2,alpha=0.8)+
      theme(legend.position = 'none')+scale_fill_manual(values= color_use)
  }
  
  if(errorbar){
    data_avg= data %>% group_by(!!variable) %>% summarise_at(vars(!!observation_sym), list(mean = mean, sd=sd))
    p1= p1+ stat_summary(mapping=aes(!!variable,!!observation_sym,colour = factor(!!variable)),data=data, fun.data = mean_se,geom='pointrange', size = 1, 
                         inherit.aes=F)+ scale_color_manual(values = color_use)
  }
  
  
  if(two_axis){ #transform the second y-axis
    data[,observation]= data[,observation] * scale_factor #this has the last step involving data, exept for stat_summary. 
    p1= p1+scale_y_continuous("Norm_expr",sec.axis = sec_axis(~ . /scale_factor, name = "Average expression"))
  }else{p1= p1+ylab('Average expression')}
  
  if(stat_summary){
    p1= p1+ stat_summary(mapping=aes(!!variable,!!observation_sym),data=data, fun = mean,geom='point', size = 8, colour = "black", shape = 95,inherit.aes=F)
  }
  
 
  if(CI){
    p1= p1+geom_ribbon(aes(ymin = y.lower, ymax = y.upper), 
                       alpha=0.1, 
                       linetype="dashed",
                       color="grey")
  }
  
  p1= p1+ theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5,size = 12)) + xlab("Zeitgeber time") +
    theme(text = element_text(size=12))+ scale_x_discrete(limits=c(1,7,13,19))+ 
    theme(legend.position= "none")+guides(colour=FALSE) +theme(aspect.ratio=1)
  
  return(list(p= p1, fit= fit))
}



# seems for psych::cosinor.plot, the fitting is not as good as cosinor::cosinor.lm() fit. Therefore, I will use cosinor for the fitting.


