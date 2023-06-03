



r.d = rbindlist(r.all, use.names = T,fill=T)

  
  getwd()
  file.exist()
  r = readRDS('rslt_1.rds')
  a = r %>% tibble() %>% unnest_wider('.',names_repair = 'unique')# unnest_auto(col='unseensmpl_part')
  a =  tibble(r=r) %>% hoist('r' , metric =c('max_criteria','metric' )) 
  
  a =  tibble(r=r, .name_repair='unique')  %>% hoist('r', cm_rates='Rate',.remove=T )  %>% 
    unnest_longer(cm_rates)
  
  
  r$unseensmpl_stats
  r$unseensmpl_part
  
a =  tibble(r=r, .name_repair='unique')  %>% hoist('r', mtr=c('metric','max f1'),.remove=T )  
a =  tibble(r=r, .name_repair='unique')  %>% hoist('r', mtr=c('cm'),.remove=T ) %>% 
  unnest_auto(mtr)  #values_to=c('metric','th','value','idx'))
                                                   
a =  tibble(r=r, .name_repair='unique')  %>% hoist('r', usp=c('unseensample_part','aucs',.remove=T )  

  
a =  tibble(r=r, .name_repair='unique')  %>% unnest_auto(col='r')
setdiff(names(a),names(r))
