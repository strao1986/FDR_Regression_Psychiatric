library(tidyverse)
library(ggplot2)
library(patchwork)
library(grid)

make_plots = function(r1, r2, func_groups, methods,anno) {
  fdr_power_long = lapply(seq_along(r1), function(i){
    df = as.data.frame(r1[[i]])
    colnames(df) = c("FDRreg_FDR","FDRreg_Power",
                      "EN_FDR","EN_Power",
                      "LASSO_FDR","LASSO_Power",
                      "MA_FDR","MA_Power")
    df_long = df %>%
      pivot_longer(cols=everything(),
                   names_sep="_",
                   names_to=c("Method","Metric"),
                   values_to="Value") %>%
      mutate(Func=paste0("f",i),
             Group=func_groups[i])
    return(df_long)
  }) %>% bind_rows()
  
  time_long = lapply(seq_along(r2), function(i){
    df = as.data.frame(r2[[i]])
    colnames(df) = methods
    df_long = df %>%
      pivot_longer(cols=everything(),
                   names_to="Method",
                   values_to="Value") %>%
      mutate(Func=paste0("f",i),
             Group=func_groups[i],
             Metric="Runtime")
    return(df_long)
  }) %>% bind_rows()
  
  all_long = bind_rows(fdr_power_long, time_long) %>%
    mutate(Method=factor(Method,levels=c("FDRreg","EN","LASSO","MA")),
    Func = factor(Func, levels = paste0("f", seq_along(r1)))
    )
  
  g1 = ggplot(filter(all_long, Metric=="FDR"),
               aes(x=Func,y=Value,fill=Method)) +
    geom_boxplot() +
    geom_hline(yintercept=0.1,color="red",linetype="dashed") +
    facet_wrap(~Group,scales="free_x",nrow=1) +
    theme_bw(base_size=24) + xlab("Function") + ylab("FDR") +
    theme(axis.text.x=element_text(angle=45,hjust=1)) + labs(title = anno[1])
  
  g2 = ggplot(filter(all_long, Metric=="Power"),
               aes(x=Func,y=Value,fill=Method)) +
    geom_boxplot() +
    facet_wrap(~Group,scales="free_x",nrow=1) +
    theme_bw(base_size=24) + xlab("Function") + ylab("Power") +
    theme(axis.text.x=element_text(angle=45,hjust=1))+ labs(title = anno[2])
  
  g3 = ggplot(filter(all_long, Metric=="Runtime"),
               aes(x=Func,y=Value,fill=Method)) +
    geom_boxplot() +
    facet_wrap(~Group,scales="free_x",nrow=1) +
    theme_bw(base_size=24) + xlab("Function") + ylab("Computational time (s)") +
    theme(axis.text.x=element_text(angle=45,hjust=1))+ labs(title = anno[3])
  
  list(FDR=g1, Power=g2, Time=g3)
}

# dataset 1: 
load("/mnt/data/zhijie/FDR/N_5000_1_para.RData")
func_groups = c(rep("Group1",4),rep("Group2",3),
                 rep("Group3",4),rep("Group4",3),
                 rep("Group5",3),rep("Group6",3))
methods = c("FDRreg","EN","LASSO","MA")
plots1 = make_plots(r1, r2, func_groups, methods, anno = c("A","B","C"))

# dataset 2: 
load("/mnt/data/zhijie/FDR/N_5000_2_para.RData")
plots2 = make_plots(r1, r2, func_groups, methods, anno = c("D","E","F"))


combined = (plots1$FDR | plots2$FDR) /
            (plots1$Power | plots2$Power) /
            (plots1$Time | plots2$Time)

final_plot = combined + 
  plot_layout(guides = "collect") &      
  theme(
    legend.position="right",             
    legend.title=element_text(size=28),
    legend.text=element_text(size=24),
    axis.title=element_text(size=26),
    axis.text=element_text(size=22),
    strip.text=element_text(size=24),
    plot.title=element_text(size=36, face="bold", hjust=0)  
  )


final_plot

ggsave("/mnt/data/zhijie/FDR/final_plot_2601.png", final_plot, dpi=300, 
        width = 100, height = 60, units = "cm", scale = 1.6, limitsize = F)
        