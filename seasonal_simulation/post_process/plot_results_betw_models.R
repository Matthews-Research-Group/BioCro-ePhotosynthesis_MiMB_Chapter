library(ggplot2)
library(dplyr)
rm(list=ls())
source('merge_rds_results.R')
years      = 2002
out_key = "ePhoto"
exp_ids = c(0,1)
folder_path = "../" #this is where the rds results' folders are

result_list = list()
#ePhoto exp 0(control) and 1(overexpress)
for (i in 1:2){
  exp_id = exp_ids[i]
  result_list[[i]] <- merge_results(years,exp_id,out_key,folder_path)
}

ephoto0 = result_list[[1]][[1]]
ephoto1 = result_list[[2]][[1]]

ephoto0_sub <- ephoto0 %>%
  mutate(case = "control") %>%
  select(fractional_doy, Grain, case,lai,canopy_assimilation_rate,sunlit_Assim_layer_0)

ephoto1_sub <- ephoto1 %>%
  mutate(case = "overexpress") %>%
  select(fractional_doy, Grain, case,lai,canopy_assimilation_rate,sunlit_Assim_layer_0)

df <- bind_rows(ephoto0_sub, ephoto1_sub)

df$case <- factor(
  df$case,
  levels = c("control", "overexpress")
)

df <- df %>%
  arrange(case, fractional_doy) %>%   # critical: correct time order
  group_by(case) %>%
  mutate(canopy_assimilation_rate_cum = 
           cumsum(canopy_assimilation_rate),
         leaf_A_cum = cumsum(sunlit_Assim_layer_0)) %>%
  ungroup()

mycb = c("#000000","#009E73")
myplot<-ggplot(df, aes(x = fractional_doy, y = lai, color = case)) +
  geom_line(linewidth = 2) +
  scale_color_manual(values=mycb)+
  labs(x = "DOY",
       #y = "LAI"
        y = "Yield (t/ha)"
       # y = bquote(Cumulative~Leaf~An~(mu*mol ~ m^-2 ~ s^-1))
       # y = bquote(Cumulative~Canopy~An~(t/ha))
  ) +
  # xlim(152,155)+
  theme_bw() +
  theme(text = element_text(size = 20),legend.position = c(0.8, 0.2))
myplot
# Save ggplot to PDF
ggsave("Yield_overexpress_2002.png", plot = myplot, width = 6, height = 6)
