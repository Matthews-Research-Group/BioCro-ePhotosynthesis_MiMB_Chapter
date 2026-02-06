library(ggplot2)
control = read.csv("output.data_Q2000_Control")
xexpress = read.csv("output.data_Q2000_Xexp")

colnames(control) = c("Q","T","Ci","An")
colnames(xexpress) = c("Q","T","Ci","An")

control$case = "control"
xexpress$case = "overexpress"

(xexpress$An_ePhoto - control$An_ePhoto)/control$An_ePhoto *100

df = rbind(control,xexpress)
myplot<-ggplot(df, aes(x = Ci, y = An_ePhoto, color = case)) +
  geom_line() + geom_point(size=3) +
  scale_color_manual(values=c("blue","green"))+
  labs(x = bquote(Ci~(ppm)),
       y = bquote(An~(mu*mol ~ m^-2 ~ s^-1))
  ) +
  theme_bw() +
  theme(text = element_text(size = 20),legend.position = c(0.8, 0.2))
# Save ggplot to PDF
ggsave("A_Ci_overexpress.png", plot = myplot, width = 6, height = 6)
