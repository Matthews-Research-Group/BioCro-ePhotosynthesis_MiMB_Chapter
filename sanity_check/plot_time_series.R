lines <- readLines("log")
num_lines <- lines[seq(3, length(lines), by = 3)]
mat <- do.call(
  rbind,
  strsplit(num_lines, "\\s+")
)
meta_time_series <- apply(mat, 2, as.numeric)
meta_names = read.table("metabolite_names")

meta_index = meta_names$V1 + 1 # c++ starts with 0. So +1 here for R indexing
meta_time_series_sub1 = meta_time_series[,meta_index]
time = meta_time_series[,100]

assim = meta_time_series[,7]
# plot(time,assim)
colnames(meta_time_series_sub1) = meta_names$V2

pdf("figs/metabolite_time_series.pdf",width=8,height=8)
for (i in 1:length(meta_index)){
  par(cex=1.5)
  ylimit = range(c(meta_time_series_sub1[,i]))
  plot(time,meta_time_series_sub1[,i],
       xlab="time (s)",ylab=paste0(meta_names$V2[i],"(mmol/L)"),
                                   ylim=ylimit,type='l',lwd=2)
  grid (NULL,NULL, lty = 6, col = "cornsilk2")
}
dev.off()

#check the last 100 iterations


var_list = c("PGCA", "ATP", "PGA","HPRc","GCA", "GCEA", "GCAc","GCEAc")
tail(meta_time_series_sub[,var_list])


write.csv(meta_time_series_sub,file = "metabolite_time_series_Q1500T20Ci300_CTL.csv",row.names = FALSE)

rates_time_series = read.table("vRates_time_series_Q1500T20Ci300_CTL",header = TRUE)
rates_names = colnames(rates_time_series)
rates_names = rates_names[-1]
timestamp = rates_time_series[,1]
pdf("vRates_time_series_Q1500T20Ci300_CTL.pdf",width=8,height=8)
for (i in 1:length(rates_names)){
  par(cex=1.0)
  plot(timestamp,rates_time_series[,i+1],
       xlab="time (s)",ylab=rates_names[i])
  grid (NULL,NULL, lty = 6, col = "cornsilk2")
}
dev.off()

plot(timestamp,rates_time_series$v123 - rates_time_series$v1in,
     xlim=c(1000,5000),
     xlab="time (s)",ylab="GCAc")
