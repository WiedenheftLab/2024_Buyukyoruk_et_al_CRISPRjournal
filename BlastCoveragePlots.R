# Load control dataframes located in provided Database files. You might need to provide to full path name.

setwd("BLAST_coverage_data")

###IA###

df <- read.delim('IAtoall.txt', header = FALSE)

colnames(df) <- c("Comparison","Coverage")

base <- ggplot(subset(df,Comparison=="IA_to_IA"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison=="IA_to_IA"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=0.25,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p1 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

base <- ggplot(subset(df,Comparison!="IA_to_IA"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="IA_to_IA"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p2 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

base <- ggplot(subset(df,Comparison!="IA_to_IA"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="IA_to_IA"), aes(x = Coverage)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p3 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

ggarrange(p1, p2, p3, nrow=3,labels=c("Correct","Uncorrect_ind","Uncorrect_all"),common.legend = TRUE, legend="top", font.label=list(color="black",size=10))

max(subset(df,Comparison!="IA_to_IA")$Coverage)

df <- read.delim('alltoIA.txt', header = FALSE)

colnames(df) <- c("Comparison","Coverage")

base <- ggplot(subset(df,Comparison!="IA_to_IA"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="IA_to_IA"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p4 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

p4

###IB###

df <- read.delim('IBtoall.txt', header = FALSE)

colnames(df) <- c("Comparison","Coverage")

base <- ggplot(subset(df,Comparison=="IB_to_IB"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison=="IB_to_IB"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=0.25,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p1 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

base <- ggplot(subset(df,Comparison!="IB_to_IB"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="IB_to_IB"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p2 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

base <- ggplot(subset(df,Comparison!="IB_to_IB"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="IB_to_IB"), aes(x = Coverage)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p3 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

ggarrange(p1, p2, p3, nrow=3,labels=c("Correct","Uncorrect_ind","Uncorrect_all"),common.legend = TRUE, legend="top", font.label=list(color="black",size=10))

max(subset(df,Comparison!="IB_to_IB")$Coverage)

df <- read.delim('alltoIB.txt', header = FALSE)

colnames(df) <- c("Comparison","Coverage")

base <- ggplot(subset(df,Comparison!="IB_to_IB"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="IB_to_IB"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p4 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

p4

###IC###

df <- read.delim('ICtoall.txt', header = FALSE)

colnames(df) <- c("Comparison","Coverage")

base <- ggplot(subset(df,Comparison=="IC_to_IC"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison=="IC_to_IC"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=0.25,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p1 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

base <- ggplot(subset(df,Comparison!="IC_to_IC"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="IC_to_IC"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p2 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

base <- ggplot(subset(df,Comparison!="IC_to_IC"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="IC_to_IC"), aes(x = Coverage)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p3 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

ggarrange(p1, p2, p3, nrow=3,labels=c("Correct","Uncorrect_ind","Uncorrect_all"),common.legend = TRUE, legend="top", font.label=list(color="black",size=10))

max(subset(df,Comparison!="IC_to_IC")$Coverage)

df <- read.delim('alltoIC.txt', header = FALSE)

colnames(df) <- c("Comparison","Coverage")

base <- ggplot(subset(df,Comparison!="IC_to_IC"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="IC_to_IC"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p4 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

p4

###ID###

df <- read.delim('IDtoall.txt', header = FALSE)

colnames(df) <- c("Comparison","Coverage")

base <- ggplot(subset(df,Comparison=="ID_to_ID"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison=="ID_to_ID"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=0.25,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p1 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

base <- ggplot(subset(df,Comparison!="ID_to_ID"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="ID_to_ID"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p2 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

base <- ggplot(subset(df,Comparison!="ID_to_ID"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="ID_to_ID"), aes(x = Coverage)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p3 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

ggarrange(p1, p2, p3, nrow=3,labels=c("Correct","Uncorrect_ind","Uncorrect_all"),common.legend = TRUE, legend="top", font.label=list(color="black",size=10))

max(subset(df,Comparison!="ID_to_ID")$Coverage)

df <- read.delim('alltoID.txt', header = FALSE)

colnames(df) <- c("Comparison","Coverage")

base <- ggplot(subset(df,Comparison!="ID_to_ID"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="ID_to_ID"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p4 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

p4

###IE###

df <- read.delim('IEtoall.txt', header = FALSE)

colnames(df) <- c("Comparison","Coverage")

base <- ggplot(subset(df,Comparison=="IE_to_IE"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison=="IE_to_IE"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=0.25,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p1 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

base <- ggplot(subset(df,Comparison!="IE_to_IE"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="IE_to_IE"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p2 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

base <- ggplot(subset(df,Comparison!="IE_to_IE"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="IE_to_IE"), aes(x = Coverage)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p3 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

ggarrange(p1, p2, p3, nrow=3,labels=c("Correct","Uncorrect_ind","Uncorrect_all"),common.legend = TRUE, legend="top", font.label=list(color="black",size=10))

max(subset(df,Comparison!="IE_to_IE")$Coverage)

df <- read.delim('alltoIE.txt', header = FALSE)

colnames(df) <- c("Comparison","Coverage")

base <- ggplot(subset(df,Comparison!="IE_to_IE"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="IE_to_IE"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p4 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

p4

###IF###

df <- read.delim('IFtoall.txt', header = FALSE)

colnames(df) <- c("Comparison","Coverage")

base <- ggplot(subset(df,Comparison=="IF_to_IF"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison=="IF_to_IF"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=0.25,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p1 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

base <- ggplot(subset(df,Comparison!="IF_to_IF"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="IF_to_IF"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p2 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

base <- ggplot(subset(df,Comparison!="IF_to_IF"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="IF_to_IF"), aes(x = Coverage)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p3 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

ggarrange(p1, p2, p3, nrow=3,labels=c("Correct","Uncorrect_ind","Uncorrect_all"),common.legend = TRUE, legend="top", font.label=list(color="black",size=10))

max(subset(df,Comparison!="IF_to_IF")$Coverage)

df <- read.delim('alltoIF.txt', header = FALSE)

colnames(df) <- c("Comparison","Coverage")

base <- ggplot(subset(df,Comparison!="IF_to_IF"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="IF_to_IF"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p4 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

p4

###IG###

df <- read.delim('IGtoall.txt', header = FALSE)

colnames(df) <- c("Comparison","Coverage")

base <- ggplot(subset(df,Comparison=="IG_to_IG"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison=="IG_to_IG"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=0.25,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p1 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

base <- ggplot(subset(df,Comparison!="IG_to_IG"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="IG_to_IG"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p2 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

base <- ggplot(subset(df,Comparison!="IG_to_IG"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="IG_to_IG"), aes(x = Coverage)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p3 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

ggarrange(p1, p2, p3, nrow=3,labels=c("Correct","Uncorrect_ind","Uncorrect_all"),common.legend = TRUE, legend="top", font.label=list(color="black",size=10))

max(subset(df,Comparison!="IG_to_IG")$Coverage)

df <- read.delim('alltoIG.txt', header = FALSE)

colnames(df) <- c("Comparison","Coverage")

base <- ggplot(subset(df,Comparison!="IG_to_IG"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="IG_to_IG"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p4 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

p4


###IIA###

df <- read.delim('IIAtoall.txt', header = FALSE)

colnames(df) <- c("Comparison","Coverage")

base <- ggplot(subset(df,Comparison=="IIA_to_IIA"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison=="IIA_to_IIA"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=0.25,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p1 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

base <- ggplot(subset(df,Comparison!="IIA_to_IIA"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="IIA_to_IIA"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p2 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

base <- ggplot(subset(df,Comparison!="IIA_to_IIA"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="IIA_to_IIA"), aes(x = Coverage)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p3 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

ggarrange(p1, p2, p3, nrow=3,labels=c("Correct","Uncorrect_ind","Uncorrect_all"),common.legend = TRUE, legend="top", font.label=list(color="black",size=10))

max(subset(df,Comparison!="IIA_to_IIA")$Coverage)

df <- read.delim('alltoIIA.txt', header = FALSE)

colnames(df) <- c("Comparison","Coverage")

base <- ggplot(subset(df,Comparison!="IIA_to_IIA"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="IIA_to_IIA"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p4 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

p4

###IIC###

df <- read.delim('IICtoall.txt', header = FALSE)

colnames(df) <- c("Comparison","Coverage")

base <- ggplot(subset(df,Comparison=="IIC_to_IIC"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison=="IIC_to_IIC"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=0.25,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p1 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

base <- ggplot(subset(df,Comparison!="IIC_to_IIC"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="IIC_to_IIC"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p2 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

base <- ggplot(subset(df,Comparison!="IIC_to_IIC"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="IIC_to_IIC"), aes(x = Coverage)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p3 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

ggarrange(p1, p2, p3, nrow=3,labels=c("Correct","Uncorrect_ind","Uncorrect_all"),common.legend = TRUE, legend="top", font.label=list(color="black",size=10))

max(subset(df,Comparison!="IIC_to_IIC")$Coverage)

df <- read.delim('alltoIIC.txt', header = FALSE)

colnames(df) <- c("Comparison","Coverage")

base <- ggplot(subset(df,Comparison!="IIC_to_IIC"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="IIC_to_IIC"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p4 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

p4

###IIIA###

df <- read.delim('IIIAtoall.txt', header = FALSE)

colnames(df) <- c("Comparison","Coverage")

base <- ggplot(subset(df,Comparison=="IIIA_to_IIIA"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison=="IIIA_to_IIIA"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=0.25,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p1 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

base <- ggplot(subset(df,Comparison!="IIIA_to_IIIA"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="IIIA_to_IIIA"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p2 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

base <- ggplot(subset(df,Comparison!="IIIA_to_IIIA"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="IIIA_to_IIIA"), aes(x = Coverage)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p3 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

ggarrange(p1, p2, p3, nrow=3,labels=c("Correct","Uncorrect_ind","Uncorrect_all"),common.legend = TRUE, legend="top", font.label=list(color="black",size=10))

max(subset(df,Comparison!="IIIA_to_IIIA")$Coverage)

df <- read.delim('alltoIIIA.txt', header = FALSE)

colnames(df) <- c("Comparison","Coverage")

base <- ggplot(subset(df,Comparison!="IIIA_to_IIIA"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="IIIA_to_IIIA"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p4 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

p4

###IIIB###

df <- read.delim('IIIBtoall.txt', header = FALSE)

colnames(df) <- c("Comparison","Coverage")

base <- ggplot(subset(df,Comparison=="IIIB_to_IIIB"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison=="IIIB_to_IIIB"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=0.25,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p1 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

base <- ggplot(subset(df,Comparison!="IIIB_to_IIIB"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="IIIB_to_IIIB"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p2 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

base <- ggplot(subset(df,Comparison!="IIIB_to_IIIB"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="IIIB_to_IIIB"), aes(x = Coverage)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p3 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

ggarrange(p1, p2, p3, nrow=3,labels=c("Correct","Uncorrect_ind","Uncorrect_all"),common.legend = TRUE, legend="top", font.label=list(color="black",size=10))

max(subset(df,Comparison!="IIIB_to_IIIB")$Coverage)

df <- read.delim('alltoIIIB.txt', header = FALSE)

colnames(df) <- c("Comparison","Coverage")

base <- ggplot(subset(df,Comparison!="IIIB_to_IIIB"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="IIIB_to_IIIB"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p4 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

p4

###IIIC###

df <- read.delim('IIICtoall.txt', header = FALSE)

colnames(df) <- c("Comparison","Coverage")

base <- ggplot(subset(df,Comparison=="IIIC_to_IIIC"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison=="IIIC_to_IIIC"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=0.25,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p1 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

base <- ggplot(subset(df,Comparison!="IIIC_to_IIIC"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="IIIC_to_IIIC"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p2 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

base <- ggplot(subset(df,Comparison!="IIIC_to_IIIC"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="IIIC_to_IIIC"), aes(x = Coverage)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p3 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

ggarrange(p1, p2, p3, nrow=3,labels=c("Correct","Uncorrect_ind","Uncorrect_all"),common.legend = TRUE, legend="top", font.label=list(color="black",size=10))

max(subset(df,Comparison!="IIIC_to_IIIC")$Coverage)

df <- read.delim('alltoIIIC.txt', header = FALSE)

colnames(df) <- c("Comparison","Coverage")

base <- ggplot(subset(df,Comparison!="IIIC_to_IIIC"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="IIIC_to_IIIC"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p4 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

p4

###IIID###

df <- read.delim('IIIDtoall.txt', header = FALSE)

colnames(df) <- c("Comparison","Coverage")

base <- ggplot(subset(df,Comparison=="IIID_to_IIID"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison=="IIID_to_IIID"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=0.25,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p1 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

base <- ggplot(subset(df,Comparison!="IIID_to_IIID"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="IIID_to_IIID"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p2 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

base <- ggplot(subset(df,Comparison!="IIID_to_IIID"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="IIID_to_IIID"), aes(x = Coverage)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p3 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

ggarrange(p1, p2, p3, nrow=3,labels=c("Correct","Uncorrect_ind","Uncorrect_all"),common.legend = TRUE, legend="top", font.label=list(color="black",size=10))

max(subset(df,Comparison!="IIID_to_IIID")$Coverage)

df <- read.delim('alltoIIID.txt', header = FALSE)

colnames(df) <- c("Comparison","Coverage")

base <- ggplot(subset(df,Comparison!="IIID_to_IIID"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="IIID_to_IIID"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p4 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

p4

###IVA###

df <- read.delim('IVAtoall.txt', header = FALSE)

colnames(df) <- c("Comparison","Coverage")

base <- ggplot(subset(df,Comparison=="IVA_to_IVA"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison=="IVA_to_IVA"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=0.25,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p1 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

base <- ggplot(subset(df,Comparison!="IVA_to_IVA"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="IVA_to_IVA"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p2 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

base <- ggplot(subset(df,Comparison!="IVA_to_IVA"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="IVA_to_IVA"), aes(x = Coverage)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p3 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

ggarrange(p1, p2, p3, nrow=3,labels=c("Correct","Uncorrect_ind","Uncorrect_all"),common.legend = TRUE, legend="top", font.label=list(color="black",size=10))

max(subset(df,Comparison!="IVA_to_IVA")$Coverage)

df <- read.delim('alltoIVA.txt', header = FALSE)

colnames(df) <- c("Comparison","Coverage")

base <- ggplot(subset(df,Comparison!="IVA_to_IVA"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="IVA_to_IVA"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p4 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

p4

###VA###

df <- read.delim('VAtoall.txt', header = FALSE)

colnames(df) <- c("Comparison","Coverage")

base <- ggplot(subset(df,Comparison=="VA_to_VA"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison=="VA_to_VA"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=0.25,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p1 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

base <- ggplot(subset(df,Comparison!="VA_to_VA"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="VA_to_VA"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p2 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

base <- ggplot(subset(df,Comparison!="VA_to_VA"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="VA_to_VA"), aes(x = Coverage)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p3 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

ggarrange(p1, p2, p3, nrow=3,labels=c("Correct","Uncorrect_ind","Uncorrect_all"),common.legend = TRUE, legend="top", font.label=list(color="black",size=10))

max(subset(df,Comparison!="VA_to_VA")$Coverage)

df <- read.delim('alltoVA.txt', header = FALSE)

colnames(df) <- c("Comparison","Coverage")

base <- ggplot(subset(df,Comparison!="VA_to_VA"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="VA_to_VA"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p4 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

p4

###VB###

df <- read.delim('VBtoall.txt', header = FALSE)

colnames(df) <- c("Comparison","Coverage")

base <- ggplot(subset(df,Comparison=="VB_to_VB"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison=="VB_to_VB"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=0.25,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p1 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

base <- ggplot(subset(df,Comparison!="VB_to_VB"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="VB_to_VB"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p2 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

base <- ggplot(subset(df,Comparison!="VB_to_VB"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="VB_to_VB"), aes(x = Coverage)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p3 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

ggarrange(p1, p2, p3, nrow=3,labels=c("Correct","Uncorrect_ind","Uncorrect_all"),common.legend = TRUE, legend="top", font.label=list(color="black",size=10))

max(subset(df,Comparison!="VB_to_VB")$Coverage)

df <- read.delim('alltoVB.txt', header = FALSE)

colnames(df) <- c("Comparison","Coverage")

base <- ggplot(subset(df,Comparison!="VB_to_VB"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="VB_to_VB"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p4 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

p4

###VK###

df <- read.delim('VKtoall.txt', header = FALSE)

colnames(df) <- c("Comparison","Coverage")

base <- ggplot(subset(df,Comparison=="VK_to_VK"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison=="VK_to_VK"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=0.25,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p1 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

base <- ggplot(subset(df,Comparison!="VK_to_VK"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="VK_to_VK"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p2 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

base <- ggplot(subset(df,Comparison!="VK_to_VK"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="VK_to_VK"), aes(x = Coverage)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p3 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

ggarrange(p1, p2, p3, nrow=3,labels=c("Correct","Uncorrect_ind","Uncorrect_all"),common.legend = TRUE, legend="top", font.label=list(color="black",size=10))

max(subset(df,Comparison!="VK_to_VK")$Coverage)

df <- read.delim('alltoVK.txt', header = FALSE)

colnames(df) <- c("Comparison","Coverage")

base <- ggplot(subset(df,Comparison!="VK_to_VK"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="VK_to_VK"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p4 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

p4

###VIA###

df <- read.delim('VIAtoall.txt', header = FALSE)

colnames(df) <- c("Comparison","Coverage")

base <- ggplot(subset(df,Comparison=="VIA_to_VIA"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison=="VIA_to_VIA"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=0.25,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p1 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

base <- ggplot(subset(df,Comparison!="VIA_to_VIA"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="VIA_to_VIA"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p2 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

base <- ggplot(subset(df,Comparison!="VIA_to_VIA"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="VIA_to_VIA"), aes(x = Coverage)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p3 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

ggarrange(p1, p2, p3, nrow=3,labels=c("Correct","Uncorrect_ind","Uncorrect_all"),common.legend = TRUE, legend="top", font.label=list(color="black",size=10))

max(subset(df,Comparison!="VIA_to_VIA")$Coverage)

df <- read.delim('alltoVIA.txt', header = FALSE)

colnames(df) <- c("Comparison","Coverage")

base <- ggplot(subset(df,Comparison!="VIA_to_VIA"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="VIA_to_VIA"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p4 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

p4

###VIB###

df <- read.delim('VIBtoall.txt', header = FALSE)

colnames(df) <- c("Comparison","Coverage")

base <- ggplot(subset(df,Comparison=="VIB_to_VIB"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison=="VIB_to_VIB"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=0.25,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p1 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

base <- ggplot(subset(df,Comparison!="VIB_to_VIB"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="VIB_to_VIB"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p2 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

base <- ggplot(subset(df,Comparison!="VIB_to_VIB"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="VIB_to_VIB"), aes(x = Coverage)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p3 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

ggarrange(p1, p2, p3, nrow=3,labels=c("Correct","Uncorrect_ind","Uncorrect_all"),common.legend = TRUE, legend="top", font.label=list(color="black",size=10))

max(subset(df,Comparison!="VIB_to_VIB")$Coverage)

df <- read.delim('alltoVIB.txt', header = FALSE)

colnames(df) <- c("Comparison","Coverage")

base <- ggplot(subset(df,Comparison!="VIB_to_VIB"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="VIB_to_VIB"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p4 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

p4

###VID###

df <- read.delim('VIDtoall.txt', header = FALSE)

colnames(df) <- c("Comparison","Coverage")

base <- ggplot(subset(df,Comparison=="VID_to_VID"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison=="VID_to_VID"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=0.25,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p1 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

base <- ggplot(subset(df,Comparison!="VID_to_VID"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="VID_to_VID"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p2 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

base <- ggplot(subset(df,Comparison!="VID_to_VID"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="VID_to_VID"), aes(x = Coverage)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p3 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

ggarrange(p1, p2, p3, nrow=3,labels=c("Correct","Uncorrect_ind","Uncorrect_all"),common.legend = TRUE, legend="top", font.label=list(color="black",size=10))

max(subset(df,Comparison!="VID_to_VID")$Coverage)

df <- read.delim('alltoVID.txt', header = FALSE)

colnames(df) <- c("Comparison","Coverage")

base <- ggplot(subset(df,Comparison!="VID_to_VID"), aes(x = Coverage)) + geom_histogram(binwidth = 1,alpha=0) +xlim(0,100)  + scale_y_continuous(name = "Count")+theme(legend.position = "none",text=element_text(size = 10,face="bold"))
dense <- ggplot(subset(df,Comparison!="VID_to_VID"), aes(x = Coverage, fill=Comparison)) + geom_density(adjust=2,alpha=0.5)+scale_y_continuous(position = "left")+xlim(0,100) +theme(panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.background = element_rect(fill='transparent'),legend.box.background = element_rect(fill='transparent'),text=element_text(size = 10,face="bold"),legend.position = "top") + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")
aligned_plots <- align_plots(base, dense, align="hv", axis="tblr")

p4 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

p4