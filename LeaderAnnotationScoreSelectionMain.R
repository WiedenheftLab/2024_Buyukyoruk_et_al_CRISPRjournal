#CURRENT currated scores
# IA_score <-0.05
# IB_score <- 0.1
# IC_score <- 0.63
# ID_score <- 0.1
# IE_score <- 0.54
# IF_score <- 0.36
# IG_score <- 0.19
# IIA_score <- 0.35
# IIB_score <- 0.05
# IIC_score <- 0.17
# IIIA_score <- 0.09
# IIIB_score <- 0.54
# IIIC_score <- 0.05
# IIID_score <- 0.1
# IVA_score <- 0.05
# IVC_score <- 0.05
# IVD_score <- 0.05
# IVE_score <- 0.05
# VA_score <- 0.05
# VB_score <- 0.05
# VF_score <- 0.05
# VIA_score <- 0.05
# VIB_score <- 0.05
# VIC_score <- 0.05
# VID_score <- 0.05
# VK_score <- 0.06

#Trial score curration
IA_score <-0
IB_score <- 0
IC_score <- 0
ID_score <- 0
IE_score <- 0
IF_score <- 0
IG_score <- 0
IIA_score <- 0
IIB_score <- 0
IIC_score <- 0
IIIA_score <- 0
IIIB_score <- 0
IIIC_score <- 0
IIID_score <- 0
IVA_score <- 0
IVC_score <- 0
IVD_score <- 0
IVE_score <- 0
VA_score <- 0
VB_score <- 0
VF_score <- 0
VIA_score <- 0
VIB_score <- 0
VIC_score <- 0
VID_score <- 0
VK_score <- 0

number_sequence <- seq(0,1,0.02)
subsubsub <- c("IB", "IC", "ID", "IE", "IF", "IG", "IIA", "IIC", "IIIA", "IIIB", "IIID", "VK")

for (l in 1:length(subsubsub)){
  rowcount <- 1
  score_df <- data.frame(matrix(nrow=51,ncol=3))
  colnames(score_df) <- c("cutoff","positive","negative")
  for (i in number_sequence){
    pos <- NA
    neg <- NA
    assign(paste0(subsubsub[l],"_score"), i)
    try(source("LeaderAnnotationScoreSelectionRun.R")) # Requires LeaderAnnotationScoreSelectionRun.R Rscript provided in the github repository.
    accuracy_CD_test
    pos <- matrix_chart[subsubsub[l],][subsubsub[l]]
    neg <- sum(data.frame(matrix_chart[subsubsub[l],])[!(row.names(data.frame(matrix_chart[subsubsub[l],])) %in% subsubsub[l]),])
    score_df$cutoff[rowcount] <- get(paste0(subsubsub[l],"_score"))
    score_df$positive[rowcount] <- pos
    score_df$negative[rowcount] <- neg
    print(subsubsub[l])
    print(get(paste0(subsubsub[l],"_score")))
    print(rowcount)
    rowcount <- rowcount + 1
  }
  assign(paste0("score_df_",subsubsub[l]),score_df)
}

CD_IB <- ggplot(score_df_IB, aes(x=cutoff)) + geom_line(aes(y=positive), color = "#56B4E9") + geom_area(aes(y=positive), fill = "#56B4E9", alpha = 1) + geom_line(aes(y=negative), color = "#E69F00") + geom_area(aes(y=negative), fill = "#E69F00", alpha = 1) + scale_x_continuous(name = "IB Score Cutoff value")+ scale_y_continuous(name = "Count")+theme(text=element_text(size = 14,face="bold")) + geom_vline (xintercept = 0.1)+ theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
CD_IC <- ggplot(score_df_IC, aes(x=cutoff)) + geom_line(aes(y=positive), color = "#56B4E9") + geom_area(aes(y=positive), fill = "#56B4E9", alpha = 1) + geom_line(aes(y=negative), color = "#E69F00") + geom_area(aes(y=negative), fill = "#E69F00", alpha = 1) + scale_x_continuous(name = "IC Score Cutoff value")+ scale_y_continuous(name = "Count")+theme(text=element_text(size = 14,face="bold")) + geom_vline (xintercept = 0.1)+ theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
CD_ID <- ggplot(score_df_ID, aes(x=cutoff)) + geom_line(aes(y=positive), color = "#56B4E9") + geom_area(aes(y=positive), fill = "#56B4E9", alpha = 1) + geom_line(aes(y=negative), color = "#E69F00") + geom_area(aes(y=negative), fill = "#E69F00", alpha = 1) + scale_x_continuous(name = "ID Score Cutoff value")+ scale_y_continuous(name = "Count")+theme(text=element_text(size = 14,face="bold")) + geom_vline (xintercept = 0.1)+ theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
CD_IE <- ggplot(score_df_IE, aes(x=cutoff)) + geom_line(aes(y=positive), color = "#56B4E9") + geom_area(aes(y=positive), fill = "#56B4E9", alpha = 1) + geom_line(aes(y=negative), color = "#E69F00") + geom_area(aes(y=negative), fill = "#E69F00", alpha = 1) + scale_x_continuous(name = "IE Score Cutoff value")+ scale_y_continuous(name = "Count")+theme(text=element_text(size = 14,face="bold")) + geom_vline (xintercept = 0.26) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
CD_IF <- ggplot(score_df_IF, aes(x=cutoff)) + geom_line(aes(y=positive), color = "#56B4E9") + geom_area(aes(y=positive), fill = "#56B4E9", alpha = 1) + geom_line(aes(y=negative), color = "#E69F00") + geom_area(aes(y=negative), fill = "#E69F00", alpha = 1) + scale_x_continuous(name = "IF Score Cutoff value")+ scale_y_continuous(name = "Count")+theme(text=element_text(size = 14,face="bold")) + geom_vline (xintercept = 0.36)+ theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
CD_IG <- ggplot(score_df_IG, aes(x=cutoff)) + geom_line(aes(y=positive), color = "#56B4E9") + geom_area(aes(y=positive), fill = "#56B4E9", alpha = 1) + geom_line(aes(y=negative), color = "#E69F00") + geom_area(aes(y=negative), fill = "#E69F00", alpha = 1) + scale_x_continuous(name = "IG Score Cutoff value")+ scale_y_continuous(name = "Count")+theme(text=element_text(size = 14,face="bold")) + geom_vline (xintercept = 0.19)+ theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
CD_IIA <- ggplot(score_df_IIA, aes(x=cutoff)) + geom_line(aes(y=positive), color = "#56B4E9") + geom_area(aes(y=positive), fill = "#56B4E9", alpha = 1) + geom_line(aes(y=negative), color = "#E69F00") + geom_area(aes(y=negative), fill = "#E69F00", alpha = 1) + scale_x_continuous(name = "IIA Score Cutoff value")+ scale_y_continuous(name = "Count")+theme(text=element_text(size = 14,face="bold")) + geom_vline (xintercept = 0.3)+ theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
CD_IIC <- ggplot(score_df_IIC, aes(x=cutoff)) + geom_line(aes(y=positive), color = "#56B4E9") + geom_area(aes(y=positive), fill = "#56B4E9", alpha = 1) + geom_line(aes(y=negative), color = "#E69F00") + geom_area(aes(y=negative), fill = "#E69F00", alpha = 1) + scale_x_continuous(name = "IIC Score Cutoff value")+ scale_y_continuous(name = "Count")+theme(text=element_text(size = 14,face="bold")) + geom_vline (xintercept = 0.04)+ theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
CD_IIIA <- ggplot(score_df_IIIA, aes(x=cutoff)) + geom_line(aes(y=positive), color = "#56B4E9") + geom_area(aes(y=positive), fill = "#56B4E9", alpha = 1) + geom_line(aes(y=negative), color = "#E69F00") + geom_area(aes(y=negative), fill = "#E69F00", alpha = 1) + scale_x_continuous(name = "IIIA Score Cutoff value")+ scale_y_continuous(name = "Count")+theme(text=element_text(size = 14,face="bold")) + geom_vline (xintercept = 0.09)+ theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
CD_IIIB <- ggplot(score_df_IIIB, aes(x=cutoff)) + geom_line(aes(y=positive), color = "#56B4E9") + geom_area(aes(y=positive), fill = "#56B4E9", alpha = 1) + geom_line(aes(y=negative), color = "#E69F00") + geom_area(aes(y=negative), fill = "#E69F00", alpha = 1) + scale_x_continuous(name = "IIIB Score Cutoff value")+ scale_y_continuous(name = "Count")+theme(text=element_text(size = 14,face="bold")) + geom_vline (xintercept = 0.54)+ theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
CD_IIID <- ggplot(score_df_IIID, aes(x=cutoff)) + geom_line(aes(y=positive), color = "#56B4E9") + geom_area(aes(y=positive), fill = "#56B4E9", alpha = 1) + geom_line(aes(y=negative), color = "#E69F00") + geom_area(aes(y=negative), fill = "#E69F00", alpha = 1) + scale_x_continuous(name = "IIID Score Cutoff value")+ scale_y_continuous(name = "Count")+theme(text=element_text(size = 14,face="bold")) + geom_vline (xintercept = 0.1)+ theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
CD_VK <- ggplot(score_df_VK, aes(x=cutoff)) + geom_line(aes(y=positive), color = "#56B4E9") + geom_area(aes(y=positive), fill = "#56B4E9", alpha = 1) + geom_line(aes(y=negative), color = "#E69F00") + geom_area(aes(y=negative), fill = "#E69F00", alpha = 1) + scale_x_continuous(name = "VK Score Cutoff value")+ scale_y_continuous(name = "Count")+theme(text=element_text(size = 14,face="bold")) + geom_vline (xintercept = 0.1)+ theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

ggplot(score_df_IB, aes(x=cutoff)) + geom_line(aes(y=positive), color = "#56B4E9") + geom_area(aes(y=positive), fill = "#56B4E9", alpha = 1) + geom_line(aes(y=negative), color = "#E69F00") + geom_area(aes(y=negative), fill = "#E69F00", alpha = 1) + scale_x_continuous(name = "IB Score Cutoff value")+ scale_y_continuous(name = "Count")+theme(text=element_text(size = 14,face="bold")) + geom_vline (xintercept = 0.1)
ggplot(score_df_IC, aes(x=cutoff)) + geom_line(aes(y=positive), color = "#56B4E9") + geom_area(aes(y=positive), fill = "#56B4E9", alpha = 1) + geom_line(aes(y=negative), color = "#E69F00") + geom_area(aes(y=negative), fill = "#E69F00", alpha = 1) + scale_x_continuous(name = "IC Score Cutoff value")+ scale_y_continuous(name = "Count")+theme(text=element_text(size = 14,face="bold")) + geom_vline (xintercept = 0.1)
ggplot(score_df_ID, aes(x=cutoff)) + geom_line(aes(y=positive), color = "#56B4E9") + geom_area(aes(y=positive), fill = "#56B4E9", alpha = 1) + geom_line(aes(y=negative), color = "#E69F00") + geom_area(aes(y=negative), fill = "#E69F00", alpha = 1) + scale_x_continuous(name = "ID Score Cutoff value")+ scale_y_continuous(name = "Count")+theme(text=element_text(size = 14,face="bold")) + geom_vline (xintercept = 0.1)
ggplot(score_df_IE, aes(x=cutoff)) + geom_line(aes(y=positive), color = "#56B4E9") + geom_area(aes(y=positive), fill = "#56B4E9", alpha = 1) + geom_line(aes(y=negative), color = "#E69F00") + geom_area(aes(y=negative), fill = "#E69F00", alpha = 1) + scale_x_continuous(name = "IE Score Cutoff value")+ scale_y_continuous(name = "Count")+theme(text=element_text(size = 14,face="bold")) + geom_vline (xintercept = 0.26)
ggplot(score_df_IF, aes(x=cutoff)) + geom_line(aes(y=positive), color = "#56B4E9") + geom_area(aes(y=positive), fill = "#56B4E9", alpha = 1) + geom_line(aes(y=negative), color = "#E69F00") + geom_area(aes(y=negative), fill = "#E69F00", alpha = 1) + scale_x_continuous(name = "IF Score Cutoff value")+ scale_y_continuous(name = "Count")+theme(text=element_text(size = 14,face="bold")) + geom_vline (xintercept = 0.36)
ggplot(score_df_IG, aes(x=cutoff)) + geom_line(aes(y=positive), color = "#56B4E9") + geom_area(aes(y=positive), fill = "#56B4E9", alpha = 1) + geom_line(aes(y=negative), color = "#E69F00") + geom_area(aes(y=negative), fill = "#E69F00", alpha = 1) + scale_x_continuous(name = "IG Score Cutoff value")+ scale_y_continuous(name = "Count")+theme(text=element_text(size = 14,face="bold")) + geom_vline (xintercept = 0.19)
ggplot(score_df_IIA, aes(x=cutoff)) + geom_line(aes(y=positive), color = "#56B4E9") + geom_area(aes(y=positive), fill = "#56B4E9", alpha = 1) + geom_line(aes(y=negative), color = "#E69F00") + geom_area(aes(y=negative), fill = "#E69F00", alpha = 1) + scale_x_continuous(name = "IIA Score Cutoff value")+ scale_y_continuous(name = "Count")+theme(text=element_text(size = 14,face="bold")) + geom_vline (xintercept = 0.3)
ggplot(score_df_IIC, aes(x=cutoff)) + geom_line(aes(y=positive), color = "#56B4E9") + geom_area(aes(y=positive), fill = "#56B4E9", alpha = 1) + geom_line(aes(y=negative), color = "#E69F00") + geom_area(aes(y=negative), fill = "#E69F00", alpha = 1) + scale_x_continuous(name = "IIC Score Cutoff value")+ scale_y_continuous(name = "Count")+theme(text=element_text(size = 14,face="bold")) + geom_vline (xintercept = 0.04)
ggplot(score_df_IIIA, aes(x=cutoff)) + geom_line(aes(y=positive), color = "#56B4E9") + geom_area(aes(y=positive), fill = "#56B4E9", alpha = 1) + geom_line(aes(y=negative), color = "#E69F00") + geom_area(aes(y=negative), fill = "#E69F00", alpha = 1) + scale_x_continuous(name = "IIIA Score Cutoff value")+ scale_y_continuous(name = "Count")+theme(text=element_text(size = 14,face="bold")) + geom_vline (xintercept = 0.09)
ggplot(score_df_IIIB, aes(x=cutoff)) + geom_line(aes(y=positive), color = "#56B4E9") + geom_area(aes(y=positive), fill = "#56B4E9", alpha = 1) + geom_line(aes(y=negative), color = "#E69F00") + geom_area(aes(y=negative), fill = "#E69F00", alpha = 1) + scale_x_continuous(name = "IIIB Score Cutoff value")+ scale_y_continuous(name = "Count")+theme(text=element_text(size = 14,face="bold")) + geom_vline (xintercept = 0.54)
ggplot(score_df_IIID, aes(x=cutoff)) + geom_line(aes(y=positive), color = "#56B4E9") + geom_area(aes(y=positive), fill = "#56B4E9", alpha = 1) + geom_line(aes(y=negative), color = "#E69F00") + geom_area(aes(y=negative), fill = "#E69F00", alpha = 1) + scale_x_continuous(name = "IIID Score Cutoff value")+ scale_y_continuous(name = "Count")+theme(text=element_text(size = 14,face="bold")) + geom_vline (xintercept = 0.1)
ggplot(score_df_VK, aes(x=cutoff)) + geom_line(aes(y=positive), color = "#56B4E9") + geom_area(aes(y=positive), fill = "#56B4E9", alpha = 1) + geom_line(aes(y=negative), color = "#E69F00") + geom_area(aes(y=negative), fill = "#E69F00", alpha = 1) + scale_x_continuous(name = "VK Score Cutoff value")+ scale_y_continuous(name = "Count")+theme(text=element_text(size = 14,face="bold")) + geom_vline (xintercept = 0.1)

library(ggpubr)
ggarrange(CD_IB,CD_IC,CD_ID,CD_IE,CD_IF,CD_IG,CD_IIA,CD_IIC,CD_IIIA,CD_IIIB,CD_IIID,CD_VK,nrow =12, heights = c(10,10,10,10,10,10,10,10,10,10,10,10),align="v")

IB_perc <- data.frame(read_csv("Score_cutoff_training_data/score_df_IB.csv"))
IC_perc <- data.frame(read_csv("Score_cutoff_training_data/score_df_IC.csv"))
ID_perc <- data.frame(read_csv("Score_cutoff_training_data/score_df_ID.csv"))
IE_perc <- data.frame(read_csv("Score_cutoff_training_data/score_df_IE.csv"))
IF_perc <- data.frame(read_csv("Score_cutoff_training_data/score_df_IF.csv"))
IG_perc <- data.frame(read_csv("Score_cutoff_training_data/score_df_IG.csv"))
IIA_perc <- data.frame(read_csv("Score_cutoff_training_data/score_df_IIA.csv"))
IIC_perc <- data.frame(read_csv("Score_cutoff_training_data/score_df_IIC.csv"))
IIIA_perc <- data.frame(read_csv("Score_cutoff_training_data/score_df_IIIA.csv"))
IIIB_perc <- data.frame(read_csv("Score_cutoff_training_data/score_df_IIIB.csv"))
IIID_perc <- data.frame(read_csv("Score_cutoff_training_data/score_df_IIID.csv"))
VK_perc <- data.frame(read_csv("Score_cutoff_training_data/score_df_VK.csv"))

CD_IB_perc <- ggplot(IB_perc, aes(x=cutoff)) + geom_line(aes(y=positive_perc), color = "#56B4E9") + geom_area(aes(y=positive_perc), fill = "#56B4E9", alpha = 1) + geom_line(aes(y=negative_perc), color = "#E69F00") + geom_area(aes(y=negative_perc), fill = "#E69F00", alpha = 1) + scale_x_continuous(name = "IB Score Cutoff value")+ scale_y_continuous(name = "Count")+theme(text=element_text(size = 14,face="bold")) + geom_vline (xintercept = 0.16)+ theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+ scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1),breaks = 10^(-1:2),limits = c(0,100))
CD_IC_perc <- ggplot(IC_perc, aes(x=cutoff)) + geom_line(aes(y=positive_perc), color = "#56B4E9") + geom_area(aes(y=positive_perc), fill = "#56B4E9", alpha = 1) + geom_line(aes(y=negative_perc), color = "#E69F00") + geom_area(aes(y=negative_perc), fill = "#E69F00", alpha = 1) + scale_x_continuous(name = "IB Score Cutoff value")+ scale_y_continuous(name = "Count")+theme(text=element_text(size = 14,face="bold")) + geom_vline (xintercept = 0.1)+ theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+ scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1),breaks = 10^(-1:2),limits = c(0,100))
CD_ID_perc <- ggplot(ID_perc, aes(x=cutoff)) + geom_line(aes(y=positive_perc), color = "#56B4E9") + geom_area(aes(y=positive_perc), fill = "#56B4E9", alpha = 1) + geom_line(aes(y=negative_perc), color = "#E69F00") + geom_area(aes(y=negative_perc), fill = "#E69F00", alpha = 1) + scale_x_continuous(name = "IB Score Cutoff value")+ scale_y_continuous(name = "Count")+theme(text=element_text(size = 14,face="bold")) + geom_vline (xintercept = 0.16)+ theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+ scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1),breaks = 10^(-1:2),limits = c(0,100))
CD_IE_perc <- ggplot(IE_perc, aes(x=cutoff)) + geom_line(aes(y=positive_perc), color = "#56B4E9") + geom_area(aes(y=positive_perc), fill = "#56B4E9", alpha = 1) + geom_line(aes(y=negative_perc), color = "#E69F00") + geom_area(aes(y=negative_perc), fill = "#E69F00", alpha = 1) + scale_x_continuous(name = "IB Score Cutoff value")+ scale_y_continuous(name = "Count")+theme(text=element_text(size = 14,face="bold")) + geom_vline (xintercept = 0.48)+ theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+ scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1),breaks = 10^(-1:2),limits = c(0,100))
CD_IF_perc <- ggplot(IF_perc, aes(x=cutoff)) + geom_line(aes(y=positive_perc), color = "#56B4E9") + geom_area(aes(y=positive_perc), fill = "#56B4E9", alpha = 1) + geom_line(aes(y=negative_perc), color = "#E69F00") + geom_area(aes(y=negative_perc), fill = "#E69F00", alpha = 1) + scale_x_continuous(name = "IB Score Cutoff value")+ scale_y_continuous(name = "Count")+theme(text=element_text(size = 14,face="bold")) + geom_vline (xintercept = 0.36)+ theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+ scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1),breaks = 10^(-1:2),limits = c(0,100))
CD_IG_perc <- ggplot(IG_perc, aes(x=cutoff)) + geom_line(aes(y=positive_perc), color = "#56B4E9") + geom_area(aes(y=positive_perc), fill = "#56B4E9", alpha = 1) + geom_line(aes(y=negative_perc), color = "#E69F00") + geom_area(aes(y=negative_perc), fill = "#E69F00", alpha = 1) + scale_x_continuous(name = "IB Score Cutoff value")+ scale_y_continuous(name = "Count")+theme(text=element_text(size = 14,face="bold")) + geom_vline (xintercept = 0.42)+ theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+ scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1),breaks = 10^(-1:2),limits = c(0,100))
CD_IIA_perc <- ggplot(IIA_perc, aes(x=cutoff)) + geom_line(aes(y=positive_perc), color = "#56B4E9") + geom_area(aes(y=positive_perc), fill = "#56B4E9", alpha = 1) + geom_line(aes(y=negative_perc), color = "#E69F00") + geom_area(aes(y=negative_perc), fill = "#E69F00", alpha = 1) + scale_x_continuous(name = "IB Score Cutoff value")+ scale_y_continuous(name = "Count")+theme(text=element_text(size = 14,face="bold")) + geom_vline (xintercept = 0.3)+ theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+ scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1),breaks = 10^(-1:2),limits = c(0,100))
CD_IIC_perc <- ggplot(IIC_perc, aes(x=cutoff)) + geom_line(aes(y=positive_perc), color = "#56B4E9") + geom_area(aes(y=positive_perc), fill = "#56B4E9", alpha = 1) + geom_line(aes(y=negative_perc), color = "#E69F00") + geom_area(aes(y=negative_perc), fill = "#E69F00", alpha = 1) + scale_x_continuous(name = "IB Score Cutoff value")+ scale_y_continuous(name = "Count")+theme(text=element_text(size = 14,face="bold")) + geom_vline (xintercept = 0.04)+ theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+ scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1),breaks = 10^(-1:2),limits = c(0,100))
CD_IIIA_perc <- ggplot(IIIA_perc, aes(x=cutoff)) + geom_line(aes(y=positive_perc), color = "#56B4E9") + geom_area(aes(y=positive_perc), fill = "#56B4E9", alpha = 1) + geom_line(aes(y=negative_perc), color = "#E69F00") + geom_area(aes(y=negative_perc), fill = "#E69F00", alpha = 1) + scale_x_continuous(name = "IB Score Cutoff value")+ scale_y_continuous(name = "Count")+theme(text=element_text(size = 14,face="bold")) + geom_vline (xintercept = 0.1)+ theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+ scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1),breaks = 10^(-1:2),limits = c(0,100))
CD_IIIB_perc <- ggplot(IIIB_perc, aes(x=cutoff)) + geom_line(aes(y=positive_perc), color = "#56B4E9") + geom_area(aes(y=positive_perc), fill = "#56B4E9", alpha = 1) + geom_line(aes(y=negative_perc), color = "#E69F00") + geom_area(aes(y=negative_perc), fill = "#E69F00", alpha = 1) + scale_x_continuous(name = "IB Score Cutoff value")+ scale_y_continuous(name = "Count")+theme(text=element_text(size = 14,face="bold")) + geom_vline (xintercept = 0.64)+ theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+ scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1),breaks = 10^(-1:2),limits = c(0,100))
CD_IIID_perc <- ggplot(IIID_perc, aes(x=cutoff)) + geom_line(aes(y=positive_perc), color = "#56B4E9") + geom_area(aes(y=positive_perc), fill = "#56B4E9", alpha = 1) + geom_line(aes(y=negative_perc), color = "#E69F00") + geom_area(aes(y=negative_perc), fill = "#E69F00", alpha = 1) + scale_x_continuous(name = "IB Score Cutoff value")+ scale_y_continuous(name = "Count")+theme(text=element_text(size = 14,face="bold")) + geom_vline (xintercept = 0.22)+ theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+ scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1),breaks = 10^(-1:2),limits = c(0,100))
CD_VK_perc <- ggplot(VK_perc, aes(x=cutoff)) + geom_line(aes(y=positive_perc), color = "#56B4E9") + geom_area(aes(y=positive_perc), fill = "#56B4E9", alpha = 1) + geom_line(aes(y=negative_perc), color = "#E69F00") + geom_area(aes(y=negative_perc), fill = "#E69F00", alpha = 1) + scale_x_continuous(name = "IB Score Cutoff value")+ scale_y_continuous(name = "Count")+theme(text=element_text(size = 14,face="bold")) + geom_vline (xintercept = 0.1)+ theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+ scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1),breaks = 10^(-1:2),limits = c(0,100))

ggarrange(CD_IB_perc,CD_IC_perc,CD_ID_perc,CD_IE_perc,CD_IF_perc,CD_IG_perc,CD_IIA_perc,CD_IIC_perc,CD_IIIA_perc,CD_IIIB_perc,CD_IIID_perc,CD_VK_perc,nrow =12, heights = c(10,10,10,10,10,10,10,10,10,10,10,10),align="v")