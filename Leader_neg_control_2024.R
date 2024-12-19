library(matrixStats)
library(ggplot2)
library(ggbreak)
library(RColorBrewer)

agg <- function(my_data) {    # Create user-defined function
  if(length(my_data$Array) == 0 ) {
    my_data
  } else {
    aggregate(Score~Array,my_data,sum)
  }
}

fun_int0 <- function(my_data) {    # Create user-defined function

  if(length(my_data) == 0 & is.integer(my_data)) {
    0
  } else {
    my_data
  }
}

# Load control dataframes located in provided Database files. You might need to provide to full path name before the Found_motifs_in_control_* base name.

for (i in 1:3) {
  dir_control <- paste(paste("Found_motifs_in_control_",i,sep=""),"/",sep="")
  setwd(dir_control)
  print(getwd())

  df_control <- paste("Control",i,sep="")

  if(i==1){
    assign(df_control,data.frame(Type="Type III",Subtype="IIIB",Taxa="Pseudomonadota_motif_2",Motif_ID="IIIB_Pseudomonadota_motif_2",y=0))
  }else{
    assign(df_control,data.frame())
  }
  for (l in 1:length(list.files())) {
    looking <- gsub(".txt","",list.files()[l])
    if (looking %in% motifs){
      line_num <- length(readLines(list.files()[l]))-1
      sub_looking <- str_split_fixed(looking,"_",2)[1]
      taxa_looking <- str_split_fixed(looking,"_",2)[2]
      type <- paste("Type",substr(sub_looking,1,nchar(sub_looking)-1))
      assign(df_control,rbind(get(df_control),(data.frame(Type=type,Subtype=sub_looking,Taxa=taxa_looking,Motif_ID=looking,y=(line_num*100/1000)))))
    }
  }
}

colnames(Control1)[ncol(Control1)] <- "Control1_perc"
colnames(Control2)[ncol(Control2)] <- "Control2_perc"
colnames(Control3)[ncol(Control3)] <- "Control3_perc"

df_control_motifs <- data.frame()
for (i in 1: 2){
  if (i==1){
    df_control_motifs <- get(paste("Control",i,sep="")) %>% left_join(get(paste("Control",i+1,sep="")), by = c("Motif_ID"="Motif_ID","Taxa"="Taxa","Type"="Type","Subtype"="Subtype"))
  }else{
    df_control_motifs <- df_control_motifs %>% left_join(get(paste("Control",i+1,sep="")),by = c("Motif_ID"="Motif_ID","Taxa"="Taxa","Type"="Type","Subtype"="Subtype"))
  }
}

df_control_motifs$Ave <- (df_control_motifs$Control1_perc+df_control_motifs$Control2_perc+df_control_motifs$Control3_perc)/3
df_control_motifs$Std = rowSds(as.matrix(df_control_motifs[,c(5,6,7)]))

false_hits<- ggplot(df_control_motifs, aes(x=Motif_ID, y=Ave)) +
  geom_bar(stat="identity", color="black",position=position_dodge()) +
  geom_errorbar(aes(ymin=Ave, ymax=Ave+Std), width=.2,position=position_dodge(.9)) + ylim(0,20) +facet_grid(cols = vars(Subtype), scales = "free_x", space = "free_x") +geom_text(aes(label=Taxa,angle=90,x=Motif_ID, y= 19),hjust =1,size=3.5) + scale_color_brewer("Set3") + labs(x="Motif ID",y="Percentage %") + theme(text = element_text(size=14,face = "bold"))

false_hits
#####1#####

# Load control dataframes located in provided Database files. You might need to provide to full path name before the Found_motifs_in_control_* base name.

if(any(grepl("package:matrixStats", search()))) detach("package:matrixStats") else message("matrixStats not loaded")
setwd("Found_motifs_in_control_1/")

for (i in 1: length(list.files())){
  assign(paste("dd",i,sep=""),read.delim(list.files()[i],header = TRUE))
}

df <- read.delim('Control_df_ids.txt', header = TRUE) # Load Control_df_ids.txt file provided in the Database files.

for (i in 1: length(list.files())){
  df <- df %>% left_join(get(paste("dd",i,sep="")), by = "Array")
}

df$Subtype_single <- NA
df$phylum <- NA
df$Subtype_CRISPRDetect <- NA

IA <- data.frame()
IB <- data.frame()
IC <- data.frame()
ID <- data.frame()
IE <- data.frame()
IF <- data.frame()
IG <- data.frame()
IIA <- data.frame()
IIB<- data.frame()
IIC<- data.frame()
IIIA<- data.frame()
IIIB<- data.frame()
IIID<- data.frame()
IVA<- data.frame()
VB<- data.frame()
VK<- data.frame()

df$Header <- df$Array

for (i in 1:length(motifs)) {
  motif_list <- (colnames(df)[grepl(gsub("-", ".",motifs[i]),colnames(df))])
  for (r in 1:(length(motif_list))){
    base <- paste(str_split_fixed(motif_list[r],"_",5)[1:4],collapse = "_")
    if (base != gsub("-", ".",motifs[i])) {
      motif_list[r] <- NA
    }
  }
  motif_list <- motif_list[!is.na(motif_list)]
  for (l in 1:(length(motif_list)/2)) {
    motif <- paste(gsub("-", "_",motifs[i]),l,sep = "_")
    motif_fetch <- paste(gsub("-", ".",motifs[i]),l,sep = "_")
    if (str_split_fixed(motif,"_",3)[1] == str_split_fixed(colnames(df)[grepl(gsub("-", ".",motifs[i]),colnames(df))][l*2-1],"_",3)[1]) {
      motif_pval <- paste(motif_fetch,"pvalue",sep = "_")
      current_motif <- (paste(paste("df_",gsub("-", "_",motifs[i]),sep=""),l,sep = "_"))
      name_motif <- (paste("df_",motifs[i],sep=""))
      current_motif_pvalue <- paste(name_motif,"pvalue",sep = "_")
      temp_df <- data.frame(Header=df$Header,Array=df$Array,name_motif=df[motif_fetch],current_motif_pvalue=df[motif_pval],Subtype_single=df$Subtype_single,Phylum=df$phylum,Subtype_CRISPRDetect=df$Subtype_CRISPRDetect)
      colnames(temp_df) <- c('Header','Array',gsub("-", "_",name_motif),gsub("-", "_",current_motif_pvalue),'Subtype_single','Phylum','Subtype_CRISPRDetect')
      assign(current_motif,temp_df)
      subtype <- str_split_fixed(motif,"_",3)[1]
      taxa <- str_split_fixed(motifs[i],"_",3)[2]
      if(l==(length(motif_list)/2)) {
        for (l in 1:(length(motif_list)/2)) {
          current_motif <- (paste(paste("df_",gsub("-", "_",motifs[i]),sep=""),l,sep = "_"))
          freq <- paste(gsub("-", "_",name_motif),"freq",sep = "_")
          updated_df <- get(current_motif) %>% left_join(get(freq), by = gsub("-", "_",name_motif))
          updated_df["Score"] <- updated_df[gsub("-", "_",current_motif_pvalue)] * updated_df$dist_score
          assign(subtype, rbind(get(subtype),subset(updated_df, Score > 0.05)[,c('Array','Score')]))
          assign(current_motif,updated_df)
        }
      }
    }
  }
}

IB_new <- agg(IB)
colnames(IB_new) <- c("Array","IB")
IC_new <- agg(IC)
colnames(IC_new) <- c("Array","IC")
ID_new <- agg(ID)
colnames(ID_new) <- c("Array","ID")
IE_new <- agg(IE)
colnames(IE_new) <- c("Array","IE")
IF_new <- agg(IF)
colnames(IF_new) <- c("Array","IF")
IG_new <- agg(IG)
colnames(IG_new) <- c("Array","IG")
IIA_new <- agg(IIA)
colnames(IIA_new) <- c("Array","IIA")
IIC_new <- agg(IIC)
colnames(IIC_new) <- c("Array","IIC")
IIIA_new <- agg(IIIA)
colnames(IIIA_new) <- c("Array","IIIA")
IIIB_new <- agg(IIIB)
colnames(IIIB_new) <- c("Array","IIIB")
IIID_new <- agg(IIID)
colnames(IIID_new) <- c("Array","IIID")
VK_new <- agg(VK)
colnames(VK_new) <- c("Array","VK")

new_df <- rbind.fill(IB_new,IC_new,ID_new,IE_new,IF_new,IG_new,IIA_new,IIC_new,IIIA_new,IIIB_new,IIID_new,VK_new)
new_df[is.na(new_df)] = 0
collapsed_new <- (aggregate(.~Array,new_df,sum))
collapsed_new$Final_prediction <- as.character(apply(collapsed_new[,-1],1,function(x) names(collapsed_new[,-1])[c(which(x>0.05),which.max(x))[duplicated(c(which(x>0.05),which.max(x)))]]))
collapsed_new[collapsed_new==0]<-NA
collapsed_new[collapsed_new=="character(0)"]<-NA

score_cutoffs <- read.delim("~/Desktop/MotifPrediction_09292023/subtype_annotation_score_cutoffs.txt",header=FALSE)
colnames(score_cutoffs) <- c("subtype","cutoff")
rownames(score_cutoffs) <- score_cutoffs$subtype

for (i in 1:nrow(collapsed_new)){
  if(collapsed_new[i,][collapsed_new$Final_prediction[i]] >= score_cutoffs[collapsed_new$Final_prediction[i],]$cutoff){
    collapsed_new$Final_prediction_after_score[i] <- collapsed_new[i,]$Final_prediction
  } else {
    collapsed_new$Final_prediction_after_score[i] <- "N/A"
  }
}

sum_count <- count(collapsed_new$Final_prediction_after_score)
colnames(sum_count) <- c('Subtype','Count')

final_prediction_Control1 <- df %>% left_join(collapsed_new,by = c("Array"))
colnames(final_prediction_Control1)[ncol(final_prediction_Control1)] <- "Final_prediction_1"

#####2#####

if(any(grepl("package:matrixStats", search()))) detach("package:matrixStats") else message("matrixStats not loaded")

setwd("Found_motifs_in_control_2")

for (i in 1: length(list.files())){
  assign(paste("dd",i,sep=""),read.delim(list.files()[i],header = TRUE))
}

df <- read.delim('Control_df_ids.txt', header = TRUE)

for (i in 1: length(list.files())){
  df <- df %>% left_join(get(paste("dd",i,sep="")), by = "Array")
}

df$Subtype_single <- NA
df$phylum <- NA
df$Subtype_CRISPRDetect <- NA

IA <- data.frame()
IB <- data.frame()
IC <- data.frame()
ID <- data.frame()
IE <- data.frame()
IF <- data.frame()
IG <- data.frame()
IIA <- data.frame()
IIB<- data.frame()
IIC<- data.frame()
IIIA<- data.frame()
IIIB<- data.frame()
IIID<- data.frame()
IVA<- data.frame()
VB<- data.frame()
VK<- data.frame()

df$Header <- df$Array

for (i in 1:length(motifs)) {
  motif_list <- (colnames(df)[grepl(gsub("-", ".",motifs[i]),colnames(df))])
  for (r in 1:(length(motif_list))){
    base <- paste(str_split_fixed(motif_list[r],"_",5)[1:4],collapse = "_")
    if (base != gsub("-", ".",motifs[i])) {
      motif_list[r] <- NA
    }
  }
  motif_list <- motif_list[!is.na(motif_list)]
  for (l in 1:(length(motif_list)/2)) {
    motif <- paste(gsub("-", "_",motifs[i]),l,sep = "_")
    motif_fetch <- paste(gsub("-", ".",motifs[i]),l,sep = "_")
    if (str_split_fixed(motif,"_",3)[1] == str_split_fixed(colnames(df)[grepl(gsub("-", ".",motifs[i]),colnames(df))][l*2-1],"_",3)[1]) {
      motif_pval <- paste(motif_fetch,"pvalue",sep = "_")
      current_motif <- (paste(paste("df_",gsub("-", "_",motifs[i]),sep=""),l,sep = "_"))
      name_motif <- (paste("df_",motifs[i],sep=""))
      current_motif_pvalue <- paste(name_motif,"pvalue",sep = "_")
      temp_df <- data.frame(Header=df$Header,Array=df$Array,name_motif=df[motif_fetch],current_motif_pvalue=df[motif_pval],Subtype_single=df$Subtype_single,Phylum=df$phylum,Subtype_CRISPRDetect=df$Subtype_CRISPRDetect)
      colnames(temp_df) <- c('Header','Array',gsub("-", "_",name_motif),gsub("-", "_",current_motif_pvalue),'Subtype_single','Phylum','Subtype_CRISPRDetect')
      assign(current_motif,temp_df)
      subtype <- str_split_fixed(motif,"_",3)[1]
      taxa <- str_split_fixed(motifs[i],"_",3)[2]
      if(l==(length(motif_list)/2)) {
        for (l in 1:(length(motif_list)/2)) {
          current_motif <- (paste(paste("df_",gsub("-", "_",motifs[i]),sep=""),l,sep = "_"))
          freq <- paste(gsub("-", "_",name_motif),"freq",sep = "_")
          updated_df <- get(current_motif) %>% left_join(get(freq), by = gsub("-", "_",name_motif))
          updated_df["Score"] <- updated_df[gsub("-", "_",current_motif_pvalue)] * updated_df$dist_score
          assign(subtype, rbind(get(subtype),subset(updated_df, Score > 0.05)[,c('Array','Score')]))
          assign(current_motif,updated_df)
        }
      }
    }
  }
}

IB_new <- agg(IB)
colnames(IB_new) <- c("Array","IB")
IC_new <- agg(IC)
colnames(IC_new) <- c("Array","IC")
ID_new <- agg(ID)
colnames(ID_new) <- c("Array","ID")
IE_new <- agg(IE)
colnames(IE_new) <- c("Array","IE")
IF_new <- agg(IF)
colnames(IF_new) <- c("Array","IF")
IG_new <- agg(IG)
colnames(IG_new) <- c("Array","IG")
IIA_new <- agg(IIA)
colnames(IIA_new) <- c("Array","IIA")
IIC_new <- agg(IIC)
colnames(IIC_new) <- c("Array","IIC")
IIIA_new <- agg(IIIA)
colnames(IIIA_new) <- c("Array","IIIA")
IIIB_new <- agg(IIIB)
colnames(IIIB_new) <- c("Array","IIIB")
IIID_new <- agg(IIID)
colnames(IIID_new) <- c("Array","IIID")
VK_new <- agg(VK)
colnames(VK_new) <- c("Array","VK")

new_df <- rbind.fill(IB_new,IC_new,ID_new,IE_new,IF_new,IG_new,IIA_new,IIC_new,IIIA_new,IIIB_new,IIID_new,VK_new)
new_df[is.na(new_df)] = 0
collapsed_new <- (aggregate(.~Array,new_df,sum))
collapsed_new$Final_prediction <- as.character(apply(collapsed_new[,-1],1,function(x) names(collapsed_new[,-1])[c(which(x>0.05),which.max(x))[duplicated(c(which(x>0.05),which.max(x)))]]))
collapsed_new[collapsed_new==0]<-NA
collapsed_new[collapsed_new=="character(0)"]<-NA

score_cutoffs <- read.delim("~/Desktop/MotifPrediction_09292023/subtype_annotation_score_cutoffs.txt",header=FALSE)
colnames(score_cutoffs) <- c("subtype","cutoff")
rownames(score_cutoffs) <- score_cutoffs$subtype

for (i in 1:nrow(collapsed_new)){
  if(collapsed_new[i,][collapsed_new$Final_prediction[i]] >= score_cutoffs[collapsed_new$Final_prediction[i],]$cutoff){
    collapsed_new$Final_prediction_after_score[i] <- collapsed_new[i,]$Final_prediction
  } else {
    collapsed_new$Final_prediction_after_score[i] <- "N/A"
  }
}

sum_count <- count(collapsed_new$Final_prediction_after_score)
colnames(sum_count) <- c('Subtype','Count')

final_prediction_Control2 <- df %>% left_join(collapsed_new,by = c("Array"))
colnames(final_prediction_Control2)[ncol(final_prediction_Control2)] <- "Final_prediction_2"

####3####

if(any(grepl("package:matrixStats", search()))) detach("package:matrixStats") else message("matrixStats not loaded")
setwd("Found_motifs_in_control_3")

for (i in 1: length(list.files())){
  assign(paste("dd",i,sep=""),read.delim(list.files()[i],header = TRUE))
}

df <- read.delim('Control_df_ids.txt', header = TRUE)

for (i in 1: length(list.files())){
  df <- df %>% left_join(get(paste("dd",i,sep="")), by = "Array")
}

df$Subtype_single <- NA
df$phylum <- NA
df$Subtype_CRISPRDetect <- NA

IA <- data.frame()
IB <- data.frame()
IC <- data.frame()
ID <- data.frame()
IE <- data.frame()
IF <- data.frame()
IG <- data.frame()
IIA <- data.frame()
IIB<- data.frame()
IIC<- data.frame()
IIIA<- data.frame()
IIIB<- data.frame()
IIID<- data.frame()
IVA<- data.frame()
VB<- data.frame()
VK<- data.frame()

df$Header <- df$Array

for (i in 1:length(motifs)) {
  motif_list <- (colnames(df)[grepl(gsub("-", ".",motifs[i]),colnames(df))])
  for (r in 1:(length(motif_list))){
    base <- paste(str_split_fixed(motif_list[r],"_",5)[1:4],collapse = "_")
    if (base != gsub("-", ".",motifs[i])) {
      motif_list[r] <- NA
    }
  }
  motif_list <- motif_list[!is.na(motif_list)]
  for (l in 1:(length(motif_list)/2)) {
    motif <- paste(gsub("-", "_",motifs[i]),l,sep = "_")
    motif_fetch <- paste(gsub("-", ".",motifs[i]),l,sep = "_")
    if (str_split_fixed(motif,"_",3)[1] == str_split_fixed(colnames(df)[grepl(gsub("-", ".",motifs[i]),colnames(df))][l*2-1],"_",3)[1]) {
      motif_pval <- paste(motif_fetch,"pvalue",sep = "_")
      current_motif <- (paste(paste("df_",gsub("-", "_",motifs[i]),sep=""),l,sep = "_"))
      name_motif <- (paste("df_",motifs[i],sep=""))
      current_motif_pvalue <- paste(name_motif,"pvalue",sep = "_")
      temp_df <- data.frame(Header=df$Header,Array=df$Array,name_motif=df[motif_fetch],current_motif_pvalue=df[motif_pval],Subtype_single=df$Subtype_single,Phylum=df$phylum,Subtype_CRISPRDetect=df$Subtype_CRISPRDetect)
      colnames(temp_df) <- c('Header','Array',gsub("-", "_",name_motif),gsub("-", "_",current_motif_pvalue),'Subtype_single','Phylum','Subtype_CRISPRDetect')
      assign(current_motif,temp_df)
      subtype <- str_split_fixed(motif,"_",3)[1]
      taxa <- str_split_fixed(motifs[i],"_",3)[2]
      if(l==(length(motif_list)/2)) {
        for (l in 1:(length(motif_list)/2)) {
          current_motif <- (paste(paste("df_",gsub("-", "_",motifs[i]),sep=""),l,sep = "_"))
          freq <- paste(gsub("-", "_",name_motif),"freq",sep = "_")
          updated_df <- get(current_motif) %>% left_join(get(freq), by = gsub("-", "_",name_motif))
          updated_df["Score"] <- updated_df[gsub("-", "_",current_motif_pvalue)] * updated_df$dist_score
          assign(subtype, rbind(get(subtype),subset(updated_df, Score > 0.05)[,c('Array','Score')]))
          assign(current_motif,updated_df)
        }
      }
    }
  }
}

# IA_new <- agg(IA)
# colnames(IA_new) <- c("Array","IA")
IB_new <- agg(IB)
colnames(IB_new) <- c("Array","IB")
IC_new <- agg(IC)
colnames(IC_new) <- c("Array","IC")
ID_new <- agg(ID)
colnames(ID_new) <- c("Array","ID")
IE_new <- agg(IE)
colnames(IE_new) <- c("Array","IE")
IF_new <- agg(IF)
colnames(IF_new) <- c("Array","IF")
IG_new <- agg(IG)
colnames(IG_new) <- c("Array","IG")
IIA_new <- agg(IIA)
colnames(IIA_new) <- c("Array","IIA")
IIC_new <- agg(IIC)
colnames(IIC_new) <- c("Array","IIC")
IIIA_new <- agg(IIIA)
colnames(IIIA_new) <- c("Array","IIIA")
IIIB_new <- agg(IIIB)
colnames(IIIB_new) <- c("Array","IIIB")
IIID_new <- agg(IIID)
colnames(IIID_new) <- c("Array","IIID")
VK_new <- agg(VK)
colnames(VK_new) <- c("Array","VK")

new_df <- rbind.fill(IB_new,IC_new,ID_new,IE_new,IF_new,IG_new,IIA_new,IIC_new,IIIA_new,IIIB_new,IIID_new,VK_new)
new_df[is.na(new_df)] = 0
collapsed_new <- (aggregate(.~Array,new_df,sum))
collapsed_new$Final_prediction <- as.character(apply(collapsed_new[,-1],1,function(x) names(collapsed_new[,-1])[c(which(x>0.05),which.max(x))[duplicated(c(which(x>0.05),which.max(x)))]]))
collapsed_new[collapsed_new==0]<-NA
collapsed_new[collapsed_new=="character(0)"]<-NA

score_cutoffs <- read.delim("~/Desktop/MotifPrediction_09292023/subtype_annotation_score_cutoffs.txt",header=FALSE)
colnames(score_cutoffs) <- c("subtype","cutoff")
rownames(score_cutoffs) <- score_cutoffs$subtype

for (i in 1:nrow(collapsed_new)){
  if(collapsed_new[i,][collapsed_new$Final_prediction[i]] >= score_cutoffs[collapsed_new$Final_prediction[i],]$cutoff){
    collapsed_new$Final_prediction_after_score[i] <- collapsed_new[i,]$Final_prediction
  } else {
    collapsed_new$Final_prediction_after_score[i] <- "N/A"
  }
}

sum_count <- count(collapsed_new$Final_prediction_after_score)
colnames(sum_count) <- c('Subtype','Count')
final_prediction_Control3 <- df %>% left_join(collapsed_new,by = c("Array"))
colnames(final_prediction_Control3)[ncol(final_prediction_Control3)] <- "Final_prediction_3"

#####Combine controls #####

if(any(grepl("package:matrixStats", search()))) detach("package:matrixStats") else message("matrixStats not loaded")

df_controls <- final_prediction_Control1 %>% left_join(final_prediction_Control2, by="Array")
df_controls <- df_controls %>% left_join(final_prediction_Control3, by="Array")

df_controls %>%select(Array,Final_prediction_1, Final_prediction_2, Final_prediction_3)

df_control_count <- data.frame(matrix(ncol = 2,nrow=12))

colnames(df_control_count) <- c("Type", "Subtype")

df_control_count$Type <- c("1Type I","1Type I","1Type I","1Type I","1Type I","1Type I","1Type III","1Type III","1Type III","2Type II","2Type II","2Type V")
df_control_count$Subtype <- c("IB","IC","ID","IE","IF","IG","IIIA","IIIB","IIID","IIA","IIC","VK")

Final_prediction_1 <- count(df_controls$Final_prediction_1)
colnames(Final_prediction_1) <- c("Subtype","False_positive1")
df_control_count <- df_control_count %>% left_join(Final_prediction_1, by="Subtype")

Final_prediction_2 <- count(df_controls$Final_prediction_2)
colnames(Final_prediction_2) <- c("Subtype","False_positive2")
df_control_count <- df_control_count %>% left_join(Final_prediction_2, by="Subtype")

Final_prediction_3 <- count(df_controls$Final_prediction_3)
colnames(Final_prediction_3) <- c("Subtype","False_positive3")
df_control_count <- df_control_count %>% left_join(Final_prediction_3, by="Subtype")

df_control_count$False_positive1 <- df_control_count$False_positive1*100/1000
df_control_count$False_positive2 <- df_control_count$False_positive2*100/1000
df_control_count$False_positive3 <- df_control_count$False_positive3*100/1000


df_control_count$Ave <- (df_control_count$False_positive1+df_control_count$False_positive2+df_control_count$False_positive3)/3
library(matrixStats)
df_control_count$Std = rowSds(as.matrix(df_control_count[,c(3,4,5)]))
if(any(grepl("package:matrixStats", search()))) detach("package:matrixStats") else message("matrixStats not loaded")

df_control_count[is.na(df_control_count)] <- 0

false_annotations <- ggplot(df_control_count, aes(x=Subtype, y=Ave)) +
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Ave, ymax=Ave+Std), width=.2,
                position=position_dodge(.9)) + ylim(0,5) +facet_grid(cols = vars(Type), scales = "free_x", space = "free_x") + scale_color_brewer("Set3")+ labs(x="Motif ID",y="Percentage %") + theme(text = element_text(size=14,face = "bold"))
false_annotations
