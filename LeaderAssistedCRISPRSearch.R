
# MSA2genemap
# 
# Author: Murat Buyukyoruk
# Affilliated Lab: Wiedenheft Lab

# Description:
# This script is developed to generate leader motif-assisted annotations, including training data, developed scoring matrix, figures, controls and enhanced annotation runs. 

### Laod packages

# Following packages will be required through out this Rscript. If you get error regarding to missing a function although the package is loaded completely, I would recommend to start over without jumping steps.

library("ape")
library("Biostrings")
library("ggplot2")
library("ggtree")
library('stringr')
library("ggnewscale")
library("tidyverse")
library("ggstance")
library("plyr")
library("scales")
library("reshape2")
library("ggseqlogo")
library("gridExtra")
library("ggpubr")
if(any(grepl("package:ggpubr", search()))) detach("package:ggpubr") else message("ggpubr not loaded")
if(any(grepl("package:matrixStats", search()))) detach("package:matrixStats") else message("matrixStats not loaded")

### Color list for phylogenetic tree

color_list <- c("#94b86c","#8b60e2","#45c558","#9b40b5","#7cc64d","#4354ca","#a5bb36","#5f7af3","#4c9d2f","#cd6ee4","#45c27f","#d540ac","#3c863b","#e13187","#42c9b8","#da402f","#458ee8","#e99929","#4267c6","#cbb43f","#6257b8","#789432","#7747a9","#9e8d28","#a375d8","#56741f","#e375cc","#346a34","#ed5286","#75c699","#d22f59","#3bbac6","#dc3748","#51a57a","#ab3f93","#3b8554","#af356f","#32856b","#e96b34","#51afdf","#b64810","#9c8fea","#e1a84c","#704993","#bc7722","#709bdc","#b03a29","#3077ab","#ea8a55","#535a99","#8f6e21","#7a7bbc","#5f6211","#bf80c6","#56642b","#e66ca2","#1a6447","#e3666d","#808948","#925595","#bdb26f","#b5a8e7","#755820","#de96c7","#8e4822","#ba6b90","#a18049","#91456a","#dca371","#aa3c50","#aa6c3e","#e08896","#ca6545","#914d5a","#f29180","#a85044","#cc7c65","#8b6055","#ae9856","#f67ead","#f6a069","#7c0b64","#eb03b7","#130d3c","#6d1704","#a046bf","#b2529f")

#### Define some functions###
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

logo_gen <- function(my_data) {
  for (m in 1:length(meme)){
    if (m==1){
      w <- 0
      anc <- 0
    }
    if(str_detect( meme[m],'ALPHABET')){
      alp <-  strsplit(str_split_fixed(meme[m]," ",2)[2], NULL)[[1]]
      meme_df <- as.numeric()
    }
    if(str_detect( meme[m],'MOTIF')){
      meme_motif_no <- (str_split_fixed(meme[m]," ",3)[2])
    }
    if (str_detect(meme[m],"letter-probability matrix")){
      wide <- str_split_fixed(str_split_fixed(meme[m],"w= ",2)[2]," ",2)[1]
      w <- 1
      anc <- m
    }
    if(m==(anc+w)){
      if(meme_motif_no==motif_no){
        arr <- (str_split_fixed((meme[m])," ",2))[2]
        meme_val <-(as.numeric(str_split_fixed(arr,"  ",4)))
        meme_df <- append(meme_df,meme_val)
      }
      if(w!=wide){
        w <- w+1
      }
    }
  }
  cs1 = make_col_scheme(chars=alp, cols=c('#CB2026', '#35459C', '#FCB316', '#0C8140'))
  custom_mat = matrix( meme_df, nrow=4, dimnames=list(alp))
  ggseqlogo( custom_mat,col_scheme=cs1 ,font="roboto_bold") + theme(text=element_text(size=18),axis.text.x = element_text(angle = 90, vjust = 0.5),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),axis.ticks = element_line(size=0.5))+scale_x_continuous(expand=c(0,0),labels=as.character(1:ncol(custom_mat)),breaks=(1:ncol(custom_mat))) +ylim(0,2)
}

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

empty <- ggplot() +
  geom_point(aes(1,1), colour="white") +
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  )

### Generate Dataframe with CRISPRDetect Default and motif data and perform a training for the N/A annotation###

setwd("Database/fimo_default_dataframe/") # set directory to the fimo_default_dataframe folder provided in the Database files.

for (i in 1: length(list.files())){
  assign(paste("dd",i,sep=""),read.delim(list.files()[i],header = TRUE))
}

df <- read.delim('CRISPRDetect_Cas_default_summary_dataframe.txt', header = TRUE) # Load CRISPRDetect_Cas_default_summary_dataframe.txt file provided in the Database files.

for (i in 1: length(list.files())){
  print(i)
  df <- df %>% left_join(get(paste("dd",i,sep="")), by = "Array")
}

df <- df[,colSums(is.na(df))<nrow(df)]

### Construct 16S phylogeny with the generated dataframe###

# Recommended to change work directory with setwd command to avoid writing new files to fimo_default_dataframe folder.

tree <- read.tree("tree_bac_arc_16S.nexus") # Load tree_bac_arc_16S.nexus file provided in the Database files.
f <- grep("NC_013849.1_Array_1",tree$tip.label) # root to archaeal genome
outgroup <- tree$tip.label[f]
tree <- root(tree, outgroup)
tree.data <- as_tibble(tree)

df_phylum <- data.frame(df$phylum)
rownames(df_phylum) <- tree$tip.label
colnames(df_phylum) <- 'Taxonomy'
df_kingdom <- data.frame(df$kingdom)
rownames(df_kingdom) <- tree$tip.label
colnames(df_kingdom) <- 'Kingdom'
df_subtype <- data.frame(df$Subtype_single)
rownames(df_subtype) <- tree$tip.label
colnames(df_subtype) <- 'Subtype'
p <- ggtree(tree, branch.length = "none", size=0.25)+theme_tree2(legend.position=c(.05, .85),legend.key.size = unit(1,'mm'),)

p <- rotate(p,69932)
p <- rotate(p,78126)
p <- rotate(p,67565)
p <- rotate(p,73616)
p <- rotate(p,67107)
p <- rotate(p,73617)
p <- rotate(p,77912)

### Visualize taxa info on 16S tree
p1 <- gheatmap(p, df_phylum, offset=10, width=0.05, hjust = 0, font.size = 3, color = NULL) + new_scale_fill()
p1 <- gheatmap(p1, df_kingdom, offset=20, width=0.05, hjust = 0, font.size = 3, color = NULL)

### Visualize CRISPR Subtype (repeat- and cas-based) info on 16S tree
p2 <- gheatmap(p, subset(df_subtype,Subtype=='IA'), offset=20, width=0.05, hjust = 0, font.size = 3, color = NULL)
p2 <- gheatmap(p2, subset(df_subtype,Subtype=='IB'), offset=30, width=0.05, hjust = 0, font.size = 3, color = NULL)
p2 <- gheatmap(p2, subset(df_subtype,Subtype=='IC'), offset=40, width=0.05, hjust = 0, font.size = 3, color = NULL)
p2 <- gheatmap(p2, subset(df_subtype,Subtype=='ID'), offset=50, width=0.05, hjust = 0, font.size = 3, color = NULL)
p2 <- gheatmap(p2, subset(df_subtype,Subtype=='IE'), offset=60, width=0.05, hjust = 0, font.size = 3, color = NULL)
p2 <- gheatmap(p2, subset(df_subtype,Subtype=='IF'), offset=70, width=0.05, hjust = 0, font.size = 3, color = NULL)

p3 <- gheatmap(p, subset(df_subtype,Subtype=='IIIA'), offset=20, width=0.05, hjust = 0, font.size = 3, color = NULL)
p3 <- gheatmap(p3, subset(df_subtype,Subtype=='IIIB'), offset=30, width=0.05, hjust = 0, font.size = 3, color = NULL)
p3 <- gheatmap(p3, subset(df_subtype,Subtype=='IIIC'), offset=40, width=0.05, hjust = 0, font.size = 3, color = NULL)
p3 <- gheatmap(p3, subset(df_subtype,Subtype=='IIID'), offset=50, width=0.05, hjust = 0, font.size = 3, color = NULL)
p3 <- gheatmap(p3, subset(df_subtype,Subtype=='IIA'), offset=60, width=0.05, hjust = 0, font.size = 3, color = NULL)
p3 <- gheatmap(p3, subset(df_subtype,Subtype=='IIC'), offset=70, width=0.05, hjust = 0, font.size = 3, color = NULL)
p3 <- gheatmap(p3, subset(df_subtype,Subtype=='VA'), offset=80, width=0.05, hjust = 0, font.size = 3, color = NULL)
p3 <- gheatmap(p3, subset(df_subtype,Subtype=='VE'), offset=90, width=0.05, hjust = 0, font.size = 3, color = NULL)

# IA <- data.frame() # no motif data available
IB <- data.frame()
IC <- data.frame()
ID <- data.frame()
IE <- data.frame()
IF <- data.frame()
IG <- data.frame()
IIA <- data.frame()
# IIB<- data.frame() # no motif data available
IIC<- data.frame()
IIIA<- data.frame()
IIIB<- data.frame()
IIID<- data.frame()
# IVA<- data.frame() # no motif data available
# VB<- data.frame() # no motif data available
VK<- data.frame()

### Following part is using motif data and known leaders for the data training, maps motifs to the 16S tree, and expands CRISPR annotation by attempting to annotate previously unannotated CRISPRs.

motifs_data <- read.delim("motif_list_w_narrow_distribution.txt",header=FALSE) # Load motif_list_w_narrow_distribution.txt file provided in the Database files.

fileConn<-file("/tmp/tree_cmd_2023.txt")
line_list <- list()
motifs <- as.character()

for (i in 1:nrow(motifs_data)){
  print(i)
  motifs <- append(motifs,as.character(motifs_data[i,1]))
}

prev_sub <- ""
p1 <- p
p4 <- p1

for (i in 1:length(motifs)) {
  temp_list <- list()
  motif_list <-  (colnames(df)[grepl(gsub("-", ".",motifs[i]),colnames(df))])
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
      temp_df <- data.frame(Array=df$Array,name_motif=df[motif_fetch],current_motif_pvalue=df[motif_pval],Subtype_single=df$Subtype_single,Phylum=df$phylum,Subtype_CRISPRDetect=df$Subtype_CRISPRDetect,Repeat_occurence = df$Repeat_occurence,Redundant=df$Redundant)
      colnames(temp_df) <- c('Array',gsub("-", "_",name_motif),gsub("-", "_",current_motif_pvalue),'Subtype_single','Phylum','Subtype_CRISPRDetect','Repeat_occurence', 'Redundant')
      assign(current_motif,temp_df)
      subtype <- str_split_fixed(motif,"_",3)[1]
      taxa <- str_split_fixed(motifs[i],"_",3)[2]
      if (nrow(subset(get(current_motif),Subtype_single==subtype & Phylum == taxa )[gsub("-", "_",name_motif)])!=0){
        temp_list <- rbind(temp_list,unlist(subset(get(current_motif),Subtype_single==subtype & Phylum == taxa & Redundant == "Non_redundant")[gsub("-", "_",name_motif)]))
      }
      if(l==(length(motif_list)/2)) {
        if(subtype!=prev_sub){
          cmd <- paste("p4 <- p4 + geom_facet(panel = 'Type ", subtype, "', data = subset(",current_motif,",Subtype_single=='", subtype ,"' ), geom = geom_col, aes(x = 200,alpha=0.25), orientation = 'y') + geom_vline(xintercept = 0)\n " , sep="")
          line_list <- append(line_list,(noquote(cmd)))
          prev_sub <- subtype
        }
        current <- (subset(get(current_motif),Subtype_single==subtype & Phylum == taxa & Redundant == "Non_redundant" ))
        size <- nrow(current)
        temp_list <- do.call(rbind.data.frame,temp_list)
        colnames(temp_list) <- "Position"
        freq <- count(temp_list,"Position")
        colnames(freq) <- c(gsub("-", "_",name_motif),"perc")
        freq[nrow(freq),][2] <- NA
        freq$dist_score <- freq$perc / max(head(freq$perc,-1))
        for (l in 1:(length(motif_list)/2)) {
          current_motif <- (paste(paste("df_",gsub("-", "_",motifs[i]),sep=""),l,sep = "_"))
          updated_df <- get(current_motif) %>% left_join(freq, by = gsub("-", "_",name_motif))
          updated_df["Score"] <- updated_df[gsub("-", "_",current_motif_pvalue)] * updated_df$dist_score
          assign(subtype, rbind(get(subtype),subset(updated_df,Subtype_single== "N/A" & Score > 0.05)[,c('Array','Score')]))
          assign(current_motif,updated_df)
        }
        for (l in 1:(length(motif_list)/2)) {
          current_motif <- (paste(paste("df_",gsub("-", "_",motifs[i]),sep=""),l,sep = "_"))
          plot_df <- get(current_motif)
          colnames(plot_df) <- c('Array',"motif_pos",gsub("-", "_",current_motif_pvalue),'Subtype_single','Phylum','Subtype_CRISPRDetect','perc','dist_score','Score')
          tryCatch({
            print(current_motif)
            cmd <- paste("p4 <- p4 + geom_facet(panel = 'Type ", subtype, "', data = subset(",current_motif,",Subtype_single=='", subtype ,"' & Phylum == '", taxa, "' ), geom = geom_point, aes(stroke=0.25,x = ",noquote(paste(gsub("-", "_",name_motif))),",color='", gsub("eobacteria","",gsub("_"," ",gsub("motif","m",str_split_fixed(paste(gsub("-", "_",name_motif)),"_",2)[2]))), "',size=Score,alpha=Score), stat = 'identity') + scale_size(range = c(0.000000000001,1)) + xlim_expand(c(200,0), panel = 'Type ", subtype,"') + theme(panel.background = element_rect(fill = '#E5E4E2',colour = '#E5E4E2'))\n" , sep="")
            line_list <- append(line_list,(noquote(cmd)))
            if (nrow(subset(plot_df,Score > 0.05 & Subtype_single=="N/A"))!=0){
              cmd <- paste("p4 <- p4 + geom_facet(panel = 'Type ", subtype, "', data = subset(",current_motif,",Score > 0.05 & Subtype_single=='N/A' ), geom = geom_point, aes(stroke=0.25,x = ",noquote(paste(gsub("-", "_",name_motif))),",color='", gsub("eobacteria","",gsub("_"," ",gsub("motif","m",str_split_fixed(paste(gsub("-", "_",name_motif)),"_",2)[2]))), "',size=Score,alpha=Score), stat = 'identity') + scale_size(range = c(0.000000000001,1)) + xlim_expand(c(200,0), panel = 'Type ", subtype,"') + theme(panel.background = element_rect(fill = '#E5E4E2',colour = '#E5E4E2'))\n" , sep="")
              line_list <- append(line_list,(noquote(cmd)))
            }
          }, error=function(e){})
        }
        current_freq <- paste(gsub("-", "_",name_motif),"freq",sep = "_")
        assign(paste("df_",gsub("-", "_",motifs[i]),sep=""),temp_list)
        assign(current_freq,freq)
        data <- get(gsub("-", "_",name_motif))
        data <- data.frame(x=data[! is.na(data)])
        colnames(data) <- "x"
        data$title <- gsub("_"," ",gsub("motif","motif",str_split_fixed(paste(gsub("-", "_",name_motif)),"_",2)[2]))
        bar <- data %>% ggplot(aes(x=x, y=..density..*max(count(data$x)[2])/max(..density..))) + geom_histogram( binwidth=1) + xlim(200, 0) + theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=0.5),text=element_text(size=18),legend.position = "none", plot.margin = margin(0, 0, 0, 1.5, "pt"))  + ylab("Count") + xlab("Distance from Leader-Repeat Junction (bp)") + geom_vline(aes(xintercept = 0),linetype="dashed") + geom_hline(aes(yintercept = 0)) + annotate("text", x = 5, y = max(count(data$x)$freq)/2 , label = "LRJ",size=5,angle=90)+scale_y_continuous(breaks= pretty_breaks())
        bar <- bar + facet_grid(.~title)
        motif_file <- paste(rev(strsplit(str_split_fixed(paste(rev(strsplit(motifs[i], NULL)[[1]]), collapse=""),"_",2)[2], NULL)[[1]]), collapse="")
        motif_no <- paste(rev(strsplit(str_split_fixed(paste(rev(strsplit(motifs[i], NULL)[[1]]), collapse=""),"_",2)[1],NULL)[[1]]),collapse="")
        logo_path <- paste(paste("meme_data/",motif_file,sep = "/"),"txt",sep = ".") # Requires meme data provided in the Database files.
        meme<-readLines(logo_path)
        logo_raw <- logo_gen(meme)
        logo <- grid.arrange(empty,logo_raw,empty,nrow=3,heights=c(0.125,1,0.125))
        assign(paste("sup",i,sep = ""),grid.arrange(bar,logo,ncol=2,widths=c(1,1.5)))
      }
    }
  }
  cmd <-  paste("p",subtype," <- p4\n", sep="")
  line_list <- append(line_list,(noquote(cmd)))
  if(motifs[i]=="VK_Cyanobacteriota_motif_3"){
    line_list <- append(line_list,"p_non_I <- p4\n")
  }
}

### OPTIONAL - for loop below is for exporting motif distance plot and PWM logos. Comment in to use and define destination to save files.

# for (i in 1: length(motifs)){
#   if(i%%7==0){
#     setEPS(); Sys.sleep(1); postscript(file=paste0(paste("~/",format(Sys.time(), "%Y%m%d_%H%M%S"),sep="/"),".eps"),width=12,height=17); (grid.arrange(get(paste("sup",i-6,sep = "")), get(paste("sup",i-5,sep = "")),get(paste("sup",i-4,sep = "")),get(paste("sup",i-3,sep = "")),get(paste("sup",i-2,sep = "")),get(paste("sup",i-1,sep = "")),get(paste("sup",i,sep = "")), nrow = 7));dev.off()
#     Sys.sleep(1)
#   }
# }
# setEPS(); Sys.sleep(1); postscript(file=paste0(paste("~/",format(Sys.time(), "%Y%m%d_%H%M%S"),sep="/"),".eps"),width=12,height=17); (grid.arrange(sup169, nrow = 7));dev.off()

### write lines generated from previous for loop to execute them to generate phylogeneti tree with PWM midpoints illustered for each subtype

writeLines(paste(unlist(line_list),collapse=""),fileConn)
source(fileConn)
close(fileConn)
p4 + scale_fill_manual(values=color_list) + scale_color_manual(values=color_list) +theme(legend.position = "none")+ scale_x_reverse()

### Assign subtype annotations to NA CRISPRs based on Leader-assisted scoring matrix

# IA_new <- agg(IA)  # no motif data available
# colnames(IA_new) <- c("Array","IA")  # no motif data available
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
# IIIC_new <- agg(IIIC)  # no motif data available
# colnames(IIIC_new) <- c("Array","IIIC")  # no motif data available
IIID_new <- agg(IIID)
colnames(IIID_new) <- c("Array","IIID")
VK_new <- agg(VK)
colnames(VK_new) <- c("Array","VK")
# VE_new <- agg(VE)  # no motif data available
# colnames(VE_new) <- c("Array","VE")  # no motif data available

new_df <- rbind.fill(IB_new,IC_new,ID_new,IE_new,IF_new,IG_new,IIA_new,IIC_new,IIIA_new,IIIB_new,IIID_new,VK_new)
new_df[is.na(new_df)] = 0
merged_new <- (aggregate(.~Array,new_df,sum))
merged_new$Final_prediction <- as.character(apply(merged_new[,-1],1,function(x) names(merged_new[,-1])[c(which(x>0),which.max(x))[duplicated(c(which(x>0),which.max(x)))]]))
merged_new[merged_new==0]<-NA
score_cutoffs <- read.delim("subtype_annotation_score_cutoffs.txt",header=FALSE) # Load subtype_annotation_score_cutoffs.txt file provided in the Database files.
colnames(score_cutoffs) <- c("subtype","cutoff")
rownames(score_cutoffs) <- score_cutoffs$subtype

for (i in 1:nrow(merged_new)){
  if(merged_new[i,][merged_new$Final_prediction[i]] >= score_cutoffs[merged_new$Final_prediction[i],]$cutoff){
    merged_new$Final_prediction_after_score[i] <- merged_new[i,]$Final_prediction
  } else {
    merged_new$Final_prediction_after_score[i] <- "N/A"
  }
}

sum_count <- count(merged_new$Final_prediction_after_score)
colnames(sum_count) <- c('Subtype','Count')

### Count comparison after applying leader-assisted CRISPR annotation

# IA_new <- fun_int0(subset(sum_count,Subtype=='IA')$Count)  # no motif data available
IB_new <- fun_int0(subset(sum_count,Subtype=='IB')$Count)
IC_new <- fun_int0(subset(sum_count,Subtype=='IC')$Count)
ID_new <- fun_int0(subset(sum_count,Subtype=='ID')$Count)
IE_new <- fun_int0(subset(sum_count,Subtype=='IE')$Count)
IF_new <- fun_int0(subset(sum_count,Subtype=='IF')$Count)
IG_new <- fun_int0(subset(sum_count,Subtype=='IG')$Count)
IIA_new <- fun_int0(subset(sum_count,Subtype=='IIA')$Count)
IIC_new <- fun_int0(subset(sum_count,Subtype=='IIC')$Count)
IIIA_new <- fun_int0(subset(sum_count,Subtype=='IIIA')$Count)
IIIB_new <- fun_int0(subset(sum_count,Subtype=='IIIB')$Count)
#IIIC_new <- fun_int0(subset(sum_count,Subtype=='IIIC')$Count)  # no motif data available
IIID_new <- fun_int0(subset(sum_count,Subtype=='IIID')$Count)
VK_new <- fun_int0(subset(sum_count,Subtype=='VK')$Count)
# VA_new <- fun_int0(subset(sum_count,Subtype=='VA')$Count)  # no motif data available
# VE_new <- fun_int0(subset(sum_count,Subtype=='VE')$Count)  # no motif data available

# IA_old <- nrow(subset(df,Subtype_CRISPRDetect=='IA'))  # no motif data available
IB_old <- nrow(subset(df,Subtype_CRISPRDetect=='IB'))
IC_old <- nrow(subset(df,Subtype_CRISPRDetect=='IC'))
ID_old <- nrow(subset(df,Subtype_CRISPRDetect=='ID'))
IE_old <- nrow(subset(df,Subtype_CRISPRDetect=='IE'))
IF_old <- nrow(subset(df,Subtype_CRISPRDetect=='IF'))
IG_old <- nrow(subset(df,Subtype_CRISPRDetect=='IG'))
IIIA_old <- nrow(subset(df,Subtype_CRISPRDetect=='IIIA'))
IIIB_old <- nrow(subset(df,Subtype_CRISPRDetect=='IIIB'))
# IIIC_old <- nrow(subset(df,Subtype_CRISPRDetect=='IIIC'))  # no motif data available
IIID_old <- nrow(subset(df,Subtype_CRISPRDetect=='IIID'))
IIA_old <- nrow(subset(df,Subtype_CRISPRDetect=='IIA'))
IIC_old <- nrow(subset(df,Subtype_CRISPRDetect=='IIC'))
VK_old <- nrow(subset(df,Subtype_CRISPRDetect=='VK'))
# VA_old <- nrow(subset(df,Subtype_CRISPRDetect=='VA'))  # no motif data available
# VE_old <- nrow(subset(df,Subtype_CRISPRDetect=='VE'))  # no motif data available

NA_old <- nrow(subset(df,Subtype_CRISPRDetect=='N/A'))

# IA_mid <- nrow(subset(df,Subtype_single=='IA'))  # no motif data available
IB_mid <- nrow(subset(df,Subtype_single=='IB'))
IC_mid <- nrow(subset(df,Subtype_single=='IC'))
ID_mid <- nrow(subset(df,Subtype_single=='ID'))
IE_mid <- nrow(subset(df,Subtype_single=='IE'))
IF_mid <- nrow(subset(df,Subtype_single=='IF'))
IG_mid <- nrow(subset(df,Subtype_single=='IG'))
IIIA_mid <- nrow(subset(df,Subtype_single=='IIIA'))
IIIB_mid <- nrow(subset(df,Subtype_single=='IIIB'))
# IIIC_mid <- nrow(subset(df,Subtype_single=='IIIC'))  # no motif data available
IIID_mid <- nrow(subset(df,Subtype_single=='IIID'))
IIA_mid <- nrow(subset(df,Subtype_single=='IIA'))
IIC_mid <- nrow(subset(df,Subtype_single=='IIC'))
VK_mid <- nrow(subset(df,Subtype_single=='VK'))
# VA_mid <- nrow(subset(df,Subtype_single=='VA'))  # no motif data available
# VE_mid <- nrow(subset(df,Subtype_single=='VE'))  # no motif data available

NA_mid <- nrow(subset(df,Subtype_single=='N/A'))

bardata <- data.frame("Subtype"=c("IB","IC","ID","IE","IF","IG","IIIA","IIIB","IIID","IIA","IIC","VK","Annotated"),"Repat-based"=c(IB_old,IC_old,ID_old,IE_old,IF_old,IG_old,IIIA_old,IIIB_old,IIID_old,IIA_old,IIC_old,VK_old,38026-NA_old),"Repeat + cas-based"=c(IB_mid,IC_mid,ID_mid,IE_mid,IF_mid,IG_mid,IIIA_mid,IIIB_mid,IIID_mid,IIA_mid,IIC_mid,VK_mid,38026-NA_mid),"Repat + cas + leader motif-based"=c(IB_mid+IB_new,IC_mid+IC_new,ID_mid+ID_new,IE_mid+IE_new,IF_mid+IF_new,IG_mid+IG_new,IIIA_mid+IIIA_new,IIIB_mid+IIIB_new,IIID_mid+IIID_new,IIA_mid+IIA_new,IIC_mid+IIC_new,VK_mid+VK_new,38026-(NA_mid-IB_new-IC_new-ID_new-IE_new-IF_new-IG_new-IIIA_new-IIIB_new-IIID_new-IIA_new-IIC_new-VK_new)),"Type"=c("1Type I","1Type I","1Type I","1Type I","1Type I","1Type I","1Type III","1Type III","1Type III","2Type II","2Type II","2Type V","1Annotated"))
bardata.long<-melt(bardata)
bardata.long %>%
arrange(desc(value)) %>%
ggplot(aes(Subtype,value,fill=variable)) + geom_bar(stat="identity",position = "dodge") +facet_grid(cols = vars(Type), scales = "free_x", space = "free_x")+geom_text(aes(label = value),position=position_dodge(width=1), vjust = -1,hjust=-0.2, size = 3.5,angle=45)+scale_fill_manual(values=c("#86677F","#97A1AB","#ECDB90")) + labs(x="Subtype",y="Count",fill="Classification Method") +theme(legend.position = "bottom",text=element_text(size = 15,face="bold"),legend.box.background = element_rect(colour = "black",size=1)) #+ scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 50),breaks = 10^(0:5))

### OPTIONAL - Line below is for exporting the  plot reporting the number of annotations between different methods. Comment in to use and define destination to save files.

# setEPS(); Sys.sleep(1); postscript(file=paste0(paste("~/",format(Sys.time(), "%Y%m%d_%H%M%S"),sep="/"),".eps"),width=8,height=5);ggplot(bardata.long,aes(Subtype,value,fill=variable)) + geom_bar(stat="identity",position="dodge") + ylim(0,10500)+facet_grid(cols = vars(Type), scales = "free_x", space = "free_x")+geom_text(aes(label = value),position=position_dodge(width=1), vjust = -1,hjust=-0.2, size = 2.5,angle=45)+scale_fill_manual(values=c("#86677F","#97A1AB","#ECDB90")) + labs(x="Subtype",y="Count",fill="Classification Method") +theme(legend.position = "bottom",text=element_text(size = 12,face="bold"),legend.box.background = element_rect(colour = "black",size=1))	 ;dev.off()

### Update dataframe with final predictions

final_prediction_default <- df %>% left_join(merged_new,by = c("Array"))
final_prediction_default$Subtype_CD_HMM_MAA <- NA

for (i in 1:nrow(final_prediction_default)) {
  if (is.na(final_prediction_default$Final_prediction_after_score[i])){
    final_prediction_default$Subtype_CD_HMM_MAA[i] <- final_prediction_default$Subtype_single[i]
  } else {
    final_prediction_default$Subtype_CD_HMM_MAA[i] <- final_prediction_default$Final_prediction_after_score[i]
  }
}

df_default <- final_prediction_default # Contains raw data, file size is too big to handle if exported and opened with Excel

### Record numbers of annotated CRISPRs after implementing leader motif-assisted method

# IA_old_edit <- IA_mid + IA_new # no motif data available
IB_old_edit <- IB_mid + IB_new
IC_old_edit <- IC_mid + IC_new
ID_old_edit <- ID_mid + ID_new
IE_old_edit <- IE_mid + IE_new
IF_old_edit <- IF_mid + IF_new
IG_old_edit <- IG_mid + IG_new
IIIA_old_edit <- IIIA_mid + IIIA_new
IIIB_old_edit <- IIIB_mid + IIIB_new
# IIIC_old_edit <- IIIC_mid + IIIC_new # no motif data available
IIID_old_edit <- IIID_mid + IIID_new
IIA_old_edit <- IIA_mid + IIA_new
IIC_old_edit <- IIC_mid + IIC_new
VK_old_edit <- VK_mid + VK_new
# VA_old_edit <- VA_mid + VA_new # no motif data available
# VE_old_edit <- VE_mid + VE_new # no motif data available

# IA_default_dr_hmm_maa <- IA_new # no motif data available
IB_default_dr_hmm_maa <- IB_new
IC_default_dr_hmm_maa <- IC_new
ID_default_dr_hmm_maa <- ID_new
IE_default_dr_hmm_maa <- IE_new
IF_default_dr_hmm_maa <- IF_new
IG_default_dr_hmm_maa <- IG_new
IIIA_default_dr_hmm_maa <- IIIA_new
IIIB_default_dr_hmm_maa <- IIIB_new
# IIIC_default_dr_hmm_maa <- IIIC_new # no motif data available
IIID_default_dr_hmm_maa <- IIID_new
IIA_default_dr_hmm_maa <- IIA_new
IIC_default_dr_hmm_maa <- IIC_new
VK_default_dr_hmm_maa <- VK_new
# VA_default_dr_hmm_maa <- VA_new # no motif data available
# VE_default_dr_hmm_maa <- VE_new # no motif data available

#########################################################################################

### Comparison of annatations (agreee/disaggree) between different annotation methods

    ### Leader motif-only annotations vs CRISPRDetect annotations

# Source Rscript file provided in the github repository. You might need to update the path to the script.
source("Leader_accuracy_CD_2024.R") # ignore error messages, it is trying all possibilities
accuracy_CD

    ### Leader motif-only annotations vs CRISPRDetect & cas-proximity (10kb) annotations

# Source Rscript file provided in the github repository. You might need to update the path to the script.
source("Leader_accuracy_CD_HMM_2024.R") # ignore error messages, it is trying all possibilities
accuracy_CD_HMM

#########################################################################################

### What if all CRISPR arrays were NA? How does leader motif-only annotation performs in terms of numbers?

# IA <- data.frame() # no motif data available
IB <- data.frame()
IC <- data.frame()
ID <- data.frame()
IE <- data.frame()
IF <- data.frame()
IG <- data.frame()
IIA <- data.frame()
# IIB<- data.frame() # no motif data available
IIC<- data.frame()
IIIA<- data.frame()
IIIB<- data.frame()
# IIIC<- data.frame() # no motif data available
IIID <- data.frame()
VK<- data.frame()
# VE<- data.frame() # no motif data available
# VIB1<- data.frame() # no motif data available

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
      temp_df <- data.frame(Array=df$Array,name_motif=df[motif_fetch],current_motif_pvalue=df[motif_pval],Phylum=df$phylum)
      colnames(temp_df) <- c('Array',gsub("-", "_",name_motif),gsub("-", "_",current_motif_pvalue),'Phylum')
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

# IA_new <- agg(IA) # no motif data available
# colnames(IA_new) <- c("Array","IA") # no motif data available
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
# IIIC_new <- agg(IIIC) # no motif data available
# colnames(IIIC_new) <- c("Array","IIIC") # no motif data available
IIID_new <- agg(IIID)
colnames(IIID_new) <- c("Array","IIID")
VK_new <- agg(VK)
colnames(VK_new) <- c("Array","VK")
# VE_new <- agg(VE) # no motif data available
# colnames(VE_new) <- c("Array","VE") # no motif data available

new_df <- rbind.fill(IB_new,IC_new,ID_new,IE_new,IF_new,IG_new,IIA_new,IIC_new,IIIA_new,IIIB_new,IIID_new,VK_new)

new_df[is.na(new_df)] = 0
merged_new <- (aggregate(.~Array,new_df,sum))
merged_new$Final_prediction <- as.character(apply(merged_new[,-1],1,function(x) names(merged_new[,-1])[c(which(x>0),which.max(x))[duplicated(c(which(x>0),which.max(x)))]]))
merged_new[merged_new==0]<-NA

for (i in 1:nrow(merged_new)){
  if(merged_new[i,][merged_new$Final_prediction[i]] >= score_cutoffs[merged_new$Final_prediction[i],]$cutoff){
    merged_new$Final_prediction_after_score[i] <- merged_new[i,]$Final_prediction
  } else {
    merged_new$Final_prediction_after_score[i] <- "N/A"
  }
}

sum_count <- count(merged_new$Final_prediction_after_score)
colnames(sum_count) <- c('Subtype','Count')

# IA_new <- fun_int0(subset(sum_count,Subtype=='IA')$Count) # no motif data available
IB_new <- fun_int0(subset(sum_count,Subtype=='IB')$Count)
IC_new <- fun_int0(subset(sum_count,Subtype=='IC')$Count)
ID_new <- fun_int0(subset(sum_count,Subtype=='ID')$Count)
IE_new <- fun_int0(subset(sum_count,Subtype=='IE')$Count)
IF_new <- fun_int0(subset(sum_count,Subtype=='IF')$Count)
IG_new <- fun_int0(subset(sum_count,Subtype=='IG')$Count)
IIA_new <- fun_int0(subset(sum_count,Subtype=='IIA')$Count)
IIC_new <- fun_int0(subset(sum_count,Subtype=='IIC')$Count)
IIIA_new <- fun_int0(subset(sum_count,Subtype=='IIIA')$Count)
IIIB_new <- fun_int0(subset(sum_count,Subtype=='IIIB')$Count)
# IIIC_new <- fun_int0(subset(sum_count,Subtype=='IIIC')$Count) # no motif data available
IIID_new <- fun_int0(subset(sum_count,Subtype=='IIID')$Count)
VK_new <- fun_int0(subset(sum_count,Subtype=='VK')$Count)
# VE_new <- fun_int0(subset(sum_count,Subtype=='VE')$Count) # no motif data available

bardata_MAA_only <- data.frame("Subtype"=c("IB","IC","ID","IE","IF","IG","IIIA","IIIB","IIID","IIA","IIC","VK","NA"),"Repeat-based"=c(IB_old,IC_old,ID_old,IE_old,IF_old,IG_old,IIIA_old,IIIB_old,IIID_old,IIA_old,IIC_old,VK_old,NA_old),"Repeat + cas-based"=c(IB_mid,IC_mid,ID_mid,IE_mid,IF_mid,IG_mid,IIIA_mid,IIIB_mid,IIID_mid,IIA_mid,IIC_mid,VK_mid,NA_mid),"Leader motif-based"=c(IB_new,IC_new,ID_new,IE_new,IF_new,IG_new,IIIA_new,IIIB_new,IIID_new,IIA_new,IIC_new,VK_new,38026-IB_new-IC_new-ID_new-IE_new-IF_new-IG_new-IIIA_new-IIIB_new-IIID_new-IIA_new-IIC_new-VK_new),"Type"=c("1 Type I","1 Type I","1 Type I","1 Type I","1 Type I","1 Type I","1 Type III","1 Type III","1 Type III","2 Type II","2 Type II","2 Type V","1 NA"))
bardata_MAA_only.long<-melt(bardata_MAA_only)
ggplot(bardata_MAA_only.long,aes(Subtype,value,fill=variable)) + geom_bar(stat="identity",position="dodge") + ylim(0,15000)+facet_grid(cols = vars(Type), scales = "free_x", space = "free_x")+geom_text(aes(label = value),position=position_dodge(width=1), vjust = -1,hjust=-0.2, size = 3.5,angle=45)+scale_fill_manual(values=c("#86677F","#97A1AB","#ECDB90")) + labs(x="Subtype",y="Count",fill="Classification Method") +theme(legend.position = "bottom",text=element_text(size = 15,face="bold"),legend.box.background = element_rect(colour = "black",size=1))

bardata_fig1_MAA_only <- data.frame("Subtype"=c("IB","IC","ID","IE","IF","IG","IIIA","IIIB","IIID","IIA","IIC","VK","NA"),"Repeat-based"=c(IB_old,IC_old,ID_old,IE_old,IF_old,IG_old,IIIA_old,IIIB_old,IIID_old,IIA_old,IIC_old,VK_old,38026-NA_old),"Repeat + cas-based"=c(IB_mid,IC_mid,ID_mid,IE_mid,IF_mid,IG_mid,IIIA_mid,IIIB_mid,IIID_mid,IIA_mid,IIC_mid,VK_mid,38026-NA_mid),"Leader motif-based"=c(IB_new,IC_new,ID_new,IE_new,IF_new,IG_new,IIIA_new,IIIB_new,IIID_new,IIA_new,IIC_new,VK_new,0-(-IB_new-IC_new-ID_new-IE_new-IF_new-IG_new-IIIA_new-IIIB_new-IIID_new-IIA_new-IIC_new-VK_new)),"Repat + cas + leader motif-based"=c(IB_mid+IB_default_dr_hmm_maa,IC_mid+IC_default_dr_hmm_maa,ID_mid+ID_default_dr_hmm_maa,IE_mid+IE_default_dr_hmm_maa,IF_mid+IF_default_dr_hmm_maa,IG_mid+IG_default_dr_hmm_maa,IIIA_mid+IIIA_default_dr_hmm_maa,IIIB_mid+IIIB_default_dr_hmm_maa,IIID_mid+IIID_default_dr_hmm_maa,IIA_mid+IIA_default_dr_hmm_maa,IIC_mid+IIC_default_dr_hmm_maa,VK_mid+VK_default_dr_hmm_maa,38026-(NA_mid-IB_default_dr_hmm_maa-IC_default_dr_hmm_maa-ID_default_dr_hmm_maa-IE_default_dr_hmm_maa-IF_default_dr_hmm_maa-IG_default_dr_hmm_maa-IIIA_default_dr_hmm_maa-IIIB_default_dr_hmm_maa-IIID_default_dr_hmm_maa-IIA_default_dr_hmm_maa-IIC_default_dr_hmm_maa-VK_default_dr_hmm_maa)),"Type"=c("1 Type I","1 Type I","1 Type I","1 Type I","1 Type I","1 Type I","1 Type III","1 Type III","1 Type III","2 Type II","2 Type II","2 Type V","1 NA"))
bardata_fig1_MAA_only.long<-melt(bardata_fig1_MAA_only)
bardata_fig1_MAA_only.long %>%
  arrange(desc(value)) %>%
  ggplot(aes(Subtype,value,fill=variable)) + geom_bar(stat="identity",position = "identity") + ylim(0,30000) +facet_grid(cols = vars(Type), scales = "free_x", space = "free_x")+geom_text(aes(label = value),position=position_dodge(width=1), vjust = -1,hjust=-0.2, size = 3.5,angle=45)+scale_fill_manual(values=c("#86677F","#97A1AB","#ECDB90","#c6763f")) + labs(x="Subtype",y="Count",fill="Classification Method") +theme(legend.position = "bottom",text=element_text(size = 15,face="bold"),legend.box.background = element_rect(colour = "black",size=1))

final_prediction_default_MAA_only <- df %>% left_join(merged_new,by = c("Array"))
final_prediction_default_MAA_only <- df %>% left_join(merged_new,by = c("Array"))
final_prediction_default_MAA_only$Subtype_CD_HMM_MAA <- NA

for (i in 1:nrow(final_prediction_default_MAA_only)) {
  if (is.na(final_prediction_default_MAA_only$Final_prediction_after_score[i])){
    final_prediction_default_MAA_only$Subtype_CD_HMM_MAA[i] <- final_prediction_default_MAA_only$Subtype_single[i]
  } else {
    final_prediction_default_MAA_only$Subtype_CD_HMM_MAA[i] <- final_prediction_default_MAA_only$Final_prediction_after_score[i]
  }
}

df_default_MAA <- final_prediction_default_MAA_only # Contains raw data, file size is too big to handle if exported and opened with Excel

### CAUTION - END of the part for the data generated using default parameter of CRISPRDetect. Some variables from previous part will be overwritten if continued.

#########################################################################################

### Positive/Negative Control scripts to test the efficiency of trained data

# Source Rscript file provided in the github repository. You might need to update the path to the script.
source("Leader_neg_control_2024.R")
false_hits
false_annotations

#########################################################################################

### CAUTION - START of the part for the data generated using relaxed parameters of CRISPRDetect. Some variables from previous part will be overwritten if continued.

setwd("Database/fimo_relaxed_dataframe/") # set directory to the fimo_relaxed_dataframe folder provided in the Database files.

for (i in 1: length(list.files())){
  assign(paste("de",i,sep=""),read.delim(list.files()[i],header = TRUE))
}

df <- read.delim('CRISPRDetect_Cas_relaxed_summary_dataframe.txt', header = TRUE) # Load CRISPRDetect_Cas_relaxed_summary_dataframe.txt file provided in the Database files.

for (i in 1: length(list.files())){
  df <- df %>% left_join(get(paste("de",i,sep="")), by = "Array")
  print(i)
}

df <- df[,colSums(is.na(df))<nrow(df)]

# Recommended to change work directory with setwd command to avoid writing new files to fimo_relaxed_dataframe folder.

# IA <- data.frame() # no motif data available
IB <- data.frame()
IC <- data.frame()
ID <- data.frame()
IE <- data.frame()
IF <- data.frame()
IG <- data.frame()
IIA <- data.frame()
# IIB<- data.frame() # no motif data available
IIC<- data.frame()
IIIA<- data.frame()
IIIB<- data.frame()
# IIIC<- data.frame() # no motif data available
IIID <- data.frame()
VK<- data.frame()
# VE<- data.frame() # no motif data available
# VIB1<- data.frame() # no motif data available

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
      temp_df <- data.frame(Array=df$Array,name_motif=df[motif_fetch],current_motif_pvalue=df[motif_pval],Subtype_single=df$Subtype_single,Phylum=df$phylum,Subtype_CRISPRDetect=df$Subtype_CRISPRDetect)
      colnames(temp_df) <- c('Array',gsub("-", "_",name_motif),gsub("-", "_",current_motif_pvalue),'Subtype_single','Phylum','Subtype_CRISPRDetect')
      assign(current_motif,temp_df)
      subtype <- str_split_fixed(motif,"_",3)[1]
      taxa <- str_split_fixed(motifs[i],"_",3)[2]
      if(l==(length(motif_list)/2)) {
        for (l in 1:(length(motif_list)/2)) {
          current_motif <- (paste(paste("df_",gsub("-", "_",motifs[i]),sep=""),l,sep = "_"))
          freq <- paste(gsub("-", "_",name_motif),"freq",sep = "_")
          updated_df <- get(current_motif) %>% left_join(get(freq), by = gsub("-", "_",name_motif))
          updated_df["Score"] <- updated_df[gsub("-", "_",current_motif_pvalue)] * updated_df$dist_score
          assign(subtype, rbind(get(subtype),subset(updated_df, Subtype_single== "N/A" & Score > 0.05)[,c('Array','Score')]))
          assign(current_motif,updated_df)
        }
      }
    }
  }
}

# IA_new <- agg(IA) # no motif data available
# colnames(IA_new) <- c("Array","IA") # no motif data available
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
# IIIC_new <- agg(IIIC) # no motif data available
# colnames(IIIC_new) <- c("Array","IIIC") # no motif data available
IIID_new <- agg(IIID)
colnames(IIID_new) <- c("Array","IIID")
VK_new <- agg(VK)
colnames(VK_new) <- c("Array","VK")
# VE_new <- agg(VE) # no motif data available
# colnames(VE_new) <- c("Array","VE") # no motif data available

new_df <- rbind.fill(IB_new,IC_new,ID_new,IE_new,IF_new,IG_new,IIA_new,IIC_new,IIIA_new,IIIB_new,IIID_new,VK_new)
new_df[is.na(new_df)] = 0
merged_new <- (aggregate(.~Array,new_df,sum))
merged_new$Final_prediction <- as.character(apply(merged_new[,-1],1,function(x) names(merged_new[,-1])[c(which(x>0),which.max(x))[duplicated(c(which(x>0),which.max(x)))]]))
merged_new[merged_new==0]<-NA
merged_new[merged_new=="character(0)"]<-NA

for (i in 1:nrow(merged_new)){
  if(merged_new[i,][merged_new$Final_prediction[i]] >= score_cutoffs[merged_new$Final_prediction[i],]$cutoff){
    merged_new$Final_prediction_after_score[i] <- merged_new[i,]$Final_prediction
  } else {
    merged_new$Final_prediction_after_score[i] <- "N/A"
  }
}

sum_count <- count(merged_new$Final_prediction_after_score)
colnames(sum_count) <- c('Subtype','Count')

df_final <- df %>% left_join(merged_new,by = c("Array"), suffix = c("", "_df2"))

# IA_new <- fun_int0(subset(sum_count,Subtype=='IA')$Count)+ nrow(subset(df,Subtype_single=='IA')) # no motif data available
IB_new <- fun_int0(subset(sum_count,Subtype=='IB')$Count)+ nrow(subset(df,Subtype_single=='IB'))
IC_new <- fun_int0(subset(sum_count,Subtype=='IC')$Count)+ nrow(subset(df,Subtype_single=='IC'))
ID_new <- fun_int0(subset(sum_count,Subtype=='ID')$Count)+ nrow(subset(df,Subtype_single=='ID'))
IE_new <- fun_int0(subset(sum_count,Subtype=='IE')$Count)+ nrow(subset(df,Subtype_single=='IE'))
IF_new <- fun_int0(subset(sum_count,Subtype=='IF')$Count)+ nrow(subset(df,Subtype_single=='IF'))
IG_new <- fun_int0(subset(sum_count,Subtype=='IG')$Count)+ nrow(subset(df,Subtype_single=='IG'))
IIA_new <- fun_int0(subset(sum_count,Subtype=='IIA')$Count)+ nrow(subset(df,Subtype_single=='IIA'))
IIC_new <- fun_int0(subset(sum_count,Subtype=='IIC')$Count)+ nrow(subset(df,Subtype_single=='IIC'))
IIIA_new <- fun_int0(subset(sum_count,Subtype=='IIIA')$Count)+ nrow(subset(df,Subtype_single=='IIIA'))
IIIB_new <- fun_int0(subset(sum_count,Subtype=='IIIB')$Count)+ nrow(subset(df,Subtype_single=='IIIB'))
# IIIC_new <- fun_int0(subset(sum_count,Subtype=='IIIC')$Count)+ nrow(subset(df,Subtype_single=='IIIC')) # no motif data available
IIID_new <- fun_int0(subset(sum_count,Subtype=='IIID')$Count)+ nrow(subset(df,Subtype_single=='IIID'))
# VA_new <- fun_int0(subset(sum_count,Subtype=='VA')$Count)+ nrow(subset(df,Subtype_single=='VA')) # no motif data available
# VE_new <- fun_int0(subset(sum_count,Subtype=='VE')$Count)+ nrow(subset(df,Subtype_single=='VE')) # no motif data available
VK_new <- fun_int0(subset(sum_count,Subtype=='VK')$Count)+ nrow(subset(df,Subtype_single=='VK'))

bardata <- data.frame("Subtype"=c("IB","IC","ID","IE","IF","IG","IIIA","IIIB","IIID","IIA","IIC","VK"),"Default Array Detection + Leader Motif-assisted CRISPR Annotation"=c(IB_old_edit,IC_old_edit,ID_old_edit,IE_old_edit,IF_old_edit,IG_old_edit,IIIA_old_edit,IIIB_old_edit,IIID_old_edit,IIA_old_edit,IIC_old_edit,VK_old_edit),"w/ Short-Degenerate Array Detection + Leader Motif-assisted CRISPR Annotation"=c(IB_new,IC_new,ID_new,IE_new,IF_new,IG_new,IIIA_new,IIIB_new,IIID_new,IIA_new,IIC_new,VK_new),"Type"=c("1Type I","1Type I","1Type I","1Type I","1Type I","1Type I","1Type III","1Type III","1Type III","2Type II","2Type II","2Type V"))

bardata.long<-melt(bardata)
ggplot(bardata.long,aes(Subtype,value,fill=variable)) + geom_bar(stat="identity",position="dodge",colour="black",size=0.25) +facet_grid(cols = vars(Type), scales = "free_x", space = "free_x")+geom_text(aes(label = value,angle = 45),position=position_dodge(width=1), vjust = -1,hjust=0.15, size = 3.5)+ labs(x="Subtype",y="Count",fill = "") +theme(legend.position = "top",text=element_text(size = 15,face="bold"),legend.box.background = element_rect(colour = "black",size=1),legend.spacing.y = unit(0.1, "cm"),legend.key.size = unit(0.4, "cm"),legend.text = element_text(size=12))+guides(fill=guide_legend(nrow=2, byrow=TRUE))

final_prediction_edit <- df %>% left_join(merged_new,by = c("Array"))
final_prediction_edit$Subtype_CD_HMM_MAA <- NA

for (i in 1:nrow(final_prediction_edit)) {
  if (is.na(final_prediction_edit$Final_prediction_after_score[i])){
    final_prediction_edit$Subtype_CD_HMM_MAA[i] <- final_prediction_edit$Subtype_single[i]
  } else {
    final_prediction_edit$Subtype_CD_HMM_MAA[i] <- final_prediction_edit$Final_prediction_after_score[i]
  }
}

df_edit <- final_prediction_edit  # Contains raw data, file size is too big to handle if exported and opened with Excel

### Properties comparison ###

df_default<-subset(df_default,CRISPR_presence=="CRISPR" & Subtype_CD_HMM_MAA != "N/A")
df_default$Type <- NA

for (i in 1:nrow(df_default)) {       # for-loop over columns
  if(df_default$Subtype_CD_HMM_MAA[i] == "IG" || df_default$Subtype_CD_HMM_MAA[i] == "IB"|| df_default$Subtype_CD_HMM_MAA[i] == "IC"|| df_default$Subtype_CD_HMM_MAA[i] == "ID"|| df_default$Subtype_CD_HMM_MAA[i] == "IE"|| df_default$Subtype_CD_HMM_MAA[i] == "IF"){
    df_default$Type[i] <- "1Type I"
  }
  if(df_default$Subtype_CD_HMM_MAA[i] == "IIA" || df_default$Subtype_CD_HMM_MAA[i] == "IIC"){
    df_default$Type[i] <- "2Type II"
  }
  if(df_default$Subtype_CD_HMM_MAA[i] == "IIIA" || df_default$Subtype_CD_HMM_MAA[i] == "IIIB"||  df_default$Subtype_CD_HMM_MAA[i] == "IIID"){
    df_default$Type[i] <- "1Type III"
  }
  if(df_default$Subtype_CD_HMM_MAA[i] == "VK" || df_default$Subtype_CD_HMM_MAA[i] == "VE"){
    df_default$Type[i] <- "2Type V"
  }
  if(is.na(df_default$Subtype_CD_HMM_MAA[i])){
    print(df_default$Array)
  }
}

df_edit<-subset(df_edit,Subtype_CD_HMM_MAA != "N/A")
df_edit$Type <- NA

for (i in 1:nrow(df_edit)) {       # for-loop over columns
  if(df_edit$Subtype_CD_HMM_MAA[i] == "IG" || df_edit$Subtype_CD_HMM_MAA[i] == "IB"|| df_edit$Subtype_CD_HMM_MAA[i] == "IC"|| df_edit$Subtype_CD_HMM_MAA[i] == "ID"|| df_edit$Subtype_CD_HMM_MAA[i] == "IE"|| df_edit$Subtype_CD_HMM_MAA[i] == "IF"){
    df_edit$Type[i] <- "1Type I"
  }
  if(df_edit$Subtype_CD_HMM_MAA[i] == "IIA" || df_edit$Subtype_CD_HMM_MAA[i] == "IIC"){
    df_edit$Type[i] <- "2Type II"
  }
  if(df_edit$Subtype_CD_HMM_MAA[i] == "IIIA" || df_edit$Subtype_CD_HMM_MAA[i] == "IIIB"||  df_edit$Subtype_CD_HMM_MAA[i] == "IIID"){
    df_edit$Type[i] <- "1Type III"
  }
  if(df_edit$Subtype_CD_HMM_MAA[i] == "VK" || df_edit$Subtype_CD_HMM_MAA[i] == "VE"){
    df_edit$Type[i] <- "2Type V"
  }
  if(is.na(df_edit$Subtype_CD_HMM_MAA[i])){
    print(df_edit$Array)
  }
}

# Calcluate spacer lengths
df_default$Spacer <- NA
df_edit$Spacer <- NA

df_default$Status <- NA
df_edit$Status <- NA

for (i in 1:nrow(df_default)) {
  df_default$Spacer[i] <- round(((df_default$Stop[i]-df_default$Start[i])-(df_default$Repeat_occurence[i]*df_default$Repeat_length[i]))/(df_default$Repeat_occurence[i]-1))
}
for (i in 1:nrow(df_edit)) {
  df_edit$Spacer[i] <- round(((df_edit$Stop[i]-as.numeric(df_edit$Start[i]))-(df_edit$Repeat_occurence[i]*nchar(df_edit$Repeat[i])))/(df_edit$Repeat_occurence[i]-1))
}
for (i in 1:nrow(df_default)) {
  df_default$Status[i] <- "Default Array Detection + Leader Motif-assisted CRISPR Annotation"
}
for (i in 1:nrow(df_edit)) {
  df_edit$Status[i] <- "w/ Short-Degenerate Array Detection + Leader Motif-assisted CRISPR Annotation"
}
for (i in 1:nrow(df_edit)) {
  df_edit$Repeat_length[i] <- nchar(df_edit$Repeat[i])
}

df_merged <- rbind(df_default[,c('Array','Subtype_CD_HMM_MAA','Repeat_occurence','Repeat_length','Spacer','Status','Type','Subtype_CRISPRDetect')],df_edit[,c('Array','Subtype_CD_HMM_MAA','Repeat_occurence','Repeat_length','Spacer',"Status",'Type','Subtype_CRISPRDetect')])

pp<- ggplot(data = subset(df_merged,Type!="NA"&Spacer>15), aes(x=Subtype_CD_HMM_MAA, y=Repeat_occurence,fill=Status))  +facet_grid(cols = vars(Type), scales = "free_x", space = "free_x") + geom_boxplot(position="dodge",outlier.size = 0.5,outlier.stroke = 0.5 ) + labs(x="Subtype",y="Array Length (Number of Repeats)",fill = "") +theme(legend.position = "top",text=element_text(size = 12,face="bold"),legend.box.background = element_rect(colour = "black",size=1),legend.spacing.y = unit(0.1, "cm"),legend.key.size = unit(0.4, "cm"),legend.text = element_text(size=12),legend.margin=margin(t = c(0.1,0.1,0.1,0.1), unit='cm'))+guides(fill=guide_legend(nrow=2, byrow=TRUE)) #+geom_jitter(position="dodge",aes(color=Status))
pp1<- ggplot(data = subset(df_merged,Type!="NA"&Spacer>15), aes(x=Subtype_CD_HMM_MAA, y=Repeat_length,fill=Status))  +facet_grid(cols = vars(Type), scales = "free_x", space = "free_x") + geom_boxplot(position="dodge",outlier.size = 0.5,outlier.stroke = 0.5) + labs(x="Subtype",y="Repeat Length (bp)",fill = "") +theme(legend.position = "top",text=element_text(size = 12,face="bold"),legend.box.background = element_rect(colour = "black",size=1),legend.spacing.y = unit(0.1, "cm"),legend.key.size = unit(0.4, "cm"),legend.text = element_text(size=12),legend.margin=margin(t = c(0.1,0.1,0.1,0.1), unit='cm'))+guides(fill=guide_legend(nrow=2, byrow=TRUE)) #+ geom_point(pch = 21, size=1 ,stroke = 0,position = position_jitterdodge())#+geom_jitter(position="dodge",aes(color=Status))
pp2<- ggplot(data = subset(df_merged,Type!="NA"&Spacer>15), aes(x=Subtype_CD_HMM_MAA, y=Spacer,fill=Status))  +facet_grid(cols = vars(Type), scales = "free_x", space = "free_x") + geom_boxplot(position="dodge",outlier.size = 0.5,outlier.stroke = 0.5) + labs(x="Subtype",y="Spacer Length (bp)",fill = "") +theme(legend.position = "top",text=element_text(size = 12,face="bold"),legend.box.background = element_rect(colour = "black",size=1),legend.spacing.y = unit(0.1, "cm"),legend.key.size = unit(0.4, "cm"),legend.text = element_text(size=12),legend.margin=margin(t = c(0.1,0.1,0.1,0.1), unit='cm'))+guides(fill=guide_legend(nrow=2, byrow=TRUE))#+ geom_point(pch = 21, size=1 ,stroke = 0,position = position_jitterdodge()) #+geom_jitter(position="dodge",aes(color=Status))

library(cowplot)
library(egg)
library(ggpubr)

phist <- ggplot(subset(df_merged,Type!="N/A"&Spacer>15), aes(x = Repeat_occurence, fill = Status ,y=100*after_stat(count)/length((subset(df_merged,Status=="w/ Short-Degenerate Array Detection + Leader Motif-assisted CRISPR Annotation")$Repeat_occurence)))) + geom_histogram(position="dodge",binwidth = 0.1,alpha=0) +scale_x_continuous(name = "Array Length (Number of Repeats)") + scale_y_continuous(name = "Percentage (%)") +theme(legend.position = "none",text=element_text(size = 14,face="bold"))
pdensity <- ggplot(subset(df_merged,Type!="N/A"&Spacer>15), aes(x = Repeat_occurence, fill = Status)) + geom_density(adjust=0.1,alpha=1)+scale_y_continuous(position = "right") + theme_half_open(11, rel_small = 1) +theme(legend.position = "none",text=element_text(size = 14,face="bold")) + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")#+rremove("legend")
aligned_plots <- align_plots(phist, pdensity, align="hv", axis="tblr")
pd1 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

phist_zoom <- ggplot(subset(df_merged,Type!="NA"&Spacer>15), aes(x = Repeat_occurence, fill = Status ,y=100*after_stat(count)/length((subset(df_merged,Status=="w/ Short-Degenerate Array Detection + Leader Motif-assisted CRISPR Annotation")$Repeat_occurence)))) + geom_histogram(position="dodge",binwidth = 0.5) +scale_x_continuous(name = "Array Length (Number of Repeats)",limits=c(0,50)) + scale_y_continuous(name = "Percentage (%)") +theme(legend.position = "none",text=element_text(size = 14,face="bold"))
pdensity_zoom <- ggplot(subset(df_merged,Type!="NA"&Spacer>15), aes(x = Repeat_occurence, fill = Status)) + geom_density(binwidth = 0.5,adjust=0.4,alpha=1)+scale_y_continuous(position = "right") + theme_half_open(11, rel_small = 1) +theme(legend.position = "none",text=element_text(size = 14,face="bold"))+ rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks") +xlim(0,50) #+rremove("legend")
aligned_plots_zoom <- align_plots(phist_zoom, pdensity_zoom, align="hv", axis="tblr")
pd1_zoom <- ggdraw(aligned_plots_zoom[[1]]) + draw_plot(aligned_plots_zoom[[2]])

pzoom <- empty + annotation_custom(ggplotGrob(pd1_zoom), xmin = 1, xmax = 4.5, ymin = 3000, ymax = 20000) +
  geom_rect(aes(xmin = 1, xmax = 4.5, ymin = 3000, ymax = 20000), color='black', linetype='dashed', alpha=0) +
  geom_rect(aes(xmin = 1, xmax = 4.75, ymin = 3000, ymax = 20000), color='black', size = 0, alpha=0) +
  geom_rect(aes(xmin = 0.1, xmax = 0.74, ymin = 1500, ymax = 20000), color='black', linetype='dashed', alpha=0) +
  geom_path(aes(x,y,group=grp),
            data=data.frame(x = c(1,0.74,1,0.74), y=c(3000,1500,20000,20000),grp=c(1,1,2,2)),
            linetype='dashed')

pd1_merge <- pd1 + draw_plot(pzoom)

phist <- ggplot(subset(df_merged,Type!="NA"&Spacer>15), aes(x = Repeat_length, fill = Status, y=100*after_stat(count)/length((subset(df_merged,Status=="w/ Short-Degenerate Array Detection + Leader Motif-assisted CRISPR Annotation")$Repeat_length)))) + geom_histogram(position="dodge",binwidth = 0.1,alpha=0)+scale_x_continuous(name = "Repeat Length (bp)",limits=c(0,60))  + scale_y_continuous(name = "Percentage (%)")+theme(legend.position = "none",text=element_text(size = 14,face="bold"))
pdensity <- ggplot(subset(df_merged,Type!="NA"&Spacer>15), aes(x = Repeat_length, fill = Status)) + geom_density(binwidth = 0.5,adjust=2,alpha=1)+scale_y_continuous(position = "right") +xlim(0,60) + theme_half_open(11, rel_small = 1) +theme(legend.position = "none",text=element_text(size = 14,face="bold"))+ rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")#+rremove("legend")
aligned_plots <- align_plots(phist, pdensity, align="hv", axis="tblr")
pd2 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

phist <- ggplot(subset(df_merged,Type!="NA"&Spacer>15), aes(x = Spacer, fill = Status ,y=100*after_stat(count)/length((subset(df_merged,Status=="w/ Short-Degenerate Array Detection + Leader Motif-assisted CRISPR Annotation")$Spacer)))) + geom_histogram(position="dodge",binwidth = 0.1,alpha=0) +scale_x_continuous(name = "Spacer Length (bp)",limits=c(0,60))+ scale_y_continuous(name = "Percentage (%)")+theme(legend.position = "none",text=element_text(size = 14,face="bold"))
pdensity <- ggplot(subset(df_merged,Type!="NA"&Spacer>15), aes(x = Spacer, fill = Status)) + geom_density(binwidth = 0.5,adjust=2,alpha=1)+scale_y_continuous(position = "right") +xlim(0,60) + theme_half_open(11, rel_small = 1) +theme(legend.position = "none",text=element_text(size = 14,face="bold"))+ rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")#+rremove("legend")
aligned_plots <- align_plots(phist, pdensity, align="hv", axis="tblr")
pd3 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

if(any(grepl("package:egg", search()))) detach("package:egg") else message("egg not loaded")

ggarrange(pd1_merge, pd2, pd3, nrow=3,labels=c("A","B","C"),common.legend = TRUE, legend="top")
ggarrange(pp, pp1, pp2, nrow=3,labels=c("A","B","C"),common.legend = TRUE, legend="top",align = "v")

# Distributions for the rest of the arrays remained NA
df_default <- final_prediction_default
df_edit <- final_prediction_edit

df_default<-subset(df_default,CRISPR_presence=="CRISPR" & Subtype_CD_HMM_MAA == "N/A")
df_default$Type <- NA
df_edit<-subset(df_edit,Subtype_CD_HMM_MAA == "N/A")
df_edit$Type <- NA

for (i in 1:nrow(df_default)) {       # for-loop over columns
  if(df_default$Subtype_CD_HMM_MAA[i] == "N/A" ){
    df_default$Type[i] <- "N/A"
  }
  if(is.na(df_default$Subtype_CD_HMM_MAA[i])){
    print(df_default$Array)
  }
}

for (i in 1:nrow(df_edit)) {       # for-loop over columns
  if(df_edit$Subtype_CD_HMM_MAA[i] == "N/A" ){
    df_edit$Type[i] <- "N/A"
  }
  if(is.na(df_edit$Subtype_CD_HMM_MAA[i])){
    print(df_edit$Array)
  }
}

df_default$Spacer <- NA
df_edit$Spacer <- NA

df_default$Status <- NA
df_edit$Status <- NA

for (i in 1:nrow(df_default)) {
  df_default$Spacer[i] <- round(((df_default$Stop[i]-df_default$Start[i])-(df_default$Repeat_occurence[i]*df_default$Repeat_length[i]))/(df_default$Repeat_occurence[i]-1))
}
for (i in 1:nrow(df_edit)) {
  df_edit$Spacer[i] <- round(((df_edit$Stop[i]-as.numeric(df_edit$Start[i]))-(df_edit$Repeat_occurence[i]*nchar(df_edit$Repeat[i])))/(df_edit$Repeat_occurence[i]-1))
}
for (i in 1:nrow(df_default)) {
  df_default$Status[i] <- "Default Array Detection + Leader Motif-assisted CRISPR Annotation"
}
for (i in 1:nrow(df_edit)) {
  df_edit$Status[i] <- "w/ Short-Degenerate Array Detection + Leader Motif-assisted CRISPR Annotation"
}
for (i in 1:nrow(df_edit)) {
  df_edit$Repeat_length[i] <- nchar(df_edit$Repeat[i])
}

df_merged <- rbind(df_default[,c('Array','Subtype_CD_HMM_MAA','Repeat_occurence','Repeat_length','Spacer','Status','Type','Subtype_CRISPRDetect')],df_edit[,c('Array','Subtype_CD_HMM_MAA','Repeat_occurence','Repeat_length','Spacer',"Status",'Type','Subtype_CRISPRDetect')])

phist <- ggplot(subset(df_merged,Type=="N/A"), aes(x = Repeat_occurence, fill = Status ,y=100*after_stat(count)/length((subset(df_merged,Status=="w/ Short-Degenerate Array Detection + Leader Motif-assisted CRISPR Annotation")$Repeat_occurence)))) + geom_histogram(position="dodge",binwidth = 0.1,alpha=0) +scale_x_continuous(name = "Array Length (Number of Repeats)") + scale_y_continuous(name = "Percentage (%)") +theme(legend.position = "none",text=element_text(size = 14,face="bold"))
pdensity <- ggplot(subset(df_merged,Type=="N/A"), aes(x = Repeat_occurence, fill = Status)) + geom_density(adjust=0.1,alpha=1)+scale_y_continuous(position = "right") + theme_half_open(11, rel_small = 1) +theme(legend.position = "none",text=element_text(size = 14,face="bold")) + rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")#+rremove("legend")
aligned_plots <- align_plots(phist, pdensity, align="hv", axis="tblr")
pd1 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

phist_zoom <- ggplot(subset(df_merged,Type=="N/A"), aes(x = Repeat_occurence, fill = Status ,y=100*after_stat(count)/length((subset(df_merged,Status=="w/ Short-Degenerate Array Detection + Leader Motif-assisted CRISPR Annotation")$Repeat_occurence)))) + geom_histogram(position="dodge",binwidth = 0.5,alpha=0) +scale_x_continuous(name = "Array Length (Number of Repeats)",limits=c(0,50)) + scale_y_continuous(name = "Percentage (%)") +theme(legend.position = "none",text=element_text(size = 14,face="bold"))
pdensity_zoom <- ggplot(subset(df_merged,Type=="N/A"), aes(x = Repeat_occurence, fill = Status)) + geom_density(adjust=1,alpha=1)+scale_y_continuous(position = "right") + theme_half_open(11, rel_small = 1) +theme(legend.position = "none",text=element_text(size = 14,face="bold"))+ rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks") +xlim(0,50) #+rremove("legend")
aligned_plots_zoom <- align_plots(phist_zoom, pdensity_zoom, align="hv", axis="tblr")
pd1_zoom <- ggdraw(aligned_plots_zoom[[1]]) + draw_plot(aligned_plots_zoom[[2]])

pzoom <- empty + annotation_custom(ggplotGrob(pd1_zoom), xmin = 1, xmax = 4.5, ymin = 3000, ymax = 20000) +
  geom_rect(aes(xmin = 1, xmax = 4.5, ymin = 3000, ymax = 20000), color='black', linetype='dashed', alpha=0) +
  geom_rect(aes(xmin = 1, xmax = 4.75, ymin = 3000, ymax = 20000), color='black', size = 0, alpha=0) +
  geom_rect(aes(xmin = 0.1, xmax = 0.74, ymin = 1500, ymax = 20000), color='black', linetype='dashed', alpha=0) +
  geom_path(aes(x,y,group=grp),
            data=data.frame(x = c(1,0.74,1,0.74), y=c(3000,1500,20000,20000),grp=c(1,1,2,2)),
            linetype='dashed')

pd1_merge <- pd1 + draw_plot(pzoom)

phist <- ggplot(subset(df_merged,Type=="N/A"), aes(x = Repeat_length, fill = Status, y=100*after_stat(count)/length((subset(df_merged,Status=="w/ Short-Degenerate Array Detection + Leader Motif-assisted CRISPR Annotation")$Repeat_length)))) + geom_histogram(position="dodge",binwidth = 0.1,alpha=0)+scale_x_continuous(name = "Repeat Length (bp)",limits=c(0,60))  + scale_y_continuous(name = "Percentage (%)")+theme(legend.position = "none",text=element_text(size = 14,face="bold"))
pdensity <- ggplot(subset(df_merged,Type=="N/A"), aes(x = Repeat_length, fill = Status)) + geom_density(adjust=0.75,alpha=1)+scale_y_continuous(position = "right") +xlim(0,60) + theme_half_open(11, rel_small = 1) +theme(legend.position = "none",text=element_text(size = 14,face="bold"))+ rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")#+rremove("legend")
aligned_plots <- align_plots(phist, pdensity, align="hv", axis="tblr")
pd2 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])


phist <- ggplot(subset(df_merged,Type=="N/A"), aes(x = Spacer, fill = Status ,y=100*after_stat(count)/length((subset(df_merged,Status=="w/ Short-Degenerate Array Detection + Leader Motif-assisted CRISPR Annotation")$Spacer)))) + geom_histogram(position="dodge",binwidth = 0.1,alpha=0) +scale_x_continuous(name = "Spacer Length (bp)",limits=c(0,60))+ scale_y_continuous(name = "Percentage (%)")+theme(legend.position = "none",text=element_text(size = 14,face="bold"))
pdensity <- ggplot(subset(df_merged,Type=="N/A"), aes(x = Spacer, fill = Status)) + geom_density(adjust=0.75,alpha=1)+scale_y_continuous(position = "right") +xlim(0,60) + theme_half_open(11, rel_small = 1) +theme(legend.position = "none",text=element_text(size = 14,face="bold"))+ rremove("x.axis")+rremove("xlab") +rremove("x.text") +rremove("x.ticks") + rremove("y.axis")+rremove("ylab") +rremove("y.text") +rremove("y.ticks")#+rremove("legend")
aligned_plots <- align_plots(phist, pdensity, align="hv", axis="tblr")
pd3 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

if(any(grepl("package:egg", search()))) detach("package:egg") else message("egg not loaded")

ggarrange(pd1_merge, pd2, pd3, nrow=3,labels=c("A","B","C"),common.legend = TRUE, legend="top")

# End for distribution plots

#########################################################################################

### CAUTION - START of the part for the data generated using relaxed parameters of CRISPRDetect on ICP1 genomes. Some variables from previous part will be overwritten if continued.

setwd("/Database/fimo_ICP1_dataframe/") # set directory to the fimo_ICP1_dataframe folder provided in the Database files.

for (i in 1: length(list.files())){
  assign(paste("de",i,sep=""),read.delim(list.files()[i],header = TRUE))
}

df <- read.delim('ICP1_CRISPRDetect3_dataframe.txt', header = TRUE) # Load ICP1_CRISPRDetect3_dataframe.txt file provided in the Database files.

for (i in 1: length(list.files())){
  df <- df %>% left_join(get(paste("de",i,sep="")), by = "Array")
}

df <- df[,colSums(is.na(df))<nrow(df)]

# Only the following motifs generated hits in this dataset
motifs <- c("IB_Campylobacterota_motif_1","IC_Pseudomonadota_motif_9","IE_Actinomycetota_motif_2","IE_Pseudomonadota_motif_7","IF_Pseudomonadota_motif_3","IF_Pseudomonadota_motif_7","IF_Pseudomonadota_motif_8","IF_Pseudomonadota_motif_9","IF_Pseudomonadota_motif_10","IIIA_Pseudomonadota_motif_1","IIIB_Pseudomonadota_motif_2","VK_Cyanobacteriota_motif_3")

# Recommended to change work directory with setwd command to avoid writing new files to fimo_ICP1_dataframe folder.

IB <- data.frame()
IC <- data.frame()
IE <- data.frame()
IF <- data.frame()
IIIA<- data.frame()
IIIB<- data.frame()
VK<- data.frame()

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

      temp_df <- data.frame(Array=df$Array,name_motif=df[motif_fetch],current_motif_pvalue=df[motif_pval],Subtype_single=df$Subtype,Subtype_CRISPRDetect=df$Subtype)
      colnames(temp_df) <- c('Array',gsub("-", "_",name_motif),gsub("-", "_",current_motif_pvalue),'Subtype_single','Subtype_CRISPRDetect')
      assign(current_motif,temp_df)

      subtype <- str_split_fixed(motif,"_",3)[1]
      taxa <- str_split_fixed(motifs[i],"_",3)[2]

      if(l==(length(motif_list)/2)) {

        for (l in 1:(length(motif_list)/2)) {
          current_motif <- (paste(paste("df_",gsub("-", "_",motifs[i]),sep=""),l,sep = "_"))
          freq <- paste(gsub("-", "_",name_motif),"freq",sep = "_")

          updated_df <- get(current_motif) %>% left_join(get(freq), by = gsub("-", "_",name_motif))
          updated_df["Score"] <- updated_df[gsub("-", "_",current_motif_pvalue)] * updated_df$dist_score

          assign(subtype, rbind(get(subtype),subset(updated_df, Subtype_single== "N/A" & Score > 0)[,c('Array','Score')]))

          assign(current_motif,updated_df)
        }
      }
    }
  }
}

# Only the following motifs generated hits in this dataset
IB_new <- agg(IB)
colnames(IB_new) <- c("Array","IB")
IC_new <- agg(IC)
colnames(IC_new) <- c("Array","IC")
IE_new <- agg(IE)
colnames(IE_new) <- c("Array","IE")
IF_new <- agg(IF)
colnames(IF_new) <- c("Array","IF")
IIIA_new <- agg(IIIA)
colnames(IIIA_new) <- c("Array","IIIA")
IIIB_new <- agg(IIIB)
colnames(IIIB_new) <- c("Array","IIIB")
VK_new <- agg(VK)
colnames(VK_new) <- c("Array","VK")

new_df <- rbind.fill(IB_new,IC_new,IE_new,IF_new,IIIA_new,IIIB_new,VK_new)
new_df[is.na(new_df)] = 0
merged_new <- (aggregate(.~Array,new_df,sum))
merged_new$Final_prediction <- as.character(apply(merged_new[,-1],1,function(x) names(merged_new[,-1])[c(which(x>0),which.max(x))[duplicated(c(which(x>0),which.max(x)))]]))
merged_new[merged_new==0]<-NA
merged_new[merged_new=="character(0)"]<-NA

for (i in 1:nrow(merged_new)){
  if(merged_new[i,][merged_new$Final_prediction[i]] >= score_cutoffs[merged_new$Final_prediction[i],]$cutoff){
    merged_new$Final_prediction_after_score[i] <- merged_new[i,]$Final_prediction
  } else {
    merged_new$Final_prediction_after_score[i] <- "N/A"
  }
}

sum_count <- count(merged_new$Final_prediction_after_score)
colnames(sum_count) <- c('Subtype','Count')

# Final annotations of ICP1 CRISPRs
df_final <- df %>% left_join(merged_new,by = c("Array"), suffix = c("", "_df2"))
