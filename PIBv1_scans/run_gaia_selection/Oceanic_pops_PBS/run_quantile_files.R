##############################################
#     DATA OCEANIA - PBS ALL CHR

library(qqman)
library(tidyverse)
library(data.table)
library(ggrepel)

temp <- list.files(pattern="*final.tsv", full.names=TRUE)

for(i in 1:length(temp)){
  df <- fread(temp[[i]], sep="\t", header=T)
  df <- na.omit(df)

  name <- gsub(pattern = ".tsv", "", basename(temp[i]))
  name_plot <- gsub(pattern = "_PIBv1_ALL_chr_pbs_sort_closest_final.tsv", "", basename(temp[i]))

  q_0.05<-quantile(df$score, probs = c(.9995)) ## 99.95th percentile = 0.05% 

  q_0.01<-quantile(df$score, probs = c(.9999)) ## 99.99th percentile = 0.01%

  q_0.001<-quantile(df$score, probs = c(.99999)) ## 99.999th percentile = 0.001% 


  write.table(df[df$score>q_0.05,], file = paste0(name, "_0.05.tsv"),
              quote=FALSE, sep="\t", row.names = FALSE)


  write.table(df[df$score>q_0.01,], file = paste0(name, "_0.01.tsv"),
              quote=FALSE, sep="\t", row.names = FALSE)

  write.table(df[df$score>q_0.001,], file = paste0(name, "_0.001.tsv"),
              quote=FALSE, sep="\t", row.names = FALSE)
}
