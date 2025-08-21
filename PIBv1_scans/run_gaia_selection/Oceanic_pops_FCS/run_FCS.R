##############################################
#     DATA OCEANIA - Fisher's combine score (FCS)

library(qqman)
library(tidyverse)
library(data.table)
library(ggrepel)

temp <- list.files(pattern="*HMP_pbs-xpehh.tsv", full.names=TRUE)

for(i in 1:length(temp)){
  df1 <- fread(temp[[i]], sep="\t", header=T)
  #df1 <- na.omit(df1)

  name <- gsub(pattern = ".tsv", "", basename(temp[i]))
  name_plot <- gsub(pattern = "_PIBv1_ALL_chr_pbs_closest_final_quantiles_window_stats_HMP_pbs-xpehh.tsv", "", basename(temp[i]))

  df1$log10_HMP <- -log10(df1$HMP)

  df1[is.infinite(df1$log10_HMP),"log10_HMP"] <- 8

  df <- df1 %>%
    rowwise() %>%
    mutate(FCS_X2=-2*sum(log(1-c(quantile_rank,max_quantile_rank_xpehh)))) %>%
    mutate(FCS=pchisq(FCS_X2, df=2*2, lower.tail=FALSE))

  df$log10_FCS <- -log10(df$FCS)

  df[is.infinite(df$log10_FCS),"log10_FCS"] <- 8

  write.table(df, file = paste0(name,"_FCS.tsv"),
              quote=FALSE, sep="\t", row.names = FALSE)

  write.table(na.omit(df[df$FCS < 0.001, ]), file = paste0(name,"_FCS_0.001.tsv"),
              quote=FALSE, sep="\t", row.names = FALSE)

  write.table(na.omit(df[df$FCS<0.0001,]), file = paste0(name,"_FCS_0.0001.tsv"),
              quote=FALSE, sep="\t", row.names = FALSE)

  write.table(na.omit(df[df$FCS<0.00001,]), file = paste0(name,"_FCS_0.00001.tsv"),
              quote=FALSE, sep="\t", row.names = FALSE)

}
