######################################################
#         PIBv1_manuscript - manhattan plot | FCS
# script adapted from:
# https://r-graph-gallery.com/101_Manhattan_plot.html
######################################################

library(qqman)
library(tidyverse)
library(data.table)
library(ggrepel)
library(ggnewscale)


pattern1 <- "FCS.tsv"
pattern2 <- "_0.01_introgressed_corehaps_uniquePOP.tsv"

# list of files
temp <- list.files(pattern = pattern1, full.names = TRUE)
temp2 <- list.files(pattern = pattern2, full.names = TRUE)

# lists to dataframes
data_list <- list()
data_list2 <- list()

# Process the data for the manhattan plot 
for (i in 1:length(temp)) {
  
  name <- gsub(pattern = ".tsv", "_manhattan", basename(temp[i]))
  name_plot <- gsub(pattern = "_PIBv1_ALL_chr_pbs_sort_closest_final_quantiles_window_stats_HMP_pbs-xpehh_FCS.tsv", "", basename(temp[i]))
  
  df <- read_tsv(temp[i]) %>%
    select(chr, start, end, population, FCS, log10_FCS, gene) %>%
    na.omit()
  
  # Merge columns chr and start with separator ":" and start and end with separator ":"
  df$coordinates <- paste(paste(df$chr, df$start, sep = ":"), df$end, sep = "-")
  
  df2 <- read_tsv(temp2[i]) %>%
    select(chr, start, end, population, test, FCS, log10_FCS, gene, TractID)
  
  data_list[[i]] <- df
  data_list2[[i]] <- df2


  # Merge the df and create the 'arch' column
  df <- df %>%
    left_join(df2 %>% select(chr, start, end, TractID), by = c("chr", "start", "end")) %>%
    mutate(arch = case_when(
      grepl("Denisovan", TractID) ~ "Denisovan",
      grepl("Neandertal", TractID) ~ "Neandertal",
      grepl("Ambiguous", TractID) ~ "Ambiguous",
      TRUE ~ NA_character_  # Set other cases to NA
    )) %>%
    select(-TractID)


  df <- df %>%
    mutate(avg_pos = (start + end) / 2)

  # Calculate cumulative position of each chromosome
  fai <- read_tsv("chr_info.tsv", col_names=c("chr","chr_len", "x1", "x2", "x3"),
                  show_col_types = FALSE) %>% 
    # Compute chromosome size
    transmute(chr=chr, tot=cumsum(chr_len)-chr_len)

  don <- df %>% 
    # Add this info to the initial dataset
    left_join(fai, by=c("chr"="chr")) %>%
    # Add a cumulative position of each SNP
    mutate(BPcum=avg_pos+tot)

  axisdf = don %>% group_by(chr) %>% summarize(center=(max(BPcum) + min(BPcum)) / 2)

  max_log_FCS <- don %>%
    group_by(arch, gene) %>%
    slice(which.max(log10_FCS)) %>%
    ungroup()

  # Select the top 10 log10_FCS values per arch
  top_max_log_FCS <-  max_log_FCS %>% #non_na_max_log_FCS %>%
    group_by(arch) %>%
    top_n(10, wt = log10_FCS) %>%
    ungroup() %>%
    filter(chr %in% unique(don$chr))

  bold_genes <- read.table("bold_genes.txt", header = FALSE, col.names = c("gene"))

  top_max_log_FCS_filtered <- top_max_log_FCS %>% 
    filter(!is.na(arch)) %>%
    mutate(bold = gene %in% bold_genes$gene)

  # Adjust size based on whether the gene is in bold or not
  top_max_log_FCS_filtered$label_size <- ifelse(top_max_log_FCS_filtered$bold, 4, 2)

  don <- don %>%
    left_join(top_max_log_FCS %>% select(arch, gene), by = c("arch", "gene")) %>%
    mutate(top_max = ifelse(gene %in% top_max_log_FCS$gene, gene, NA_character_)) 

  ######################################################
    # Manhattan plot 
  ######################################################
  xpehh_plot <- ggplot(don, aes(x = BPcum, y = log10_FCS)) +
    geom_point(aes(color = as.factor(chr)), 
               size = 0.2, 
               alpha = 0.5,
               show.legend = FALSE) +
    scale_color_manual(values = rep(c("#7D7D7D", "#C5C5C5"), 22)) +
  
    # Override the colors 
    new_scale_color() +  # Create a new color scale
    geom_point(data = don %>% filter(arch %in% c("Denisovan", "Ambiguous", "Neandertal")),
               aes(x = BPcum, y = log10_FCS, color = arch), size = 0.4, alpha = 0.4, 
               inherit.aes = FALSE) +
    scale_color_manual(values = c(NA, "#7D7D7D", "#C5C5C5", 
                                  Denisovan = "#41BBEC", 
                                  Ambiguous = "#B588AF", 
                                  Neandertal = "#E84B42",
                                  guide = "none")) +
  
    # Custom X axis:
    scale_x_continuous(label = axisdf$chr, 
                       breaks = axisdf$center, 
                       expand=c(0,0)) +
  
    # Custom Y axis:
    ylim(0, 15) +
    ylab(expression(paste("-log"[10], " (", italic("p")["FCS"], ")")))+
    xlab("Chromosomes") +
    ggtitle(paste0(name_plot)) +
    theme_classic()+
    theme(
      text = element_text(size = 10)) +
  
    # threshold
    geom_hline(yintercept = 2, linetype = "dashed", color = "red", linewidth = 0.15) +
    geom_hline(yintercept = 3, linetype = "dashed", color = "red", linewidth = 0.15) +
    geom_hline(yintercept = 4, linetype = "dashed", color = "red", linewidth = 0.15) +
  
    # Adding labels 
    geom_label_repel(
      data = top_max_log_FCS_filtered,
     # aes(label = gene, fontface = 'bold', fill = arch),
      aes(label = gene, fontface = ifelse(bold, 'bold', 'plain'), fill = arch, size = label_size),
      alpha = 0.8, 
      box.padding = 0.2,
      hjust = 0.6, 
      size = 2.5,
      segment.size = 0.35,
      segment.color = '#272026',
      segment.linetype = 1, 
      segment.alpha = 0.6,    
      max.overlaps = Inf,
      max.iter = 50000,
      min.segment.length = 0.1, 
      nudge_y=4
    )+
  
    # Set the fill colors
    scale_fill_manual(values = c(Denisovan = "#41BBEC", 
                                 Ambiguous = "#B588AF", 
                                 Neandertal = "#E84B42",
                                 guide = "none"))+
    guides(alpha = "none")+
    theme(legend.position = "none")  # Remove the legend
  

  ggsave(paste0(name,".png"), xpehh_plot, width = 10, height = 3, dpi = 300, limitsize = FALSE)

  gc()

}


dev.off()

