##################################################
#    Heatmap - PIBv1_manuscript
##################################################

library(tidyverse)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)


# Select files that contain the introgressed core haplotypes per individual population 
temp <- list.files(pattern = "*uniquePOP.tsv", full.names = TRUE)

data_list <- list()

# Choose the Archaic: Ambiguous, Neandertal or Denisovan 
arch <- "Denisovan"

# Subset the necessary columns
for (i in 1:length(temp)) {
  name <- gsub(pattern = ".tsv", "", basename(temp[i]))
  df <- read_tsv(temp[[i]])
  df <- df %>%
    select(chr, start, end, population, FCS, log10_FCS, gene, TractID)
  df <- na.omit(df)
  df <- df %>%
    filter(str_detect(TractID, arch))
  data_list[[i]] <- df
}

# Merge windows within 10kb 
merge_windows_within_10kb <- function(df) {
  df %>%
    arrange(chr, start) %>%
    group_by(chr) %>%
    mutate(
      overlapping = start < lag(end, default = first(end)),
      distance = start - lag(end, default = first(start)),
      group_id = cumsum(start > lag(end, default = first(start)) + 10000)
    ) %>%
    group_by(chr, group_id) %>%
    summarize(
      new_start = min(start),
      new_end = max(end),
      gene = paste(unique(unlist(strsplit(gene, ",", fixed =TRUE))),
                   sep=",", collapse = ","),
      population = first(population),
      FCS = first(FCS),
      log10_FCS = max(log10_FCS)
    ) %>%
    ungroup() %>%
    select(-group_id)
}

overlapping_10kb <- lapply(data_list, merge_windows_within_10kb)

# 10kb large windows per population
merged_overlapping_10kb <- bind_rows(overlapping_10kb)

##################################################
  # Top windows per population 
##################################################

#  Subset the top windows with the highest log10_FCS for each population
top_windows_10kb_list <- lapply(overlapping_10kb, function(df) {
  df %>%
    arrange(desc(log10_FCS)) %>%
    slice_head(n = 10) %>%
    arrange(chr, new_start, new_end) %>%
    ungroup()
})

# combine all populations top windows per population
merged_top_windows_10kb <- bind_rows(top_windows_10kb_list)


##################################################
# Subset key genes
# In this case from the Interferon Gamma pathway 
##################################################
keygene_list <- readLines("interferon_gamma_genes.txt", warn = FALSE)

match_genes <- function(gene_line) {
  genes <- strsplit(gene_line, ",")[[1]]
  any(genes %in% keygene_list)
}

# Subset keygene_list from overlapping_10kb
filtered_key_genes_10kb_list <- lapply(overlapping_10kb, function(df) {
  df %>%
    filter(sapply(gene, match_genes))
})

filtered_key_genes_10kb_list


##################################################
 	 # IMPORTANT 
# Get the respective merged-windows for the other 
# populations directly from the FCS output
##################################################

temp2 <- list.files(pattern = "*FCS.tsv", full.names = TRUE)

process_window <- function(full_df_path, interval_df) {
  
  full_df <- read_tsv(full_df_path, col_select = c(chr, start, end, population, FCS, log10_FCS, gene))
  
  # Create a new df
  merged_df <- bind_rows(interval_df)
  
  # Set p0.001 to TRUE in the merged_df data frame
  merged_df$p0.01 <- TRUE
  
  # Create an empty data frame to store the merged rows and intervals
  merged_rows_df <- data.frame(chr = integer(),
                               new_start = integer(),
                               new_end = integer(),
                               population = character(),
                               gene = character(),
                               log10_FCS = numeric(),
                               p0.01 = logical())
  
  # Filter and merge the rows from full_df that fall within the intervals specified in interval_df
  for (i in seq_len(nrow(interval_df))) {
    interval <- interval_df[i, ]
    merged_rows <- full_df %>%
      filter(chr == interval$chr, start >= interval$new_start, end <= interval$new_end) %>%
      mutate(gene = paste(unique(unlist(strsplit(gene, ",", fixed =TRUE))),
                          sep=",", collapse = ",")) %>%
      select(-c(start, end)) 
    
    # Add the interval values to merged rows
    merged_rows$new_start <- interval$new_start
    merged_rows$new_end <- interval$new_end
    
    # Append the current interval's merged rows to the result data frame
    merged_rows_df <- rbind(merged_rows_df, merged_rows)
  }
  
  # Combine the rows from interval_df and the merged rows with intervals
  merged_df <- bind_rows(merged_df, merged_rows_df)
  
  # Group and summarize
  merged_df <- merged_df %>%
    group_by(chr, new_start, new_end, population, gene) %>%
    summarize(
      log10_FCS = max(log10_FCS, na.rm = TRUE),
      p0.01 = all(p0.01),
    ) %>%
    ungroup()
  
  merged_df$p0.01[is.na(merged_df$p0.01)] <- FALSE  
  
  rm(full_df);
  gc();
  
  return(merged_df)
}

  ##################################################
 # Apply function TO: 

  #overlapping_10kb
  #top_windows_10kb_list
  #filtered_key_genes_10kb_list
##################################################

windows_list_10kb <- lapply(temp2, process_window, 
                            interval_df= bind_rows(top_windows_10kb_list)) 

merged_df <- bind_rows(windows_list_10kb)

##################################################
# OPTION 1: For ALL selective sweeps 
##################################################

#merged_df <- bind_rows(windows_list_10kb)

#merged_df <- merged_df[merged_df$p0.01 == FALSE, ]

#merged_df

##################################################
# OPTION 2: Arch selective sweeps 
##################################################

merged_df <- bind_rows(windows_list_10kb)

merged_df <- merged_df %>%
  group_by(chr, new_start, new_end, population, gene, log10_FCS) %>%
  arrange(desc(p0.01)) %>%
  slice(1) %>%
  ungroup()

##################################################
# Apply the second windowed to remove redundancy 
# merge overlapping or adjacent genomic windows into larger, non-redundant regions
# Within each merged region (group of overlapping windows), take the maximum FCS value
##################################################

process_window2 <- function(df) {
  df %>%
    arrange(chr, new_start) %>%
    group_by(chr) %>%
    mutate(
      overlapping = new_start <= lag(new_end, default = first(new_start)),
      distance = new_start - lag(new_end, default = first(new_start)),
      group_id = cumsum(new_start >= lag(new_end, default = first(new_start)) & !overlapping)
    ) %>%
    group_by(chr,group_id, population) %>%
    summarize(
      new_start2 = min(new_start),
      new_end2 = max(new_end[length(new_end)]),
      log10_FCS = max(log10_FCS),
      gene = paste(unique(unlist(strsplit(gene, ",", fixed =TRUE))),
                   sep=",", collapse = ","),
      p0.01_true = any(p0.01),          
    ) %>%
    ungroup() %>%
    mutate(
      new_coordinates2 = paste(chr, new_start2, sep = ":"),
      new_coordinates2 = paste(new_coordinates2, new_end2, sep = "-")
    ) %>%
    select(-group_id)
}

merged_df2 <- process_window2(merged_df)


##################################################
  	# 10KB windows Heatmap "merged-windows_10kb" 
##################################################



# Color categories based on log10_FCS values
merged_df2$p_val <- cut(merged_df2$log10_FCS,
                        breaks = c(-Inf, 2, 3, 4, Inf),
                        labels = c("<2",">2",">3",">4"),
                        right = FALSE)


merged_df2$population <- factor(merged_df2$population, levels = c("Sepik-Goroka","Kove","Nakanai-Mangseng","Mamusi","Ata","Melamela","Baining-Kagat","Baining-Mali","Lavongai-Mussau","Nailik-Notsi-Tigak","Saposa","Nasioi","Vella-Lavella","Malaita","Bellona-Rennell","Santa-Cruz","Tikopia"))

#Separate by \n each 8 genes (just if needed)
#merged_df2$gene <- gsub("((?:\\w+,\\s*){3}\\w+),", "\\1,\n", merged_df2$gene)

# Sort 
merged_df2 <- merged_df2 %>% 
  arrange(chr, new_start2, new_end2)


# Create a separate dataframe for organizing gene labels on the right
gene_labels_df <- merged_df2 %>%
  filter(population == unique(merged_df2$population)[1]) %>%
  arrange(new_start2)


# Heatmap plot 

heatmap_plot <- ggplot(merged_df2, aes(
  x = population, 
  y = factor(new_coordinates2, levels = rev(unique(new_coordinates2))), 
  fill = ifelse(!!rlang::sym("p0.01_true"), as.character(p_val), "non_arch")
)) +
  geom_tile(color = "black", size = 0.3) +
  scale_fill_manual(
    values = c(
      "non_arch" = "#BEBEBE",
      ">2" = "#ffd3aa",
      ">3" = "#ffa755",
      ">4" = "#ff7b00"
    ),
    labels = c(
      #"<2" = ">0.01",
      ">2" = "<0.01",
      ">3" = "<0.001",
      ">4" = "<0.0001"
    )
  ) +
  labs(fill = expression(italic("p")[FCS])) +
  labs(x = "", y = "", title = "Ambiguous AI regions - IFNy") +
  theme_void() +
  theme(
    plot.title = element_text(size = 6),
    axis.text.y = element_text(hjust = 0, size = 4, vjust = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size =5),
    axis.text.x.bottom = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.text.x.top = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y.left = element_text(angle = 0, vjust = 0.5, hjust = 1),  # Rotate gene labels vertically
    plot.margin = margin(0.5, 0.5, 0.5, 1, "cm")
  ) 

# Organize the gene labels on the right side
heatmap_plot <- heatmap_plot +
  geom_tile(color = "black", size = 0.1) +  
  geom_text(
    data = gene_labels_df,
    aes(x = population, y = factor(new_coordinates2, levels = unique(new_coordinates2)), label = gene),
    nudge_x = 28,
    show.legend = FALSE,
    hjust = 1,
    vjust = 0.5, 
    size = 1,
  )

heatmap_plot

ggsave("merged-windows2_heatmap_Denisovan_AllPops_arch.pdf", 
       heatmap_plot, 
       width = 5.5, 
       height = 2, 
       dpi = 300, 
       limitsize = FALSE)
