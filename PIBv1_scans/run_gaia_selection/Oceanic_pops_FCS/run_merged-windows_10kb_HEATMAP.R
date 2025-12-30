##########################################################
          # 10kb windows and heatmap 
##########################################################
library(tidyverse)
library(ggrepel)

run_heatmap_pipeline <- function(arch) {
  
  # Select files that contain the introgressed core haplotypes per individual population 
  temp <- list.files(pattern = "*uniquePOP.tsv", full.names = TRUE)
  temp2 <- list.files(pattern = "*FCS.tsv", full.names = TRUE)
  
  cat("Reading in adaptive introgression candidate windows of", arch, "origin\n")

  data_list <- list()
  
  # Subset the necessary columns
  for (i in 1:length(temp)) {
    name <- gsub(pattern = ".tsv", "", basename(temp[i]))
    df <- read_tsv(temp[[i]],
                   col_select=c(chr, start, end, population, FCS, log10_FCS, gene, TractID)) %>%
      na.omit() %>%
      filter(str_detect(TractID, arch))
    data_list[[i]] <- df
  }
  
  # Merge windows within 10kb 
  merge_windows_within_10kb <- function(df) {
    df %>%
      arrange(chr, start) %>%
      group_by(chr) %>%
      mutate(group_id = cumsum(start > lag(cummax(end), default = first(start)) + 10000)) %>%
      group_by(chr, group_id) %>%
      summarize(new_start=min(start),
                new_end=max(end),
                gene=paste(unique(unlist(strsplit(gene, ",", fixed =TRUE))),
                           sep=",",
                           collapse = ","),
                population=first(population),
                FCS=min(FCS),
                log10_FCS=max(log10_FCS)) %>%
      ungroup() %>%
      select(-group_id)
  }
  
  cat("Merging windows of", arch, "origin within 10 kbp of each other within population\n")
  merged_10kb_windows <- lapply(data_list, merge_windows_within_10kb)
  
  interval_df <- merged_10kb_windows %>%
    bind_rows()
  # all 10kb windowns
  outfn <- paste0("merged_10kb_windows_", arch, ".tsv")
  cat("Saving 10 kbp within-pop merged", arch, "regions to", outfn, "\n")
  interval_df %>%
    write.table(file=outfn,
                sep="\t",
                quote=FALSE,
                row.names=FALSE)
  # Top windows per population
  top_windows_10kb_df <- interval_df %>%
    arrange(desc(log10_FCS)) %>%
    group_by(population) %>%
    slice_head(n=10) %>%
    arrange(population, chr, new_start, new_end) %>%
    ungroup()
  outfn <- paste0("merged_10kb_windows_", arch, "_TOP10.tsv")
  cat("Saving top 10 per population 10 kbp within-pop merged", arch, "regions to", outfn, "\n")
  top_windows_10kb_df %>%
    write.table(file=outfn,
                sep="\t",
                quote=FALSE,
                row.names=FALSE)
  # Subset key genes - IFNy pathway
  keygene_list <- readLines("interferon_gamma_genes.txt", warn = FALSE)
  match_genes <- function(gene_line) {
    genes <- strsplit(gene_line, ",")[[1]]
    any(genes %in% keygene_list)
  }
  filtered_key_genes_10kb_list <- interval_df %>%
    filter(map_lgl(gene, match_genes))
  outfn <- paste0("merged_10kb_windows_", arch, "_IFNy.tsv")
  cat("Saving key genes 10 kbp within-pop merged", arch, "regions to", outfn, "\n")
  filtered_key_genes_10kb_list %>%
    write.table(file=outfn,
                sep="\t",
                quote=FALSE,
                row.names=FALSE)
  #Retrieve regions from a pop:
  retrieve_FCS_regions <- function(fcs_path, target_regions) {
    #Read in the TSV of selection scan windows for the current population:
    fcs_df <- read_tsv(fcs_path,
                       col_select=c(chr, start, end, population, FCS, log10_FCS, gene))
    #Function to extract windows within a given interval and merge together:
    extract_region <- function(chr_a, start_a, end_a, pop_a, window_df) {
      window_df %>%
        filter(chr == chr_a, start >= start_a, end <= end_a) %>%
        summarize(chr=first(chr),
                  new_start=min(start),
                  new_end=max(end),
                  gene=paste(unique(unlist(strsplit(gene, ",", fixed=TRUE))),
                             sep=",",
                             collapse=","),
                  population=first(population),
                  FCS=min(FCS, na.rm=TRUE),
                  log10_FCS=max(log10_FCS, na.rm=TRUE),
                  p0.01=case_when(population == pop_a ~ TRUE,
                                  TRUE ~ FALSE))
    }
    #Extract windows within each target region and merge windows so that
    # the merged region is the same as the target region:
    target_regions_annot <- target_regions %>%
      mutate(p0.01=TRUE)
    pop_merged_regions <- target_regions %>%
      dplyr::select(chr, new_start, new_end, population) %>%
      pmap(\(chr, new_start, new_end, population) extract_region(chr, new_start, new_end, population, window_df=fcs_df)) %>%
      bind_rows()
    #Now return:
    pop_merged_regions
  }
  
  # Apply function 
  cat("Extracting FCS windows for", arch, "adaptive introgression regions for all populations\n")
#  windows_list_10kb <- lapply(temp2, process_window, interval_df = interval_df)
  windows_list_10kb <- temp2 %>%
    map(\(fcs_path) retrieve_FCS_regions(fcs_path, target_regions=interval_df)) %>%
    bind_rows()

  # Arch selective sweeps 
  merged_df <- windows_list_10kb %>%
    group_by(chr, new_start, new_end, population, gene, log10_FCS) %>%
    arrange(desc(p0.01)) %>%
    slice_head(n=1) %>%
    ungroup()

  outfn <- paste0("merged_df_", arch, "_PFRtest.tsv")
  cat("Saving deduplicated FCS windows for", arch, "regions from all populations to", outfn, "\n")
  merged_df %>%
    write.table(file=outfn,
                sep="\t",
                quote=FALSE,
                row.names=FALSE)

  top10_merged_df <- merged_df %>%
    inner_join(top_windows_10kb_df %>%
                 dplyr::select(chr, new_start, new_end) %>%
                 distinct(),
               by=c("chr", "new_start", "new_end"))

  outfn <- paste0("merged_df_", arch, "_TOP10_PFRtest.tsv")
  cat("Saving deduplicated FCS windows for top 10", arch, "regions from all populations to", outfn, "\n")
  top10_merged_df %>%
    write.table(file=outfn,
                sep="\t",
                quote=FALSE,
                row.names=FALSE)

  key_genes_merged_df <- merged_df %>%
    inner_join(filtered_key_genes_10kb_list %>%
                 dplyr::select(chr, new_start, new_end) %>%
                 distinct(),
               by=c("chr", "new_start", "new_end"))

  outfn <- paste0("merged_df_", arch, "_IFNy_PFRtest.tsv")
  cat("Saving deduplicated FCS windows for key genes", arch, "regions from all populations to", outfn, "\n")
  key_genes_merged_df %>%
    write.table(file=outfn,
                sep="\t",
                quote=FALSE,
                row.names=FALSE)

  #Merge adaptive introgression regions across populations for heatmap:
  merge_windows_across_pops <- function(df) {
    df %>%
      arrange(chr, new_start) %>%
      group_by(chr) %>%
      mutate(group_id=cumsum(new_start >= lag(cummax(new_end), default=first(new_start)))) %>%
      group_by(chr, group_id, population) %>%
      summarize(new_start=min(new_start),
                new_end=max(new_end),
                log10_FCS=max(log10_FCS, na.rm=TRUE),
                gene=paste(unique(unlist(strsplit(gene, ",", fixed=TRUE))),
                           sep=",",
                           collapse=","),
                p0.01=any(p0.01)) %>%
      ungroup() %>%
      mutate(coords=str_c(chr, ":", new_start, "-", new_end)) %>%
      dplyr::select(-group_id)
  }

  pop_order <- c("Sepik-Goroka", "Kove", "Nakanai-Mangseng",
                 "Mamusi", "Ata", "Melamela",
                 "Baining-Kagat", "Baining-Mali", "Lavongai-Mussau",
                 "Nailik-Notsi-Tigak", "Saposa", "Nasioi",
                 "Vella-Lavella", "Malaita", "Bellona-Rennell",
                 "Santa-Cruz", "Tikopia")

  #Merge, filter, output tables, and plot heatmaps for the three plot sets:
  # (all regions, top 10 regions, and key genes subset)
  analysis_types <- c("all", "top10", "key_genes")
  title_suffixes <- c("all", "TOP10", "IFNy")
  heights <- c(30, 8, 2)
  merged_dfs <- list()
  merged_dfs[[analysis_types[1]]] <- merged_df
  merged_dfs[[analysis_types[2]]] <- top10_merged_df
  merged_dfs[[analysis_types[3]]] <- key_genes_merged_df

  for (i in seq_along(analysis_types)) {
    cat("Merging", arch, "re-extracted regions across populations\n")
    across_pop_merged_df <- merged_dfs[[i]] %>%
      merge_windows_across_pops() %>%
      mutate(population=factor(population,
                               levels=pop_order),
             AdIntHit=case_when(p0.01 & log10_FCS > 4 ~ ">4",
                                p0.01 & log10_FCS > 3 ~ ">3",
                                p0.01 & log10_FCS > 2 ~ ">2",
                                TRUE ~ "non_arch"),
             coords=factor(coords, levels=rev(unique(coords))))

    outfn <- paste0("merged_df2_", arch, "_", title_suffixes[i], "_PFRtest.tsv")
    cat("Saving across-population merged", arch, "regions to", outfn, "\n")
    across_pop_merged_df %>%
      write.table(file=outfn,
                  sep="\t",
                  quote=FALSE,
                  row.names=FALSE)

  # 10KB windows Heatmap "merged-windows_10kb"
    #Prep the labels for the y axes:
    cat("Establishing heatmap y-axis labels of both merged", title_suffixes[i], arch, "region coordinates and gene lists\n")
    coords_genes_df <- across_pop_merged_df %>%
      dplyr::select(coords, gene) %>%
      distinct()
    coords_genes_breaks <- 1:length(levels(coords_genes_df$coords))
    coords_labels <- levels(coords_genes_df$coords)
    gene_labels <- rev(coords_genes_df$gene)

    #Now generate the heatmap:
    cat("Generating heatmap for merged", title_suffixes[i], arch, "regions\n")
    heatmap_plot <- across_pop_merged_df %>%
      ggplot(aes(x=population,
                 y=as.numeric(coords),
                 fill=AdIntHit)) +
        geom_tile(color="black",
                  linewidth=0.3) +
        scale_y_continuous(expand=c(0, 0.6),
                           breaks=coords_genes_breaks,
                           labels=coords_labels,
                           sec.axis=sec_axis(~.,
                                             breaks=coords_genes_breaks,
                                             labels=gene_labels)) +
        scale_fill_manual(values=c("nonarch"="#BEBEBE",
                                   ">2"="#ffd3aa",
                                   ">3"="#ffa755",
                                   ">4"="#ff7b00"),
                          labels=c(">2"="<0.01",
                                   ">3"="<0.001",
                                   ">4"="<0.0001")) +
        labs(x="",
             y="",
             fill=expression(italic("p")[FCS]),
             title=paste(arch, "adaptive introgression regions -", title_suffixes[i])) +
        theme_void() +
        theme(plot.title=element_text(size=6),
              axis.text.x=element_text(angle=90, size=5, hjust=1, vjust=0.5),
              axis.text.x.bottom=element_text(angle=90, hjust=1, vjust=0.5),
              axis.text.x.top=element_text(angle=90, hjust=1, vjust=0.5),
              axis.text.y=element_text(size=4, hjust=0, vjust=0.5),
              axis.text.y.left=element_text(angle=0, hjust=1, vjust=0.5),
              plot.margin=margin(0.5, 0.5, 0.5, 1, "cm"),
              legend.position=c(2.0, 0.9))
  
    # Save plot
    outfn <- paste0("merged_10kb_windows_heatmap_", arch, "_", title_suffixes[i], ".pdf")
    cat("Saving heatmap to", outfn, "\n")
    ggsave(outfn,
           heatmap_plot, 
           width = 5.5, 
           height = heights[i],
           dpi = 300, 
           limitsize = FALSE)

    # Save RData
    outfn <- paste0("selection_scans_10kb_windows_", title_suffixes[i], "_", arch, ".RData")
    cat("Saving R data.frames and objects to", outfn, "\n")
    save(
      interval_df,
      top_windows_10kb_df,
      filtered_key_genes_10kb_list,
      windows_list_10kb,
      merged_df,
      across_pop_merged_df,
      coords_genes_df,
      heatmap_plot,
      file=outfn
    )
  }
} 

# Run for all combinations
archaics <- c("Ambiguous", "Neandertal", "Denisovan")

for (arch in archaics) {
  run_heatmap_pipeline(arch)
}
