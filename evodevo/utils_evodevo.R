library(tidyverse)


# function to load normalized gene expression data and time information
load_evodevo_data <- function(datadir, metadir, countdir, save = FALSE){
  if(!file.exists(file.path(datadir, "expression_all_species.RData")) |
     !file.exists(file.path(datadir, "time_points.RData")) ) {
    
    species <- c("Human", "Mouse", "Opossum", "Rabbit", "Rat")
    organs <- c("Brain", "Cerebellum", "Heart", "Kidney", "Liver", "Ovary", "Testis")
    germ_layers <- data.frame(Tissue = organs, 
                              Germ_layer = c("Ectoderm", "Ectoderm",
                                             "Mesoderm", "Mesoderm", 
                                             "Endoderm", "Mesoderm", "Mesoderm"))
    
    # normalized gene expression data
    expression_all_species <- read.delim(paste0(datadir,
                                                "normalized_orthologues_all_species.txt"), 
                                         sep = "\t", header = TRUE)
    
    # gene correspondence according to Nature paper
    orthologues <- read.table(paste0(metadir, "orthologues_species.txt"),
                              header = TRUE)
    
    # time correspondence according to Nature paper (DTW)
    time_points <- read.delim(paste0(metadir, "corresponding_time.txt"), sep = "\t", header = TRUE)
    colnames(time_points)[1] <- "sample"
    
    # add numeric time ranging from 1 to N_g per group (unaligned times)
    time_points$sample <- factor(time_points$sample, levels = unique(time_points$sample))
    time_points <- time_points %>% group_by(specie) %>%
      mutate(time_numeric = as.numeric(factor(sample))) %>% ungroup()
    
    # match gene ids to corresponding human gene
    get_matches <- function(specie) {
      match(filter(expression_all_species, species == specie)$gene, 
            orthologues[, paste0(specie, "_ID")])
    }
    corresponding_genes <- lapply(species, function(specie) {
      orthologues[get_matches(specie), "Human_ID"]
    })
    expression_all_species %<>% mutate(cor_h_gene = unlist(corresponding_genes))

    # add column with sample name (time_specie)
    expression_all_species %<>% mutate(sample = paste(time, species, sep = "_"))
    
    # expand to include all tissues for a given sample (to ensure same sample name and time for missing tissues when expanding times below)
    expression_all_species <- complete(expression_all_species,
                                       expand(expression_all_species, sample,
                                              tissue, cor_h_gene)) # complete samples across views and features
    
    # add matched times
    expression_all_species %<>% left_join(time_points, by= "sample", suffix = c("_id", "_matched"))
    expression_all_species %<>% select(-species) # remove incomplete specie column (species column from time point is complete and kept)
    
    # make feature names uniuq
    expression_all_species %<>% mutate(cor_h_gene = paste(cor_h_gene, tissue, sep = "_"))
    
    if(save) {
      save(expression_all_species, file = file.path(datadir, "expression_all_species.RData"))
      save(time_points, file = file.path(datadir, "time_points.RData"))
    }
  
  } else {
      load(file.path(datadir, "expression_all_species.RData"))
      load(file.path(datadir, "time_points.RData"))
    
  }
  
  return(list(expression_all_species = expression_all_species, time_points = time_points))
}

# function to obtain dataframe for input to MOFA
get_data4MOFA <- function(matched_times, species2include, views2include){
  
  if(!matched_times) {
    df_mofa <- subset(expression_all_species,
                      select = c("specie", "tissue", "sample",
                                 "cor_h_gene", "replicate_mean", "time_numeric"))
  } else {
    df_mofa <- subset(expression_all_species,
                      select = c("specie", "tissue", "sample",
                                 "cor_h_gene", "replicate_mean", "time_matched"))
  }
  colnames(df_mofa) <- c("group", "view", "sample", "feature", "value", "time")
  
  if(species2include != "all"){
    df_mofa %<>% filter(group %in% species2include)
  }
  
  if(views2include != "all"){
    df_mofa %<>% filter(view %in% views2include)
  }
  

  return(df_mofa)
}

# function to write data to csv for training in python
export_data <- function(model, dir){
  if(!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  dat <- get_data(model)
  for(g in groups_names(model)) {
    for (m in views_names(model)) {
      write.csv(dat[[m]][[g]],
                file = paste0(dir, "view_",m,"_group_",g, ".csv"))
      }
      write.csv(get_covariates(model)[[g]],
                file = paste0(dir, "times_group_",g,".csv"))
  }
}

# function to plot top weights per factor
plot_top_features <- function(object, factor, views = "all", n_views =5, n_genes = 3, gene_description){
  
  # get weights and extract ensemble id from feature name (ensid_view)
  df_weights <- get_weights(object, factors = factor, as.data.frame = T) %>%
    mutate(ens_id = gsub("_.*","",feature)) %>% group_by(ens_id, view, feature)
  
  if(views == "all"){
    views <- views_names(object)
  }
  
  factornm <- MOFA2:::.check_and_get_factors(object, factor)
  # get views with highest R2 by factor
  topviews <- get_variance_explained(object, as.data.frame = T)$r2_per_factor %>%
    filter(view %in% views) %>%
    group_by(factor, view) %>%
    summarise(mean = mean(abs(value))) %>%
    ungroup()%>%
    filter(factor == factornm) %>% arrange(-mean) %>%
    head(n_views) %>% select(view) %>% unlist()
  
  # subset to top n features per view (by absolute value)
  dftopn <- df_weights %>%
    filter(view %in% topviews) %>%
    group_by(view) %>%
    slice_max(abs(value), n = n_genes) %>%
    ungroup()
  
  dftopn <- left_join(dftopn, gene_description, by = "ens_id")
  
  data <- get_data(object, as.data.frame = TRUE) %>%
    rename(expression_value = value) %>%
    filter(feature %in% dftopn$feature)
  covari <- get_covariates(object, as.data.frame = TRUE, warped = TRUE) %>%
    rename(covariate_value = value)
  
  df_data <- left_join(data, covari, by ="sample")
  df_data %<>% filter(feature %in% dftopn$feature)
  dftopn <- left_join(df_data, dftopn, by = c("feature", "view"))
  dftopn %<>% mutate(symbol = ifelse(is.na(symbol), ens_id, symbol))
  dftopn %<>% mutate(symbol_view = paste(view, symbol, sep="-"))
  dftopn %<>% rename(species = group)
  
  ggplot(dftopn, aes(x=covariate_value, y = expression_value, col = species)) +
    geom_point() + stat_summary(fun = "mean", geom = "line") +
    facet_wrap(~symbol_view, ncol = n_genes, scales = "free_y") + theme_bw() +
    xlab("Developmental time") + ylab("Normalized expression") +
    theme(legend.position = "top")
}

# nicer heatmap
my_heatmap <- function(data, divergent_colors = TRUE, midpoint =0,
                       ncols = 100, x_label="", y_label="", min_c = NULL, max_c = NULL,...){
  if(divergent_colors){
    if(is.null(min_c))
      min_dat <- min(min(data, na.rm = TRUE), midpoint)
    else
      min_dat <- min_c
    if(is.null(max_c))
      max_dat <- max(max(data, na.rm = TRUE), midpoint)
    else
      max_dat <- max_c
    seq_breaks <- c(seq(min_dat, midpoint, (midpoint-min_dat)/ncols * 2), seq(midpoint, max_dat, (max_dat-midpoint)/ncols * 2)[-1])
    pheatmap::pheatmap(data, color = rev(colorRampPalette((RColorBrewer::brewer.pal(n = 7, name ="RdBu")))(ncols)), breaks = seq_breaks, ...)
  } else {
    pheatmap::pheatmap(data, color = colorRampPalette((RColorBrewer::brewer.pal(n = 7, name ="Blues")))(ncols), ...)
  }
}



get_paper_times <- function(metadir){
  dev_corres <- read.xlsx(file.path(metadir, "DevCorresp.xlsx"))[1:5,]
  dev_corres <- mutate_all(dev_corres, .funs = function(x) gsub("\\)", "", gsub(".*\\(", "", x))) # take unique match
  dev_corres[3,] <- ifelse(dev_corres[3,] == "tod", "toddler",
                           ifelse(dev_corres[3,] == "ya", "youngAdult",
                                  ifelse(dev_corres[3,] == "yma", "youngMidAge", dev_corres[3,] )))
  colnames(dev_corres) <- c("group", dev_corres[1,-1])
  df_paper_times <- gather(dev_corres[-1,], key = "id_mouse", value = "id_paper", -group)
  df_paper_times <- gather(dev_corres[-1,], key = "id_mouse", value = "id_paper", -group)
  # add intermediate times
  df_paper_times <- rbind(df_paper_times,
                          c(group = "rat", id_mouse = "e11.5", id_paper = "e13"),
                          c(group = "rat", id_mouse = "e14.5", id_paper = "e17"),
                          c(group = "rat", id_mouse = "P3", id_paper = "P7"),
                          c(group = "human", id_mouse = "e11.5", id_paper = "6w"),
                          c(group = "human", id_mouse = "e14.5", id_paper = "9w"),
                          c(group = "human", id_mouse = "e14.5", id_paper = "10w"),
                          c(group = "human", id_mouse = "e17.5", id_paper = "16w"),
                          c(group = "human", id_mouse = "e17.5", id_paper = "18w"),
                          c(group = "human", id_mouse = "P3", id_paper = "newborn"),
                          c(group = "human", id_mouse = "P14", id_paper = "infant"),
                          c(group = "human", id_mouse = "P14", id_paper = "school"),
                          c(group = "human", id_mouse = "P28", id_paper = "teenager"),
                          c(group = "human", id_mouse = "P63", id_paper = "olderMidAge"),
                          c(group = "human", id_mouse = "P63", id_paper = "senior"),
                          c(group = "opossum", id_mouse = "e16.5", id_paper = "P10"),
                          c(group = "opossum", id_mouse = "e16.5", id_paper = "P10"),
                          c(group = "opossum", id_mouse = "P3", id_paper = "P28"),
                          c(group = "opossum", id_mouse = "P14", id_paper = "P42"),
                          c(group = "opossum", id_mouse = "P63", id_paper = "P150"),
                          c(group = "opossum", id_mouse = "P63", id_paper = "P180"),
                          c(group = "rabbit", id_mouse = "e11.5", id_paper = "e14"),
                          c(group = "rabbit", id_mouse = "e12.5", id_paper = "e16.5"),
                          c(group = "rabbit", id_mouse = "e17.5", id_paper = "e24"))
  
  
  
  # capitalize species and rename
  df_paper_times$group <- paste(toupper(substr(df_paper_times$group, 1, 1)), substr(df_paper_times$group, 2, nchar(df_paper_times$group)), sep="")
  df_paper_times %<>% mutate(id_paper = gsub("w$", "wpc", id_paper))
  df_paper_times %<>% mutate(id_paper = gsub("e19.5 ", "e19.5", id_paper))
  
  return(df_paper_times)
}
    

plot_gsea_across_organs <- function(object, df, fac,
                                    min_r2 = 5, alpha = 0.05, Nmin = 4, show_bar = TRUE){
  
  combis <- object@cache$variance_explained$r2_per_factor %>%
    Reduce(pmax,.) %>% melt(.) %>% filter(value >= min_r2) %>%
    mutate(view_factor = paste(Var2, Var1, sep="-"))
  
  df %<>% mutate(view_factor = paste(view, factor, sep="-"))
  df %<>% filter(view_factor %in% combis$view_factor)
  
  
  # Show all organs with p-values below 5% FDR on these pathways
  df4plot <- df %>%
    filter(factor == fac) %>%
    filter(padj < alpha)
  
  # Show # organs enriched
  dfn <- df4plot %>% group_by(pathway) %>%
    summarise(Nview = n()) %>% ungroup()
  
  df4plot %<>% rename(organ = view)
  
  # arrange by number of organ
  dfn %<>% arrange(Nview)
  dfn %<>% filter(Nview >= Nmin)
  
  df4plot %<>% left_join(dfn, by = c("pathway"))
  df4plot %<>% mutate(pathway = factor(pathway, levels = dfn$pathway))
  df4plot %<>% filter(Nview >= Nmin)
  
  gg <- ggplot(df4plot, aes(x=pathway, y = -log(padj), col = organ)) +
    geom_point() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1)) +
    coord_flip()  + theme_bw()
    
    if(show_bar){
      gg <- gg + geom_bar(data = dfn, aes(x=pathway, Nview*10),
                       col = "gray", alpha = 0.3, stat = "identity") + 
        geom_hline(yintercept = 10 * object@dimensions$M,
                   lty = "dashed", col = "gray") +
        scale_y_continuous(sec.axis = sec_axis(~./10,
                                               name = "Number of organs with enrichment",
                                               breaks = c(1, 3, 5)))
    }
  gg
}


make_fig_weight_factor <- function(object, views, fac){
  top_list <- lapply(views, function(m){
    top2_brain <- plot_top_weights(object,
                                   view = m,
                                   factors = fac) +
      theme(strip.text = element_blank()) + ggtitle(m) +
      geom_point(col = col4organs[m])+
      theme(plot.title = element_text(colour = col4organs[m]))
  })
  
  top <- cowplot::plot_grid(plotlist = top_list, nrow =1)
  
  temp_list <- lapply(views, function(m){
    temp <- plot_data_vs_cov(MOFAobject_symbol,
                             factor = fac, features = 3, 
                             view = m, dot_size = 1, color_by = "species") +
      stat_summary(geom = "line", aes(col = color_by), fun = "mean") + 
      guides(col = FALSE) +
      theme(legend.position = "top",
            strip.text  = element_text(colour = col4organs[m]),
            axis.ticks.x = element_blank(),
            text = element_text(size = 13)) +
      scale_x_continuous(breaks = c(1.5,13.5), labels = c("early", "late")) +
      xlab("Developmental time")
    
    if(length(views) > 1){
      temp <- temp + facet_wrap(~feature, ncol = 1) 
    if(m != views[2]){
      temp <- temp + xlab("")
    }
    } else{
      temp <- temp + facet_wrap(~feature, nrow = 1) 
    }
    temp
  })
 
  if(length(views) == 1){
    cowplot::plot_grid(top,
                       temp_list[[1]] ,
                       ncol =1, rel_heights = c(1,1.2),
                       labels = c("A","B",""))
  } else {
    plot_grid(top,
              plot_grid(plotlist = lapply(temp_list, function(gg) gg + guides(fill = FALSE)),
                                 nrow =1), #axis = "l", align = "v",
                       get_legend(temp_list[[1]]),
                       ncol =1, rel_heights = c(1,2,0.1),
                       labels = c("A","B",""))
  }
  
}

plot_data_overview_custom <- function (object, covariate = 1, colors = NULL, show_covariate = FALSE, 
          show_dimensions = TRUE) 
{
  if (!is(object, "MOFA")) 
    stop("'object' has to be an instance of MOFA")
  data <- object@data
  if (any(object@dimensions[["C"]] < 1, is.null(object@covariates))) 
    covariate <- NULL
  if (!is.null(covariate)) {
    if (is.numeric(covariate)) {
      if (covariate > object@dimensions[["C"]]) 
        stop("Covariate index out of range")
      covariate <- covariates_names(object)[covariate]
    }
    if (!is.character(covariate) | !covariate %in% covariates_names(object)) 
      stop("Covariate mispecified. Please read the documentation")
    covari <- MOFA2:::.set_xax(object, covariate)
  }
  M <- get_dimensions(object)[["M"]]
  G <- get_dimensions(object)[["G"]]
  if (M == 1 & G == 1) 
    warning("This function is not useful when there is just one view and one group")
  if (is.null(dim(data[[1]][[1]]))) 
    stop("Data not found")
  if (is.null(colors)) {
    palette <- c("#FF7F50", "#D95F02", "#377EB8", "#E6AB02", 
                 "#31A354", "#7570B3", "#E7298A", "#66A61E", "#A6761D", 
                 "#666666", "#E41A1C", "#4DAF4A", "#984EA3", "#FF7F00", 
                 "#FFFF33", "#A65628", "#F781BF", "#1B9E77")
    if (M < 18) 
      colors <- palette[seq_len(M)]
    else colors <- rainbow(M)
    names(colors) <- views_names(object)
  }
  else {
    if (length(colors) != M) 
      stop("Length of 'colors' does not match the number of views")
    if (is.null(names(colors))) {
      names(colors) <- views_names(object)
    }
    else {
      stopifnot(sort(names(colors)) == sort(views_names(object)))
    }
  }
  tmp <- lapply(data, function(m) sapply(m, function(g) apply(g, 
                                                              2, function(x) !all(is.na(x)))))
  ovw <- do.call(cbind, lapply(seq_len(M), function(m) {
    do.call(rbind, lapply(tmp[[m]], as.data.frame))
  }))
  rownames(ovw) <- object@samples_metadata$sample
  colnames(ovw) <- views_names(object)
  ovw$sample <- object@samples_metadata$sample
  ovw$group <- object@samples_metadata$group
  to.plot <- melt(ovw, id.vars = c("sample", "group"), var = c("view"))
  if (!is.null(covariate)) {
    to.plot <- left_join(to.plot, covari, by = "sample")
    to.plot$sample <- factor(to.plot$sample, levels = unique(to.plot$sample[order(to.plot$covariate_value)]))
  }
  else {
    to.plot$sample <- factor(to.plot$sample, levels = rownames(ovw))
  }
  n <- length(unique(to.plot$sample))
  to.plot$combi <- ifelse(to.plot$value, as.character(to.plot$view), 
                          "missing")
  if (show_dimensions) {
    to.plot$ntotal <- paste("N=", sapply(data[[1]], function(e) ncol(e))[as.character(to.plot$group)], 
                            sep = "")
    to.plot$ptotal <- paste("D=", sapply(data, function(e) nrow(e[[1]]))[as.character(to.plot$view)], 
                            sep = "")
    if (length(unique(to.plot$group)) == 1) {
      to.plot <- mutate(to.plot, view_label = paste(view, 
                                                    ptotal, sep = "\n"), group_label = ntotal)
    }
    else {
      to.plot <- mutate(to.plot, view_label = paste(view, 
                                                    ptotal, sep = "\n"), group_label = paste(group, 
                                                                                             ntotal, sep = "\n"))
    }
  }
  else {
    to.plot <- mutate(to.plot, view_label = view, group_label = group)
  }
  p <- ggplot(to.plot, aes(x = as.factor(covariate_value), y = view_label, 
                                  fill = combi)) + geom_tile() +
    scale_fill_manual(values = c(missing = "grey", colors)) +
    guides(fill = FALSE) + facet_wrap(~group_label, nrow = length(unique(to.plot$view_label))) + 
    theme(panel.background = element_rect(fill = "grey"), 
          text = element_text(size = 14), axis.line = element_blank(), 
          axis.title.y = element_blank(), 
          axis.text.y = element_text(color = "black"), 
          axis.ticks.y = element_blank(),
          strip.background = element_blank(), panel.grid = element_blank()) +
    xlab("Unmatched time points") 
  if (show_covariate) {
    p2 <- ggplot(to.plot, aes_string(x = "sample", y = "covariate_value")) + 
      geom_point(size = 0.5) + theme_bw() + theme(text = element_text(size = 10), 
                                                  axis.ticks.x = element_blank(), axis.title.x = element_blank(), 
                                                  axis.text.x = element_blank(), strip.background = element_blank(), 
                                                  strip.text = element_blank()) + ylab(covariate) + 
      facet_wrap(~group_label, ncol = 1)
    if (object@dimensions["G"] == 1) {
      p <- cowplot::plot_grid(p, p2, align = "v", ncol = 1, 
                              rel_heights = c(1, 0.2))
    }
    else {
      p <- cowplot::plot_grid(p, p2, align = "h", nrow = 1, 
                              rel_widths = c(1, 1))
    }
  }
  return(p)
}
