make_genera_plot <- function(df_weights, factor_name, type = "delivery", ngenera = 5){
  # filter out genera with missing genera information and shorten name
  df_weights %<>% filter(factor == factor_name)
  df_weights %<>% filter(!is.na(genus), genus != "g__")
  df_weights %<>% mutate(genus = gsub("g__","", genus))
  df_weights %<>% mutate(genus = ifelse(genus == " Pseudoramibacter_Eubacterium",
                                         "Pseudoramibacter_Eub.", genus))
  df_weights %<>% mutate(genus = gsub("\\[","", genus))
  df_weights %<>% mutate(genus = gsub("\\]","", genus))
  # summarize by mean weights across all species in the genus and 
  # filter to top 10 positive and negative ones
  df_top <- df_weights %>% group_by(genus) %>% 
    summarize(mean_weight = mean(value), n_spec = n()) %>%
    ungroup() %>% filter(n_spec > 2) 
    # cesarean is linked to higher latent values
    if(type == "delivery"){
      df_top <- mutate(df_top, type = ifelse(mean_weight > 0, "Cesarean", "Vaginal"))
      cols <- cols4delivery
    } else if (type == "diet"){
      df_top <- mutate(df_top, type = ifelse(mean_weight > 0, "bd", "fd")) 
      cols <- cols4diet
    }
  df_top <- df_top %>%
    group_by(type) %>%
    slice_max(abs(mean_weight), n= ngenera) %>% ungroup() %>% arrange(mean_weight) %>%
    mutate(genus = factor(genus, levels = .$genus))
  
  # plot
  gg_w1 <-  df_top %>% 
    ggplot(aes(x= genus, y = mean_weight, fill = type)) + geom_bar(stat="identity") +
    coord_flip() + theme_bw()  + scale_fill_manual(values = cols) +
    theme(legend.position = "top") + xlab("") + guides(fill = guide_legend(title="")) +
    geom_point(data = filter(df_weights, genus %in% df_top$genus),
               aes(x = genus, y = value), inherit.aes = FALSE, alpha = 0.3)  +
    ylab(paste("Weight", factor_name))
}