plot_sharedness <- function(df, methods ="MEFISTO") {
  
  df <- df %>% filter(method %in% methods,
                      !(param == "Number of groups" & value ==1),
                      ! is.na(true_sharedness_0)) 
  # inferred sharedness
  df1 <-  df %>%
    select(param, value, starts_with("inferred_sharedness"), method, seed) %>% 
    gather(starts_with("inferred_sharedness"), key = "factor",
           value = "inferred_sharedness") %>%
    mutate(factor = gsub("inferred_sharedness", "Factor" , factor))
  
  # simulated sharedness
  df2 <- df %>% 
    select(param, value, starts_with("true_sharedness"), method, seed) %>%
    gather(starts_with("true_sharedness"), key = "factor", value = "true_sharedness") %>%
    mutate(factor = gsub("true_sharedness", "Factor",factor))
  
  # join and rename 
  dfjoint <- inner_join(df1, df2, by =c("factor", "param",
                                        "value", "method", "seed")) 
  dfjoint %<>% mutate(factor = ifelse(true_sharedness == 0, "non-shared", 
                                      ifelse(true_sharedness == 1, "shared", "partially shared")))
  
  # set colors
  cols <- RColorBrewer::brewer.pal(7, "Dark2")[c(3,5,7)]
  names(cols) <- c("non-shared", "shared", "partially shared")
  
  # plot
  gg <- dfjoint %>%
    ggplot(aes(x= value, y= inferred_sharedness,
               col = factor, lty = method))
  
  if(length(methods) > 1){
    gg <- gg +  stat_summary(fun.data = "mean_se", geom = "line")
  } else {
    gg <- gg +  stat_summary(fun.data = "mean_se", geom = "line", lty = 1) + guides(lty = FALSE)
  }
  if(length(unique(dfjoint$param)) > 1){
    gg <- gg +
      facet_wrap(~ param, scales = "free_x") + xlab("")
  } else {
    stopifnot(unique(dfjoint$param) == "N")
    gg <- gg +
      xlab("Number of timepoints")
  }
  
  gg <- gg +
    stat_summary(fun.data = "mean_se", size = 0.1) +
    stat_summary(aes( y = true_sharedness),
                 geom = "line",lty = "dashed")  + theme_bw() +
    ylab("Sharedness") +
    scale_color_manual(values = cols) +
    scale_x_continuous(breaks = scales::pretty_breaks(4), limits = c(NA, NA))
}


plot_smoothness <- function(df, methods ="MEFISTO") {
  
  # get inferred scales
  df_learnt <- df %>% filter(method %in% methods) %>%
    select(param, value, starts_with("scales_learnt"), method, seed) %>%
    gather(starts_with("scales_learnt"), key = "factor", value = "scale_learnt") %>%
    mutate(factor = gsub("scales_learnt", "Factor" , factor))
  
  # get simulated scales
  df_sim <- df %>% filter(method %in% methods) %>%
    select(-starts_with("scales_learnt")) %>%
    select(param, value, starts_with("scales"), method, seed) %>%
    gather(starts_with("scales"), key = "factor", value = "scale") %>%
    mutate(factor = gsub("scales", "Factor",factor))
  
  # join and rename factors
  df_joint <- inner_join(df_learnt, df_sim,
                         by =c("factor", "param", "value", "method", "seed")) %>%
    mutate(factor = ifelse(factor == "Factor_0", "smooth", 
                           ifelse(factor == "Factor_1", "partially smooth", "non-smooth"))) 
  
  # plot
  gg <- ggplot(df_joint, aes(x= value, y= scale_learnt,
                             col = factor, lty = method)) 
  if(length(methods) > 1){
    gg <- gg +  stat_summary(fun.data = "mean_se", geom = "line")
  } else {
    gg <- gg +  stat_summary(fun.data = "mean_se", geom = "line", lty = 1) + guides(lty = FALSE)
  }
  if(length(unique(df_joint$param)) > 1){
    gg <- gg +
      facet_wrap(~ param, scales = "free_x") + xlab("")
  } else {
    stopifnot(unique(df_joint$param) == "N")
    gg <- gg +
      xlab("Number of timepoints")
  }
  
  # set colors
  cols <- RColorBrewer::brewer.pal(7, "Dark2")[1:3]
  names(cols) <- c("smooth", "partially smooth", "non-smooth")
  
  gg <- gg +
    stat_summary(fun.data = "mean_se", size = 0.1) +
    geom_hline(aes(yintercept = scale, col = factor), lty = "dashed") +
    theme_bw() + scale_color_manual(values = cols) +
    ylab("Smoothness") +
    scale_x_continuous(breaks = scales::pretty_breaks(4), limits = c(NA, NA))
}



# helper function to plot hyperparameters (learnt vs true)
plot_hyper <- function(df, par, methods ="MEFISTO", lty = TRUE) {
  par_learnt  <- paste(par, "learnt", sep = "_")
  
  df1 <- df %>% filter(method %in% methods) %>%
    select(param, value, starts_with(par_learnt), method) %>%
    gather(starts_with(par_learnt), key = "factor",
           value = !!par_learnt) %>%
    mutate(factor = gsub(par_learnt, "Factor" , factor))
  
  df2 <- df %>% filter(method %in% methods) %>%
    select(-starts_with(par_learnt)) %>%
    select(param, value, starts_with(par), method) %>%
    gather(starts_with(par), key = "factor", value = !!par) %>%
    mutate(factor = gsub(par, "Factor",factor))
  
  gg <- inner_join(df1, df2, by =c("factor", "param",
                                   "value", "method")) %>%
    mutate(factor = ifelse(factor == "Factor_0", "smooth", 
                           ifelse(factor == "Factor_1", "partially smooth", "non-smooth"))) %>%
    ggplot(aes_string(x= "value", y= par_learnt,
                      col = "factor", lty = "method")) 
  if(lty){
    gg <- gg +  stat_summary(fun.data = "mean_se", geom = "line")
  } else {
    gg <- gg +  stat_summary(fun.data = "mean_se", geom = "line", lty = 1) +
      guides(lty = FALSE)
  }
  gg <- gg +
    facet_wrap(~ param, scales = "free_x") +
    stat_summary(fun.data = "mean_se", size = 0.1) +
    geom_hline(aes_string(yintercept = par, col = "factor"),
               lty = "dashed") +
    theme_bw() + scale_color_brewer(palette="Dark2") +
    scale_x_continuous(breaks = scales::pretty_breaks(4), limits = c(NA, NA))
}


make_nice_names <- function(df){
  stopifnot("param" %in% colnames(df))
  df <- mutate(df, param = ifelse(param == "N", "Number of timepoints",
                        ifelse(param == "G", "Number of groups",
                               ifelse(param == "missing", "% missing samples", 
                                      ifelse(param == "noise", "Noise variance", param)))))
  df <- mutate(df, param = factor(param, levels = c("Number of timepoints",
                                                        "Number of groups",
                                                        "% missing samples",
                                                        "Noise variance")))
  return(df)
}