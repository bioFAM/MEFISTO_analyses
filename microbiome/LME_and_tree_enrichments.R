# ##############################################################################
#
## Explore associations with metadata
#
# ##############################################################################

library('tidyverse')
library("lmerTest")
library("progress")
library("cowplot")
library("ggrepel")

datadir <- "data/processed_data/"

# sample data
meta <- read_csv(paste0(datadir, 'samples_metadata.csv'))

# factors
factors.ctf <- read_csv(paste0(datadir, 'CTF_factors.csv')) %>% 
  full_join(meta %>% transmute(group=child_id, month, abx_all_sources))
factors.mefisto <- read_csv(paste0(datadir, 'MEFISTO_factors.csv')) %>% 
  full_join(meta %>% transmute(group=child_id, month, abx_all_sources))

# weights
weights.mefisto <- read_csv(paste0(datadir, 'MEFISTO_weights.csv'))
weights.ctf <- read_csv(paste0(datadir, 'CTF_weights.csv'))


# sOTU table
feat.mat.rclr <- read_csv(paste0(datadir, 'rclr_mat_filtered.csv'))
feat.mat.rclr <- as.data.frame(feat.mat.rclr)
rownames(feat.mat.rclr) <- feat.mat.rclr$SampleID
feat.mat.rclr$SampleID <- NULL
feat.mat.rclr <- as.matrix(feat.mat.rclr)

# ##############################################################################
# compute associations on ASV level

# test ASVs for associations with metadata
df.res.assoc <- list()
meta.test <- meta %>% 
  select(SampleID, month, child_id, abx_all_sources, 
         delivery, diet, sex) %>% 
  as.data.frame()

rownames(meta.test) <- meta.test$SampleID
for (var in setdiff(colnames(meta.test), c('month', 'child_id', 'SampleID'))){
  message(var)
  
  pb <- progress::progress_bar$new(total=nrow(feat.mat.rclr))
  
  # rclr abundances
  results <- vapply(rownames(feat.mat.rclr), FUN = function(x){
    feat <- feat.mat.rclr[x,rownames(meta.test)]
    df.test <- cbind(meta.test, feat)
    if (var %in% c('early_late', 'year')){
      .f <- paste0(paste0('feat~', var, '+(1|child_id)'))
      df.test <- df.test[!is.na(df.test[[var]]),]
    } else {
      .f <- paste0(paste0('feat~', var, '+month+(1|child_id)'))
    }
    if (mean(df.test$feat!=0) < 0.05){
      pb$tick()
      return(c('estimate'=NA_real_, 'p.value'=NA_real_))
    }
    fit <- suppressMessages(suppressWarnings(
      lmer(.f, data=df.test)))
    res <- coefficients(summary(fit))
    pb$tick()
    return(c('estimate'=res[2,1], 'p.value'=res[2,5]))
  }, FUN.VALUE = double(2))
  results <- t(results) %>% 
    as_tibble(rownames = 'ASV') %>% 
    mutate(variable=var)
  df.res.assoc[[(length(df.res.assoc)+1)]] <- results
}

df.res.assoc <- bind_rows(df.res.assoc)
df.res.assoc %>% 
  ggplot(aes(x=estimate, y=-log10(p.value))) + 
    geom_point() + 
    facet_grid(~variable)

# ##############################################################################
# correlate associations with factor loading

# compare to factor weights
df.plot.mefisto <- df.res.assoc %>% 
  transmute(feature=ASV, estimate, p.value, variable) %>% 
  full_join(weights.mefisto, by='feature') %>% 
  left_join(enframe(rowMeans(feat.mat.rclr!=0), name='feature', 
                    value='prevalence'))
df.plot.ctf <- df.res.assoc %>% 
  transmute(feature=ASV, estimate, p.value, variable) %>% 
  full_join(weights.ctf, by='feature') %>% 
  left_join(enframe(rowMeans(feat.mat.rclr!=0), name='feature', 
                    value='prevalence'))

# correlation
cor.ctf <- df.plot.ctf %>% 
  filter(p.value < 0.1) %>%
  filter(prevalence>0.05) %>% 
  group_by(variable) %>% 
  mutate(estimate=-estimate) %>% 
  summarise(cor1=cor(estimate, Factor1, method='spearman', 
                     use='pairwise.complete.obs'),
            cor2=cor(estimate, Factor2, method='spearman', 
                     use='pairwise.complete.obs'))
cor.mefisto <- df.plot.mefisto %>% 
  filter(p.value < 0.1) %>%
  filter(prevalence>0.05) %>% 
  group_by(variable) %>% 
  mutate(estimate=-estimate) %>% 
  summarise(cor1=cor(estimate, Factor1, method='spearman', 
                     use='pairwise.complete.obs'),
            cor2=cor(estimate, Factor2, method='spearman', 
                     use='pairwise.complete.obs'))

g1 <- bind_rows(cor.ctf %>% mutate(data='ctf'), 
          cor.mefisto %>% mutate(data='mefisto')) %>% 
  pivot_longer(cols=c(cor1, cor2)) %>% 
  mutate(name=case_when(name=='cor1'~'Factor1', 
                        name=='cor2'~'Factor2')) %>% 
  mutate(variable=case_when(variable=='abx_all_sources'~'antibiotics',
                            TRUE~variable)) %>% 
  mutate(variable=str_to_title(variable)) %>% 
  mutate(data=toupper(data)) %>% 
  ggplot(aes(x=name, y=value, fill=variable)) + 
    geom_bar(stat='identity', position = position_dodge()) +
    facet_grid(data~.) + 
    theme_bw() + 
    xlab('') + ylab("Spearman correlation") + 
    ggembl::theme_publication(panel.grid = 'major_y') + 
    scale_fill_manual(values=c('#009F4D', '#307FE2', '#E40046', '#FFA300'), 
                      name='') +
    scale_y_continuous(limits=c(-0.8, 0.8))

# ##############################################################################
# example plots

# Factor 1
g2 <- df.plot.mefisto %>% 
  mutate(type=case_when(p.value < 0.1 & estimate > 0~'cesarean',
                        p.value < 0.1 & estimate < 0~'vaginal',
                        TRUE~'other')) %>% 
  mutate(type=factor(type, levels = c('cesarean', 'vaginal', 
                                      'other'))) %>% 
  filter(prevalence>0.05) %>% 
  filter(variable=='delivery') %>% 
  arrange(type) %>% 
  mutate(label=case_when(abs(estimate) > 0.48 ~ genus,
                         TRUE~'')) %>% 
  ggplot(aes(x=-estimate, y=Factor1, col=type)) + 
  geom_point() + 
  theme_bw() + 
  xlab('LME model coefficient') + 
  ylab('Factor 1 weight') + 
  theme(panel.grid.minor = element_blank()) + 
  scale_colour_manual(values=c(c("#d95f02", "#e6ab02"), 
                               alpha('#707372', alpha=0.5)),
                      name='associated with') +
  ggrepel::geom_text_repel(aes(label=label))

# Factor 2
g3 <- df.plot.mefisto %>% 
  mutate(type=case_when(p.value < 0.1 & estimate > 0~'fd',
                        p.value < 0.1 & estimate < 0~'bd',
                        TRUE~'other')) %>% 
  mutate(type=factor(type, levels = c('fd', 'bd', 
                                      'other'))) %>% 
  filter(prevalence>0.05) %>% 
  filter(variable=='diet') %>% 
  arrange(type) %>% 
  mutate(label=case_when(abs(estimate) > 0.48 ~ genus,
                         TRUE~'')) %>% 
  ggplot(aes(x=-estimate, y=Factor2, col=type)) + 
    geom_point() + 
    theme_bw() + 
    xlab('LME model coefficient') + 
    ylab('Factor 2 weight') + 
    theme(panel.grid.minor = element_blank()) + 
    scale_colour_manual(values=c(c("#1f78b4", "#b2df8a"), 
                                 alpha('#707372', alpha=0.5)),
                        name='associated with') +
    ggrepel::geom_text_repel(aes(label=label))

g.all <- plot_grid(g1, g2, g3, nrow = 1)
ggsave(g.all, filename = 'compare_mefisto_supervised.pdf',
       width = 19, height = 6)

# same for CTF?
g1.ctf <- df.plot.ctf %>% 
  mutate(type=case_when(p.value < 0.1 & estimate > 0~'cesarean',
                        p.value < 0.1 & estimate < 0~'vaginal',
                        TRUE~'other')) %>% 
  mutate(type=factor(type, levels = c('cesarean', 'vaginal', 
                                      'other'))) %>% 
  filter(prevalence>0.05) %>% 
  filter(variable=='delivery') %>% 
  arrange(type) %>% 
  mutate(label=case_when(abs(estimate) > 0.48 ~ genus,
                         TRUE~'')) %>% 
  ggplot(aes(x=-estimate, y=Factor1, col=type)) + 
  geom_point() + 
  theme_bw() + 
  xlab('LME model coefficient') + 
  ylab('Factor 1 weight') + 
  theme(panel.grid.minor = element_blank()) + 
  scale_colour_manual(values=c(c("#d95f02", "#e6ab02"), 
                               alpha('#707372', alpha=0.5)),
                      name='associated with') +
  ggrepel::geom_text_repel(aes(label=label))

g2.ctf <- df.plot.ctf %>% 
  mutate(type=case_when(p.value < 0.1 & estimate > 0~'fd',
                        p.value < 0.1 & estimate < 0~'bd',
                        TRUE~'other')) %>% 
  mutate(type=factor(type, levels = c('fd', 'bd', 
                                      'other'))) %>% 
  filter(prevalence>0.05) %>% 
  filter(variable=='diet') %>% 
  arrange(type) %>% 
  mutate(label=case_when(abs(estimate) > 0.48 ~ genus,
                         TRUE~'')) %>% 
  ggplot(aes(x=-estimate, y=Factor2, col=type)) + 
  geom_point() + 
  theme_bw() + 
  xlab('LME model coefficient') + 
  ylab('Factor 2 weight') + 
  theme(panel.grid.minor = element_blank()) + 
  scale_colour_manual(values=c(c("#1f78b4", "#b2df8a"), 
                               alpha('#707372', alpha=0.5)),
                      name='associated with') +
  ggrepel::geom_text_repel(aes(label=label))

g.ctf <- cowplot::plot_grid(g1.ctf, g2.ctf)
ggsave(g.ctf, filename = 'ctf_scatters.pdf', 
       width = 11, height = 4)

# ##############################################################################
# viz tree in iTOL

# select strong weights for MEFISTO and CTF
# export mean weights as tsv
# which weights are significant?
# put everything together in iTOL

# calculate enrichments for genera using one-sided Wilcoxon
.f_genus_enrichment <- function(weights){
  weights <- weights %>% 
    filter(!is.na(genus))
  df.enrich <- tibble(genus=character(0), direction=character(0), 
                      factor=double(0),
                      effect_size=double(0), p_value=double(0))
  for (i in unique(weights$genus)){
    # positive, factor 1
    temp <- weights %>% 
      filter(Factor1>0) %>% 
      mutate(test=genus==i)
    if (sum(temp$test) >= 3){
      t <- wilcox.test(temp$Factor1~temp$test, alternative='less')
      df.enrich <- df.enrich %>% 
        add_row(genus=i, direction='up',
                effect_size=mean(temp %>% filter(test) %>% pull(Factor1)),
                factor=1, p_value=t$p.value)
    }
    # negative, factor 1
    temp <- weights %>% 
      filter(Factor1<0) %>% 
      mutate(test=genus==i)
    if (sum(temp$test) >= 3){
      t <- wilcox.test(temp$Factor1~temp$test, alternative='greater')
      df.enrich <- df.enrich %>% 
        add_row(genus=i, direction='down',
                effect_size=mean(temp %>% filter(test) %>% pull(Factor1)),
                factor=1, p_value=t$p.value)
    }
    # positive, factor 2
    temp <- weights %>% 
      filter(Factor2>0) %>% 
      mutate(test=genus==i)
    if (sum(temp$test) >= 3){
      t <- wilcox.test(temp$Factor2~temp$test, alternative='less')
      df.enrich <- df.enrich %>% 
        add_row(genus=i, direction='up',
                effect_size=mean(temp %>% filter(test) %>% pull(Factor2)),
                factor=2, p_value=t$p.value)
    }
    # negative, factor 2
    temp <- weights %>% 
      filter(Factor2<0) %>% 
      mutate(test=genus==i)
    if (sum(temp$test) >= 3){
      t <- wilcox.test(temp$Factor2~temp$test, alternative='greater')
      df.enrich <- df.enrich %>% 
        add_row(genus=i, direction='down',
                effect_size=mean(temp %>% filter(test) %>% pull(Factor2)),
                factor=2, p_value=t$p.value)
    }
  }
  
  return(df.enrich)
}

enrichment.mefisto <- .f_genus_enrichment(weights.mefisto)
enrichment.ctf <- .f_genus_enrichment(weights.ctf)

# mefisto main tree
mefisto.red <- enrichment.mefisto %>% 
  select(-p_value) %>% 
  group_by(genus) %>% 
  # effect size filter
  filter(sum(abs(effect_size) > 0.2) >= 1) %>% 
  pivot_wider(names_from = c(direction, factor), values_from = effect_size,
              values_fill = 0) %>% 
  mutate(genus=str_remove(genus, 'g__')) %>% 
  filter(genus!='')
write_csv(mefisto.red, file = 'tree_mefisto_weights.csv')

# adj.p.values
enrichment.mefisto %>% 
  group_by(genus) %>%
  filter(sum(abs(effect_size) > 0.2) >= 1) %>%
  mutate(p.adj=p.adjust(p_value, method='fdr')) %>% 
  arrange(p.adj) %>% 
  filter(p.adj < 0.05)

# compare mefisto/CTF tree
ctf.red <- enrichment.ctf %>% 
  select(-p_value) %>% 
  filter(genus %in% mefisto.red$genus) %>% 
  pivot_wider(names_from = c(direction, factor), values_from = effect_size,
              values_fill = 0) %>% 
  mutate(genus=str_remove(genus, 'g__')) %>% 
  filter(genus!='')
write_csv(ctf.red, file = 'tree_ctf_weights.csv')

# adj.p.values
enrichment.ctf %>% 
  filter(genus %in% mefisto.red$genus) %>% 
  mutate(p.adj=p.adjust(p_value, method='fdr')) %>% 
  arrange(p.adj) %>% 
  filter(p.adj < 0.05)
