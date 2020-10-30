library(tidyverse)
library(magrittr)
library(DESeq2)

# path pointing to the count data from https://www.nature.com/articles/s41586-019-1338-5
path_input <- "data/counts/"

# path pointing to suppl information from https://www.nature.com/articles/s41586-019-1338-5 etc
path_meta <- "data/meta/"

# path to save normalized counts to
path_output <- "data/normalized_counts/"
if(!dir.exists(path_output)) dir.create(path_output)

species <- c("Human", "Mouse", "Opossum", "Rabbit", "Rat")

# file containing the orthologue information between species (from Ensembl BioMart (E77))
orthologues <- read.table(paste0(path_meta, "orthologues_species.txt"), header = TRUE)

result <- lapply(species, function(s) {
  print(paste("Processing counts for", s))
  
  # csv files pointing to count .txt-file for each sample
  df_specie <- read.table(paste0(path_input, s, ".csv"), header = TRUE, sep = ",")
  
  sampleTable <- data.frame(sampleName = df_specie$sample, fileName = df_specie$files, 
                            species = df_specie$species, tissue = df_specie$tissue, time = df_specie$time)
  
  # read counts into DESeq2
  DESeq2Table <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, 
                                            directory = paste0(path_input, s, "_counts"), 
                                            design = ~ 1)
  
  # Library size correction 
  DESeq2Table <- estimateSizeFactors(DESeq2Table)
  
  # Normalization
  normalized <- varianceStabilizingTransformation(DESeq2Table)
  specie_normalized <- rownames_to_column(as.data.frame(assay(normalized)))
  
  specie_values <- gather(specie_normalized, key = "sample", value = "expression", -rowname)
  specie_values %<>% mutate(gene = rowname)
  specie_values %<>% select(-rowname)
  
  specie_values <- specie_values %>% left_join(sampleTable, by = c("sample" = "sampleName"))

  # rename human time points that are inconsistent between tissues
  if (s == "Human"){
    specie_values$time[specie_values$time == "Senior"] <- "senior"
    specie_values$time[specie_values$time == "youngTeenager"] <- "teenager"
    specie_values$time[specie_values$time == "oldTeenager"] <- "teenager"
  }
  
  #filter orthologues
  specie_values %<>% filter(gene %in% orthologues[ , paste0(s, "_ID")])
  
  # Average replicate samples from same time point and tissue
  mean_replicates <- specie_values %>% group_by(gene, tissue, time) %>% dplyr::summarise(
    species = s,
    number_of_replicates = n(),
    replicate_mean = mean(expression),
    variance = var(expression)) %>% ungroup()
  
}) %>% bind_rows()

# save results
write.table(result, paste0(path_output, "normalized_orthologues_all_species.txt"), 
            quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)


