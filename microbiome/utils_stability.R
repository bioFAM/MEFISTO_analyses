load_files <- function(outdir, pattern, exl_pattern = NULL){
  idxshift <- as.numeric(grepl("_", pattern)) # shift index for strsplit
  files2load <- list.files(outdir)

  files2load <- files2load[grepl(pattern, files2load)]
  if(!is.null(exl_pattern)){
    files2load <- files2load[!grepl(exl_pattern, files2load)]
  }
  model_list <- lapply(files2load, function(fnm) {
    load_model(file.path(outdir,fnm))
  })
  Z <- lapply(seq_along(files2load), function(i) 
    Reduce(rbind, get_factors(model_list[[i]])) %>% as.data.frame() %>%
      rownames_to_column("sample_id") %>%
      mutate(seed = gsub("\\.hdf5","",strsplit(files2load[i], "_")[[1]][3 +idxshift]),
             miss_times =  strsplit(files2load[i], "_")[[1]][2 + idxshift])
  ) %>% bind_rows()
  return(Z)
}

load_files_CTF <- function(outdir, patternZ = "CTF_straj"){ 
  files2load <- list.files(outdir)
  Zfiles <- files2load[grepl(patternZ, files2load)]

  Z <- lapply(Zfiles, function(fnm) {
    df <- read.csv(file.path(outdir,fnm))
    df <- mutate(df, sample_id = paste(subject_id, month, sep= "_")) %>%
      select(sample_id, Factor1 = PC2, Factor2 = PC2)
    df %>% mutate(seed = strsplit(fnm, "_")[[1]][3],
                  miss_times =  gsub("\\.csv","",strsplit(fnm, "_")[[1]][2]))
  }) %>% bind_rows()
  
  return(Z)
}

load_files_others <- function(outdir, patternZ){
  files2load <- list.files(outdir)
  method <- strsplit(patternZ, "_")[[1]][1]
  Zfiles <- files2load[grepl(patternZ,files2load)]
  Z <- lapply(Zfiles, function(fnm) {
    df <- read.csv(file.path(outdir,fnm))
    colnames(df) <- c("sample_id","Factor1", "Factor2")
    df %>% mutate(seed = strsplit(fnm, "_")[[1]][3],
                  miss_times =  gsub("\\.csv","",strsplit(fnm, "_")[[1]][2]),
                  method = method)
  })%>% bind_rows()
  return(Z)
}