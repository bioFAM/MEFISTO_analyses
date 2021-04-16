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