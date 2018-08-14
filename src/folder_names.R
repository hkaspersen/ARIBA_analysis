# Identifies folder names for extracting report.tsv
folder_names <- function(filepath) {
  folders <- list.files(path = filepath, pattern = "_ariba")
  return(folders)
}