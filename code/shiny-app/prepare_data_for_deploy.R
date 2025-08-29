#!/usr/bin/env Rscript
# Prepare minimal data bundle for deploying the Shiny app to shinyapps.io

message("Preparing data subset for deployment...")

project_root <- normalizePath(file.path(dirname(sys.frame(1)$ofile %||% "./"), "..", ".."))
src_data <- file.path(project_root, "data")
dest_data <- file.path(dirname(sys.frame(1)$ofile %||% getwd()), "data")

if(!dir.exists(src_data)) stop("Source data directory not found: ", src_data)
if(!dir.exists(dest_data)) dir.create(dest_data, recursive = TRUE)

copy_file <- function(rel){
  from <- file.path(src_data, rel)
  to <- file.path(dest_data, rel)
  dir.create(dirname(to), showWarnings = FALSE, recursive = TRUE)
  file.copy(from, to, overwrite = TRUE)
  message("Copied ", rel)
}

# Core metadata files
core_files <- c(
  "trial_metadata.xlsx",
  file.path("size", "trial_size.xlsx")
)
lapply(core_files, copy_file)

# Plate files (exclude date 20250630 explicitly to shrink bundle)
plate_dir <- file.path(src_data, "plate-files")
dates <- list.dirs(plate_dir, full.names = FALSE, recursive = FALSE)
dates <- dates[dates != "20250630"]

for(d in dates){
  xls <- list.files(file.path(plate_dir, d), pattern = "[.]xlsx$", full.names = FALSE)
  for(f in xls){
    rel <- file.path("plate-files", d, f)
    copy_file(rel)
  }
}

message("Done. Data copied to: ", dest_data)
message("Files copied: ")
print(list.files(dest_data, recursive = TRUE))

invisible(TRUE)
