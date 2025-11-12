#!/usr/bin/env Rscript
library(GENESPACE)

# Function to parse command-line arguments
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  # Default values
  params <- list(
    working_dir = "missing",
    mcscanx_path = "missing",
    ploidy = 1,
    ref_genome = "hap1",
    sameChr = FALSE,
    threads = 1,
    output = "test.tsv"
  )

  # Parse arguments
  i <- 1
  while (i <= length(args)) {
    if (args[i] %in% c("-w", "--working_dir")) {
      params$working_dir <- args[i + 1]
      i <- i + 2
    } else if (args[i] %in% c("-m", "--mcscanx_path")) {
      params$mcscanx_path <- args[i + 1]
      i <- i + 2
    } else if (args[i] %in% c("-p", "--ploidy")) {
      params$ploidy <- as.integer(args[i + 1])
      i <- i + 2
    } else if (args[i] %in% c("-r", "--ref_genome")) {
      params$ref_genome <- args[i + 1]
      i <- i + 2
    } else if (args[i] %in% c("-c", "--sameChr")) {
      params$sameChr <- as.logical(args[i + 1])
      i <- i + 2
    } else if (args[i] %in% c("-t", "--threads")) {
      params$threads <- as.integer(args[i + 1])
      i <- i + 2
    } else if (args[i] %in% c("-o", "--output")) {
      params$output <- args[i + 1]
      i <- i + 2
    } else {
      i <- i + 1
    }
  }

  return(params)
}

# Parse command-line arguments
opt <- parse_args()

# Initialize GENESPACE parameters
gpar <- init_genespace(
  wd = opt$working_dir,
  path2mcscanx = opt$mcscanx_path,
  ploidy = opt$ploidy,
  onlySameChrs = opt$sameChr,
  nCores = opt$threads
)

# Run GENESPACE
out <- run_genespace(gpar, overwrite = FALSE)

# Query pangenes
pangenes <- query_pangenes(
  gsParam = out, refGenome = opt$ref_genome
)

# Convert the result to a data frame
pangenes <- as.data.frame(pangenes)

# Select the columns starting with "hap" (without dplyr)
hap_cols <- grep("^hap", names(pangenes), value = TRUE)
pangenes <- pangenes[, hap_cols, drop = FALSE]

# Convert lists to comma-separated strings
pangenes[] <- lapply(pangenes, function(col) {
  if (is.list(col)) {
    sapply(col, toString)
  } else {
    col
  }
})

# Fill empty cells with NA
pangenes[pangenes == ""] <- NA

# Write the data frame to a file
write.table(pangenes, file = opt$output, sep = "\t", quote = FALSE, row.names = FALSE)

# Print the result
print(head(pangenes))
