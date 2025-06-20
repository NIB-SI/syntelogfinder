library(GENESPACE)
library(dplyr)
library(optparse)

# Define command-line arguments
option_list <- list(
  make_option(c("-w", "--working_dir"), type = "character",
              default = "/scratch/nadjafn/potato-allelic-orthogroups/results/input_genespace/",
              help = "Working directory [default= %default]", metavar = "character"),
  make_option(c("-m", "--mcscanx_path"), type = "character",
              default = "/scratch/nadjafn/MCScanX",
              help = "Path to MCScanX [default= %default]", metavar = "character"),
  make_option(c("-p", "--ploidy"), type = "integer",
              default = 1,
              help = "Ploidy level [default= %default]", metavar = "integer"),
  make_option(c("-r", "--ref_genome"), type = "character",
              default = "hap1",
              help = "Reference genome for query_pangenes [default= %default]", metavar = "character"),
  make_option(c("-c", "--sameChr"), type = "logical",
              default = TRUE,
              help = "Use only same chromosomes to find pangenes [default= %default]", metavar = "logical"),
  make_option(c("-t", "--threads"), type = "integer",
              default = 1,
              help = "Number of cores for running genespace", metavar = "integer"),
  make_option(c("-o", "--output"), type = "character",
              default = "test.tsv",
              help = "Output of pangenes [default= %default]", metavar = "character")
)

# Parse command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

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

# Select the hap1, hap2, hap3, and hap4 columns
pangenes <- pangenes[, c("hap1", "hap2", "hap3", "hap4")]

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
write.table(pangenes, file = opt$output , sep = "\t", quote = FALSE, row.names = FALSE)

# Print the result
print(head(pangenes))

