# panther_prep.R
args <- commandArgs(trailingOnly = TRUE)

# Check if input arguments are provided
if (length(args) < 2) {
  stop("Usage: Rscript panther_prep.R input.txt Label")
}

input_file <- args[1]
label <- args[2]

# Read the tab-delimited input file
df <- read.delim(input_file, stringsAsFactors = FALSE, header = TRUE)

# Add a new column with the label, without naming it
df <- cbind(df, label)

# Construct output filename by appending "_fixed.txt"
output_file <- sub("\\.txt$", "_fixed.txt", input_file)

# Write to output file, no quotes, tab-separated, no row names, and no column names
write.table(df, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

cat("Output written to:", output_file, "\n")