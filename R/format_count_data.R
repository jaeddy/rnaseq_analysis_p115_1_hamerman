library(dplyr)
library(readr)
library(stringr)
library(biomaRt)

# Specify flowcell & project
fcDir <- "/Volumes/genomics/Illumina/150707_D00565_0090_AHMMF2ADXX"
processedDir <- "Project_P115-1Processed"

# Get file path for counts file
countsDir <- file.path(fcDir, processedDir, "counts")

countsFile <- list.files(countsDir) %>% 
    str_extract(".*combined_counts.csv") %>% 
    na.omit() %>% 
    file.path(countsDir, .)

# Read in combined counts file
counts <- read_csv(countsFile)

# Get Entrez Gene IDs corresponding to Ensembl Gene IDs
mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
ens2Gene <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"), 
                  filters = "ensembl_gene_id", 
                  values = counts$geneName, mart = mart)
ens2Gene = ens2Gene[match(counts$geneName, ens2Gene$ensembl_gene_id), ]

# Insert Entrez Gene IDs for genes in counts data
countsMgi <- counts %>% 
    mutate(mgiSymbol = ens2Gene$mgi_symbol)
colOrder <- c(names(countsMgi)[1], names(countsMgi)[length(countsMgi)],
              names(countsMgi)[2:(length(countsMgi)-1)])
countsMgi <- countsMgi %>% 
    dplyr::select(one_of(colOrder)) %>% 
    mutate(mgiSymbol = ifelse(is.na(mgiSymbol), "NA", mgiSymbol))

# Write the new file
countsFileNew <- str_replace(countsFile, ".csv", "_wMgiSymbol.csv")
write_csv(countsMgi, countsFileNew)
