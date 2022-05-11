library(tidyverse)

setwd('~/data/atacseq/')
basic <- '~/data/atacseq/'

## normalize count

peak <- '~/data/atacseq/new_result/'
scale <- '~/data/atacseq/new_result/scale_factor/'

file_list <-list.files(path = peak, pattern = '*_new.txt')


datafile <- lapply(paste0(peak, file_list), FUN = read.table, col.names = c('merged', 'counts'))
names(datafile) <- file_list


data_count <- as.data.frame(datafile) %>%
  dplyr::select(ends_with('counts')) %>%
  as.matrix()


rownames(data_count) <- datafile[[1]]$merged

colnames(data_count) <- c('si673_12h_1','si673_12h_2','si673_12h_3',
                          'si673_24h_1', 'si673_24h_2','si673_24h_3', 
                          'si673_4h_1', 'si673_4h_2', 'si673_4h_3', 
                          'siNT_1', 'siNT_2', 'siNT_3')

NormFactor <- edgeR::calcNormFactors(object = data_count, method = "TMM")

## if you prefer to use the DESeq2 strategy use method="RLE" instead

## raw library size:
LibSize <- colSums(data_count)

## calculate size factors:
SizeFactors <- NormFactor * LibSize / 1000000

## Reciprocal, please read section below:   
SizeFactors.Reciprocal <- 1/SizeFactors




n <- 1
for (d in SizeFactors.Reciprocal){
  write.table(d, file = paste0(scale, 'scale_factor_',names(datafile)[n]),
              row.names = FALSE, col.names = FALSE)
  n = n + 1
}


write.table(data_count, file = paste0(peak, 'raw_count_new.txt'))
