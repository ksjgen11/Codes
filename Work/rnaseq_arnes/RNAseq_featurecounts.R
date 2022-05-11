library(tidyverse)
library(Rsamtools)
library(Rsubread)
llibrary(Mus.musculus)



file_storage <- '~/data/rnaseq_arnes/files/star_mapping/'
gtf_file <- '~/data/rnaseq_arnes/index/Mus_musculus.GRCm38.102.gtf.gz'
kl3 <- '~/data/rnaseq_arnes/files/KL003.bam'

file_list <- list.files(path = file_storage, pattern = '*.bam$')
filenames <- file.path(path = file_storage, file_list)
bamfiles <- BamFileList(filenames, yieldSize = 2000000)


# make reference 

# Get gtf files from ENSENBL ftp download site
# http://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz

gtffile <- gtf_file


# Using Rsubread package

fc <- featureCounts(files=filenames,
                    annot.ext=gtffile,
                    isGTFAnnotationFile=TRUE,
                    isPairedEnd=TRUE)


fc_kl3 <- featureCounts(files=kl3, 
                    annot.ext=gtffile,
                    isGTFAnnotationFile=TRUE,
                    isPairedEnd=FALSE)

count <- cbind(fc$counts, fc_kl3$counts)

write.table(count, file = paste0(dataset, 'raw_count.txt'))

