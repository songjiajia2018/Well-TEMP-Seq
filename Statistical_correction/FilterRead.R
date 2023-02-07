args <- commandArgs(T)
samplename <- args[1]

info <- read.table(paste0(samplename,'_MAPQ20_barcode_and_gene.txt'),header = F)
info_freq <- table(info$V2)

discard_barcode <- info_freq[which(info_freq==1)]

count_file <- read.table(paste0(samplename,'_count.txt'),header = F)
count_filter <-  count_file[-which(count_file$V1 %in% names(discard_barcode)),]
write.table(count_filter,paste0(samplename,'_filter_count.txt'),col.names = F,row.names = F)