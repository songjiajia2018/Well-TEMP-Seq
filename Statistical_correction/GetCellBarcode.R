args <- commandArgs(T)
rds_file <- readRDS(args[1])
mtx <- as.matrix(rds_file)
cell_filter <- colnames(mtx)
prefix <- strsplit(args[1],'_')[[1]][1]
outname <- paste0(prefix,'_cell_filter_name.txt')
write.table(cell_filter,file = outname, quote = F, col.names = F, row.names = F)