#!/usr/bin/Rscript
args <- commandArgs(T)
a=read.table(args[1], header=F, stringsAsFactors=F)
x=cumsum(a$V1)
x=x/max(x)
pdf_name <- paste0(args[2], sep = '_', 'cell_selection.pdf')
pdf(pdf_name)
plot(1:length(x), x, type='l', col="blue", xlab="cell barcodes sorted by number of reads [descending]",ylab="cumulative fraction of reads", xlim=c(1,10000))
dev.off()
