matrix.list<-generateMatrix(gwasFile = "genfiles\\hcm.fix.2024.format.shiny.tsv.gz", configFile = "hcmr.config.txt", pos.split = 3E6)

write.table(matrix.list$melt, "3mbp.HCMFIX1.melt.txt", col.names = T, row.names = F, sep = "\t", quote = FALSE)
write.table(matrix.list$annot, "3mbp.HCMFIX1.annot.txt", col.names = T, row.names = F, sep = "\t", quote = FALSE)
write.table(matrix.list$chr.len, "3mbp.HCMFIX1.chrlen.txt", col.names = T, row.names = F, sep = "\t", quote = FALSE)
write.table(matrix.list$annot.melt, "3mbp.HCMFIX1.annotmelt.txt", col.names = T, row.names = F, sep = "\t", quote = FALSE)
write.table(matrix.list$pos.df, "3mbp.HCMFIX1.positions.txt", col.names = T, row.names = F, sep = "\t", quote = FALSE)

# for(i in 1:22){
#   matrix.list<-generateMatrix(gwasFile = paste0("hcm.fix.2024.format.shiny.chr",i,".tsv.gz"), 
#                               configFile = "..//hcmr.config.txt", pos.split = 1E5)
#   
#   ##write.table(matrix.list$matrix, paste0("100kbp.HCMFIX1.chr",i,".matrix.txt"), col.names = F, row.names = F, sep = "\t", quote = FALSE)
#   write.table(matrix.list$melt, paste0("100kbp.HCMFIX1.chr",i,".melt.txt"), col.names = T, row.names = F, sep = "\t", quote = FALSE)
#   write.table(matrix.list$annot, paste0("100kbp.HCMFIX1.chr",i,".annot.txt"), col.names = T, row.names = F, sep = "\t", quote = FALSE)
#   write.table(matrix.list$chr.len, paste0("100kbp.HCMFIX1.chr",i,".chrlen.txt"), col.names = T, row.names = F, sep = "\t", quote = FALSE)
#   write.table(matrix.list$annot.melt, paste0("100kbp.HCMFIX1.chr",i,".annotmelt.txt"), col.names = T, row.names = F, sep = "\t", quote = FALSE)
#   write.table(matrix.list$pos.df, paste0("100kbp.HCMFIX1.chr",i,".positions.txt"), col.names = T, row.names = F, sep = "\t", quote = FALSE)
# }
