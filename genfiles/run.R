matrix.list<-generateMatrix(gwasFile = "genfiles\\cad.add.160614_manhformat_rsid_geneinfo_pvalmatch.txt.gz", 
                            configFile = "cad.config.txt", 
                            pos.split = 3E6)

write.table(matrix.list$melt, "3mbp.cad.melt.txt", 
            col.names = T, row.names = F, sep = "\t", quote = FALSE)
write.table(matrix.list$annot, "3mbp.cad.annot.txt", 
            col.names = T, row.names = F, sep = "\t", quote = FALSE)
write.table(matrix.list$chr.len, "3mbp.cad.chrlen.txt", 
            col.names = T, row.names = F, sep = "\t", quote = FALSE)
write.table(matrix.list$annot.melt, "3mbp.cad.annotmelt.txt", 
            col.names = T, row.names = F, sep = "\t", quote = FALSE)
write.table(matrix.list$pos.df, "3mbp.cad.positions.txt", 
            col.names = T, row.names = F, sep = "\t", quote = FALSE)

for(i in 1:22){
  in.df<-read.table(file = "genfiles\\cad.add.160614_manhformat_rsid_geneinfo_pvalmatch.txt.gz", sep="\t", header = T)
  in.df<-in.df[in.df$chr==i,]
  write.table(x = in.df, file = "gwas.tmp", append = F, quote = F, sep = "\t", row.names = F)
  
  matrix.list<-generateMatrix(gwasFile = "gwas.tmp",
                              configFile = "cad.config.txt", pos.split = 1E5)

  write.table(matrix.list$melt, paste0("100kbp.cad.chr",i,".melt.txt"), 
              col.names = T, row.names = F, sep = "\t", quote = FALSE)
  write.table(matrix.list$annot, paste0("100kbp.cad.chr",i,".annot.txt"), 
              col.names = T, row.names = F, sep = "\t", quote = FALSE)
  write.table(matrix.list$chr.len, paste0("100kbp.cad.chr",i,".chrlen.txt"), 
              col.names = T, row.names = F, sep = "\t", quote = FALSE)
  write.table(matrix.list$annot.melt, paste0("100kbp.cad.chr",i,".annotmelt.txt"), 
              col.names = T, row.names = F, sep = "\t", quote = FALSE)
  write.table(matrix.list$pos.df, paste0("100kbp.cad.chr",i,".positions.txt"), 
              col.names = T, row.names = F, sep = "\t", quote = FALSE)
}
