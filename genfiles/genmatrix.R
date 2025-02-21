library(reshape2)

generateMatrix<-function(gwasFile="cad.add.160614_manhformat_rsid_geneinfo_pvalmatch.txt.gz", configFile = "config.txt",
                         FDR=1E-3, MAF=0.05,
                         pos.split=3E6,pval.split=0.125,max.pval=20){
  message("Reading the GWAS file (",gwasFile,")...")
  d<-read.table(gwasFile, header=T)
  d$maf[d$maf > 0.5]<-(1-(d$maf[d$maf>0.5]))
  
  message("Reading the config file (",configFile,")...")
  config<-read.table(configFile,sep="\t", header =T,stringsAsFactors = F, skip=10)
  
  ## generate the pvalue bins (using the pval.split parameter)
  pvals<-seq(from=0, to=max.pval, by=pval.split)
  pvals.cells.index<-data.frame(id=1:length(pvals),LP=pvals,UP=c(pvals[2:length(pvals)],max.pval))
  
  ## create an empty matrix
  final<-matrix(0, nrow = length(pvals), ncol = 0)
  max.cellcount<-0
  
  ## create an empty matrix
  annotation.matrix<-matrix("", nrow = length(pvals), ncol = 0)
  
  chr.matrix.len<-as.data.frame(matrix(nrow=length(unique(d$chr)),ncol=4))
  names(chr.matrix.len)<-c("chr","length", "cumm", "mid")
  len.idx<-1
  
  pos.idx<-1
  pos.df<-data.frame(idx=c(-1), chr=c(-1), st=c(-1), en=c(-1))
  
  for(chr in min(d$chr):max(d$chr)){
    message(paste("Processing chromosome ",chr,"\r", sep=""),appendLF = F)
    
    ## extract chromosome specific data from GWAS input
    chr.slice<-d[d$chr==chr,]
    
    ## Generate the chromosome position bins (using the pos.split parameter)
    chunks<-seq(from=min(chr.slice$pos), to=max(chr.slice$pos), by=pos.split)
    chunks[length(chunks)+1]<-max(chr.slice$pos)
    
    ## create the matrix for this chromosome (pvals * position)
    mdat<-matrix(0, nrow = length(pvals), ncol = length(chunks)-1)
    annot.mdat<-matrix(0, nrow = length(pvals), ncol = length(chunks)-1)
    
    for (i in 1:(length(chunks)-1)){
      slice<-chr.slice[chr.slice$pos >= chunks[i] & chr.slice$pos < chunks[i+1],]
      
      st<-chunks[i] 
      en<-chunks[i+1]
      
      pos.df<-rbind(pos.df, data.frame(idx=c(pos.idx), chr=c(chr), st=c(st), en=c(en)))
      pos.idx<-pos.idx+1
      
      for (j in 1:length(pvals)){
        ## get slice of gwas input for each of the pvalue bins (and position bins)
        if(j == length(pvals)){ ## last element in pval array - include all variants gt then this
          p.val.slice<-slice[-log10(slice$pvalue) >= pvals[j],]  
        }
        else{ ## otherwise slice between the current and next values in the pval array
          p.val.slice<-slice[-log10(slice$pvalue) >= pvals[j] & -log10(slice$pvalue) < pvals[j+1],]
        }      
        
        ## determine the number of variants in the slice
        len<-dim(p.val.slice)[1]
        idx<-0
        
        ## determine the number of variants in the slice, which have the HIGH consequence flag set to 1
        conseq.len<-dim(p.val.slice[p.val.slice$conseq==1,])[1] ## test if any high consequences
        ## determine the number of variants in the slice, which have MAF less than 5%
        maf.len<-dim(p.val.slice[p.val.slice$maf<=MAF,])[1] ## test if any MAF < 5%
        
        annot.mdat[j,i]<-""
        
        if(len == 0){
          if(pvals[j] <= -log10(max(d$pvalue))){
            ## determine whether the chromosome is odd and even
            if((chr %% 2) != 0){idx<-config$idx[config$type=="oddchr"]}
            else{idx<-config$idx[config$type=="evenchr"]}
          }
          else{
            idx<-0
          }
        } ## the case with no variants - blank cell.
        else{
          ## this defines the region which should be greyed out
          if(pvals[j] <= -log10(FDR)){
            ## determine whether the chromosome is odd and even
            if((chr %% 2) != 0){idx<-config$idx[config$type=="oddchr"]}
            else{idx<-config$idx[config$type=="evenchr"]}
          }
          else{
            for(k in 1:length(config$idx)){
              if(config$type[k] == "val"){ ## Check the config is a valid type.
                len.chk<-FALSE
                conseq.chk<-FALSE
                maf.chk<-FALSE
                
                ## check counts
                if(len >= config$min.count[k]){ ## check that the current cell length is gt or eq than the config min count
                  if(is.na(config$max.count[k])){ ## if the config max is NA then accept length condition
                    len.chk<-TRUE
                  }
                  else if(len <= config$max.count[k]){ ## if the current cell length is lt the config max count then accept the length condition
                    len.chk<-TRUE
                  }
                }
                
                ## check HIGH impact
                if(config$conseq[k] == TRUE && conseq.len > 0){conseq.chk<-TRUE} ## accept if high impact is active in config and there are 1 or more high impact variants within the cell
                else if(config$conseq[k] == FALSE && conseq.len == 0){conseq.chk<-TRUE} ## accept if high impact is inactive and there are 0 high impact variants within the cell.
                
                ## check MAF
                if(config$maf[k] == TRUE && maf.len > 0){maf.chk<-TRUE} ## accept if MAF 5% active and there is one or more variant with MAF lt 5% in the cell
                else if(config$maf[k] == FALSE && maf.len == 0){maf.chk<-TRUE} ## accept if MAF 5% active and there are 0 variants with MAF lt 5% in the cell
                
                if(len.chk==TRUE && conseq.chk==TRUE && maf.chk==TRUE){ ## if all three clauses are correct then accept the idx for the config and exit the for loop
                  idx<-config$idx[k]
                  break ## found config which fulfils cell criteria accept and exit group.
                }
                
              }
            }
            
            if(idx == 0){ ## if no idx was found in the annotations, generate a remaining idx (max idx in annotations + 1)
              idx<-max(config$idx)+1; ## assign to others
            }
            if(len > max.cellcount){max.cellcount<-len} ## assign the maximum cell count if length of current cell is gt than current
          
            ## get the best SNP and information
            bestpval<-min(p.val.slice$pvalue)
            best.all<-p.val.slice[p.val.slice$pvalue == bestpval,]
            bestsnp<-as.character(best.all$snp[1])
            bestgenes<-as.character(best.all$genes[1])
            bestchr<-as.character(best.all$chr[1])
            bestpos<-as.character(best.all$pos[1])
           
            conseq.str<-""
            if(conseq.len > 0){
              conseq.slice<-p.val.slice[p.val.slice$conseq==1,]
              conseq.pval<-min(conseq.slice$pvalue)
              conseq.snp<-conseq.slice$snp[conseq.slice$pvalue == conseq.pval][1]
              conseq.str<-paste("</br><b>High impact:</b> ",
                                conseq.snp,
                                sep="")
              
              ##print(conseq.str)
            }

            maf.str<-""
            if(maf.len > 0){
              maf.slice<-p.val.slice[p.val.slice$maf<=MAF,]
              maf.maf<-min(maf.slice$maf)
              maf.snp<-maf.slice$snp[maf.slice$maf == maf.maf][1]
              maf.str<-paste("</br><b>Lowest MAF:</b> ",
                                maf.snp,
                                " (MAF=",formatC(maf.maf, format = "E", digits = 2),")",
                                sep="")
              ##print(maf.str)
              
            }
                         
            # annot.mdat[j,i]<-paste("<b>'#'Variants:</b> ", len,
            #                        "</br>","<b>Best Variant:</b> ",bestsnp,
            #                        "</br>","<b>Best p-value:</b> ",bestpval,
            #                        "</br>","<b>chr:bp:</b> ", chr,":",bestpos,
            #                        "</br>","<b>Prox. gene(s):</b> ",bestgenes, 
            #                        conseq.str,maf.str,sep="")  
            annot.mdat[j,i]<-paste0(len, "</br>",
                                   bestsnp, "</br>",
                                   bestpval, "</br>",
                                   chr,":",bestpos, "</br>",
                                   bestgenes, "</br>",
                                   conseq.str, "</br>",
                                   maf.str)  
          }
        }
        
        mdat[j,i]<-idx
      }
    }
    
    ## cell sizes of the chromosomes
    chr.matrix.len$chr[len.idx]<-chr
    chr.matrix.len$length[len.idx]<-dim(mdat)[2]
    chr.matrix.len$cumm[len.idx]<-dim(mdat)[2]+ dim(final)[2]
    chr.matrix.len$mid[len.idx]<-chr.matrix.len$cumm[len.idx] - (chr.matrix.len$length[len.idx] / 2)
    
    len.idx<-(len.idx+1)
    
    ## bind to the final matrix
    final<-cbind(final,mdat)
    annotation.matrix<-cbind(annotation.matrix,annot.mdat)
  }
  
  message("\nMelting matrix...")
  m<-melt(final)
  annot.m<-melt(annotation.matrix)
  
  message("Generation cell annotations...")
  cell.annot<-data.frame(idx=config$idx, name=rep("", length(config$idx)), stringsAsFactors = F)
  
  for(j in 1:length(config$idx)){
    if(config$type[j]=="val"){
      max<-config$max.count[j]
      
      if(is.na(max)){max<-ceiling(max.cellcount/100)*100}
      
      if(config$min.count[j] == max){
        text<-paste(config$min.count[j], sep="")
      }
      else{
        text<-paste(config$min.count[j], "-",max, sep="")
      }
      
      if(config$maf[j]== TRUE && config$conseq[j] == TRUE){text<-paste(text, "(BOTH)")}
      else if(config$maf[j]  == TRUE){text<-paste(text, "(MAF)")}
      else if(config$conseq[j] == TRUE){text<-paste(text, "(HIGH impact)")}
    
      cell.annot$name[j]<-text
      }
    else{
      cell.annot$name[j]<-config$type[j]
    }
  }

  message(paste("Generated matrix: ",dim(final)[1],"x",dim(final)[2],"(pvals x chr:bp)"))
  
  ## excluding empty cells and chromosome cells.
  cells<-m[m$value!=0,]
  cells.1<-cells[cells$value!=config$idx[config$type=="oddchr"],]
  cells.2<-cells.1[cells.1$value!=config$idx[config$type=="evenchr"],]
  
  message(paste("Melted matrix: ", dim(m)[1],"data points", "excluding empty cells: ", dim(cells)[1], "excluding chromosomes: ", dim(cells.2)[1]))

  message(paste(length(config$idx), "states in config file"))
  return.list<-list(matrix=final, melt=m, annot = cell.annot, chr.len = chr.matrix.len, annot.melt=annot.m, pos.df=pos.df)
}
  