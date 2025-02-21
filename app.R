##options(repos = BiocManager::repositories())

library(shiny)
library(ggplot2)
library(dplyr)
library(reshape2)
library(shinyjs)
##library("BiocManager")
##library("biomaRt")
##library(pryr)

### input files
## the list of studies to show
gwas.list<-"gwas_table.txt"
gwas.df<-read.table(gwas.list,header=T, sep="\t", stringsAsFactors = F)
gwas.list<-list()
for(i in 1:dim(gwas.df)[1]){
  gwas.list[[gwas.df$name[i]]]<-gwas.df$prefix[i]
}

## the configuration file
createConfig<-function(dataset="cad"){
  config.in<-paste0(dataset,".config.txt") 
  config<-read.table(config.in,sep="\t", header =T,stringsAsFactors = F, skip=10)  
  return (config)
}

## list of genes for search.
genes.in<-"gene_master_list_v3.txt"
genes.df<-read.table(genes.in, sep="\t", header = T)

## creates list of datasources for a specified resolution and chromosome
# bins - the resolution
createList<-function(bins = "500kbp", dataset="cad"){
  ## the melted matrix with config states as values - generates the heatmap
  melt.in<-paste0(bins,".",dataset,".melt.txt")
  m<-read.table(melt.in, header=T)
  
  ## Description for each of the config states for the legendf
  annot.in<-paste0(bins,".",dataset,".annot.txt")
  annot<-read.table(annot.in,sep="\t", header=T)
  
  ## chromosome information, to generate the chr axis breaks / text
  chrlen.in<-paste0(bins,".",dataset,".chrlen.txt")
  chr.matrix.len<-read.table(chrlen.in, header=T)
  
  ## the melted matrix with description of each cell as values - floating text over heatmap cells
  annotmelt.in<-paste0(bins,".",dataset,".annotmelt.txt")
  t<-read.table(annotmelt.in, header=T, sep="\t", stringsAsFactors = F)
  
  # maps of positions to heatmap cell
  positions.in<-paste0(bins,".",dataset,".positions.txt")
  pos.df<-read.table(positions.in, header = T)
  
  list.out<-list(m = m,
                 annot = annot,
                 chr.matrix.len = chr.matrix.len,
                 t = t,
                 pos.df = pos.df)
  
  return(list.out)
}

## creates list of datasources for a specified resolution and chromosome
# chr - the chromosome
# bins - the resolution
createChrList<-function(chr, bins = "500kbp", dataset="cad"){
  ## the melted matrix with config states as values - generates the heatmap
  melt.in<-paste0(bins,".",dataset,".chr",chr,".melt.txt")
  m<-read.table(melt.in, header=T)
  
  ## Description for each of the config states for the legendf
  annot.in<-paste0(bins,".",dataset,".chr",chr,".annot.txt")
  annot<-read.table(annot.in,sep="\t", header=T)
  
  ## chromosome information, to generate the chr axis breaks / text
  chrlen.in<-paste0(bins,".",dataset,".chr",chr,".chrlen.txt")
  chr.matrix.len<-read.table(chrlen.in, header=T)
  
  ## the melted matrix with description of each cell as values - floating text over heatmap cells
  annotmelt.in<-paste0(bins,".",dataset,".chr",chr,".annotmelt.txt")
  t<-read.table(annotmelt.in, header=T, sep="\t", stringsAsFactors = F)
  
  # maps of positions to heatmap cell
  positions.in<-paste0(bins,".",dataset,".chr",chr,".positions.txt")
  pos.df<-read.table(positions.in, header = T)
  
  list.out<-list(m = m,
                    annot = annot,
                    chr.matrix.len = chr.matrix.len,
                    t = t,
                    pos.df = pos.df)
  
  return(list.out)
}

## url to the locuszoom.js cgi script
weburl<-"https://euclid.well.ox.ac.uk/cgi-bin/manh/dyn.LZ_with_betas.cgi"

## instructions when cell information pane is not frozen
instructions<-"<b>Click on a cell to freeze</b></br>"
## instructions when cell information pane is frozen
instructions2<-"<b>Click anywhere to unfreeze</b></br>"

## the default text in the cell information freeze pane
defaultText<-paste(instructions,"<b>#Variants:</b></br><b>Best Variant:</b></br><b>Best p-value:</b></br><b>chr:start-end:</b></br><b>Prox. gene(s):</b>")

## create the UI object
ui <- fluidPage(theme="manhplot.css", title = "Interactive Manhattan++", useShinyjs(),
   ## titlePanel("Interactive Manhattan++"),
    sidebarPanel(style = "height: 100%; overflow-y: auto; position:fixed; width: 30%;",
      h2("Interactive Manhattan++"),selectInput("study", "Study:", gwas.list),
      uiOutput("title_info"),
      uiOutput("hover_info"), ## handles hover over heatmap functionality
      textInput("position", "Zoom to gene, SNP or position \n(eg. BAG3, rs72840788, chr10:121415685, chr10:121165685-121665685)"),
      actionButton("do", "Zoom"), ## search button
      ##verbatimTextOutput("searchval"), ## text shows the results of the search
      tags$br(),HTML("<b>Association plot Link:</b>"),
      wellPanel(uiOutput("Link")), ## link to the locuszoom plot (if available)
      # HTML("<h5 id = \"myP\">Instructions:</br>
      #   Zoom in: Select the region using the cursor and double click.</br>
      #   Zoom out: Double click without a cursor selection.</br>
      #   Cell information: Displayed when the cursor hovers over the cell. If cells contain p-values greater than 0.001 then cell information is not available (Cells are greyed out).Click on cell to freeze information pane.</br>
      #   Loci association plot: Link to an association plot displayed when query (position, gene or SNP) is less than 1 Mbp</br>
      #   </br>References:</br> 
      #   1. <a href=\"https://www.ncbi.nlm.nih.gov/pubmed/31775616\" target=\"_blank\">Grace, C. <i>et al</i>. Manhattan++: displaying genome-wide association summary statistics with multiple annotation layers. BMC Bioinformatics 20, 610 (2019).</a></br>
      #   2. <a href=\"https://www.ncbi.nlm.nih.gov/pubmed/20634204\"  target=\"_blank\">Pruim, R.J. <i>et al</i>. LocusZoom: regional visualization of genome-wide association scan results. Bioinformatics 26, 2336-7 (2010).</a></br></br>
      #   Contact: <a href=\"mailto:cgrace@well.ox.ac.uk\"  target=\"_blank\">cgrace@well.ox.ac.uk</a></br></br>
      #   Tested using Firefox, Safari & Chrome browsers in Windows 10, Mac OS 10.15.4.</h5>")
      ),
    mainPanel(tabsetPanel(type="tabs",
      tabPanel("Plot",verticalLayout(
        verbatimTextOutput("zoomval"), ## text to display region zoomed into on heatmap
        splitLayout(
          uiOutput("zoominui", inline = T), ## zoom in button
          uiOutput("zoomoutui", inline = T), ## zoom out button
          cellWidths = 105
        ), 
        ## CIRCLE ELEMENT
        uiOutput("circle"), ## the selection circle on the heatmap
        plotOutput("distPlot", ## the manhattan++ heatmap
          hover = hoverOpts(id ="plot_hover", delay = 150, delayType = "debounce"),
          dblclick = "plot1_dblclick",
          click="plot1_click",
          brush = brushOpts(
            id = "plot1_brush",
            resetOnNew = TRUE)
        ),
        tags$form(id = "myform",method="post", action=weburl , target="_blank", ## invisible form which contains POSTing functionality
          tags$input(type="hidden", name="chr", id="chrid" , value=""),
          tags$input(type="hidden", name="st", id = "stid", value=""),
          tags$input(type="hidden", name="en", id="enid", value=""),
          tags$input(type="hidden", name="snp", id = "snpid", value=""),
          tags$input(type="hidden", name="study", id="studyid", value=""),
          tags$input(type="hidden", name="prefix", id="prefixid", value="cad"),
          tags$input(type="hidden", name="PPA", id="PPAid", value=1),
          tags$input(type="delay", value="dummy", class="link-button")),
        tags$script("
          // Fix for form bug in shiny (https://community.rstudio.com/t/issuing-post-from-within-a-shiny-app/60413)
          var wait = function() {
            if (Shiny.setInputValue) {
              console.log('done!');
              $('input[type=\"delay\"]').attr('type', 'submit');
              document.getElementById('myform').style.visibility = 'hidden';
            } else {
              console.log('waiting');
              setTimeout(wait, 1);
            }
          }
        
          wait();
    
          // implement catching enter key on search box
          $(document).keyup(function(event) {
            if(event.key == 'Enter'){
              if(document.getElementById('position').matches(':focus')){
                document.getElementById('do').click();
              }
          }
        });
            
      ")
  )
  ), ## <img src=\"searchgene.png\" ></br></br>         Loci association plot: Link to an association plot displayed when query (position, gene or SNP) is less than 1 Mbp</br>
  tabPanel("Help", 
        HTML("<h4 id=\"top\"><b>Contents:</b></h4>
            <h4><a href=\"#search\">1. Search box</a></br>
             <a href=\"#myP\">2. Zoom functionality</a></br>
             <a href=\"#ref\">3. References</a></br>
             <a href=\"#ref\">4. Contact</a></br></br></h4>"),
        HTML("<h4 id=\"search\"><b>1. Search box </b></h4>
        The plot can be searched for genes, SNPs and coordinates using the search box and the corresponding zoom button.
        </br><b>Search gene</b></br>
        Searches for gene name and identifies cells which overlap with start and positions for gene. A link is generated for the locuszoom plot with 250kbp flanks around the gene coordinates.
        <img src=\"searchgene.png\" ></br></br>
        <b>Search SNP</b></br>
        Searches for rsid or chr:pos and then identifies and zooms to cells which overlap with this position. A link is generated for the locuszoom plot with 250kbp flanks around the SNP position
        <img src=\"searchsnp.png\" ></br></br>
        <b>Search coordinates</b></br>
        Will zoom into the coordinates specified in the search box. If the region between the coordinates is less than 1 Mbp a locuszoom plot link is generated for the length of the coordinates.
        <img src=\"searchcoord.png\"></br>"), 
        HTML("<h5 id = \"myP\"><h4><b>2. Zoom functionality</b></h4>
        Zoom in: Select the region using the cursor and double click.</br>
        <img src=\"zoomin.png\" ></br></br>
        Zoom out: Double click without a cursor selection.</br>
        <img src=\"zoomout.png\" ></br></br>
        Cell information: Displayed when the cursor hovers over the cell. Click on cell to freeze information pane.</br>
        <img src=\"circle.png\" ></br></br>
        "),
        HTML("</br><h4 id=\"ref\"><b>3. References</b></h4> 
        1. <a href=\"https://www.ncbi.nlm.nih.gov/pubmed/31775616\" target=\"_blank\">Grace, C. <i>et al</i>. Manhattan++: displaying genome-wide association summary statistics with multiple annotation layers. BMC Bioinformatics 20, 610 (2019).</a></br>
        2. <a href=\"https://www.ncbi.nlm.nih.gov/pubmed/20634204\"  target=\"_blank\">Pruim, R.J. <i>et al</i>. LocusZoom: regional visualization of genome-wide association scan results. Bioinformatics 26, 2336-7 (2010).</a></br>
        </br><h4 id=\"contact\"><b>4. Contact</b></h4>
        <a href=\"mailto:cgrace@well.ox.ac.uk\"  target=\"_blank\">cgrace@well.ox.ac.uk</a></br></br>
        Tested using Firefox, Safari & Chrome browsers in Windows 10, Mac OS 10.15.4.</h5>
        <a href=\"https://github.com/cgrace1978/manhplot/wiki\"  target=\"_blank\">Manhattan++ wiki</a></br></br>
        </br></br></br>"))))
)

## server object
server <- function(input, output, session) {
  ## output of the search:
  gene<-reactiveValues(name=NULL)
  ## text to display in the zoom box (heatmap zoom coordinates)
  zoomtext<-reactiveValues(data=NULL)
  ## the heatmap options
  plotopts<-reactiveValues(height = 2000, width=800, 
                           x=NULL, y = NULL, 
                           pval.units = 5, 
                           list=NULL,
                           config=NULL,
                           study=NULL,
                           zoom = FALSE)
  ## current study
  study<-reactiveValues(prefix=NULL,name=NULL, url=NULL, PPA=NULL)
  ## data to both to the LZ cgi script
  LZpos<-reactiveValues(chr=7, st=149887626, en=151190176, show=F, text="NOS3", snp=NULL)
  ## data to configure freezing cell information panel on mouse click
  panelFreeze<-reactiveValues(active=F, text="", instructions=instructions)
  ## data to configure display of circle on freeze
  circle.dat<-reactiveValues(active=F,left=300, top=275, initialSearchState=F)
  ## data to display the SNP marker
  snpmark.dat<-reactiveValues(active=F, x= 1, y = 1, pval=NULL)
  ## data for the zoom functionality
  zoom.dat<-reactiveValues(show = F, chr=NULL, st=NULL, en=NULL, out.max=FALSE, in.max=FALSE)
  
  ## display text under the heatmap (region displayed)
  output$zoomval<- renderText({zoomtext$data})
  
  ## display the results of the search
  ##output$searchval<- renderText({gene$name})
 
  ## display the zoom in button 
  output$zoominui <- renderUI(
    {
      if(zoom.dat$show == T && zoom.dat$in.max == F){
        actionButton(inputId = "zoomin", label = "Zoom in",icon = icon("search-plus"))
      }
      else{
        actionButton(inputId = "zoomin", label = "Zoom in",icon = icon("search-plus"), disabled=T)
      }
    }
  )
  
  ## display the zoom out button
  output$zoomoutui <- renderUI(
    {
      if(zoom.dat$show == T && zoom.dat$out.max == F){
        actionButton(inputId = "zoomout", label = "Zoom out",icon = icon("search-minus"))
      }
      else{
        actionButton(inputId = "zoomout", label = "Zoom out",icon = icon("search-minus"), disabled=T)
      }
    }
  )
  
  observeEvent(input$study, {
    message(input$study)
    ##print(gwas.df)
    
    study$prefix<-input$study
    study$name<-gwas.df[gwas.df$prefix == input$study,]$name
    study$url<-gwas.df[gwas.df$prefix == input$study,]$url
    study$PPA<-gwas.df[gwas.df$prefix == input$study,]$PPA
    
    plotopts$study<-input$study
    plotopts$list<-createList(bins = "3mbp", dataset=input$study) ## place datasource list here see createChrList function
    plotopts$config<-createConfig(dataset = input$study)
    
    plotopts$height=2000
    plotopts$x = NULL
    plotopts$y = NULL
    plotopts$zoom = FALSE
    
    zoomtext$data=NULL
    gene$name=NULL
    
    plotopts$pval.units=5
    
    LZpos$show=F
    LZpos$snp=NULL
    
    zoom.dat$show=F
    
    panelFreeze$active=F
    panelFreeze$instructions=instructions
    
    circle.dat$active=F
    circle.dat$initialSearchState=F
    snpmark.dat$active=F
  })
  
  
  ## draw the link to the LocusZoom server
  output$Link <- renderUI({
    if(LZpos$show==T && study$PPA != 1){
      snpval=""
      if(!is.null(LZpos$snp)){snpval = LZpos$snp}
      
      ## submit a POST request using the invsible form 'myform'
      tags$a(style="cursor:pointer", 
            onclick=
            paste0("document.getElementById('studyid').value='",study$name,"';
            document.getElementById('snpid').value='",snpval,"';
            document.getElementById('chrid').value='",LZpos$chr,"';
            document.getElementById('stid').value='",LZpos$st,"';
            document.getElementById('enid').value='",LZpos$en,"';
            document.getElementById('prefixid').value='",study$prefix,"';
            document.getElementById('PPAid').value='",study$PPA,"';
            document.getElementById('myform').submit();"), 
            HTML(LZpos$text,"<img src=\"icon.png\" width=15 height=15>")) ## ,"<img src=\"icon.png\" width=15 height=15>"
    }
    else{
      tags$p("Not available")
    }
  })
  
  widthFunct<-function(){
    return(plotopts$width)
  }
  
  heightFunct<-function(){
    return(plotopts$height)
  }
  
  ## Draw the manhattan++ heatmap
  output$distPlot <- renderPlot(width = reactive(widthFunct()), height = reactive(heightFunct()),expr = {

    pval.split=0.125
    max.pval<-20
    pval.seq<-seq(from=plotopts$pval.units,to=max.pval,by=plotopts$pval.units)
    
    pvals<-seq(from=0, to=max.pval, by=pval.split)
    pvals.cells.index<-data.frame(id=1:length(pvals),LP=pvals,UP=c(pvals[2:length(pvals)],max.pval))
    
    p.val.cell<-function(p.val,pval.chunks){
      gws<--1*log10(p.val) ## convert the p-value to -log10
      ##print(gws)
      idx<-log10.cell(gws) ## call the log10.index function.
      return(idx)
    }
    
    log10.cell<-function(log, pval.chunks){
      ##print(pvals.cells.index)
      if(log >= max(pvals.cells.index$UP)){
        idx<-max(pvals.cells.index$id)
      }
      else{
        idx<-pvals.cells.index$id[log >= pvals.cells.index$LP & log < pvals.cells.index$UP]
      }
      return (idx)
    }
    
    ### function returns y index on the heatmap that maps to a specific -log10(p-value)
    ## log - log10 pval to convert to index
    log10.index<-function(log){
      idx<-(log / pval.split) + 0.5
      return(idx)
    }
    
    x.labels<-c(0,pval.seq)
    x.breaks<-c(0.5,log10.index(pval.seq))
    
    # print(x.labels)
    # print(x.breaks)
    # print(plotopts$x)
    
    y.labels<-plotopts$list$chr.matrix.len$chr
    y.breaks<-plotopts$list$chr.matrix.len$mid
    
    chromosome.label<-"chromosome"
    
    if(length(y.labels) == 1){
      chromosome.label<-paste0("chromosome ", y.labels[1])
      y.labels<-NULL
    }
    
    grad.breaks<-seq(from=0,to=max(plotopts$config$idx))
    grad.labels<-c("",as.character(plotopts$list$annot$name))
    
    p <-ggplot(data=plotopts$list$m,aes(Var1, Var2, fill = value)) + geom_tile() + #colour="black", stroke=1E-6) +
      scale_x_continuous("-log(pval)",
                         breaks=x.breaks, 
                         labels=x.labels,
                         trans="reverse", 
                         position = "top",
                         expand=c(0,0)) +
      scale_y_continuous(chromosome.label,
                         breaks=y.breaks, 
                         labels=y.labels,
                         trans="reverse", 
                         position = "right",
                         expand=c(0,0)) +
      scale_fill_gradientn(colours = c("white",plotopts$config$col), 
                        guide="legend", 
                        breaks=grad.breaks, 
                        labels=grad.labels,
                        name = "Key\n#Variants (Consequence)", 
                        limits = c(0,max(plotopts$config$idx))) +
      coord_cartesian(xlim = plotopts$x, 
                        ylim = plotopts$y, 
                        expand = TRUE)  +
      theme(axis.line.x = element_line(color="black", size = 0.5),
                        axis.line.y = element_line(color="black", size = 0.5),
                        axis.text=element_text(size=10),
                        axis.title=element_text(size=12,face="bold"),
                        legend.position="left",
                        legend.key.size=unit(1.0,"line"),
                        legend.title=element_text(size=12, face="bold"),
                        legend.text=element_text(size=12),
                        panel.grid.major =element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank()) + 
      geom_vline(xintercept=p.val.cell(5E-8)+0.5, linetype="dashed", colour="red")
      
    if(snpmark.dat$active == T){
     ## print("mark")
      pval.x<-p.val.cell(snpmark.dat$pval)
     ## print(snpmark.dat$pval)
     ## print(pval.x)
      p<-p+annotate("point", x = pval.x, y = snpmark.dat$y, 
                    colour="red", size=10, shape = 7)  
    }
    
    ##p<-p + expand_limits(x = 200)
    
    return(p)
  })

  ## Search for HCNC gene symbol in local data.frame
  searchGene<-function(geneName){
    gene$name<-geneName ## get the chromosome
    
    filt.df<-genes.df[toupper(as.character(genes.df$HGNC)) == toupper(gene$name),]
    
    if(dim(filt.df)[1] > 0){
      gene$name<-paste(filt.df$HGNC[1], filt.df$Chromosome, 
                       format(filt.df$st, big.mark = ",",scientific = F), 
                       format(filt.df$en, big.mark = ",", scientific = F))
    }
    else{
      return(FALSE)
    }
    
    if(filt.df$Chromosome < 1 || filt.df$Chromosome > 22){return()}
    
    chrList<-createChrList(filt.df$Chromosome, bins = "100kbp", dataset = plotopts$study)
    
    ## subset the dataframe
    chr.sub<-chrList$pos.df[chrList$pos.df$chr == filt.df$Chromosome,]
    bp.sub<-chr.sub[chr.sub$en > filt.df$st & chr.sub$st < filt.df$en,]
    
    st.txt<-chr.sub$st[chr.sub$idx == min(chr.sub$idx)]
    en.txt<-chr.sub$en[chr.sub$idx == max(chr.sub$idx)]
    
    height<-500
    
    if(dim(bp.sub)[1] == 0){ ## if the position does not intersect, plot the whole chromosome
      st<-min(chr.sub$idx)
      en<-max(chr.sub$idx)
    }
    else{
      st<-min(bp.sub$idx)
      en<-max(bp.sub$idx)
      
      st.txt<-bp.sub$st[bp.sub$idx == min(bp.sub$idx)]
      en.txt<-bp.sub$en[bp.sub$idx == max(bp.sub$idx)]
      
      height<-275
    }
    
    zoomtext$data<-paste(paste0(filt.df$HGNC[1], ": chr", filt.df$Chromosome, ":" ,
                               format(filt.df$st, big.mark = ",",scientific = F), "-", 
                               format(filt.df$en, big.mark = ",", scientific = F)),"\nViewing: chr",filt.df$Chromosome,":",
                         format(st.txt, big.mark = ",",scientific = F),"-",
                         format(en.txt, big.mark = ",",scientific = F), sep="")
    
    LZpos$chr=filt.df$Chromosome
    LZpos$st=(filt.df$st-2.5E5)
    LZpos$en=(filt.df$en+2.5E5)
    LZpos$text=paste("LocusZoom plot for ",filt.df$HGNC[1]," coordinates: chr",LZpos$chr,":",LZpos$st,"-",LZpos$en, sep="")
    LZpos$show=T
    LZpos$snp=NULL
    
    plotopts$x<-c(165,0) ## all the pvalues
    plotopts$y<-c((en+0.5),(st-0.5)) ## start and end of the region
    plotopts$zoom = TRUE
    
    plotopts$list = NULL
    gc()
    plotopts$list = chrList
    
    dist<-en-st
    
    plotopts$height<-height
    
    plotopts$pval.units=5
    
    snpmark.dat$active = F
    
    zoom.dat$show=T
    zoom.dat$chr = filt.df$Chromosome
    zoom.dat$st = st.txt
    zoom.dat$en = en.txt
    
    return(TRUE)
  }
  
  ## Search for specific coordinates in the genome (chr##:##-##)
  searchCoord<-function(position){
    chr<-as.numeric(gsub(x=strsplit(position, ":")[[1]][1], pattern = "chr", replacement = ""))
    stpos<-as.numeric(gsub(x=strsplit(strsplit(position, ":")[[1]][2],"-")[[1]][1], pattern = ",", replacement = ""))
    enpos<-as.numeric(gsub(x=strsplit(strsplit(position, ":")[[1]][2],"-")[[1]][2], pattern = ",", replacement = ""))
      
    if(is.na(chr) || chr < 1 || chr > 22){return(FALSE)}
    
    ## subset the dataframe
    chrList<-createChrList(chr, bins = "100kbp", dataset = plotopts$study)
    
    chr.sub<-chrList$pos.df[chrList$pos.df$chr == chr,]
    
    st.txt<-chr.sub$st[chr.sub$idx == min(chr.sub$idx)]
    en.txt<-chr.sub$en[chr.sub$idx == max(chr.sub$idx)]
    
    height<-500
    
    if(is.na(enpos) && !is.na(stpos)){
      print("HERE 3")
      bp.sub<-chr.sub[chr.sub$en > (as.numeric(stpos)-5E5) & chr.sub$st < (as.numeric(stpos)+5E5),]
      if(dim(bp.sub)[1] == 0){ ## if the position does not intersect, plot the whole chromosome
        st<-min(chr.sub$idx)
        en<-max(chr.sub$idx)
      }
      else{
        st<-min(bp.sub$idx)
        en<-max(bp.sub$idx)
        
        st.txt<-bp.sub$st[bp.sub$idx == min(bp.sub$idx)]
        en.txt<-bp.sub$en[bp.sub$idx == max(bp.sub$idx)]
        
        height<-275
      }
    } else if(is.na(stpos)){
      st<-min(chr.sub$idx)
      en<-max(chr.sub$idx)
    }
    else{
      bp.sub<-chr.sub[chr.sub$en > as.numeric(stpos) & chr.sub$st < as.numeric(enpos),]
      if(dim(bp.sub)[1] == 0){ ## if the position does not intersect, plot the whole chromosome
        st<-min(chr.sub$idx)
        en<-max(chr.sub$idx)
      }
      else{
        st<-min(bp.sub$idx)
        en<-max(bp.sub$idx)
        
        st.txt<-bp.sub$st[bp.sub$idx == min(bp.sub$idx)]
        en.txt<-bp.sub$en[bp.sub$idx == max(bp.sub$idx)]
        
        height<-275
      }
    }
    
    max.text<-""
    if(st == min(chr.sub$idx) && en == max(chr.sub$idx)){
      zoom.dat$out.max=T
      max.text<-" (Maximum zoom out reached)"
    }
    else{
      zoom.dat$out.max=F
    }
    
    gene$name = position 
    zoomtext$data<-paste(position,"\nViewing: chr",chr,":", 
                         format(st.txt, big.mark = ",", scientific = F),"-",
                         format(en.txt, big.mark = ",", scientific = F), 
                         max.text,sep="")
    
    prefix<-""
    
    if(is.na(enpos) && !is.na(stpos)){ ## if only one coordinate
      tmp<-stpos
      stpos<-(as.numeric(stpos)-5E5)
      enpos<-(as.numeric(tmp)+5E5)
      ##gene$name=paste0(position,": ", "variant not found")
      zoomtext$data<-paste(position,": ", "variant not found - displaying coordinates","\nViewing: chr",chr,":", 
                           format(st.txt, big.mark = ",", scientific = F),"-",
                           format(en.txt, big.mark = ",", scientific = F), 
                           max.text,sep="")
      ##prefix<-"variant not found:\n"
    }

    region<-as.numeric(enpos)-as.numeric(stpos)

    print(paste(region,enpos,stpos))
    print(gene$name)
    
    if(region <= 1E6 && !is.na(stpos) && !is.na(enpos)){
      LZpos$chr=chr
      LZpos$st=format(stpos,scientific = F)
      LZpos$en=format(enpos,scientific = F)
      LZpos$text=paste(prefix,"LocusZoom plot for coordinates: chr",LZpos$chr,":",LZpos$st,"-",LZpos$en, sep="")
      LZpos$show=T
      LZpos$snp=NULL
    }
    else{
      print("HERE2")
      LZpos$show=F
      LZpos$snp=NULL
    }
    
    plotopts$x<-c(165,0) ## all the pvalues
    plotopts$y<-c((en+0.5),(st-0.5)) ## start and end of the region
    plotopts$zoom = TRUE
    
    plotopts$list = NULL
    gc()
    plotopts$list = chrList
    
    dist<-(en-st)
    
    plotopts$height<-height
    
    plotopts$pval.units=5
    
    snpmark.dat$active = F
    
    zoom.dat$show=T
    zoom.dat$chr = chr
    zoom.dat$st = st.txt
    zoom.dat$en = en.txt
    
    return(TRUE)
  }
  
  ## function to search for SNPs in file using grep
  searchSNP<-function(snpName){
    if(length(grep("rs", snpName))!=0){
      
      grep.out<-system(paste0("grep -w ",snpName," ",study$prefix,".add.snpinfo.txt"), intern=T)
      
      out<-strsplit(grep.out, "\t")

      if(length(out) > 0){
        LZpos$chr=out[[1]][2]
        LZpos$st=(as.numeric(out[[1]][3])-2.5E5)
        LZpos$en=(as.numeric(out[[1]][3])+2.5E5)
        LZpos$text=paste("LocusZoom plot for ",out[[1]][1]," coordinates: chr",LZpos$chr,":",LZpos$st,"-",LZpos$en, sep="")
        LZpos$show=T
        LZpos$snp = out[[1]][1]
        
        chrList<-createChrList(LZpos$chr, bins = "100kbp", dataset = plotopts$study)
        ## subset the dataframe
        chr.sub<-chrList$pos.df[chrList$pos.df$chr == LZpos$chr,]
        bp.sub<-chr.sub[chr.sub$en >= as.numeric(out[[1]][3]) & chr.sub$st < as.numeric(out[[1]][3]),]
 
        st<-min(bp.sub$idx)
        en<-max(bp.sub$idx)
        
        st.txt<-bp.sub$st[bp.sub$idx == min(bp.sub$idx)]
        en.txt<-bp.sub$en[bp.sub$idx == max(bp.sub$idx)]
        
        plotopts$x<-c(165,0) ## all the pvalues
        plotopts$y<-c((en+0.5),(st-0.5)) ## start and end of the region
        plotopts$zoom = TRUE
        
        plotopts$list = NULL
        gc()
        plotopts$list = chrList
        
        plotopts$height<-275
    
        gene$name<-paste(out[[1]][1]," chr",out[[1]][2],
                         ":" ,format(as.numeric(out[[1]][3]), big.mark = ",", scientific = F), 
                         sep="") ## get the chromosome
        
        zoomtext$data<-paste(paste(out[[1]][1]," chr",out[[1]][2],
                                   ":" ,format(as.numeric(out[[1]][3]), big.mark = ",", scientific = F), 
                                   sep=""),"\nViewing: chr",LZpos$chr,":", 
                             format(st.txt, big.mark = ",", scientific = F),"-",
                             format(en.txt, big.mark = ",", scientific = F), sep="")
            
        plotopts$pval.units=5
        
        snpmark.dat$active = T
        snpmark.dat$x = 100
        snpmark.dat$y = st##((st+en) / 2)
        snpmark.dat$pval = as.numeric(out[[1]][4])
        
        zoom.dat$show=T
        zoom.dat$chr = LZpos$chr
        zoom.dat$st = st.txt
        zoom.dat$en = en.txt
        
        return(TRUE)
      }
      else{
        print("HERE 5")
        gene$name<-paste("Cannot find ",snpName, sep="") ## get the chromosome
        
        return(FALSE)
      }
    }
    else{
      chr<-as.numeric(gsub(x=strsplit(snpName, ":")[[1]][1], 
                           pattern = "chr", replacement = ""))
      pos<-as.numeric(gsub(strsplit(strsplit(snpName,":")[[1]][2],"_")[[1]][1],
                           pattern=",",replacement=""))
      
      if(!is.na(chr) && !is.na(pos)){
        out.str<-character(0)
        out.grep<-system(paste0("grep -w ",pos," ",study$prefix,".add.snpinfo.txt"), intern = T)
        print(out.grep)
        
        if(length(out.grep) > 0){
          out<-strsplit(out.grep,"\t")
          
          for(j in 1:length(out)){
            if(out[[j]][2] == chr){
              out.str<-out[[j]]
              break
              }
          }
        }
      
        if(length(out.str) == 0){
          print("HERE")
          gene$name = paste(snpName,"not found")
          return(FALSE)
          }
          
        LZpos$chr=chr
        LZpos$st=(pos-2.5E5)
        LZpos$en=(pos+2.5E5)
        LZpos$text=paste("LocusZoom plot for ",snpName," coordinates: chr",LZpos$chr,":",LZpos$st,"-",LZpos$en, sep="")
        LZpos$show=T
        LZpos$snp=gsub(snpName, pattern="chr", replacement="")
        
        gene$name = snpName
        
        ## subset the dataframe
        chrList<-createChrList(chr, bins = "100kbp", dataset = plotopts$study)
        chr.sub<-chrList$pos.df[chrList$pos.df$chr == chr,]
        bp.sub<-chr.sub[chr.sub$en >= pos & chr.sub$st < pos,]
      
        st<-min(bp.sub$idx)
        en<-max(bp.sub$idx)
        
        st.txt<-bp.sub$st[bp.sub$idx == min(bp.sub$idx)]
        en.txt<-bp.sub$en[bp.sub$idx == max(bp.sub$idx)]
    
        plotopts$x<-c(165,0) ## all the pvalues
        plotopts$y<-c((en+0.5),(st-0.5)) ## start and end of the region
        plotopts$zoom = TRUE
        
        plotopts$list = NULL
        gc()
        plotopts$list = chrList
        
        plotopts$height<-275
        
        zoomtext$data<-paste(snpName,"\nViewing: chr",LZpos$chr,":", 
                             format(st.txt, big.mark = ",", scientific = F),"-",
                             format(en.txt, big.mark = ",", scientific = F), sep="")
        
        plotopts$pval.units=5
    
        snpmark.dat$active = T
        snpmark.dat$x = 100
        snpmark.dat$y = st
        snpmark.dat$pval = as.numeric(out.str[4])
        
        return(TRUE)  
      }
      else{return(FALSE)}
    }
    plotopts$x <- NULL
    plotopts$y <- NULL
    plotopts$height<-2000
    
    zoomtext$data<-NULL
    gene$name<-NULL
    
    plotopts$pval.units=5  
    return(TRUE)
  }
  
  output$title_info<-renderUI(
    {HTML(paste("<h1><a href=",
                study$url," target=\"_blank\">",
                study$name,"</a></h1>",sep=""))}
  )
  
  ## handle hover events on the heatmap: Handles cell information
  output$hover_info<-renderUI({
      if(panelFreeze$active == FALSE){
        if(!is.null(input$plot_hover)){
          hover=input$plot_hover
            
          if(round(hover$x,0) == 0 || round(hover$y,0) == 0 ||
             round(hover$x,0) > 161 || round(hover$y,0) > 931){
              if(round(hover$x,0) > 161){
                  return(wellPanel(p(HTML(paste0(defaultText)))))
              }
            
            wellPanel(
              p(HTML(paste0(defaultText))))
          }
            
          panelText<-as.character(plotopts$list$t$value[plotopts$list$t$Var1==round(hover$x,0)&plotopts$list$t$Var2==round(hover$y,0)])
          
          if(nchar(panelText) > 0){
            out<-strsplit(panelText, split="</br>")
            
            nvar<-out[[1]][1]
            snp<-out[[1]][2]
            pval<-out[[1]][3]
            chrbp<-out[[1]][4]
            gene<-out[[1]][5]
            
            mafconseqstr<-""
            
            if(length(out[[1]]) == 9){
              mafstr<-out[[1]][7]
              conseqstr<-out[[1]][9]
              
              mafconseqstr<-paste0(mafstr,"</br>",conseqstr)
            }
            else if(length(out[[1]]) == 8){
              mafconseqstr<-out[[1]][8]
            }
            
            out2<-strsplit(chrbp, split=":")
            chr<-as.numeric(out2[[1]][1])
            bp<-as.numeric(out2[[1]][2])
                    
            ## submit a POST request using the invsible form 'myform'
            chrbplink<-tags$a(style="cursor:pointer",
                   onclick=
                     paste0("document.getElementById('studyid').value='",study$name,"';
                      document.getElementById('snpid').value='","","';
                      document.getElementById('chrid').value='",chr,"';
                      document.getElementById('stid').value='",(bp - 1E5),"';
                      document.getElementById('enid').value='",(bp + 1E5),"';
                      document.getElementById('prefixid').value='",study$prefix,"';
                      document.getElementById('PPAid').value='",study$PPA,"';
                      document.getElementById('myform').submit();"),
                   paste0("chr",chr,":",(bp - 1E5),"-",(bp + 1E5)))
            
            snplink<-tags$a(style="cursor:pointer",
                              onclick=
                                paste0("document.getElementById('studyid').value='",study$name,"';
                      document.getElementById('snpid').value='",snp,"';
                      document.getElementById('chrid').value='",chr,"';
                      document.getElementById('stid').value='",(bp - 1E5),"';
                      document.getElementById('enid').value='",(bp + 1E5),"';
                      document.getElementById('prefixid').value='",study$prefix,"';
                      document.getElementById('PPAid').value='",study$PPA,"';
                      document.getElementById('myform').submit();"),
                              snp)
            
            
            ##print(study$PPA)
            
            if(study$PPA == 1){ 
              panelText<-paste0("<b>#Variants:</b> ",nvar,"</br>",
                                "<b>Best Variant:</b> ",snp,"</br>",
                                "<b>Best p-value: </b>",pval,"</br>",
                                "<b>chr:start-end:</b> ",chrbp, "</br>",
                                "<b>Prox. gene(s):</b> ",gene,
                                "</br>",mafconseqstr)
            }
            else{
              panelText<-paste0("<b>#Variants:</b> ",nvar,"</br>",
                                "<b>Best Variant:</b> ",snplink,"&nbsp;<img src=\"icon.png\" width=15 height=15></br>",
                                "<b>Best p-value: </b>",pval,"</br>",
                                "<b>chr:start-end:</b> ",chrbplink, "&nbsp;<img src=\"icon.png\" width=15 height=15></br>",
                                "<b>Prox. gene(s):</b> ",gene,
                                "</br>",mafconseqstr)
            }
         }
          
          panelFreeze$text<-panelText
          
          if((!identical(panelFreeze$text, character(0))) && nchar(panelFreeze$text) > 0){ 
            wellPanel(
              p(HTML(paste0(panelFreeze$instructions,panelFreeze$text))))
          }
          else{
            wellPanel(
              p(HTML(paste0(defaultText))))
          }
        }
        else{
          wellPanel(
            p(HTML(paste0(defaultText))))
        }
      }
    else{
      wellPanel(
        p(HTML(paste0(panelFreeze$instructions,panelFreeze$text))))
    }
  })
  
  ## draw the selection circle on the heatmap
  output$circle <- renderUI({
    if(circle.dat$active==T){
      ## CIRCLE POSITION
      ## draws a svg of a circle on the heatmap using the coordinates from the click event
      # HTML(paste0("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"50px\" height=\"50px\" id=\"svgid\"
      #          style=\"position:fixed;z-index:10;left:",circle.dat$left,"px;top:",circle.dat$top,"px;\">
      #            <circle
      #            cx=\"15px\" cy=\"15px\" r=\"15\" fill=\"red\" 
      #            fill-opacity=\"0.4\" onclick=\"Shiny.setInputValue('plot1_click','foo');\"/>
      #            </svg>"))
    }
    else{
      return(cat(""))
    }
    
    })
  
  ## observe the zoom out event
  observeEvent(input$zoomout, {
    ## subset the dataframe
    chr.sub<-plotopts$list$pos.df[plotopts$list$pos.df$chr == zoom.dat$chr,]
    bp.sub<-chr.sub[chr.sub$en > (zoom.dat$st - 1E5) & chr.sub$st < (zoom.dat$en + 1E5),]
    
    st<-min(bp.sub$idx)
    en<-max(bp.sub$idx)
      
    plotopts$x<-c(165,0) ## all the pvalues
    plotopts$y<-c((en+0.5),(st-0.5)) ## start and end of the region
    plotopts$zoom = TRUE
    
    max.text<-""
    
    if(st == min(chr.sub$idx) && en == max(chr.sub$idx)){
      zoom.dat$out.max=T
      max.text<-" (Maximum zoom out reached)"
    }
    else if(st == min(chr.sub$idx)){
      max.text<-" (Lower bounds reached)"
      zoom.dat$en = zoom.dat$en+1E5
    }
    else if(en == max(chr.sub$idx)){
      max.text<-" (Upper bounds reached)"
      zoom.dat$st = zoom.dat$st-1E5
    }
    else{
      zoom.dat$st = zoom.dat$st-1E5
      zoom.dat$en = zoom.dat$en+1E5
    }
    
    zoom.dat$in.max=F
    
    circle.dat$active=F
    panelFreeze$active=F
    circle.dat$initialSearchState=F
    
    zoomtext$data<-paste(paste0("Viewing: chr",zoom.dat$chr,":",
                                format(zoom.dat$st, big.mark = ",",scientific = F),"-",
                                format(zoom.dat$en, big.mark = ",",scientific = F), 
                                max.text, 
                                sep=""))
  })
  
  ## observe a zoom in event
  observeEvent(input$zoomin, {
    ## subset the dataframe
    chr.sub<-plotopts$list$pos.df[plotopts$list$pos.df$chr == zoom.dat$chr,]
    bp.sub<-chr.sub[chr.sub$en > (zoom.dat$st + 1E5) & chr.sub$st < (zoom.dat$en - 1E5),]
    
    st<-min(bp.sub$idx)
    en<-max(bp.sub$idx)
    
    plotopts$x<-c(165,0) ## all the pvalues
    plotopts$y<-c((en+0.5),(st-0.5)) ## start and end of the region
    plotopts$zoom = TRUE
    
    zoom.dat$st = zoom.dat$st+1E5
    zoom.dat$en = zoom.dat$en-1E5
    zoom.dat$out.max=F
    
    min.text<-""
    ## if not possible to zoom in any further
    if((zoom.dat$st+1E5) >= (zoom.dat$en - 1E5)){
      zoom.dat$in.max=T
      min.text<-" (Minimum zoom in reached)"
    }
    
    circle.dat$active=F
    panelFreeze$active=F
    circle.dat$initialSearchState=F
    
    zoomtext$data<-paste(paste0("Viewing: chr",zoom.dat$chr,":",
                         format(zoom.dat$st, big.mark = ",",scientific = F),"-",
                         format(zoom.dat$en, big.mark = ",",scientific = F), 
                         min.text,
                         sep=""))
    
    })
  
  ## observe the action button: activates search functionality
  observeEvent(input$do, {
    circle.dat$active=F
    panelFreeze$active=F
    
    result<-searchGene(input$position)
    
    if(result==TRUE){
      zoom.dat$show=T 
      zoom.dat$in.max=T
      zoom.dat$out.max=F
      circle.dat$initialSearchState=T
      zoomtext$data=paste0(zoomtext$data, " (minimum zoom in reached)")
      return()
      }
    
    result<-searchSNP(input$position)
    
    if(result == TRUE){
      zoom.dat$show=T 
      zoom.dat$in.max=T
      zoom.dat$out.max=F
      circle.dat$initialSearchState=T
      zoomtext$data=paste0(zoomtext$data, " (minimum zoom in reached)")
      return()
      }
    
    result<-searchCoord(input$position)
    
    if(result == TRUE){
      zoom.dat$show=T 
      zoom.dat$in.max=F
      #zoom.dat$out.max=F
      
      min.text<-""
      ## if not possible to zoom in any further
      if((zoom.dat$st+1E5) >= (zoom.dat$en - 1E5)){
        zoom.dat$in.max=T
        min.text<-" (Minimum zoom in reached)"
      }
      
      zoomtext$data=paste0(zoomtext$data, min.text)
      
      circle.dat$initialSearchState=T
      
      return()
    }
    
    alert(paste0("query not found: ", input$position))
  })
  
  ## handle double click events on the heatmap, activates zoom functionality
  observeEvent(input$plot1_dblclick, {
    brush <- input$plot1_brush
    
    if (!is.null(brush)) {
      plotopts$x <- c(brush$xmax, brush$xmin)
      plotopts$y <- c(brush$ymax, brush$ymin)
      
      plotopts$height<-500
      
      ymin.val<-ceiling(plotopts$y[1])
      if(plotopts$y[1]%%1 < 0.5){ymin.val<-floor(plotopts$y[1])}
      
      ymax.val<-ceiling(plotopts$y[2])
      if(plotopts$y[2]%%1 < 0.5){ymax.val<-floor(plotopts$y[2])}
      
      if(ymin.val > max(plotopts$list$pos.df$idx)){
        ymin.val<-max(plotopts$listpos.df$idx)
      }
      
      if(ymax.val > max(plotopts$list$pos.df$idx)){
        ymax.val<-max(plotopts$list$pos.df$idx)
      }
      
      st.chr.txt<-plotopts$list$pos.df$chr[plotopts$list$pos.df$idx==ymin.val]
      st.txt<-plotopts$list$pos.df$st[plotopts$list$pos.df$idx==ymin.val]
      en.chr.txt<-plotopts$list$pos.df$chr[plotopts$list$pos.df$idx==ymax.val]
      en.txt<-plotopts$list$pos.df$en[plotopts$list$pos.df$idx==ymax.val]
      
      if(en.chr.txt == st.chr.txt){
        zoomtext$data<-paste("Viewing: chr",st.chr.txt,":",
                             format(st.txt, big.mark = ",",scientific = F),"-",
                             format(en.txt, big.mark = ",",scientific = F), sep="")
      }
      else{
        zoomtext$data<-paste("Viewing: chr",st.chr.txt,":",
                             format(st.txt, big.mark = ",",scientific = F),"-",
                             "chr",en.chr.txt,":",format(en.txt, big.mark = ",",scientific = F), sep="")
      }
      
      height<-500
      
      if(st.chr.txt==en.chr.txt){
        dist<-en.txt - st.txt
        
        if(dist <= 9E6){
          height<-275
        }  
      }
      
      plotopts$height<-height
      
      x.sel<-brush$xmax - brush$xmin
      
      if(x.sel < 10){
        plotopts$pval.units = 0.125
      }
      else if(x.sel< 20){
        plotopts$pval.units = 0.5
      }
      else if(x.sel < 50){
        plotopts$pval.units = 1
      }
      else if(x.sel < 100){
        plotopts$pval.units = 2.5
      }
      else if(x.sel < 2){
        plotopts$pval.units=0.0625
      }
      
      plotopts$zoom = TRUE
      zoom.dat$show=F
      
    } else {
      plotopts$height=2000
      plotopts$x = NULL
      plotopts$y = NULL
      plotopts$zoom = FALSE
      
      plotopts$list = NULL
      gc()
      plotopts$list = createList(bins = "3mbp", dataset = plotopts$study)
      
      zoomtext$data=NULL
      gene$name=NULL
      
      plotopts$pval.units=5
      
      LZpos$show=F
      LZpos$snp=NULL
      
      zoom.dat$show=F
    }
    
    panelFreeze$active=F
    panelFreeze$instructions=instructions
    
    circle.dat$active=F
    circle.dat$initialSearchState=F
    snpmark.dat$active=F
  })
  
  ## handle click events on the heatmap: activates freeze functionality
  observeEvent(input$plot1_click, {
    if(panelFreeze$active == FALSE && 
       nchar(panelFreeze$text) > 0){ ##only active when cell information is in panel
      panelFreeze$active = TRUE
      panelFreeze$instructions=instructions2

      left_pct <- (input$plot1_click$x - input$plot1_click$domain$left) / (input$plot1_click$domain$right - input$plot1_click$domain$left)
      top_pct <- (input$plot1_click$domain$top - input$plot1_click$y) / (input$plot1_click$domain$top - input$plot1_click$domain$bottom)
      
      # calculate distance from left and bottom side of the picture in pixels
      left_px <- input$plot1_click$range$left + (left_pct * (input$plot1_click$range$right - input$plot1_click$range$left))
      top_px <- input$plot1_click$range$top + (top_pct * (input$plot1_click$range$bottom - input$plot1_click$range$top))
      
      factor_top<-1
      factor_left<-1
      message(paste0("top: ",top_px))
      message(paste0("left:",left_px))
      
      circle.dat$active=T
      
      nudge_left<-175
      nudge_top<-15
      
      ## CIRCLE COORDS
      ## have had to nudge these coordinates
      if(is.null(plotopts$x)){
        circle.dat$left=left_px#nudge_left#left_pct*factor_left#left_px # - 2.5
        circle.dat$top=top_px#nudge_top#top_pct*factor_top # - nudge# + 17.5
      }
      else{
        circle.dat$left=left_px#nudge_left#left_pct*factor_left
        
        if(circle.dat$initialSearchState == T){
          circle.dat$top=top_px#nudge_top#top_px - nudge# + 85.0
        }
        else{
          circle.dat$top=top_px#nudge_top#top_px - nudge# + 70.0
        }
      }
    }
    else if(panelFreeze$active == TRUE){
      panelFreeze$active = FALSE
      panelFreeze$instructions=instructions
      
      circle.dat$active=F
    }
  
  })
  
  observe({
    query <- parseQueryString(session$clientData$url_search)
    if (!is.null(query[['dataset']])) {
      ##alert(query[['dataset']])
      
      dataset<-query[['dataset']]
      
     ## print(gwas.df)
      
      study$prefix<-dataset
      study$name<-gwas.df[gwas.df$prefix == dataset,]$name
      study$url<-gwas.df[gwas.df$prefix == dataset,]$url
      study$PPA<-gwas.df[gwas.df$prefix == dataset,]$PPA
      
      plotopts$study<-dataset
      plotopts$list<-createList(bins = "3mbp", dataset=dataset) ## place datasource list here see createChrList function
      plotopts$config<-createConfig(dataset = dataset)
      
      plotopts$height=2000
      plotopts$x = NULL
      plotopts$y = NULL
      plotopts$zoom = FALSE
      
      zoomtext$data=NULL
      gene$name=NULL
      
      plotopts$pval.units=5
      
      LZpos$show=F
      LZpos$snp=NULL
      
      zoom.dat$show=F
      
      panelFreeze$active=F
      panelFreeze$instructions=instructions
      
      circle.dat$active=F
      circle.dat$initialSearchState=F
      snpmark.dat$active=F
      
      updateSelectInput(session = session, inputId = "study", selected = dataset)
    }
  })
}

## run the shiny app
shinyApp(ui = ui,server = server)