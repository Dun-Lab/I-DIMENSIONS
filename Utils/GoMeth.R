library(knitr)
library(limma)
library(minfi)
library(IlluminaHumanMethylationEPICv2manifest)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(RColorBrewer)
library(missMethyl)
library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)
library(org.Hs.eg.db)
library(tidyverse) 
library(ggrepel)
library(sva)
library(reshape2)
library(ggsignif)
library(ggpubr)


.getFlatAnnotation <- function(array.type=c("450K","EPICv2"), anno=NULL)
  # flatten 450k or EPIC array annotation
  # Jovana Maksimovic
  # 18 September 2018
  # Updated 18 September 2018
  # Modified version of Belida Phipson's .flattenAnn code
{
  library(minfi)
  library(org.Hs.eg.db)
  library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
  
  if(is.null(anno)){
    if(array.type=="450K"){
      anno <- minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19::IlluminaHumanMethylation450kanno.ilmn12.hg19)
    } 
    else if (array.type=="EPICv2"){
      anno <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
    }
    else {
      anno <- minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19::IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
    }
  }
  
  # get rid of the non-CpG sites
  ann.keep<-anno[grepl("^cg",anno$Name),]
  
  # get rid of CpGs that are not annotated
  missing<-ann.keep$UCSC_RefGene_Name==""
  ann.keep<-ann.keep[!missing,]
  
  # get individual gene names for each CpG
  geneslist<-strsplit(ann.keep$UCSC_RefGene_Name,split=";")
  names(geneslist)<-rownames(ann.keep)
  
  grouplist<-strsplit(ann.keep$UCSC_RefGene_Group,split=";")
  names(grouplist)<-rownames(ann.keep)
  
  flat<-data.frame(symbol=unlist(geneslist),group=unlist(grouplist))
  flat$symbol<-as.character(flat$symbol)
  flat$group <- as.character(flat$group)
  
  flat$cpg<- substr(rownames(flat),1,10)
  
  #flat$cpg <- rownames(flat)
  flat$alias <- suppressWarnings(limma::alias2SymbolTable(flat$symbol))
  
  #eg <- toTable(org.Hs.egSYMBOL2EG)
  eg <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, 
                                keys=AnnotationDbi::keys(org.Hs.eg.db), 
                                columns=c("ENTREZID","SYMBOL"), 
                                keytype="ENTREZID"))
  colnames(eg) <- c("gene_id","symbol")
  
  m <- match(flat$alias,eg$symbol)
  flat$entrezid <- eg$gene_id[m]
  flat <- flat[!is.na(flat$entrezid),]
  
  # keep unique cpg by gene name annotation
  id<-paste(flat$cpg,flat$entrezid,sep=".")
  d <- duplicated(id)
  flat.u <- flat[!d,]
  flat.u
}

getMappedEntrezIDs <- function(sig.cpg, all.cpg=NULL, 
                               array.type=c("450K","EPICv2"), anno=NULL, 
                               genomic.features = c("ALL", "TSS200","TSS1500",
                                                    "Body","1stExon","3'UTR",
                                                    "5'UTR","ExonBnd"))
    # From a list of CpG sites, obtain the Entrez Gene IDs that are used for 
    # testing pathway enrichment
    # Belinda Phipson & Jovana Maksimovic
    # 10 February 2016
    # Updated 12 May 2020 to allow restricting sig.cpg by genomic features
{
    # check input
    sig.cpg <- as.character(sig.cpg)
    sig.cpg <- sig.cpg[!is.na(sig.cpg)]
    
    # array.type <- match.arg(toupper(array.type), c("450K","EPICv2"))  
    array.type <- "EPICv2"  
    genomic.features <- match.arg(genomic.features, c("ALL", "TSS200","TSS1500",
                                                      "Body", "1stExon","3'UTR",
                                                      "5'UTR","ExonBnd"), 
                                  several.ok = TRUE)
    
    if(length(genomic.features) > 1 & any(grepl("ALL", genomic.features))){
        message("All input CpGs are used for testing.") 
        genomic.features <- "ALL"   
    } 
    
    # Get annotaton in appropriate format
    if(is.null(anno)){
        flat.u <- .getFlatAnnotation(array.type)
    } else {
        flat.u <- .getFlatAnnotation(array.type,anno)
    }
    
    if(is.null(all.cpg)) {
        all.cpg <- unique(flat.u)
    } else {
        all.cpg <- as.character(all.cpg)
        all.cpg <- all.cpg[!is.na(all.cpg)]
        all.cpg <- unique(all.cpg)
    }
    
    # remove CpGs from annotation that are not in all.cpg
    m_all <- match(rownames(flat.u), all.cpg)
    flat.u = flat.u[!is.na(m_all),]
    
    # map CpG sites to entrez gene id's
    sig.cpg <- unique(sig.cpg)
    
    m1 <- match(rownames(flat.u), sig.cpg)

    # eg.sig <- flat.u$entrezid[!is.na(m1)]
    if(any(grepl("ALL", genomic.features))){
        eg.sig <- flat.u$entrezid[!is.na(m1)]
    } else {
        # select only genes with sig. CpGs mapping to certain genomic features
        eg.sig <- flat.u$entrezid[!is.na(m1) & flat.u$group %in% genomic.features]
    }
    
    eg.sig <- unique(eg.sig)
    if(length(eg.sig)==0) {
        stop("There are no genes annotated to the significant CpGs")
    }
    
    m2 <- match(rownames(flat.u),all.cpg)
    eg.all <- flat.u$entrezid[!is.na(m2)]
    
    freq_genes <- table(eg.all)
    eg.universe <- names(freq_genes)
    test.de <- as.integer(eg.universe %in% eg.sig)
    
    sorted.eg.sig <- eg.universe[test.de==1]
    
    multimap <- data.frame(table(rownames(flat.u)))
    multimap$Var1 <- as.character(multimap$Var1)
    m3 <- match(rownames(flat.u), multimap$Var1)

    flat.u$multimap <- multimap$Freq[m3]
    flat.u$inv.multimap <- 1/flat.u$multimap
    
    equivN <- tapply(flat.u$inv.multimap,flat.u$entrezid,sum)
    mm <- match(eg.universe,names(equivN))
    equivN <- equivN[mm]
    
    sig.flat <- flat.u[!is.na(m1),]

    if(any(grepl("ALL", genomic.features))){
        sig.flat <- flat.u[!is.na(m1),]
    } else {
        # select only CpGs that map to certain genomic features
        sig.flat <- flat.u[!is.na(m1) & flat.u$group %in% genomic.features, ]
    }
    
    fract <- data.frame(weight=pmin(tapply(1/sig.flat$multimap,
                                           sig.flat$entrezid,sum),
                                    1))
    
    m4 <- match(sorted.eg.sig,rownames(fract))
    fract.counts <- fract$weight[m4]
    
    out <- list(sig.eg = sorted.eg.sig, universe = eg.universe, 
                freq = freq_genes, equiv =  equivN, de = test.de, 
                fract.counts = data.frame(sigid=sorted.eg.sig,frac=fract.counts))
    out
}
gsameth <- function(sig.cpg, all.cpg=NULL, collection, 
                    array.type = c("450K","EPIC", "EPICv2"), plot.bias=FALSE, 
                    prior.prob=TRUE, anno=NULL, equiv.cpg = TRUE,
                    fract.counts = TRUE, 
                    genomic.features = c("ALL", "TSS200","TSS1500","Body",
                                         "1stExon","3'UTR","5'UTR","ExonBnd"),
                    sig.genes = FALSE)
  # Generalised version of gometh with user-specified gene sets 
  # Gene sets collections must be Entrez Gene ID
  # Can take into account probability of differential methylation 
  # based on numbers of probes on array per gene.
  # Belinda Phipson
  # 10 February 2016
  # Updated 21 March 2019
{
  
  if(!is.vector(sig.cpg))
    stop("Input CpG list is not a character vector")
#   array.type <- match.arg(toupper(array.type),c("450K","EPIC"))
  array.type <- "EPICv2"
  genomic.features <- match.arg(genomic.features, c("ALL", "TSS200","TSS1500",
                                                    "Body", "1stExon","3'UTR",
                                                    "5'UTR","ExonBnd"), 
                                several.ok = TRUE)
  
  if(length(genomic.features) > 1 & any(grepl("ALL", genomic.features))){
    # message("All input CpGs are used for testing.") 
    genomic.features <- "ALL"   
  } 
  
  if(array.type == "450K" & any(("ExonBnd" %in% genomic.features))){
      stop("'ExonBnd' is not an annotated feature on 450K arrays,
           please remove it from your genomic.feature parameter
           specification.") 
  }
  
  # Get mapped entrez gene IDs from CpG probe names
  if(!is.null(anno)){
    out <- getMappedEntrezIDs(sig.cpg=sig.cpg,all.cpg=all.cpg,
                              array.type=array.type, anno, 
                              genomic.features = genomic.features)
  } else {
    out <- getMappedEntrezIDs(sig.cpg=sig.cpg,all.cpg=all.cpg,
                              array.type=array.type, 
                              genomic.features = genomic.features)
  }
  sorted.eg.sig <- out$sig.eg
  eg.universe <- out$universe
  freq_genes <- out$freq
  test.de <- out$de
  frac <- out$fract.counts
  equiv <- out$equiv
  
  # Check collection is a list with character vectors
  if(!is.list(collection))
    collection <- list(collection=collection)
  collection <- lapply(collection, as.character)
  # Make sure gene set collections don't have any NAs
  collection <- lapply(collection, function(x) x[!is.na(x)])
  # Remove genes that are NOT in the universe from collections
  collection <- lapply(collection, function(x) x[x %in% eg.universe])
  # Remove collections with no genes left after universe filter
  inUniv <- sapply(collection, function(x) length(x) > 0)
  collection <- collection[inUniv]

  # Estimate prior probabilities
  if(prior.prob){
    if(equiv.cpg){ 
        # use "equivalent" no. of cpgs in odds calculation
        pwf <- .estimatePWF(D=test.de,bias=as.vector(equiv))
        if(plot.bias)
            .plotBias(D=test.de,bias=as.vector(equiv))
    } else {
        pwf <- .estimatePWF(D=test.de,bias=as.vector(freq_genes))
        if(plot.bias)
            .plotBias(D=test.de,bias=as.vector(freq_genes))
    }
  }
  
  results <- matrix(NA,ncol=4,nrow=length(collection))
  colnames(results) <- c("N","DE","P.DE","FDR")
  rownames(results) <- names(collection)
  results[,"N"] <- unlist(lapply(collection,length))
  if(sig.genes) SigGenesInSet <- rep(NA,length(collection))
  
  if(prior.prob & fract.counts){ 
      # use fractional counting to account for cpgs that map to multiple genes
      results[,"DE"] <- unlist(lapply(collection, function(x) 
          sum((sorted.eg.sig %in% x) * frac$frac)))
  } else {
      results[,"DE"] <- unlist(lapply(collection, function(x) 
          sum((sorted.eg.sig %in% x))))
  } 
  
  Nuniverse <- length(eg.universe)
  m <- length(sorted.eg.sig)
  
  # Hypergeometric test with prior probabilities
  if(prior.prob){
    for(i in 1:length(collection)){
      InSet <- eg.universe %in% collection[[i]]
      pw.red <- sum(pwf[InSet])/results[i,"N"]
      pw.white <- sum(pwf[!InSet])/(Nuniverse-results[i,"N"])
      odds <- pw.red/pw.white
      results[i,"P.DE"] <- BiasedUrn::pWNCHypergeo(results[i,"DE"],
                                                   results[i,"N"],
                                                   Nuniverse-results[i,"N"],
                                                   m,odds,lower.tail=FALSE) + 
          BiasedUrn::dWNCHypergeo(results[i,"DE"],
                                  results[i,"N"],
                                  Nuniverse-results[i,"N"],
                                  m,odds)
      if(sig.genes){
        # Get gene symbols of significant genes
        SigGenesEntrezID <- sorted.eg.sig[sorted.eg.sig %in% collection[[i]]]
        SigGenesSymbol <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, 
                                             keys = SigGenesEntrezID,
                                             columns = "SYMBOL"))
        SigGenesInSet[i] <- paste(SigGenesSymbol$SYMBOL,collapse=",")
      }
      if(results[i,"P.DE"]==0){
        message("Pvalue of exactly zero detected. Performing hypergeometric 
                test for gene set ", rownames(results)[i])
        results[i,"P.DE"] <- stats::phyper(q=results[i,"DE"]-0.5,m=m,
                                           n=Nuniverse-m,k=results[i,"N"],
                                           lower.tail=FALSE)
      }
    }
  }
  # Hypergeometric test without prior probabilities
  else{
    for(i in 1:length(collection)){
      results[i,"P.DE"] <- stats::phyper(q=results[i,"DE"]-0.5,m=m,
                                         n=Nuniverse-m,k=results[i,"N"],
                                         lower.tail=FALSE)
      if(sig.genes){
        # Get gene symbols of significant genes
        SigGenesEntrezID <- sorted.eg.sig[sorted.eg.sig %in% collection[[i]]]
        SigGenesSymbol <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, 
                                                               keys = SigGenesEntrezID,
                                                               columns = "SYMBOL"))
        SigGenesInSet[i] <- paste(SigGenesSymbol$SYMBOL,collapse=",")
      }
    }
  }
  results[,"FDR"] <- stats::p.adjust(results[,"P.DE"],method="BH")
  results[,"DE"] <- floor(results[,"DE"])
  if(sig.genes) data.frame(results, SigGenesInSet)
  else data.frame(results)
}

gometh <- function(sig.cpg, all.cpg=NULL, collection=c("GO","KEGG"), 
                   array.type = c("450K","EPIC", "EPICv2"), plot.bias=FALSE, 
                   prior.prob=TRUE, anno=NULL, equiv.cpg = TRUE, 
                   fract.counts = TRUE, 
                   genomic.features = c("ALL", "TSS200","TSS1500","Body",
                                       "1stExon","3'UTR","5'UTR","ExonBnd"),
                   sig.genes = FALSE)
  # Gene ontology testing or KEGG pathway analysis for Illumina methylation 
  # arrays based on goseq
  # Takes into account probability of differential methylation based on
  # numbers of probes on array per gene
  # Belinda Phipson
  # 28 January 2015. Last updated 1 September 2020.
  # EPIC functionality contributed by Andrew Y.F. Li Yim
{
  # array.type <- match.arg(toupper(array.type), c("450K","EPIC")) 
  array.type <- "EPICv2"    
  collection <- match.arg(toupper(collection), c("GO","KEGG"))
  genomic.features <- match.arg(genomic.features, c("ALL", "TSS200","TSS1500",
                                                    "Body", "1stExon","3'UTR",
                                                    "5'UTR","ExonBnd"), 
                                several.ok = TRUE)
  
  if(length(genomic.features) > 1 & any(grepl("ALL", genomic.features))){
    message("All input CpGs are used for testing.") 
    genomic.features <- "ALL"
  } 
  
  if(array.type == "450K" & any(grepl("ExonBnd", genomic.features))){
      stop("'ExonBnd' is not an annotated feature on 450K arrays,\n
           please remove it from your genomic.feature parameter\n
           specification.") 
  }
   
  if(collection == "GO"){
    go <- .getGO()
    result <- gsameth(sig.cpg=sig.cpg, all.cpg=all.cpg, collection=go$idList, 
                      array.type=array.type, plot.bias=plot.bias, 
                      prior.prob=prior.prob, anno=anno, equiv.cpg=equiv.cpg,
                      fract.counts=fract.counts, 
                      genomic.features = genomic.features,
                      sig.genes = sig.genes)
    result <- merge(go$idTable,result,by.x="GOID",by.y="row.names")
    rownames(result) <- result$GOID

  } else if(collection == "KEGG"){
    kegg <- .getKEGG()
    result <- gsameth(sig.cpg=sig.cpg, all.cpg=all.cpg, collection=kegg$idList, 
                      array.type=array.type, plot.bias=plot.bias, 
                      prior.prob=prior.prob, anno=anno, equiv.cpg=equiv.cpg,
                      fract.counts=fract.counts, 
                      genomic.features = genomic.features,
                      sig.genes = sig.genes)
    result <- merge(kegg$idTable,result,by.x="PathwayID",by.y="row.names")
    rownames(result) <- result$PathwayID
  }
  
  result[,-1]
}  

.getGO <- function(){
  if(!requireNamespace("org.Hs.eg.db", quietly = TRUE))
    stop("org.Hs.eg.db package required but not installed.")
  egGO2ALLEGS <- utils::getFromNamespace("org.Hs.egGO2ALLEGS", "org.Hs.eg.db")
  GeneID.PathID <- AnnotationDbi::toTable(egGO2ALLEGS)[,c("gene_id", "go_id", "Ontology")]
  d <- !duplicated(GeneID.PathID[, c("gene_id", "go_id")])
  GeneID.PathID <- GeneID.PathID[d, ]
  GOID.TERM <- suppressMessages(AnnotationDbi::select(GO.db::GO.db, 
                                                      keys=unique(GeneID.PathID$go_id), 
                                                      columns=c("GOID","ONTOLOGY","TERM"), 
                                                      keytype="GOID"))
  go <- tapply(GeneID.PathID$gene_id, GeneID.PathID$go_id, list)
    
  list(idList=go, idTable=GOID.TERM)
}

.getKEGG <- function(){
  GeneID.PathID <- limma::getGeneKEGGLinks(species.KEGG = "hsa", convert = TRUE)
  GeneID.PathID$PathwayID <- gsub("path:", "", GeneID.PathID$PathwayID)
  isna <- rowSums(is.na(GeneID.PathID[, 1:2])) > 0.5
  GeneID.PathID <- GeneID.PathID[!isna, ]
  ID.ID <- paste(GeneID.PathID[, 1], GeneID.PathID[, 2], sep = ".")
  d <- !duplicated(ID.ID)
  GeneID.PathID <- GeneID.PathID[d, ]
  PathID.PathName <- limma::getKEGGPathwayNames(species.KEGG = "hsa", 
                                         remove.qualifier = TRUE)
  #PathID.PathName$PathwayID <- paste0("path:", PathID.PathName$PathwayID)
  GeneID.PathID <- merge(GeneID.PathID, PathID.PathName, by="PathwayID")
  kegg <- tapply(GeneID.PathID$GeneID, GeneID.PathID$PathwayID, list)
  
  list(idList = kegg, idTable = PathID.PathName)
}  

.plotBias <- function(D,bias)
  # Plotting function to show gene level CpG density bias
  # Belinda Phipson
  # 5 March 2015
{
  o <- order(bias)
  splitf <- rep(1:100,each=200)[1:length(bias)]
  avgbias <- tapply(bias[o],factor(splitf),mean)
  sumDM <- tapply(D[o],factor(splitf),sum)
  propDM <- sumDM/table(splitf)
  graphics::par(mar=c(5,5,2,2))
  graphics::plot(avgbias,as.vector(propDM),
                 xlab="Number of CpGs per gene (binned)",
                 ylab="Proportion Differential Methylation",cex.lab=1.5,
                 cex.axis=1.2)
  graphics::lines(stats::lowess(avgbias,propDM),col=4,lwd=2)
}

.estimatePWF <- function(D,bias)
  # An alternative to goseq function nullp, which is transformation invariant
  # Belinda Phipson and Gordon Smyth
  # 6 March 2015
{
  prior.prob <- bias
  o <- order(bias)
  prior.prob[o] <- limma::tricubeMovingAverage(D[o],span=0.5)
  prior.prob
}

.getFlatAnnotation <- function(array.type=c("450K","EPIC"),anno=NULL)
  # flatten 450k or EPIC array annotation
  # Jovana Maksimovic
  # 18 September 2018
  # Updated 18 September 2018
  # Modified version of Belida Phipson's .flattenAnn code
{
  if(is.null(anno)){
    if(array.type=="450K"){
      anno <- minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19::IlluminaHumanMethylation450kanno.ilmn12.hg19)
    } else {
      anno <- minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19::IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
    }
  }
  
  # get rid of the non-CpG sites
  ann.keep<-anno[grepl("^cg",anno$Name),]
  
  # get rid of CpGs that are not annotated
  missing<-ann.keep$UCSC_RefGene_Name==""
  ann.keep<-ann.keep[!missing,]
  
  # get individual gene names for each CpG
  geneslist<-strsplit(ann.keep$UCSC_RefGene_Name,split=";")
  names(geneslist)<-rownames(ann.keep)
  
  grouplist<-strsplit(ann.keep$UCSC_RefGene_Group,split=";")
  names(grouplist)<-rownames(ann.keep)
  
  flat<-data.frame(symbol=unlist(geneslist),group=unlist(grouplist))
  flat$symbol<-as.character(flat$symbol)
  flat$group <- as.character(flat$group)
  
  flat$cpg<- substr(rownames(flat),1,10)
  
  #flat$cpg <- rownames(flat)
  flat$alias <- suppressWarnings(limma::alias2SymbolTable(flat$symbol))
  
  #eg <- toTable(org.Hs.egSYMBOL2EG)
  eg <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, 
                                keys=AnnotationDbi::keys(org.Hs.eg.db), 
                                columns=c("ENTREZID","SYMBOL"), 
                                keytype="ENTREZID"))
  colnames(eg) <- c("gene_id","symbol")
  
  m <- match(flat$alias,eg$symbol)
  flat$entrezid <- eg$gene_id[m]
  flat <- flat[!is.na(flat$entrezid),]
  
  # keep unique cpg by gene name annotation
  id<-paste(flat$cpg,flat$entrezid,sep=".")
  d <- duplicated(id)
  flat.u <- flat[!d,]
  flat.u
  # This randomly samples only 1 gene ID for multimapping CpGs
  #.reduceMultiMap(flat.u)
}

dotplot.enrichResult <- function(df, y = NULL, font.size = 18, title = "GO term analysis") {
  
  # Convert df to data frame if it's not already
  if (!is.data.frame(df)) {
    df <- as.data.frame(df)
  }
  
  # Calculate -log10(P.DE) and store it in a new column
  df$log10P <- -log10(df$P.DE)
  
  # Use ggplot to create the dot plot
  ggplot(df, aes(x = GeneRatio, y = .data[[y]], size = DE, color = log10P)) +
    geom_point() +
    scale_color_continuous(low = "red", high = "blue", name = "log P-value",
                           guide = guide_colorbar(reverse = TRUE)) +
    ggtitle(title) +
    labs(size = "Gene Count") +
    scale_size(range = c(3, 8)) +
    theme(plot.title = element_text(hjust = 0.5))  # Center-align the plot title

}

# gometh_plot <- function(sigCpGs=NULL, regulated=NULL, n=10, filename=NULL) {
#   # Identify go terms by given CpGs
#   gst_go <- gometh(sig.cpg=sigCpGs, all.cpg=all, plot.bias=TRUE, 
#                 collection = "GO", array.type="EPICv2", anno=annEPICv2)

#   # Top n GO categories
#   go_res <- topGSA(gst_go, number=n)
#   go_res$GeneRatio <- go_res$DE / go_res$N # Calculate the ratio of DE genes to total genes within the GO category

#   # Identify pathway by given CpGs
#   gst_kegg <- gometh(sig.cpg=sigCpGs, all.cpg=all, plot.bias=TRUE, 
#                 collection = "KEGG", array.type="EPICv2", anno=annEPICv2)

#   # Top n KEGG pathways
#   kegg_res <- topGSA(gst_kegg, number=n)
#   kegg_res$GeneRatio <- kegg_res$DE / kegg_res$N # Calculate the ratio of DE genes to total genes within the GO category

#   # Identify pathway by given CpGs
#   gst_kegg <- gometh(sig.cpg=sigCpGs, all.cpg=all, plot.bias=TRUE, 
#                 collection = GeneSet, array.type="EPICv2", anno=annEPICv2)

#   # Top n KEGG pathways
#   kegg_res <- topGSA(gst_kegg, number=n)
#   kegg_res$GeneRatio <- kegg_res$DE / kegg_res$N # Calculate the ratio of DE genes to total genes within the GO category

#   go_p <- dotplot.enrichResult(go_res, y = "TERM", title = "GO term analysis")
#   kegg_p <- dotplot.enrichResult(kegg_res, y = "Description", title = "KEGG pathway analysis")
#   return (list(go_p, kegg_p))
# }

gometh_plot <- function(sigCpGs=NULL, n=10, pathwaylist=NULL ) {

    # Identify go terms by given CpGs
    gst_go <- gometh(sig.cpg=sigCpGs, all.cpg=all, plot.bias=FALSE, 
                collection = "GO", array.type="EPICv2", anno=annEPICv2)

    # Top n GO categories
    go_res <- topGSA(gst_go, number=n)
    go_res$GeneRatio <- go_res$DE / go_res$N # Calculate the ratio of DE genes to total genes within the GO category

    # Identify pathway by given CpGs
    gst_kegg <- gometh(sig.cpg=sigCpGs, all.cpg=all, plot.bias=FALSE, 
                collection = "KEGG", array.type="EPICv2", anno=annEPICv2)

    # Top n KEGG pathways
    kegg_res <- topGSA(gst_kegg, number=n)
    kegg_res$GeneRatio <- kegg_res$DE / kegg_res$N # Calculate the ratio of DE genes to total genes within the GO category

    go_p <- dotplot.enrichResult(go_res, y = "TERM", title = "GO term analysis")
    kegg_p <- dotplot.enrichResult(kegg_res, y = "Description", title = "KEGG pathway analysis")
    print(go_p)
    print(kegg_p)

    # Identify other pathway by given CpGs
    res <- list()
    res.plot <- list()
    for (GeneSet in names(pathwaylist)){
        gmt_file <- pathwaylist[[GeneSet]]
        gmt <- gmtPathways(gmt_file)  # Load GMT file
        tmp <- gsameth(sig.cpg = sigCpGs, all.cpg = all, collection = gmt, anno= annEPICv2,
                plot.bias = FALSE, prior.prob = TRUE, sig.genes=TRUE, equiv.cpg=FALSE)
        # Top n pathways in other gene sets
        res[[GeneSet]] <- topGSA(tmp, number=n)
        res[[GeneSet]]$GeneRatio <- res[[GeneSet]]$DE / res[[GeneSet]]$N # Calculate the ratio of DE genes to total genes within the GO category
        res[[GeneSet]]$Description <- rownames(res[[GeneSet]])
        res.plot[[GeneSet]] <- dotplot.enrichResult(res[[GeneSet]], y = "Description", title = paste(GeneSet, "gene set pathway analysis", sep = " "))
        print(res.plot[[GeneSet]])
    }

    return (c(go_p, kegg_p, res.plot))
}

stat_signif_f <- function(tmp){
    # Separate the data into two groups
    group1 <- tmp[tmp$pheno == unique(tmp$pheno)[1],]$Beta_val
    group2 <- tmp[tmp$pheno == unique(tmp$pheno)[2],]$Beta_val

    # Perform Wilcoxon rank-sum test (Mann-Whitney U test)
    wilcox_test_result <- wilcox.test(group1, group2)

    # Extract p-value
    p_value <- wilcox_test_result$p.value

    # Convert p-value to star annotation
    star_annotation <- ifelse(p_value < 0.001, "***",
                    ifelse(p_value < 0.01, "**",
                    ifelse(p_value < 0.05, "*",
                    "NS")))

    return(star_annotation)
}