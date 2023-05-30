#' @title Pathways dictionary.
#' @description Define 3 targets:
#' 1. get pathways from GO
#' 2. get Pathways from Kegg
#' 3. merge data
#' @return A list of target objects.
#' @export
#' @param KEGG boolean, TRUE = include KEGG pathways 
#' @param GO boolean, TRUE = include GO pathways


pathways_dict <- function(KEGG=T,GO=T ) {
  sym_KEGG = as.symbol(KEGG)
  sym_GO = as.symbol(GO)
  command_GO <- substitute(getGO(DO=GO),env = list(GO = GO))
  command_KEGG <- substitute(getKEGG(DO=KEGG),env = list(KEGG = KEGG))
  list(
    tar_target_raw("GO_dict",command_GO , deployment = "main"),
    tar_target_raw("KEGG_dict", command_KEGG, deployment = "main"),
    tar_target_raw("p_dict", quote(rbind(GO_dict,KEGG_dict)),
                   deployment = "main") 
  )
}

#' @title Pathways target factory.
#' @description Define 4 targets:
#' 1. Track the user-supplied data file.
#' 2. Read the data using `read_data()` (defined elsewhere).
#' 3. Strip category away.
#' 4. Generate category data.table.
#' @return A list of target objects.
#' @export
#' @param ann Character, vector with cg probe names.
#' @param dmrs gRanges object with DMRs
#' @param n Number of pathways to show on the summary
#' @param type hyper/hypo/all
pathways <- function(ann=F, dmrs=F, n=Inf, type=c("Hyper","Hypo","All"),KEGG=T,GO=T ) {
  list(
    tar_target_raw("gsa",quote(data.table::fwrite(p_dict,"path_dict.csv")))
  )
}
#' @title Go pathways dictionary.
#' @description returns ENTREZID genes from GO pathways
#' @return A data.table with obs=ENTREZID & cols=[(GO)ID,ENTREZID,TERM,ONTOLOGY]
#' @export
#' @param ev.filter Character, vector with evidence filters to remove from query
#' @param DO boolean if FALSE returns an empty data.table
getGO <- function(ev.filter=NULL,DO=T){
  GO <- data.table::data.table(ID=character(0),ENTREZID=character(0),TERM=character(0),ONTOLOGY=character(0))
  if(DO==TRUE){
      
    if(!requireNamespace("org.Hs.eg.db", quietly = TRUE))
      stop("org.Hs.eg.db package required but not installed.")
    egGO2ALLEGS <- utils::getFromNamespace("org.Hs.egGO2ALLEGS", "org.Hs.eg.db")
    GeneID.PathID <- data.table::data.table(
      AnnotationDbi::toTable(egGO2ALLEGS))
    names(GeneID.PathID)<-c("GeneID","GOID","Evidence","ONTOLOGY")
    GeneID.PathID<-GeneID.PathID[!Evidence %in% ev.filter,list(GeneID,GOID,ONTOLOGY)]
    GOID.TERM <- suppressMessages(AnnotationDbi::select(GO.db::GO.db, 
                                                        keys=unique(GeneID.PathID$GOID), 
                                                        columns=c("GOID","ONTOLOGY","TERM"), 
                                                        keytype="GOID"))
    dt_GO_dict<-merge(unique(GeneID.PathID),GOID.TERM)
    GO <- dt_GO_dict[,.(ID=GOID,ENTREZID=GeneID,TERM=TERM,ONTOLOGY=ONTOLOGY)]
  }
  return(GO)
}
#' @title KEGG pathways dictionary.
#' @description returns ENTREZID genes from KEGG pathways
#' @return A data.table with obs=ENTREZID & cols=[(KEGG)ID,ENTREZID,TERM,ONTOLOGY]
#' @export
#' @param DO boolean if FALSE returns an empty data.table
getKEGG <- function(DO=T){
  KEGG <- data.table::data.table(ID=character(0),ENTREZID=character(0),TERM=character(0),ONTOLOGY=character(0))
  if(DO==TRUE){
    if(!requireNamespace("data.table", quietly = TRUE))
      stop("data.table package required but not installed.")
    dt_gene_pathway <- unique(na.omit(data.table::data.table(limma::getGeneKEGGLinks(
      species.KEGG = "hsa", convert = TRUE))))
    dt_pathway_term <- unique(na.omit(data.table::data.table(limma::getKEGGPathwayNames(
      species.KEGG = "Hsa", remove.qualifier = TRUE))))
    dt_KEGG_dict <- merge(dt_gene_pathway,dt_pathway_term)
    KEGG<-dt_KEGG_dict[,.(ID=PathwayID,ENTREZID=GeneID,TERM=Description,ONTOLOGY="KEGG")]
  }
  return(KEGG)
}



prepare_dmrs<-function(dmrs){
  require(GenomicRanges)
  object<-GRanges(object)
  cont<-unique(object$Contrast)
}
