# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline 

# Load packages required to define the pipeline:
library(targets)
library(RcppParallel)
library(renv)
library(data.table)
library(tarchetypes) 


ncores=RcppParallel::defaultNumThreads() #Not nedded actually since tar_make is called from run.* file


# Set target options:
tar_option_set(
  packages = c("tibble","foreach","S4Vectors","ggfortify","ggrepel","gplots","ggplot2"), # packages that your targets need to run
  #imports = "cnv.methyl",
  format = "qs" # default storage format ="rds", more @ https://docs.ropensci.org/targets/reference/tar_target.html#storage-formats
  # Set other options as needed.
)

# Don't modify anything here:
# tar_make_clustermq() configuration (okay to leave alone):
options(clustermq.scheduler = "multicore"# multiprocess or multicore, LSF, SGE, Slurm etc.
)
# tar_make_future() configuration (okay to leave alone):
future::plan(future.callr::callr)


# Load the R scripts with your custom functions:
#lapply(list.files("R", full.names = TRUE, recursive = TRUE), source)

# source("R/functions.R") # Calls default function @meth_functions and adds/modifies
                        # them as needed for this project.
# source("R/deps.R") # Makes sure all needed packages are installed.
source("../meth_functions.R")


# produce_data <- function() {
#   expand.grid(samplesheet = c("a", "b"), model = c("c", "d"), normalization = c(1, 2, 3))
# }
# list(
#   tar_group_by(data, produce_data(), samplesheet, model),
#   tar_target(group, data, pattern = map(data))
# )

# VARs:
results_folder = "./results/"
analysis_folder = "./analysis/"
samplesheets <- c(
  All="data/SS.rds",
  noCereb = "data/ss_NoCereb.rds",
  noCereb_sub = "data/ss_NoCereb_sub.rds",
  adults = "data/ss_Adults.rds",
  adults_sub = "data/ss_Adults_sub.rds"

)


values <- tibble::tibble( # Use all possible combinations of input settings.
  method_function = rlang::syms(c("noob")),#, "pq","funn","noob_pq","Em")),
  data_paths = samplesheets,
  data_names = names(samplesheets)
  # Contrasts = paste(c("Sample_GroupControl_No-Sample_GroupExon1_No","Sample_GroupControl_No-Sample_GroupTAT_No" , "Sample_GroupExon1_No-Sample_GroupTAT_No",
                # "Sample_GroupControl_No-Sample_GroupControl_Dox","Sample_GroupExon1_No-Sample_GroupExon1_Dox" , "Sample_GroupTAT_No-Sample_GroupTAT_Dox"),
                # collapse = ";")
  # model=c("Sample_Group",)
)

# Pipline:
top_betas_N_target <- tar_target(top_betas_N,c(100,1000,5000,10000))
targets <- tarchetypes::tar_map(
  values = values,
  names = data_names, #"data_source", # Select columns from `values` for target names.
  tar_target(fparams,                                                            # Save paths & parameters
             make_results_dirs(subf=data_names, results_folder = results_folder,
                               analysis_folder = analysis_folder)),
  
  tar_target(samplesheet_path, data_paths, format = "file"),                     # Checks samplesheet for changes
  tar_target(samplesheet, readRDS(samplesheet_path)),                            # Reads sample sheet
  tar_target(ss,samplesheet[-1,]),                                               # Strips category tags
  tar_target(category,
             data.table::as.data.table(t(samplesheet[1,]),keep.rownames = T)),   # Category tags dict
  tar_target(nrgSet, cnv.methyl::read.metharray.exp.par(
    ss,                       # Read idats (cnv.methyl)
    arraytype="EPIC",force = T,ncores=ncores,folder="analysis/intermediate/"),     # force=TRUE for different nrows
    deployment = "worker", resources = tar_resources(                              # deploymernt = worker for parallel 
      future = tar_resources_future(resources = list(n_cores = 32)))),
  # tar_target(nrgSet,                                                             # Reads idats
  #            cnv.methyl::read.metharray.exp.par(
  #              folder = paste(data_names,"/analysis/intermediate/"), targets = ss,
  #              extended = T, force = T)),
  tar_target(rgSet,name_rgset2(nrgSet,ss)),                                       # Makes rgSet rownames == ss colnames
  # Qc report:
  tar_target(QC_plots, qc(                                                       # Makes qc plots:qcReport, density, beanplot, mean_qc 
    rgSet,sampGroups = "Sample_Group",sampNames="barcode",
    qc_folder = fparams[["qc_folder"]]),error = "continue",
    packages=c("S4Vectors","Biostrings","Biobase","minfi" )),

  # Calculate purity:
  tar_target(purity, cnv.methyl::purify(myLoad=rgSet)),                          # Calculates tumor purity as fraction of overall sample [0-1].
  tar_target(filtered, filter(                                                   # 1.- Filters: -probes: pval<0.01 
    targets=ss, rgSet=rgSet,sampGroups="Sample_Group",                           #              -samples: 10% probes fail
    qc_folder = fparams[["qc_folder"]])),                                        # 2.- Plots: Sample p.values barplot (colMeans)

  tar_target(normalize, method_function(filtered)),                              # Apply normalization; default --> Noob
  tar_target(clean, prep(normalize)),                                            # Preprocess: remove snps, remove Xreactive, Sex pred & removal
  
  tar_target(ss_clean,addcol(clean,purity)),                                     # Add tumor purity to sample sheet
  tar_target(save_ss_clean,                                                      # Save updated sample sheet
             write.table(ss_clean,paste0(fparams[["ss_clean_path"]],"/",
                                         "ss_clean.csv"),quote = F,sep = ",")),

  tar_target(plotvars,                                                           # Variables to plot on correlation matrix
             c(data.table::last(colnames(ss_clean)),"predictedSex",
               category[V1 %in% c("covs","batch"),rn])),

  tar_target(ann, minfi::getAnnotation(clean)),                                  # Annotate rgSet
  tar_target(betas, minfi::getBeta(clean),memory = "persistent"),                # Calculate beta values
  tar_target(top,top_beta(betas,n=top_betas_N),pattern=map(top_betas_N)),                                        # subsample for PCA & ploting
  tar_target(pca, pca_res(top),pattern=map(top)),                                                 # PCA

  tar_target(pca_corrplot,corpca(beta_top100 = top,                              # Correlation plot
                                 metadata=ss_clean,
                                 path=paste0(fparams[["corrplot_folder"]],"/",NROW(top)),
                                 filename=paste0(data_names,"_pca_corrplot",NROW(top),".png"),
                                 title=paste0("PC1-6 correlations with ",data_names," clinical vars(top ",NROW(top)," )")
                                 ),pattern=map(top)
  ),


  tar_target(bplots, bplot(pca,                                                  # Bi plots for PCA components
                           ss=ss_clean,
                           colgroup=category[V1 %in% c("covs","batch"),rn],
                           s="Type",
                           cols=NULL,
                           folder = paste0(fparams$bplots_folder,"/",NROW(pca$rotation))),
             packages = c("ggfortify","ggrepel","gplots","ggplot2"),
             pattern=map(pca)
  ),

  tar_target(model, mod(object = betas, group_var = "Sample_Group",              # Model with limma for diff meth 
                        contrasts = NULL,singular=T,
                        covs= category[V1 %in% c("batch"),rn],
                        metadata = ss_clean)
  ),

  tar_target(dmps_mod1, cnv.methyl::DMPextr(fit = model,                         # Toptable & stats
                                            ContrastsDM = colnames(model$contrasts),
                                            beta_normalized = betas,
                                            p.value = 0.95,
                                            mDiff = 0.01,
                                            ann = ann,
                                            writeOut = F
  )),
  tar_target(dmp_battery,priority = 1,                                           # DMPs distribution along params.
             apply_filter_dmps(
               dmps = dmps_mod1,path=paste0(fparams$dmp_folder,data_names)),error ="continue"),
  tar_target(dmps_f , filter_dmps(dmps_mod1, p.value = 0.05, mDiff = 0.2),error ="continue"),      # Choose filter for DMPs
  tar_target(save_dmps, data.table::fwrite(
    dmps_f,paste0(fparams$dmp_folder,data_names,"_dmps.txt")),error ="continue"),                  # Save DMPs
  tar_target(betas_DIF,betasdmps(betas,dmps_f,rgSet),error="continue"),                           # 
  tar_target(betas_DIF_full,betasdmps(betas,ann,rgSet),error="continue"),                         #
  # tar_target(write_betas_DIF, writexl::write_xlsx(betas_DIF,paste0(fparams$dmp_folder,"betas_DIF_",data_names,".xlsx"))),
  # tar_target(write_betas_DIF_full, writexl::write_xlsx(betas_DIF_full,paste0(fparams$dmp_folder,"betas_DIF_full_",data_names,".xlsx"))),


  tar_target(dmps_summary,                                                       # Summary statistics for DMPs
             summary_dmps(dmps_mod1, dir = fparams$dmp_folder,name=data_names),error ="continue"),
  tar_target(dmpplot_mod1, plotDMP(dmps_f,path=fparams[["dmpplots_folder"]]),error ="continue"),   # Barplots hipo/hyper, genomic region, CpG islands.

  tar_target(dmrs1,                                                              # Finds DMRs with dmrcate
             find_dmrs(betas,model,
                       fdr = 0.25, p.value = 0.1,betacutoff = 0.1, min.cpg=3)),
  tar_target(dmrs_battery,  priority = 1, error ="continue",                                       # DMRs distribution along params
             apply_filter_dmrs(
               dmrs = dmrs1,path=paste0(fparams$dmrs_folder,data_names))),
  tar_target(dmrs, filter_dmrs(dmrs1,p.value = 0.1, mDiff = 0.25, min.cpg=3)),   # Filters DMRs
  tar_target(save_dmrs,                                                          # Saves DMRs
             writexl::write_xlsx(
               dmrs, paste0(fparams$dmrs_folder,"_",data_names,".xlsx"))),
  tar_target(top4dmrs, {top4<-dmrs[,head(.SD,4),by=Contrast,
                                   .SDcols=c("Contrast","seqnames","start","no.cpgs","HMFDR","maxdiff","meandiff","overlapping.genes")]
  
    data.table::fwrite(top4,paste0(fparams$dmrs_folder,"_top4",data_names,".csv"))
    }),
  tar_target(dmrs_summary,                                                       # Summary stats for DMRs
             summary_dmrs(
               dmrs,path=paste0(fparams$dmrs_folder,"full_dmrs_summary",data_names,".csv"))),
  tar_target(vennDiags,  error ="continue",                                      # venn plot for genes comparing groupvar(contrasts?)
             venns(
               dmrs,groupvar="Contrast",res=paste0(fparams$dmrs_folder,"/VENNS/"))),
  
  tar_target(allpathways, gopath(dmrs,all.cpg=rownames(betas),n=Inf,ann=ann)),   # Pathways with all
  tar_target(hyperpathways,                                                      # Pathways hyper
             gopath(dmrs[meandiff>0,],
             all.cpg=rownames(betas),n=Inf,ann=ann,savepath=paste0(fparams$pathway_folder,"pathways_hyper.csv"))),
  tar_target(hypopathways,                                                       # Pathways with hypo
             gopath(dmrs[meandiff<0,],
                    all.cpg=rownames(betas),n=Inf,ann=ann,savepath=paste0(fparams$pathway_folder,"pathways_hypo.csv"))),
  tar_target(save_gopathall,
             writexl::write_xlsx(allpathways[FDR<1,],
                                 paste0(fparams$pathway_folder,data_names,"_allPathways.xlsx"))),
  tar_target(save_gopathhypo,
             writexl::write_xlsx(hypopathways[FDR<1,],
                                 paste0(fparams$pathway_folder,data_names,"_hypoPathways.xlsx"))),
  tar_target(save_gopathhyper,
             writexl::write_xlsx(hyperpathways[FDR<1,],
                                 paste0(fparams$pathway_folder,data_names,"_hyperPathways.xlsx")))
  )

AMN_values <- expand.grid( # Use all possible combinations of input settings.
  method_function = rlang::syms(c( "pq","funn","noob_pq","Em")),
  data_paths = samplesheets)

AMN_target<-tarchetypes::tar_map(
  values = AMN_values,
  names = method_function, #"data_source", # Select columns from `values` for target names.
  tar_target(fparams,                                                            # Save paths & parameters
             make_results_dirs(subf="AMN", results_folder = results_folder,
                               analysis_folder = analysis_folder)),
  
  tar_target(samplesheet_path, "data/ss_AMN.rds", format = "file"),                     # Checks samplesheet for changes
  tar_target(samplesheet, readRDS(samplesheet_path)),                            # Reads sample sheet
  tar_target(ss,samplesheet[-1,]),                                               # Strips category tags
  tar_target(category,
             data.table::as.data.table(t(samplesheet[1,]),keep.rownames = T)),   # Category tags dict
  tar_target(rgSet_input_file, "data/rgset_AMN.qs",format = "file"),
  tar_target(nrgSet, qs::qread(rgSet_input_file,
                               nthreads = RcppParallel::defaultNumThreads()),     # force=TRUE for different nrows
    deployment = "worker", resources = tar_resources(                              # deploymernt = worker for parallel 
      future = tar_resources_future(resources = list(n_cores = ncores)))),
  # tar_target(nrgSet,                                                             # Reads idats
  #            cnv.methyl::read.metharray.exp.par(
  #              folder = paste("AMN","/analysis/intermediate/"), targets = ss,
  #              extended = T, force = T)),
  tar_target(rgSet,name_rgset(nrgSet,ss)),                                       # Makes rgSet rownames == ss colnames
  # Qc report:
  tar_target(QC_plots, qc(                                                       # Makes qc plots:qcReport, density, beanplot, mean_qc 
    rgSet,sampGroups = "Sample_Group",sampNames="barcode",
    qc_folder = fparams[["qc_folder"]]),error = "continue" ),
  
  # Calculate purity:
  tar_target(purity, cnv.methyl::purify(myLoad=rgSet)),                          # Calculates tumor purity as fraction of overall sample [0-1].
  tar_target(filtered, filter(                                                   # 1.- Filters: -probes: pval<0.01 
    targets=ss, rgSet=rgSet,sampGroups="Sample_Group",                           #              -samples: 10% probes fail
    qc_folder = fparams[["qc_folder"]])),                                        # 2.- Plots: Sample p.values barplot (colMeans)
  
  tar_target(normalize, noob(filtered)),                              # Apply normalization; default --> Noob
  tar_target(clean, prep(normalize)),                                            # Preprocess: remove snps, remove Xreactive, Sex pred & removal
  tar_target(ss_clean,addcol(clean,purity)),                                     # Add tumor purity to sample sheet
  tar_target(save_ss_clean,                                                      # Save updated sample sheet
             write.table(ss_clean,paste0(fparams[["ss_clean_path"]],"/",
                                         "ss_clean.csv"),quote = F,sep = ",")),
  
  tar_target(plotvars,                                                           # Variables to plot on correlation matrix
             c(data.table::last(colnames(ss_clean)),"predictedSex",
               category[V1 %in% c("covs","batch"),rn])),
  
  tar_target(ann, minfi::getAnnotation(clean)),                                  # Annotate rgSet
  tar_target(betas, minfi::getBeta(clean),memory = "persistent"),                # Calculate beta values
  tar_target(top,top_beta(betas,n=top_betas_N),pattern=map(top_betas_N)),                                        # subsample for PCA & ploting
  tar_target(pca, pca_res(top),pattern=map(top)),                                                 # PCA
  
  tar_target(pca_corrplot,corpca(beta_top100 = top,                              # Correlation plot
                                 metadata=ss_clean,
                                 path=paste0(fparams[["corrplot_folder"]],"/",NROW(top)),
                                 filename=paste0("AMN","_pca_corrplot",NROW(top),".png"),
                                 title=paste0("PC1-6 correlations with ","AMN"," clinical vars(top ",NROW(top)," )")
  ),pattern=map(top)
  ),
  
  tar_target(ss_clean2,ss_clean[,sentrixID:=sapply(barcode,function(x) unlist(strsplit(x,"_"))[1])]),
  tar_target(bplots, bplot(pca,                                                  # Bi plots for PCA components
                           ss=ss_clean2,
                           colgroup=c(category[V1 %in% c("covs","batch"),rn],sentrixID),
                           s="Type",
                           cols=NULL,
                           folder = paste0(fparams$bplots_folder,"/",NROW(pca$rotation),"/")),
             packages = c("ggfortify","ggrepel","gplots","ggplot2"),
             pattern=map(pca)
  ),
  
  tar_target(model, mod(object = betas, group_var = "Sample_Group",              # Model with limma for diff meth 
                        contrasts = NULL,singular=T,
                        covs= category[V1 %in% c("batch"),rn],
                        metadata = ss_clean)
  ),
  
  tar_target(dmps_mod1, cnv.methyl::DMPextr(fit = model,                         # Toptable & stats
                                            ContrastsDM = colnames(model$contrasts),
                                            beta_normalized = betas,
                                            p.value = 0.95,
                                            mDiff = 0.01,
                                            ann = ann,
                                            writeOut = F
  )),
  tar_target(dmp_battery,priority = 1,                                           # DMPs distribution along params.
             apply_filter_dmps(
               dmps = dmps_mod1,path=paste0(fparams$dmp_folder,"AMN")),error ="continue"),
  tar_target(dmps_f , filter_dmps(dmps_mod1, p.value = 0.05, mDiff = 0.2),error ="continue"),      # Choose filter for DMPs
  tar_target(save_dmps, data.table::fwrite(
    dmps_f,paste0(fparams$dmp_folder,"AMN","_dmps.txt")),error ="continue"),                  # Save DMPs
  tar_target(betas_DIF,betasdmps(betas,dmps_f,rgSet),error="continue"),                           # 
  tar_target(betas_DIF_full,betasdmps(betas,ann,rgSet),error="continue"),                         #
  # tar_target(write_betas_DIF, writexl::write_xlsx(betas_DIF,paste0(fparams$dmp_folder,"betas_DIF_","AMN",".xlsx"))),
  # tar_target(write_betas_DIF_full, writexl::write_xlsx(betas_DIF_full,paste0(fparams$dmp_folder,"betas_DIF_full_","AMN",".xlsx"))),
  
  
  tar_target(dmps_summary,                                                       # Summary statistics for DMPs
             summary_dmps(dmps_mod1, dir = fparams$dmp_folder,name="AMN"),error ="continue"),
  tar_target(dmpplot_mod1, plotDMP(dmps_f,path=fparams[["dmpplots_folder"]]),error ="continue"),   # Barplots hipo/hyper, genomic region, CpG islands.
  
  tar_target(dmrs1,                                                              # Finds DMRs with dmrcate
             find_dmrs(betas,model,
                       fdr = 0.25, p.value = 0.1,betacutoff = 0.1, min.cpg=3)),
  tar_target(dmrs_battery,  priority = 1, error ="continue",                                       # DMRs distribution along params
             apply_filter_dmrs(
               dmrs = dmrs1,path=paste0(fparams$dmrs_folder,"AMN"))),
  tar_target(dmrs, filter_dmrs(dmrs1,p.value = 0.1, mDiff = 0.25, min.cpg=3)),   # Filters DMRs
  tar_target(save_dmrs,                                                          # Saves DMRs
             writexl::write_xlsx(
               dmrs, paste0(fparams$dmrs_folder,"_","AMN",".xlsx"))),
  tar_target(top4dmrs, {top4<-dmrs[,head(.SD,4),by=Contrast,
                                   .SDcols=c("Contrast","seqnames","start","no.cpgs","HMFDR","maxdiff","meandiff","overlapping.genes")]
  
  data.table::fwrite(top4,paste0(fparams$dmrs_folder,"_top4","AMN",".csv"))
  }),
  tar_target(dmrs_summary,                                                       # Summary stats for DMRs
             summary_dmrs(
               dmrs,path=paste0(fparams$dmrs_folder,"full_dmrs_summary","AMN",".csv"))),
  tar_target(vennDiags,  error ="continue",                                      # venn plot for genes comparing groupvar(contrasts?)
             venns(
               dmrs,groupvar="Contrast",res=paste0(fparams$dmrs_folder,"/VENNS/"))),
  
  tar_target(allpathways, gopath(dmrs,all.cpg=rownames(betas),n=Inf,ann=ann)),   # Pathways with all
  tar_target(hyperpathways,                                                      # Pathways hyper
             gopath(dmrs[meandiff>0,],
                    all.cpg=rownames(betas),n=Inf,ann=ann,savepath=paste0(fparams$pathway_folder,"pathways_hyper.csv"))),
  tar_target(hypopathways,                                                       # Pathways with hypo
             gopath(dmrs[meandiff<0,],
                    all.cpg=rownames(betas),n=Inf,ann=ann,savepath=paste0(fparams$pathway_folder,"pathways_hypo.csv"))),
  tar_target(save_gopathall,
             writexl::write_xlsx(allpathways[FDR<0.05,],
                                 paste0(fparams$pathway_folder,"AMN","_allPathways.xlsx"))),
  tar_target(save_gopathhypo,
             writexl::write_xlsx(hypopathways[FDR<0.05,],
                                 paste0(fparams$pathway_folder,"AMN","_hypoPathways.xlsx"))),
  tar_target(save_gopathhyper,
             writexl::write_xlsx(hyperpathways[FDR<0.05,],
                                 paste0(fparams$pathway_folder,"AMN","_hyperPathways.xlsx")))
  
)

list(top_betas_N_target,AMN_target)
#
  
  # tar_target(save_gopath,
  #            writexl::write_xlsx(pathways[FDR<1,],
  #                                paste0(fparams$pathway_folder,data_names,"_Pathways.xlsx"))),
  # NULL
# )
# # list(
# #   tar_target(samplesheet_path, "data/ss_H358.rds", format = "file"),
# #   tar_target(samplesheet, readRDS(samplesheet_path)),
# #   # tar_target(ss,samplesheet),
# #   tar_target(ss,samplesheet[-1,]),
# #   tar_target(category,samplesheet[1,]),
# #   tar_target(rgSet, cnv.methyl::read.metharray.exp.par(ss,folder="ana",extended=T)),
# #   # Qc report:
# #   tar_target(QC_plots, cnv.methyl::qc(rgSet,sampGroups = "condition")),
# #
# #   # Calculate purity:
# #   tar_target(purity, cnv.methyl::purify(myLoad=rgSet)),
# #   tar_target(filtered, filter(targets=ss, rgSet=rgSet,sampGroups="condition")),
# #   targets,
# #   # tar_target( alldmps,
# #   #             list(apply(tidyr::expand_grid(c("dmps_ANA","dmps_SANDRA","dmps_WT"),c("noob", "pq","funn","noob_pq","Em")),1,function(x)paste(x,collapse = "_")))
# #   # )
# #
# #   #tar_combine(combined_gopath_ANA,
# #   #list(apply(tidyr::expand_grid(c("gopath_ANA"),c("noob", "pq","funn","noob_pq","Em")),1,function(x)paste(x,collapse = "_")))
# #   NULL
# #   )
# #
# # subset(ss,!is.na(ss$ANA_dom1))
# # technical_features <- c("Sample_Name", "organism", "Basename", "barcode")
# #
# #
# #  mval<- getM(npq[,!is.na(ss$ANA_dom1)])
# #  pheno<-ss[!is.na(ANA_dom1)]
# #  mod <- model.matrix(
# #    formula(paste(" ~ ANA_dom1 +", paste(batch,sep="+",collapse="+")))
# #    ,data = pheno
# #  )
# #  mod0 <- model.matrix(
# #    formula(paste(" ~", paste(batch,sep="+",collapse="+")))
# #    ,data = pheno
# #  )
# #  sva.results <- sva(mval, mod, mod0)
# # #
# # design <- model.matrix(
# #   formula(paste(" ~", paste(covs,sep="+",collapse="+")))
# #   ,data = ss
# #   )
# # n.sv = sva::num.sv(betas,design,method="leek")
# #
# # svobj = sva(rgSet,mod,mod0,n.sv=n.sv)
# #
# # 
# # # subset ss
# # BiocManager::install(c("Biobase", "conumee", "GenomicRanges", "IlluminaHumanMethylation450kanno.ilmn12.hg19", "IlluminaHumanMethylationEPICanno.ilm10b4.hg19", "impute", "IRanges", "limma", "maxprobes", "minfi", "SummarizedExperiment"))
# 
