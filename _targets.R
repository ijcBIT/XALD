# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline 

# Load packages required to define the pipeline:
library(targets)
library(RcppParallel)
# library(renv)
library(data.table)
library(tarchetypes) 
library(future)
library(future.batchtools)
NCORES=RcppParallel::defaultNumThreads() #Not nedded actually since tar_make is called from run.* file
Sys.setenv(R_LIBS_USER="/mnt/beegfs/idevillasante/Projects/XALD/renv/library/R-4.2/x86_64-pc-linux-gnu")

# Set target options:
tar_option_set(
  # packages = c("tibble","foreach","S4Vectors","ggfortify","ggrepel","gplots","ggplot2"), # packages that your targets need to run
  packages = c("minfi"),
  #imports = "cnv.methyl",
  format = "qs" # default storage format ="rds", more @ https://docs.ropensci.org/targets/reference/tar_target.html#storage-formats
  # Set other options as needed.
)

# Don't modify anything here:
# tar_make_clustermq() configuration (okay to leave alone):
options(clustermq.scheduler = "multicore"# multiprocess or multicore, LSF, SGE, Slurm etc.
)


# SSH:
# options(clustermq.scheduler = "ssh",
#         clustermq.ssh.host = "idevillasante@minastirith", # set this up in your local ~/.ssh/config
#         clustermq.ssh.log = "~/ssh_proxy.log", # log file on your HPC
#         clustermq.ssh.timeout = 30, # if changing the default connection timeout
#         clustermq.template = "ssh.template" # if using your own template
# )

# # SLURM:
minastirith<-"10.50.2.100"
# options(
#   clustermq.scheduler = "slurm",
#   clustermq.template = "~/Projects/XALD/clustermq.template", #"~/.clustermq.template", # if using your own template
#   clustermq.host = "10.50.2.100",
#   clustermq.worker.timeout=6000
#   )
# 
# # )


# tar_make_future() configuration (okay to leave alone):
future::plan(future.callr::callr,workers=NCORES)
# future::plan(batchtools_slurm,
#              template = "batchtools.slurm.tmpl",
#              resources = list(nodes = "1:ppn=12", vmem = "5gb",ntasks=100,walltime="2-0:0:0",memory="32G"))

# Load the R scripts with your custom functions:
#lapply(list.files("R", full.names = TRUE, recursive = TRUE), source)

# source("R/functions.R") # Calls default function @meth_functions and adds/modifies
                        # them as needed for this project.
# source("R/deps.R") # Makes sure all needed packages are installed.
source("~/Projects/meth_functions.R")

projectdir<-"~/Projects/XALD/"

# VARs:
results_folder = "~/Projects/XALD/results/batch/"
analysis_folder = "~/Projects/XALD/analysis/batch/"
samplesheets <- c(
  inhouse = "data/SS.rds",
  cereb = "data/ss_Cereb.rds",
  noCereb = "data/ss_NoCereb.rds",
  noCereb_sub = "data/ss_NoCereb_sub.rds",
  adults = "data/ss_Adults.rds",
  adults_sub = "data/ss_Adults_sub.rds"
  # GEO = "data/ss_AMN.rds"
)

# samplesheets<-samplesheets[2]
# 
# [model 1: Adults] meth ~ Sample_Group/Condition (XALD,cALD,ctl.Adult) adults=ss_Adults.rds
# [model 2: Cereb] meth ~ Sample_Group/Condition (child + Adult cereb)
# [model 3: Age] meth ~ Type + Age (no sample group only age and disease/control , same samples as model2)
# 
values <- tibble::tibble(
  method_function = c(rep(rlang::syms(c("noob")),6),rlang::syms("swan"),rlang::syms("swan")),
  data_paths = c(m1 = "data/ss_Adults.rds",
                 m2a = "data/ss_NoCereb.rds",
                 m2 = "data/ss_NoCereb.rds",
                 m3 = "data/ss_NoCereb.rds",
                 m4 = "data/SS.rds",
                 m5 = "data/ss_NoCereb.rds",
                 mparams = "data/SS.rds",
                 mparams_full = "data/ss_Adults.rds"
                 ),
  covs.formula = c(
    m1 = "~0 + Condition",
    m2a = "~0 + Condition",
    m2 = "~0 + Condition + Age",
    m3 = "~0 + Type + Age",
    m4 = "~0 + Condition + Age",
    m5 = "~0 + Condition + Age + purity",
    mparams = "~0 + Condition ",
    mparams_full = "~0 + Condition "
  ),
  data_names = c("m1","m2a","m2","m3","m4","m5","mparams", "mparams_full"),
  group = c("Condition","Condition","Condition","Type","Condition","Condition", "Condition","Condition")
  
  
)

# values <- tibble::tibble( # Use all possible combinations of input settings.
#   method_function = rlang::syms(c("noob")),#, "pq","funn","noob_pq","Em")),
#   data_paths = samplesheets,
#   data_names = names(samplesheets)
# )
values$ncores<-sapply(values$data_paths,function(x)min(RcppParallel::defaultNumThreads(),NROW(readRDS(x))))

idcol<-"Sample_Name"
# Pipeline:
# method_function_target <- tar_target(method_function,rlang::syms(c("noob", "pq","funn","noob_pq")))#,"Em"))
top_betas_N_target <- tar_target(top_betas_N,c(100,1000,5000,10000),deployment = "main")


# Pipline:

targets <- tarchetypes::tar_map(
  # !!!!! tar_target(dmps_f , filter_dmps(dmpsod1, p.value = 0.01, mDiff = 0.05) ¡¡¡¡ & minfi::preprocessswan
  values = values[7,],

  names = data_names, #"data_source", # Select columns from `values` for target names.
  tar_target(fparams,                                                            # Save paths & parameters
             make_results_dirs(subf=data_names, results_folder = results_folder,
                               analysis_folder = analysis_folder),deployment="main"),

  tar_target(samplesheet_path, data_paths, format = "file"),                     # Checks samplesheet for changes
  tar_target(ss, readRDS(samplesheet_path),deployment="main"),                                     # Reads sample sheet
  tar_target(var_category, as.data.table(attributes(ss)[c("category","names")]),    # type of variable 
             deployment="main"),      
  # tar_target(nrgSet, cnv.methyl::read.metharray.exp.par(
  #   ss,                       # Read idats (cnv.methyl)
  #   arraytype="EPIC",force = T,ncores=ncores,folder="analysis/intermediate/"),     # force=TRUE for different nrows
  #   deployment = "worker", resources = tar_resources(                              # deploymernt = worker for parallel
  #     future = tar_resources_future(resources = list(n_cores = ncores)))),
  # tar_target(nrgSet,                                                             # Reads idats
  #            cnv.methyl::read.metharray.exp.par(
  #              folder = paste(data_names,"/analysis/intermediate/"), targets = ss,
  #              extended = T, force = T,ncores=ncores),
  #            deployment = "worker",
  #            # resources = tar_resources(
  #            #   future = tar_resources_future(resources = list(n_cores = ncores))
  #            # )
  #            ),
  tar_target(nrgSet, minfi::read.metharray.exp(targets = ss)),

#   tar_target(QC_plots, qc(                                                       # Makes qc plots:qcReport, density, beanplot, mean_qc 
#     rgSet,sampGroups = "Sample_Group",sampNames="barcode",
#     qc_folder = fparams[["qc_folder"]]),error = "continue",
#     packages=c("S4Vectors","Biostrings","Biobase","minfi" )),
# 
#   # Calculate purity:
#   tar_target(purity, cnv.methyl::purify(myLoad=rgSet)),                          # Calculates tumor purity as fraction of overall sample [0-1].
#   tar_target(filtered, filter(                                                   # 1.- Filters: -probes: pval<0.01 
#     targets=ss, rgSet=rgSet,sampGroups="Sample_Group",                           #              -samples: 10% probes fail
#     qc_folder = fparams[["qc_folder"]])),    
# 
# 
# 
# 
# 
# 
# AMN_target<-list(
#   tar_target(fparams,                                                            # Save paths & parameters
#              make_results_dirs(subf="AMN_reduced", results_folder = results_folder,
#                                analysis_folder = analysis_folder),deployment = "main"),
#   tar_target(samplesheet_path, c("~/Projects/XALD/data/SS.rds")),#,"data/ss_AMN.rds")),                     # Checks samplesheet for changes
#   tar_target(files, samplesheet_path, format = "file"), 
#              # pattern = map(samplesheet_path)),
#   tar_target(samplesheet, readRDS(samplesheet_path),deployment = "main"),
#              # pattern = map(files)),       # Reads sample sheet
#   tar_target(ss,samplesheet[-1,Age:=as.numeric(Age)],deployment = "main"),
#              # pattern = map(samplesheet)),                                               # Strips category tags
#   tar_target(category,
#              data.table::as.data.table(t(samplesheet[1,]),keep.rownames = T),
#              deployment = "main"),
#              # pattern = map(ss)),   # Category tags dict
#   
#   tar_target(rgSet_input_file, "~/Projects/XALD/data/rgset_AMN.qs",format = "file",
#   deployment = "main"),
#   tar_target(nrgSet, qs::qread(rgSet_input_file,
#                     nthreads = RcppParallel::defaultNumThreads()
#                     ),
#   deployment = "main",memory="transient",
#   resources = tar_resources(
#     future = tar_resources_future(
#       # plan = tweak(
#       #   batchtools_slurm,
#       #   template = "batchtools.slurm.tmpl",
#       #   resources = list(num_cores = ncores)
#       # )
#       resources=list(n_cores=ncores)
#       )
#     )
#   ),
#   # tar_target(nrgSet,                                                             # Reads idats
#   #            cnv.methyl::read.metharray.exp.par(
#   #              , targets = ss,
#   #              extended = T, force = T)),
  tar_target(rgSet,
             name_rgset(nrgSet,ss,exclude="206702460034_R05C01",newname=idcol),
             deployment="worker",
             memory="persistent"),                                               # Makes rgSet rownames == ss colnames
  # Qc report:
  tar_target(QC_plots, qc(                                                       # Makes qc plots:qcReport, density, beanplot, mean_qc 
    rgSet,sampGroups = "Sample_Group",sampNames=idcol,idcol=idcol,
    qc_folder = fparams[["qc_folder"]]),packages = "minfi",error = "continue",deployment = "worker" ),
  
  # Calculate purity:
  tar_target(purity, cnv.methyl::purify(myLoad=rgSet),error="continue"),         # Calculates tumor purity as fraction of overall sample [0-1].
  tar_target(filtered, filter(                                                   # 1.- Filters: -probes: pval<0.01 
    targets=ss, rgSet=rgSet,sampGroups="Sample_Group",                           #              -samples: 10% probes fail
    qc_folder = fparams[["qc_folder"]]),                                        # 2.- Plots: Sample p.values barplot (colMeans)
    deployment = "worker",memory = "transient"  ),
  # Until here only loaded and cleaned the data now try different normalizations: 


  tar_target(normalize, method_function(filtered),
             deployment = "worker", packages = "minfi"),                              # Apply normalization; default --> Noob
  tar_target(clean, prep(normalize, remove_sex=F, arraytype="EPIC"),
             deployment = "worker"),
  tar_target(clean2, createbatch(clean), deployment = "worker"),
  tar_target(ss_clean,addcol(clean2,newcol=purity,cname = "purity"),
             deployment = "main"),                                               # Add tumor purity to sample sheet
  tar_target(save_ss_clean,                                                      # Save updated sample sheet
             write.table(ss_clean,paste0(fparams[["ss_clean_path"]],"/",
                                         "ss_clean.csv"),quote = F,sep = ","),
             deployment = "main"),

  tar_target(plotvars,                                                           # Variables to plot on correlation matrix
             c(data.table::last(colnames(ss_clean)),"predictedSex",
               var_category[category %in% c("covs","batch"),names]),deployment = "main"),

  tar_target(ann, minfi::getAnnotation(clean)),                                  # Annotate rgSet
  # tar_target(betas, minfi::getBeta(clean),memory = "persistent"),                # Calculate beta values
  tar_target(betas,getBetas(clean2)),
  tar_target(top,top_beta(betas,n=top_betas_N),pattern=map(top_betas_N)),                                        # subsample for PCA & ploting
  tar_target(pca, pca_res(top),pattern=map(top)),                                                 # PCA

  tar_target(pca_corrplot,corpca(beta_top100 = top,                              # Correlation plot
                                 metadata=ss_clean,
                                 idcol=idcol,
                                 path=paste0(fparams[["corrplot_folder"]],"/",NROW(top)),
                                 filename=paste0(data_names,"_pca_corrplot",NROW(top),".png"),
                                 title=paste0("PC1-6 correlations with ",data_names," clinical vars(top ",NROW(top)," )")
  ),pattern=map(top)
  ),
# 
#   tar_target(bplots, bplot(pca,                                                  # Bi plots for PCA components
#                            ss=ss_clean,
#                            colgroup=c(var_category[category %in% c("covs","batch"),names],"predictedSex","purity"),
#                            s="Type",
#                            cols=NULL,
#                            overlap=1,
#                            alfa=0.8,
#                            labs=T,
#                            idcol=idcol,
#                            folder = paste0(fparams$bplots_folder,"/",NROW(pca$rotation),"/")),
#              packages = c("ggfortify","ggrepel","gplots","ggplot2"),
#              pattern=map(pca),priority = 1
#   ),

  ## -- (added  to rgSet betas & ss_clean from here to end)

 tar_target(model, mod(object = betas, group_var = group,              # Model with limma for diff meth
                        singular=F,
                        covs.formula = covs.formula,
                        covs= "Condition",#var_category[category %in% c("batch"),names],
                        metadata = ss_clean,
                        idcol=idcol),
             deployment = "worker"
             ),

  tar_target(dmpsod1, cnv.methyl::DMPextr(fit = model,                         # Toptable & stats
                                            ContrastsDM = colnames(model$contrasts),
                                            beta_normalized = betas,
                                            p.value = 0.95,
                                            mDiff = 0.01,
                                            ann = ann,
                                            writeOut = F),
             deployment = "worker"
             ),
  tar_target(save_dmps1,write.table(dmpsod1,paste0(fparams$dmp_folder,as.character(quote(method.function)),"dmps1.txt"))),
  tar_target(dmp_battery,priority = 1,                                           # DMPs distribution along params.
             apply_filter_dmps(
               dmps = dmpsod1,path=paste0(fparams$dmp_folder,data_names)),
             error ="continue",deployment = "worker",memory = "transient"),
  # tar_target(dmps_f , filter_dmps(dmpsod1, p.value = 0.01, mDiff = 0.05),
  #            error ="continue"),      # Choose filter for DMPs
  tar_target(dmps_f ,dmpsod1[abs(dmpsod1$diff_meanMeth) > 0.05 & dmpsod1$P.Value > 0.01 & abs(dmpsod1$t)>2,] ,
             error ="continue"),      # Choose filter for DMPs


  tar_target(save_dmps, data.table::fwrite(
    dmps_f,paste0(fparams$dmp_folder,as.character(quote(method.function)),"dmps.txt")),
    error ="continue",deployment = "main"),                  # Save DMPs
  # tar_target(betas_DIF,betasdmps(betas,dmps_f,rgSet),error="continue",priority = 0.1),                           #
  # tar_target(betas_DIF_full,betasdmps(betas,ann,rgSet),error="continue",priority = 0.1),                         #
  # tar_target(write_betas_DIF, writexl::write_xlsx(betas_DIF,paste0(fparams$dmp_folder,"betas_DIF_",data_names,".xlsx"))),
  # tar_target(write_betas_DIF_full, writexl::write_xlsx(betas_DIF_full,paste0(fparams$dmp_folder,"betas_DIF_full_",data_names,".xlsx"))),


  tar_target(dmps_summary,                                                       # Summary statistics for DMPs
             summary_dmps(dmpsod1, dir = fparams$dmp_folder,name=data_names),error ="continue"),
  tar_target(dmpplotod1, plotDMP(dmps_f,path=fparams[["dmpplots_folder"]]),error ="continue"),   # Barplots hipo/hyper, genomic region, CpG islands.

  tar_target(dmrs1,                                                              # Finds DMRs with dmrcate can be relaxed here and filter by HMFDR later
             find_dmrs(betas,model,
                       fdr = 0.01,betacutoff = 0.05, min.cpg=3),
             deployment = "worker"),
  tar_target(dmrs_battery,  priority = 1, error ="continue",                                       # DMRs distribution along params
             apply_filter_dmrs(
               dmrs = dmrs1,path=paste0(fparams$dmrs_folder,data_names))),
  tar_target(dmrs, filter_dmrs(dmrs1,p.value = "FDR", mDiff = 0.05, min.cpg=3)),   # Filters DMRs, p.value=HMFDR
  tar_target(save_dmrs,                                                          # Saves DMRs
             writexl::write_xlsx(
               dmrs, paste0(fparams$dmrs_folder,"_",data_names,".xlsx"))),
  tar_target(top4dmrs, {top4<-dmrs[,head(.SD,4),by=Contrast,
                                   .SDcols=c("Contrast","seqnames","start","no.cpgs","HMFDR","maxdiff","meandiff","overlapping.genes")]

  data.table::fwrite(top4,paste0(fparams$dmrs_folder,"_top4",data_names,".csv"))
  }),
  tar_target(dmrs_summary,                                                       # Summary stats for DMRs
             summary_dmrs(
               dmrs,path=paste0(fparams$dmrs_folder,"full_dmrs_summary",data_names,".csv")),
             error = "continue"),
  tar_target(vennDiags,  error ="continue",                                      # venn plot for genes comparing groupvar(contrasts?)
             venns(
               dmrs,groupvar="Contrast",res=paste0(fparams$dmrs_folder,"/VENNS/"))),

  tar_target(hyperpathways,                                                      # Pathways hyper
             gopath(dmrs[meandiff>0,],
                    all.cpg=rownames(betas),n=Inf,ann=ann,
                    savepath=paste0(fparams$pathway_folder,"pathways_hyper.csv")),
             resources = tar_resources(
               future = tar_resources_future(
                 # plan = tweak(
                 #   batchtools_slurm,
                 #   template = "batchtools.slurm.tmpl",
                 #   resources = list(num_cores = ncores)
                 # )
                 resources=list(n_cores=ceiling(NCORES/4))
               )
             )
             ),
  tar_target(hypopathways,                                                       # Pathways with hypo
             gopath(dmrs[meandiff<0,],
                    all.cpg=rownames(betas),n=Inf,ann=ann,
                    savepath=paste0(fparams$pathway_folder,"pathways_hypo.csv")),
             resources = tar_resources(
               future = tar_resources_future(
                 # plan = tweak(
                 #   batchtools_slurm,
                 #   template = "batchtools.slurm.tmpl",
                 #   resources = list(num_cores = ncores)
                 # )
                 resources=list(n_cores=ceiling(NCORES/4))
               )
             )
             ),
  tar_target(results_gopathhypo,
             path_results(pathway=hypopathways,topN=50,group="method",pval=0.05,path=paste0(fparams$pathway_folder,data_names,"_hypoPathways.csv"))),
  tar_target(results_gopathhyper,
             path_results(pathway=hyperpathways,topN=50,group="method",pval=0.05,path=paste0(fparams$pathway_folder,data_names,"_hyperPathways.csv"))),

  tar_target(save_gopathhypo,
             writexl::write_xlsx(hypopathways[FDR<0.05,],
                                 paste0(fparams$pathway_folder,"hypoPathways.xlsx"))),
  tar_target(save_gopathhyper,
             data.table::fwrite(hyperpathways[FDR<0.05,],
                                paste0(fparams$pathway_folder,"hyperPathways.csv"))),
  NULL
  )

#
#   ## -- Discard females:
#   tar_target(ss_clean_M, ss_clean2[predictedSex=="M",]),
#   tar_target(betas_M,betas[,ss_clean_M$barcode]),
#   tar_target(rgSet_M,{
#     library(minfi)
#     library(Biostrings)
#     rg<-rgSet[,ss_clean$predictedSex=="M"] #rgSet[,ss_clean_M$barcode]
#     rg@colData<-as(ss_clean_M,"DataFrame")
#     rg
#   }),
#   # Repeat pca:
#   tar_target(top_M,top_beta(betas_M,n=top_betas_N),pattern=map(top_betas_N)),                                        # subsample for PCA & ploting
#   tar_target(pca_M, pca_res(top_M),pattern=map(top_M)),                                                 # PCA
#   tar_target(bplots_M, bplot(pca_M,                                                  # Bi plots for PCA components
#                            ss=ss_clean_M,
#                            colgroup=c(category[V1 %in% c("covs","batch"),rn],"predictedSex","sentrixID"),
#                            s="Condition",
#                            cols=NULL,
#                            folder = paste0(fparams$bplots_folder,"/MALES/",NROW(pca_M$rotation))),
#              # packages = c("ggfortify","ggrepel","gplots","ggplot2"),
#              pattern=map(pca_M)
#   ),
#
#   ## -- (added _M to rgSet betas & ss_clean from here to end)
#   tar_target(model, mod(object = betas_M, group_var = "Sample_Group",              # Model with limma for diff meth
#                         contrasts = NULL,singular=T,
#                         covs= category[V1 %in% c("batch"),rn],
#                         metadata = ss_clean_M)
#   ),
#
#   tar_target(dmps_mod1, cnv.methyl::DMPextr(fit = model,                         # Toptable & stats
#                                             ContrastsDM = colnames(model$contrasts),
#                                             beta_normalized = betas_M,
#                                             p.value = 0.95,
#                                             mDiff = 0.01,
#                                             ann = ann,
#                                             writeOut = F
#   )),
#   tar_target(dmp_battery,priority = 1,                                           # DMPs distribution along params.
#              apply_filter_dmps(
#                dmps = dmps_mod1,path=paste0(fparams$dmp_folder,data_names)),error ="continue"),
#   tar_target(dmps_f , filter_dmps(dmps_mod1, p.value = 0.05, mDiff = 0.2),error ="continue"),      # Choose filter for DMPs
#   tar_target(save_dmps, data.table::fwrite(
#     dmps_f,paste0(fparams$dmp_folder,data_names,"_dmps.txt")),error ="continue"),                  # Save DMPs
#   tar_target(betas_DIF,betasdmps(betas_M,dmps_f,rgSet),error="continue"),                           #
#   tar_target(betas_DIF_full,betasdmps(betas_M,ann,rgSet),error="continue"),                         #
#   # tar_target(write_betas_DIF, writexl::write_xlsx(betas_DIF,paste0(fparams$dmp_folder,"betas_DIF_",data_names,".xlsx"))),
#   # tar_target(write_betas_DIF_full, writexl::write_xlsx(betas_DIF_full,paste0(fparams$dmp_folder,"betas_DIF_full_",data_names,".xlsx"))),
#
#
#   tar_target(dmps_summary,                                                       # Summary statistics for DMPs
#              summary_dmps(dmps_mod1, dir = fparams$dmp_folder,name=data_names),error ="continue"),
#   tar_target(dmpplot_mod1, plotDMP(dmps_f,path=fparams[["dmpplots_folder"]]),error ="continue"),   # Barplots hipo/hyper, genomic region, CpG islands.
#
#   tar_target(dmrs1,                                                              # Finds DMRs with dmrcate
#              find_dmrs(betas_M,model,
#                        fdr = 0.25, p.value = 0.1,betacutoff = 0.1, min.cpg=3)),
#   tar_target(dmrs_battery,  priority = 1, error ="continue",                                       # DMRs distribution along params
#              apply_filter_dmrs(
#                dmrs = dmrs1,path=paste0(fparams$dmrs_folder,data_names))),
#   tar_target(dmrs, filter_dmrs(dmrs1,p.value = 0.1, mDiff = 0.25, min.cpg=3)),   # Filters DMRs
#   tar_target(save_dmrs,                                                          # Saves DMRs
#              writexl::write_xlsx(
#                dmrs, paste0(fparams$dmrs_folder,"_",data_names,".xlsx"))),
#   tar_target(top4dmrs, {top4<-dmrs[,head(.SD,4),by=Contrast,
#                                    .SDcols=c("Contrast","seqnames","start","no.cpgs","HMFDR","maxdiff","meandiff","overlapping.genes")]
#
#   data.table::fwrite(top4,paste0(fparams$dmrs_folder,"_top4",data_names,".csv"))
#   }),
#   tar_target(dmrs_summary,                                                       # Summary stats for DMRs
#              summary_dmrs(
#                dmrs,path=paste0(fparams$dmrs_folder,"full_dmrs_summary",data_names,".csv"))),
#   tar_target(vennDiags,  error ="continue",                                      # venn plot for genes comparing groupvar(contrasts?)
#              venns(
#                dmrs,groupvar="Contrast",res=paste0(fparams$dmrs_folder,"/VENNS/"))),
#
#   tar_target(allpathways, gopath(dmrs,all.cpg=rownames(betas_M),n=Inf,ann=ann)),   # Pathways with all
#   tar_target(hyperpathways,                                                      # Pathways hyper
#              gopath(dmrs[meandiff>0,],
#                     all.cpg=rownames(betas_M),n=Inf,ann=ann,savepath=paste0(fparams$pathway_folder,"pathways_hyper.csv"))),
#   tar_target(hypopathways,                                                       # Pathways with hypo
#              gopath(dmrs[meandiff<0,],
#                     all.cpg=rownames(betas_M),n=Inf,ann=ann,savepath=paste0(fparams$pathway_folder,"pathways_hypo.csv"))),
#   tar_target(save_gopathall,
#              writexl::write_xlsx(allpathways[FDR<0.05,],
#                                  paste0(fparams$pathway_folder,data_names,"_allPathways.xlsx"))),
#   tar_target(save_gopathhypo,
#              writexl::write_xlsx(hypopathways[FDR<0.05,],
#                                  paste0(fparams$pathway_folder,data_names,"_hypoPathways.xlsx"))),
#   tar_target(save_gopathhyper,
#              writexl::write_xlsx(hyperpathways[FDR<0.05,],
#                                  paste0(fparams$pathway_folder,data_names,"_hyperPathways.xlsx")))
#   )
# )

list(#method_function_target,
     top_betas_N_target,
     targets)
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
