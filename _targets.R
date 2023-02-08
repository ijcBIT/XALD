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
  packages = c("tibble","foreach","S4Vectors"), # packages that your targets need to run
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

source("R/functions.R") # Calls default function @meth_functions and adds/modifies
                        # them as needed for this project.
# source("R/deps.R") # Makes sure all needed packages are installed.
# source("../meth_functions.R")


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
  sample_Sheet1="data/ss.rds"

)

values <- tibble::tibble( # Use all possible combinations of input settings.
  method_function = rlang::syms(c("noob")),#, "pq","funn","noob_pq","Em")),
  data_paths = samplesheets,
  data_names = c("full"),
  Contrasts = paste(c("Sample_GroupControl_No-Sample_GroupExon1_No","Sample_GroupControl_No-Sample_GroupTAT_No" , "Sample_GroupExon1_No-Sample_GroupTAT_No",
                "Sample_GroupControl_No-Sample_GroupControl_Dox","Sample_GroupExon1_No-Sample_GroupExon1_Dox" , "Sample_GroupTAT_No-Sample_GroupTAT_Dox"),
                collapse = ";")
  # model=c("Sample_Group",)
)


# Pipline:
targets <- tar_map(
  values = values,
  names = data_names, #"data_source", # Select columns from `values` for target names.
  tar_target(fparams,make_results_dirs(subf=data_names, results_folder = results_folder, analysis_folder = analysis_folder)),
  tar_target(samplesheet_path, data_paths, format = "file"),
  tar_target(samplesheet, readRDS(samplesheet_path)),
  tar_target(ss,samplesheet[-1,]),
  tar_target(category,samplesheet[1,]),
  tar_target(nrgSet, minfi::read.metharray.exp(base = NULL, targets = ss, extended = T,
                                       recursive = FALSE, verbose = FALSE, force = T)),
  tar_target(rgSet,name_rgset(nrgSet,ss)),
  # Qc report:
  tar_target(QC_plots, qc(rgSet,sampGroups = "Sample_Group",sampNames="barcode", qc_folder = fparams[["qc_folder"]])),

  # Calculate purity:
  tar_target(purity, cnv.methyl::purify(myLoad=rgSet)),
  tar_target(filtered, filter(targets=ss, rgSet=rgSet,sampGroups="Sample_Group", qc_folder = fparams[["qc_folder"]])),

  tar_target(normalize, method_function(filtered)),
  tar_target(clean, prep(normalize)),
  tar_target(ss_clean, droplevels.data.frame( cbind(clean@colData,purity))),
  tar_target(save_ss_clean,write.table(ss_clean,paste0(fparams[["ss_clean_path"]],"/","ss_clean.csv"),quote = F,sep = ",") ),

  tar_target(plotvars, c(data.table::last(colnames(ss_clean)),"predictedSex",colnames(as.matrix(category[1,]))[as.matrix(category[1,])%in%c("batch","covs")])),

  tar_target(ann, minfi::getAnnotation(clean)),
  tar_target(betas, minfi::getBeta(clean),memory = "persistent"),
  tar_target(top,top_beta(betas,n=1000)),
  tar_target(pca, pca_res(top)),

  tar_target(pca_corrplot,corpca(beta_top100 = top,
                                 metadata=ss_clean[,plotvars],
                                 title=paste0("PC1-6 correlations with ",data_names," clinical vars"))
  ),
  tar_target(save_pca_corrplot,save_plot(object=pca_corrplot,path=fparams[["corrplot_folder"]],filename=paste0(data_names,"_pca_corrplot.png"))),

  tar_target(bplots, bplot(pca,
                           ss=ss_clean,
                           colgroup=colnames(as.matrix(category[1,]))[as.matrix(category[1,])%in%c("batch","covs")],
                           s="Type",
                           cols=NULL,
                           folder = fparams$bplots_folder)
  ),

  tar_target(model, mod(object = betas, group_var = "Sample_Group", contrasts = Contrasts,
                        covs= colnames(as.matrix(category[1,]))[as.matrix(category[1,])%in%c("batch")],
                        metadata = ss_clean)
  ),

  tar_target(dmps_mod1, cnv.methyl::DMPextr(fit = model,
                                            ContrastsDM = colnames(model$contrasts),
                                            beta_normalized = betas,
                                            p.value = 0.95,
                                            mDiff = 0.01,
                                            ann = ann,
                                            writeOut = F
  )),
  tar_target(dmp_battery,apply_filter_dmps(dmps = dmps_mod1,path=paste0(fparams$dmp_folder,data_names))),
  tar_target(dmps_f , filter_dmps(dmps_mod1, p.value = 0.01, mDiff = 0.4)),
  # tar_target(save_dmps, writexl::write_xlsx(dmps_mod1,paste0(fparams$dmp_folder,data_names,"_dmps.xlsx"))),
  tar_target(save_dmps, data.table::fwrite(dmps_f,paste0(fparams$dmp_folder,data_names,"_dmps.txt"))),
  tar_target(betas_DIF,betasdmps(betas,dmps_f,rgSet)),
  tar_target(betas_DIF_full,betasdmps(betas,ann,rgSet)),
  # tar_target(write_betas_DIF, writexl::write_xlsx(betas_DIF,paste0(fparams$dmp_folder,"betas_DIF_",data_names,".xlsx"))),
  # tar_target(write_betas_DIF_full, writexl::write_xlsx(betas_DIF_full,paste0(fparams$dmp_folder,"betas_DIF_full_",data_names,".xlsx"))),


  tar_target(dmps_summary,summary_dmps(dmps_mod1, dir = fparams$dmp_folder,name=data_names)),
  tar_target(dmpplot_mod1, plotDMP(dmps_f,path=fparams[["dmpplots_folder"]])),

  # tar_target(dmrs, find_dmrs(betas,model,pcutoff = 0.01, betacutoff = 0.25, min.cpg=5)),
  tar_target(dmrs1, find_dmrs(betas,model,fdr = 0.25, p.value = 0.1,betacutoff = 0.1, min.cpg=3)),
  tar_target(dmrs_battery,apply_filter_dmrs(dmrs = dmrs1,path=paste0(fparams$dmrs_folder,data_names))),
  tar_target(dmrs_par_plot,apply_filter_dmrs),
  tar_target(dmrs, filter_dmrs(dmrs1,p.value = 0.1, mDiff = 0.25, min.cpg=3)),
  
  tar_target(save_dmrs,
             writexl::write_xlsx(as.data.frame(dmrs),
                                 paste0(fparams$dmrs_folder,"_",data_names,".xlsx"))),
  tar_target(top4dmrs, {top4<-dmrs[,head(.SD,4),by=Contrast,
                                   .SDcols=c("Contrast","seqnames","start","no.cpgs","HMFDR","maxdiff","meandiff","overlapping.genes")]
  
    data.table::fwrite(top4,paste0(fparams$dmrs_folder,"_top4",data_names,".csv"))
    }),
  tar_target(vennDiags,venns(dmrs,groupvar="Contrast",res=paste0(fparams$dmrs_folder,"/VENNS/"))),
  tar_target(dmrs_summary, summary_dmrs(dmrs,path=paste0(fparams$dmrs_folder,"full_dmrs_summary",data_names,".csv"))),
  tar_target(allpathways, gopath(dmrs,all.cpg=rownames(betas),n=Inf,ann=ann)),
  tar_target(hyperpathways, gopath(dmrs[meandiff>0,],all.cpg=rownames(betas),n=Inf,ann=ann,savepath=paste0(fparams$pathway_folder,"pathways_hyper.csv"))),
  tar_target(hypopathways, gopath(dmrs[meandiff<0,],all.cpg=rownames(betas),n=Inf,ann=ann,savepath=paste0(fparams$pathway_folder,"pathways_hypo.csv"))),
  #
  # tar_target(save_gopath,
  #            writexl::write_xlsx(pathways[FDR<1,],
  #                                paste0(fparams$pathway_folder,data_names,"_Pathways.xlsx"))),
  NULL
)
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
