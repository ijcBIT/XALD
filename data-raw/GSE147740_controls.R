## Download Whole Blood samples.

if (!requireNamespace("BiocManager", quietly = TRUE))
 install.packages("BiocManager")
BiocManager::install("GEOquery")
install.packages(c("data.table","foreach","itertools","doParallel","RcppParallel"))
library(data.table)
setDTthreads(0L)
library(GEOquery)
library(foreach)
library(itertools)

n=ncores=RcppParallel::defaultNumThreads()
dt<-data.table::fread("data-raw/GSE147740_samples_LPS.csv")
gsms<-dt$id
cdir="/mnt/beegfs/idevillasante/Projects/XALD/data-raw/AMN_PBS_controls"
if(!dir.exists(cdir))dir.create(cdir)

#Fetch all 48 Females Age >= 18
cl<- parallel::makePSOCKcluster(n, outfile="")
parallel::clusterEvalQ(cl,{
  library("GEOquery")
  cdir="/mnt/beegfs/idevillasante/Projects/XALD/data-raw/AMN_PBS_controls"
})
doParallel::registerDoParallel(cl)

res<-foreach::foreach(k=isplitIndices(length(gsms),chunks=ncores),
                      .combine='c',
                      .multicombine = F,
                      .inorder=F,
                      .packages = "GEOquery",
                      # .export="cdir",
                      .errorhandling = "pass"
                      )%dopar%{

                        sapply(as.list(k),function(i)GEOquery::getGEOSuppFiles(
                          gsms[i],fetch_files = T,baseDir = cdir))
                      }

parallel::stopCluster(cl)

