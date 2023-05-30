require(readxl)
require(data.table)

bp_stopifnot = getFromNamespace("stopifnot", "backports")


# Some paths to find the idats:
iscan_folder<-"/home/idevillasante/shares/BDESTELLER/ISCAN/"
idats_dir<-"/mnt/beegfs/idevillasante/Projects/XALD/data-raw/idats/"
PBS_idats <- "/mnt/beegfs/idevillasante/Projects/XALD/data-raw/AMN_PBS_controls/"
#Import excel samplesheets from metadata:
# Samples_on_Array_last <- readxl::read_excel("/mnt/beegfs/idevillasante/Projects/XALD/metadata/Samples_on_Array_last.xlsx", 
                                            # skip = 7)
Sample_sheet<-data.table::fread("/mnt/beegfs/idevillasante/Projects/XALD/data-raw/Samplesheet_XALD.csv",colClasses =  "character")
paths<-Sample_sheet[,.(path_from=paste0(iscan_folder,Sentrix_ID,"/",Barcode),
                path_to=paste0(idats_dir,Sentrix_ID,"/",Barcode))]

# paths[,sapply(1:.N, function(x){
#   dir.create( dirname(paths$path_to[x]),recursive = T)
#   system(paste0("scp ijc20250:/",paths$path_from[x],"*.idat"," ",dirname(paths$path_to[x])))
#   
#   
# })
# ]


# # Make a standard sample sheet for the samples. 
# ##It must contain the following columns:
# - Sample_Name
# - Basename
# 
# ##Other recommended columns:
# - Project
# - Pool_ID
# - Sample_Plate
# - Sample_Well
# - Sample_Group
# - Sentrix_ID
# - Sentrix_Position
# 
# ##Pheno columns:
# - Gender
# - Type
# - Condition

Basename<-paste0(idats_dir,Sample_sheet$Sentrix_ID,"/",Sample_sheet$Sentrix_ID,"_",Sample_sheet$Sentrix_Position)
bp_stopifnot("Malformed idats path " = Basename == paths$path_to)
bp_stopifnot("Missing idats " = file.exists(paste0(Basename,"_Grn.idat")) & file.exists(paste0(Basename,"_Red.idat")))

# stri_sub_all(Sample_sheet$Sample_Name, stri_locate_all_regex(Sample_sheet$Sample_Name, ' ', omit_no_match=TRUE)) <- '_'
ss<-data.table(Sample_Name=Sample_sheet$Sample_Name,
               Basename=Basename,
               barcode=basename(Basename),
               Sample_Plate=Sample_sheet$Sample_Plate,
               Age=as.numeric(Sample_sheet$Age),
               Sample_Well =Sample_sheet$Sample_Well,
               Sample_Group=paste0(Sample_sheet$Sample_Group),
               Sentrix_ID = Sample_sheet$Sentrix_ID,
               Sentrix_Position = Sample_sheet$Sentrix_Position,
               Type = as.factor(Sample_sheet$Type),
               Condition = as.factor(make.names(Sample_sheet$Condition))
)

# Fix NAs in Age:
ss[,avAge:=round(mean(Age,na.rm=T)),by=Sample_Group]
ss[is.na(Age),Age:=avAge]

 
# # ids:
vars<-list(
  ids = c("Sample_Name",  # Unique identifier for study
          "barcode",      # Unique idats identifiersentrix_ID+Sentrix_position
          "Basename"),     # Path to idats/raw files.
  batch = c("Sentrix_ID","Sentrix_Position","Sample_Well","Sample_Plate"), # Sentrix chip used
  covs = c(
    "Type",           # Induction: No= no Doxicycline, Dox= induced with Doxi.
    "Condition",      # GENE: Exon1: only exon1, TAT: Whole gene, Control
    "Age",
    NULL),
  mgroups =(
    "Sample_Group"    # condition_type
  )
)

if(length(colnames(ss))>length(unlist(vars)))vars$pheno = names(ss)[!names(ss) %in% unlist(vars)] #Empty
category<-unlist(lapply(seq_along(vars), function(y, n, i) {
  rep(n[[i]],length(y[[i]]))
},
y=vars,
n=names(vars)
))

ss_clean <- ss[,.SD,.SDcols=unlist(vars)] # Redundant
setattr(ss_clean,"category",category)
out<-ss_clean
saveRDS(out,"data/SS.rds")


# Controls for AMN from database 
samp_controls <- data.table::fread("/mnt/beegfs/idevillasante/Projects/XALD/data-raw/GSE147740_samples_LPS.csv")
samp_controls<-samp_controls[Sex=="M",]
samp_controls[,Basename:=paste0(PBS_idats,id,"/",id,"_",geo_accession)]
sapply(samp_controls$Basename, function(x){
  print(x)
  bp_stopifnot("Missing idats " = file.exists(paste0(x,"_Grn.idat.gz")) & file.exists(paste0(x,"_Red.idat.gz")))
  
})

ss_ctl<- data.table(
  Sample_Name=samp_controls$id,
  barcode=samp_controls$geo_accession,
  Sentrix_ID=sapply(samp_controls$geo_accession,function(x) unlist(strsplit(x,"_"))[1]),
  Basename=samp_controls$Basename,
  Type="Control",
  Condition="LPS_ctl",
  Sample_Group="Control",
  Age=as.numeric(samp_controls$Age)
)

ss_AMN <-rbind(ss_ctl, out[Age>18,.SD,.SDcols=names(ss_ctl)])
ss_AMN[,Sample_Group:=Type]
ss_AMN[Type=="Disease" & Condition=="None",Condition:="AMN"]
setattr(ss_AMN,"category",attributes(out)$category[match(names(ss_AMN),names(out))])
saveRDS(ss_AMN,"data/ss_AMN.rds")

# Model1: No cereb(condition2)
ss_NoCereb<-out[!endsWith(Sample_Group,"Cereb"),]
setattr(ss_NoCereb,"category",attributes(out)$category[match(names(ss_NoCereb),names(out))])

saveRDS(ss_NoCereb,"data/ss_NoCereb.rds")

# Model2: Only Adults
ss_Adults<-out[Age>18,]
setattr(ss_Adults,"category",attributes(out)$category[match(names(ss_Adults),names(out))])
saveRDS(ss_Adults,"data/ss_Adults.rds")

# Model3.1: Balanced design:
# Subsample AMN (Sample_group == Adult_Disease_None):
ss_sub<-rbind(out[Sample_Group!="Adult_Disease_None"],
      out[Sample_Group=="Adult_Disease_None"][sample(.N,12)]
)
setattr(ss_sub,"category",attributes(out)$category[match(names(ss_sub),names(out))])
saveRDS(ss_sub,"data/ss_sub.rds")
# Model3.2:
saveRDS(ss_sub[Condition!="Cereb",],"data/ss_NoCereb_sub.rds")
# Model3.3:
saveRDS(ss_sub[Age>18,],"data/ss_Adults_sub.rds")
# usethis::use_data(out, overwrite = TRUE)
