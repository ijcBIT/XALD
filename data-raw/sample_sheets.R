require(readxl)
require(data.table)

bp_stopifnot = getFromNamespace("stopifnot", "backports")


# Some paths to find the idats:
iscan_folder<-"/home/idevillasante/shares/BDESTELLER/ISCAN/"
idats_dir<-"/mnt/beegfs/idevillasante/Projects/SG0002/data-raw/idats/"

#Import excel samplesheets from metadata:
Samples_on_Array_last <- readxl::read_excel("/mnt/beegfs/idevillasante/Projects/SG0002/metadata/Samples_on_Array_last.xlsx", 
                                            skip = 7)
dt<-data.table::as.data.table(Samples_on_Array_last)
ss2<-dt[Investigator=="Nuria Climent (IDIBAPs/PCB)"]
ss2$id<-ss2$Sample_Name
ss2$Sample_Name<-NULL
ss1<-data.table::as.data.table(readxl:: read_excel("/mnt/beegfs/idevillasante/Projects/SG0002/metadata/SampleSheet_TAT.xlsx"))
ss1[,id:=sapply(strsplit(Sample_Name," "),"[",2)]
ss1[,Sample_Name:=paste0(sapply(strsplit(Sample_Name," "),"[",1),"_",id)]#gsub()
Sample_sheet<-merge(ss1,ss2, by="id")

paths<-Sample_sheet[,.(path_from=paste0(iscan_folder,Sentrix_ID,"/",barcode),
                path_to=paste0(idats_dir,Sentrix_ID,"/",barcode))]

paths[,sapply(1:.N, function(x){
  dir.create( dirname(paths$path_to[x]),recursive = T)
  system(paste0("scp ijc20250:/",paths$path_from[x],"*.idat"," ",dirname(paths$path_to[x])))
  
  
})
]


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
               Pool_ID=Sample_sheet$Pool_ID,
               Sample_Well =Sample_sheet$Sample_Well,
               Sample_Group=paste0(Sample_sheet$Condition,"_",Sample_sheet$Treatment),
               Sentrix_ID = Sample_sheet$Sentrix_ID,
               Sentrix_Position = Sample_sheet$Sentrix_Position,
               Type = as.factor(Sample_sheet$Treatment),
               Condition = as.factor(Sample_sheet$Condition)
)

# 
# 
# # ids:
ids <- c("Sample_Name",  # Unique identifier for study
         "barcode",      # Unique idats identifiersentrix_ID+Sentrix_position
         "Basename")     # Path to idats/raw files.


batch <- c("Sentrix_ID") # Sentrix chip used

covs <- c(
  "Type",           # Induction: No= no Doxicycline, Dox= induced with Doxi.
  "Condition",      # GENE: Exon1: only exon1, TAT: Whole gene, Control
  NULL)

mgroups <-c(
  "Sample_Group"    # condition_type

)
# # # Pheno: Not aplicable

ss_clean <- ss[,.SD,.SDcols=c(ids,batch,covs,mgroups)]
category <- as.list(c(rep("ids",length(ids)),rep("batch",length(batch)),rep("covs",length(covs)),rep("mgroups",length(mgroups))))
out <- rbindlist(list(category,ss_clean))
saveRDS(out,"data/ss.rds")
# usethis::use_data(out, overwrite = TRUE)
