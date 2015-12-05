#' # 1-assemble_data
#' 
#' - Scott Funkhouser and Ting Shen -
#' - STT855 LUAD Final Project -
#' 
#' Required package for use of TCGA-Assembler
# install.packages("HGNChelper")

#' Location of TCGA-Assembler on my HD
setwd("/Users/sfunkhouser/Desktop/MSU/Ernst_Lab_Notebook/PROTOCOLS:NOTES/LUAD/")

#' Load TCGA-Assembler modules.
source("./Module_A.r")
source("./Module_B.r")

#' Use of functions within TCGA-Assemlber to download platform specific data for all available samples.
RNASeqRawData <- 
  DownloadRNASeqData(traverseResultFile = "./DirectoryTraverseResult_Jul-08-2014.rda",
                     saveFolderName = "./data/",
                     cancerType = "LUAD",
                     assayPlatform = "RNASeqV2",
                     dataType = "rsem.genes.normalized_results")

CNVRawData <-
  DownloadCNAData(traverseResultFile = "./DirectoryTraverseResult_Jul-08-2014.rda",
                  saveFolderName = "./data/",
                  cancerType = "LUAD")

miRNASeqRawData <- 
  DownloadmiRNASeqData(traverseResultFile = "./DirectoryTraverseResult_Jul-08-2014.rda", 
                       saveFolderName = "./data/",
                       cancerType = "LUAD", 
                       assayPlatform = "miRNASeq")

#' Save raw data
save(miRNASeqRawData, CNVRawData, RNASeqRawData, file = "./data/raw_data.RData")

#' TCGA-Assemlber also includes functions to process data.
#' 
#' *Need to understand what processing is done*
miRNASeqData_ga <- 
  ProcessmiRNASeqData(inputFilePath = "./data/LUAD__bcgsc.ca__illuminaga_mirnaseq__GRCh37__Jul-08-2014.txt", 
                      outputFileName = "LUAD__illuminaga_mirnaseq",
                      outputFileFolder = "./data/processed/")

miRNASeqData_hiseq <-
  ProcessmiRNASeqData(inputFilePath = "./data/LUAD__bcgsc.ca__illuminahiseq_mirnaseq__GRCh37__Jul-08-2014.txt", 
                      outputFileName = "LUAD__illuminahiseq_mirnaseq",
                      outputFileFolder = "./data/processed/")

cnvData_nocnv19 <-
  ProcessCNAData(inputFilePath = "./data/LUAD__broad.mit.edu__genome_wide_snp_6__nocnv_hg19__Jul-08-2014.txt", 
                 outputFileName = "LUAD_cnv_nocnv_19",
                 outputFileFolder = "./data/processed/",
                 refGenomeFile = "./SupportingFiles/Hg19GenePosition.txt")

cnvData_cnv19 <-
  ProcessCNAData(inputFilePath = "./data/LUAD__broad.mit.edu__genome_wide_snp_6__hg19__Jul-08-2014.txt", 
                 outputFileName = "LUAD_cnv_19",
                 outputFileFolder = "./data/processed/",
                 refGenomeFile = "./SupportingFiles/Hg19GenePosition.txt")

mRNAexpression <-
  ProcessRNASeqData(inputFilePath = "./data/LUAD__unc.edu__illuminahiseq_rnaseqv2__rsem.genes.normalized_results__Jul-08-2014.txt", 
                    outputFileName = "LUAD__illuminahiseq_rnaseqv2__GeneExp",
                    outputFileFolder = "./data/processed/",
                    dataType = "GeneExp",
                    verType = "RNASeqV2")

#' Save all processed data.
save(miRNASeqData_ga, miRNASeqData_hiseq, cnvData_nocnv19, mRNAexpression, file = "./data/processed/processed_data.RData")

#' We are mostly interested in working with the CNV and/or mRNAexpression data. Modify these objects
#' so that they are in a more friendly n x p matrix format, with the outcome (tumor sample/normal sample)
#' as the first column and sample IDs as the rownames.
mrna <- mRNAexpression$Data
mrna <- t(mrna)
colnames(mrna) <- mRNAexpression$Des[, "EntrezID"]

cnv <- cnvData_nocnv19$Data
cnv <- t(cnv)
colnames(cnv) <- cnvData_nocnv19$Des[, "GeneSymbol"]

cnv_all <- cnvData_cnv19$Data
cnv_all <- t(cnv_all)
colnames(cnv_all) <- cnvData_cnv19$Des[, "GeneSymbol"]

#' Create "Tumor" column. 1 For tumor, 0 for normal.
tumor <- ifelse(substr(rownames(mrna), 14, 15) >= 10, 0, 1)
mrna <- cbind(tumor, mrna)

tumor <- ifelse(substr(rownames(cnv), 14, 15) >= 10, 0, 1)
cnv <- cbind(tumor, cnv)

tumor <- ifelse(substr(rownames(cnv_all), 14, 15) >= 10, 0, 1)
cnv_all <- cbind(tumor, cnv_all)

save(cnv, mrna, file = "./data/processed/for_analysis/data_for_analysis.RData")

#' Lastly, download and format drug and clinical data for potential use.
DownloadClinicalData(traverseResultFile = "./DirectoryTraverseResult_Jul-08-2014.rda",
                     saveFolderName = "./data/",
                     cancerType = "LUAD",
                     clinicalDataType = c("patient"))

DownloadClinicalData(traverseResultFile = "./DirectoryTraverseResult_Jul-08-2014.rda",
                     saveFolderName = "./data/",
                     cancerType = "LUAD",
                     clinicalDataType = c("drug"))
drug <-
  read.csv("/Users/sfunkhouser/Desktop/MSU/Ernst_Lab_Notebook/PROTOCOLS:NOTES/LUAD/data/nationwidechildrens.org_clinical_drug_luad.csv",
           header = TRUE)

drug <- drug[-1, ]
drug_mod <- drug[, c("bcr_patient_barcode", "drug_name", "therapy_type", "days_to_drug_therapy_start",
                     "days_to_drug_therapy_end", "measure_of_response")]

patient <-
  read.csv("/Users/sfunkhouser/Desktop/MSU/Ernst_Lab_Notebook/PROTOCOLS:NOTES/LUAD/data/nationwidechildrens.org_clinical_patient_luad.csv",
           header = TRUE)

patient <- patient[-1, ]
patient_mod <- patient[, c("bcr_patient_barcode",
                           "gender",
                           "history_other_malignancy",
                           "tobacco_smoking_history_indicator",
                           "tobacco_smoking_pack_years_smoked",
                           "vital_status",
                           "death_days_to",
                           "birth_days_to")]

#' Combine with tumor CNV and mRNA data
# idx <- substr(rownames(cnv), 1, 12) %in% drug_mod[, 1]
# drug_mod_cnv <- cnv[idx, ]
# drug_mod_cnv_tumor <- drug_mod_cnv[drug_mod_cnv[, 1] == 1, ]
# idx <- duplicated(substr(rownames(drug_mod_cnv_tumor), 1, 12))
# drug_mod_cnv_tumor <- drug_mod_cnv_tumor[!idx, ]
# sum(substr(rownames(drug_mod_cnv_tumor), 1, 12) %in% drug_mod[, 1])

#' Combine "death" data with tumor cnv and mrna data
patient_life <- patient_mod[patient_mod$birth_days_to != "[Not Available]", ]
patient_life$birth_days_to <- as.character(patient_life$birth_days_to) %>% 
                                as.numeric()


