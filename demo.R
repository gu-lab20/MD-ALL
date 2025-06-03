
#library ----
library(dplyr)
library(stringr)
library(Rphenograph)
library(MDALL)

# For single sample
file_count="/scratch/zuhu/project/rnaseq/out_raw/COH000333_D1.15M/TRANSCRIPTOME/HTSeq/COH000333_D1.15M.DUX4patched.HTSeq"
file_vcf="/scratch/zuhu/project/rnaseq/out_raw/COH000333_D1.15M/TRANSCRIPTOME/Mutation/COH000333_D1.15M.HaplotypeCaller.vcf"
file_fusioncatcher="/scratch/zuhu/project/rnaseq/out_raw/COH000333_D1.15M/TRANSCRIPTOME/FusionCatcher/final-list_candidate-fusion-genes.txt"
file_cicero=""



df_out_testOne=run_one_sample(sample_id = "TestId",file_count = file_count,
                              file_vcf = file_vcf,
                              file_fusioncatcher = file_fusioncatcher,
                              file_cicero = file_cicero,
                              featureN_PG = c(100))

# For multiple samples
setwd("/home/zgu_labs/bin/R/shinyApp/MDALL")
df_listing=read.table("/home/zgu_labs/bin/R/shinyApp/MDALL/test/file_list.tsv",sep  = "\t",header = T)


out_testMul=run_multiple_samples(file_listing = "test/file_list.tsv",featureN_PG = c(100,1058))

out_mul=out_testMul$df_sums



#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# df_out_testOne=run_one_sample(sample_id = "COH000333_D1.15M",file_count = file_count,
#                               file_vcf = file_vcf,
#                               file_fusioncatcher = file_fusioncatcher,
#                               file_cicero = file_cicero,
#                               featureN_PG = c(100))


# out_testMul=run_multiple_samples(file_listing = "test/file_list.tsv",featureN_PG = c(100))
# out_mul=out_testMul$df_sums

