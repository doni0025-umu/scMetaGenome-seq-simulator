library(tidyverse)

# Table S3A, related to Figure 3. Estimation of inter- and intra-individual varability of the 384 microbial species. ICC is the intraclass correlation coefficient calculated as the ratio between the intra-individual and the total variance. Prevalence of microbial species is calculated for all 300 samples (column H) or for all individuals (in at least one of the 4 samples, column I). Reference intervals for the relative abundance (normal ranges) are calculated as the central 95% of the populations for the 300 samples. 

df_S3A <- read.csv("/husky/douglas/sim/TableS3_VariabilityPatterns_384_species.csv",
                   header = T)
colnames(df_S3A)[1] <- "Organism.Name"

# Dataset gathered from NCBI based on parameters specified in the link.
# https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=201174,1224,976,1239&assembly_level=1:3
df_NCBI_fourbig_scaffold_up <- read.delim("/husky/douglas/sim/scaffold_and_up_429Kgens_NCBI.tsv")

# Merging these two dfs so that all data present in df_S3A that has a RefSeq gets added to the new df.
# Joining adata with metadata and adjusting attributes in adata for UMAP coloring ----
rel_NCBI_asmblys <- merge(df_S3A, df_NCBI_fourbig_scaffold_up) %>% filter(startsWith(Assembly.Paired.Assembly.Accession, "GCF")) %>% filter(ANI.Check.status == "OK")
rel_NCBI_asmblys$Assembly.Level <- factor(rel_NCBI_asmblys$Assembly.Level, 
                                          levels = c("Complete Genome", "Chromosome", "Scaffold")) 

# Filtering so unique species only are present.
spec_in_dataset <- rel_NCBI_asmblys$Organism.Name %>% unique()
final_unique_species_df <- data.frame()
for (strain_name in spec_in_dataset) {
  temp_df <- arrange(filter(rel_NCBI_asmblys, Organism.Name == strain_name), Assembly.Level)
  final_unique_species_df <- rbind(final_unique_species_df, temp_df[1,])
}

final_unique_species_df <- final_unique_species_df %>% filter(Relative.abundance..0.1.>0.0002)

if(length(commandArgs(trailingOnly = T)) != 0) {
  num_of_SPCs <- commandArgs(trailingOnly = T) %>% as.numeric()
} else {
  num_of_SPCs <- "100" %>% as.character()
}


final_unique_species_df$bacteria_per_species <- c(ceiling(10^(final_unique_species_df$log10.relative.abundance)/sum(10^(final_unique_species_df$log10.relative.abundance))*num_of_SPCs))

###### Select my bacteria ######
final_unique_species_df <- final_unique_species_df %>% filter(Organism.Name == "Bacteroides uniformis")
################################


json_string <- "{"
for (i in 1:nrow(final_unique_species_df)) {
  if (i == nrow(final_unique_species_df)) {
    json_string <- paste0(json_string,
                          "\n\t\"",
                          final_unique_species_df[i,]$Assembly.Paired.Assembly.Accession,
                          "\": \"",
                          final_unique_species_df[i,]$bacteria_per_species,
                          "\"")
  } else {
    json_string <- paste0(json_string,
                          "\n\t\"",
                          final_unique_species_df[i,]$Assembly.Paired.Assembly.Accession,
                          "\": \"",
                          final_unique_species_df[i,]$bacteria_per_species,
                          "\",")
  }
  
}
descript_setupjson <- "ish_SPC"
json_string <- paste0(json_string,"\n}")
write(json_string, sprintf("~/scMetaGenome-seq-simulator/active_run_params/%d%s-run-composition.json", 
                           num_of_SPCs, 
                           descript_setupjson))
write(final_unique_species_df$Assembly.Paired.Assembly.Accession,
      "~/scMetaGenome-seq-simulator/active_run_params/accessions.txt",
      ncolumns = 1)
