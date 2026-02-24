
library(tidyverse)
library(ggplot2)

if(length(commandArgs(trailingOnly = T)) == 0) {
  run_identifier <- "cpynm1_B_uniformis_1vNoplas_spectratest_260209__09_52"
} else {
  run_identifier <- commandArgs(trailingOnly = T)
}

mean_cos_theta <- read.csv(sprintf("/husky/douglas/sim/%s/costheta.csv", run_identifier),
                           header = F)
metafile <- read.delim(sprintf("/husky/douglas/sim/%s/meta_and_setup/meta_%s.1.tsv", run_identifier, run_identifier),
                       header = F,
                       sep = "\t")
mean_cos_theta$cell_id <- read.csv(sprintf("/husky/douglas/sim/%s/cellids.csv", run_identifier),
                                   header = T)
mean_cos_theta$cell_id <- factor(mean_cos_theta$cell_id$adata.meta.data.cellids)
mean_cos_theta <- arrange(mean_cos_theta, V1) %>% select(cell_id, everything()) 


frequency_plot <- ggplot(mean_cos_theta, aes(cell_id, V1)) + 
                  geom_point()

ggsave(filename = sprintf("/husky/douglas/sim/%s/freq_plot_%s.jpg", run_identifier, run_identifier),
       height = 7, 
       width = 12,
       plot = frequency_plot,
       quality = 50)


