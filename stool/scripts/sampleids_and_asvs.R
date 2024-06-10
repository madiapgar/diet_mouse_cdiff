## 5-24-23

## putting together total sum scaled sample ids and asvs so I can
## match the taxonomic information with my core diversity

## reading in tss meta contrib file 
tss_contrib <- read_tsv(
  "~/gut_microbiome_metabolomics/total_sum_scaled/butyrate_bile_pred_outputs/untrustable_meta_contrib.tsv")

tss_contrib %>% 
  select(sample, taxon) -> tss_sample_asv

names(tss_sample_asv)[names(tss_sample_asv) == 'sample'] <- 'sampleid'
names(tss_sample_asv)[names(tss_sample_asv) == 'taxon'] <- 'asv'

write_tsv(tss_sample_asv, "~/gut_microbiome_metabolomics/total_sum_scaled/tss_sample_asv.tsv" )



