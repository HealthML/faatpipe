library(ukbtools)
dat_rel <- read.table('./relatedness.tsv', header=TRUE) # UK Biobank Resource 531
eid_exome <- read.table('./data/raw/iid_exome.txt', header=FALSE) # IIDs of participants for which exome sequencing is available
eid_exome <- eid_exome$V1
idx_unrel <- ukb_gen_samples_to_remove(dat_rel, ukb_with_data = eid_exome, cutoff = 0.0884)
write.table(data.frame(idx_unrel), col.names=F, row.names=F, file='./iid_related.txt')
