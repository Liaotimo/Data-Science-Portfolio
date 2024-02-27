# Q1: What is the mean expression across all features for sample 5 in the ALL dataset (from the ALL package)?

library(ALL)
data(ALL)
mean(exprs(ALL[,5]))

# 5.6296

# Q2 .
# We will use the biomaRt package to annotate an Affymetrix microarray. We want our results in the hg19 build of the human genome and we therefore need to connect to Ensembl 75 which is the latest release on this genome version. How to connect to older versions of Ensembl is described in the biomaRt package vignette; it can be achived with the command 
# mart <- useMart(host='feb2014.archive.ensembl.org', biomart = "ENSEMBL_MART_ENSEMBL")
library(biomaRt)
mart <- useMart(host="feb2014.archive.ensembl.org", biomart = "ENSEMBL_MART_ENSEMBL")
ensembl<-useDataset("hsapiens_gene_ensembl", mart)
affy_probe_ids <- featureNames(ALL) # values to annotate
annotation_ALL <- getBM(attributes=c("ensembl_gene_id","affy_hg_u95av2"),filters="affy_hg_u95av2",values=affy_probe_ids,mart=ensembl)
summary(annotation_ALL)
# use table function to compute frequency of features (number of occurences in probe set)
sum(table(annotation_ALL[,2])>1)

# Q3 How many probesets (Affymetrix IDs) are annotated with one or more genes on the autosomes (chromosomes 1 to 22).

annotation_ALL <- getBM(attributes=c("chromosome_name","affy_hg_u95av2"),filters=c("affy_hg_u95av2", "chromosome_name"),values=list(affy_probe_ids,c(1:22)),mart=ensembl)

annotation_ALL <- subset(annotation_ALL, chromosome_name %in% 1:22)
sum(table(table(annotation_ALL[,2])))

# Q4 What is the mean value of the Methylation channel across the features for sample “5723646052_R04C01”?

meth <- getMeth(MsetEx)
sample_meth <- meth[, "5723646052_R04C01"]
mean(sample_meth)

# Q5: Access the processed data from NCBI GEO Accession number GSE788, what is the mean expression level of sample GSM9024?
  
library(GEOquery)
geo_lst <- getGEO("GSE788")
geo_series <- geo_lst[[1]]
mean(exprs(geo_series)[,"GSM9024"])


# Q6: We are using the airway dataset. What is the average of the average length across the samples in the experiment?
library(airway)
library(GenomicRanges)
data("airway")
mean(airway$avgLength)

# Q7: We are using the airway dataset from the airway package. The features in this dataset are Ensembl genes.
# Question: What is the number of Ensembl genes which have a count of 1 read or more in sample SRR1039512?
sum(assay(airway)[,"SRR1039512"]>=1)


# Q8: The airway dataset contains more than 64k features. How many of these features overlaps with transcripts on the autosomes (chromosomes 1-22) as represented by the TxDb.Hsapiens.UCSC.hg19.knownGene package?
# Clarification: A feature has to overlap the actual transcript, not the intron of a transcript. So you will need to make sure that the transcript representation does not contain introns.


# find transcripts on the autosomes
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
txdb<-TxDb.Hsapiens.UCSC.hg19.knownGene
exons_hg19 <- exons(txdb)
autosomes <- c(paste("chr", 1:22, sep=""))
autosome_exons <- keepSeqlevels(exons_hg19, autosomes, pruning.mode="coarse")
airway_gene_ranges <- rowRanges(airway)

seqlevelsStyle(autosome_exons) <- "NCBI"
seqlevelsStyle(airway_gene_ranges)<-"NCBI"
sum(countOverlaps(airway_gene_ranges, autosome_exons)>=1)


# Q9:The expression measures of the airway dataset are the number of reads mapping to each feature. In the previous question we have established that many of these features do not overlap autosomal transcripts from the TxDb.Hsapiens.UCSC.hg19.knownGene. But how many reads map to features which overlaps these transcripts?
# 
# Question: For sample SRR1039508, how big a percentage (expressed as a number between 0 and 1) of the total reads in the airway dataset for that sample, are part of a feature which overlaps an autosomal TxDb.Hsapiens.UCSC.hg19.knownGene transcript?

sample_data <- assay(airway)[, "SRR1039508"]
counts <- countOverlaps(airway_gene_ranges, autosome_exons)

# find genes/features that overlap autosomal exons
features_overlapping_autosome_tx <- labels(counts[counts>=1])

# count reads from probes that map to genes that overlap with autosomal exons
sample_data_overlap <- sample_data[labels(sample_data)%in%features_overlapping_autosome_tx]
# get total number of those overlapping reads and divide by total reads for sample
sum(sample_data_overlap)/ sum(sample_data)



# # Question 10
# Consider sample SRR1039508 and only consider features which overlap autosomal transcripts from TxDb.Hsapiens.UCSC.hg19.knownGene. We should be able to very roughly divide these transcripts into expressed and non expressed transcript. Expressed transcripts should be marked by H3K4me3 at their promoter. The airway dataset have assayed “airway smooth muscle cells”. In the Roadmap Epigenomics data set, the E096 is supposed to be “lung”. Obtain the H3K4me3 narrowPeaks from the E096 sample using the AnnotationHub package.
# 
# Question: What is the median number of counts per feature (for sample SRR1039508) containing a H3K4me narrowPeak in their promoter (only features which overlap autosomal transcripts from TxDb.Hsapiens.UCSC.hg19.knownGene are considered)?
#   
#   Clarification: We are using the standard 2.2kb default Bioconductor promotor setting.
# 
# Conclusion Compare this to the median number of counts for features without a H3K4me3 peak. Note that this short analysis has not taken transcript lengths into account and it compares different genomic regions to each other; this is highly suscepticle to bias such as sequence bias.

# get the methylation data
library(AnnotationHub)
ah <- AnnotationHub()
ah_query <- query(ah, c("H3K4me3","narrowPeak", "E096"))
meth_gr <- ah_query[["AH30596"]]
meth_gr

# get genes (feature labels) using same sample_data_overlap from Q9
overlap_genes <- labels(sample_data_overlap)

# get promoters of genes

hg19_genes_ah_query <- query(ah, c("TxDb.Hsapiens.UCSC.hg19.knownGene"))  
hg19_genes_db <- hg19_genes_ah_query[["AH52258"]]
hg19_transcripts <- transcripts(hg19_genes_db)
hg19_promoters <- promoters(hg19_transcripts) # using default 2200 bases promoter
overlap_ranges <- rowRanges(airway)[overlap_genes]

# need to normalize seqlevel style
seqlevelsStyle(overlap_ranges) <- "UCSC"
seqlevelsStyle(hg19_promoters) <- "UCSC"
# need to remove extraneous chromosomes because they have mismatch between comparison sets
hg19_promoters <- keepSeqlevels(hg19_promoters, paste("chr", 1:22, sep=""), pruning.mode="coarse")
overlap_ranges <- keepSeqlevels(overlap_ranges, paste("chr", 1:22, sep=""), pruning.mode="coarse")


overlap_gr <- overlap_ranges[queryHits(findOverlaps(overlap_ranges, hg19_promoters))])
relevant_genes <- unique(names(overlap_gr))

sample_data <- assay(airway)[, "SRR1039508"]
median(sample_data[relevant_genes])
median(sample_data[!labels(sample_data)%in%relevant_genes])

# upon review, could have just found all overlap genes across the genome and then overlapped with probes
