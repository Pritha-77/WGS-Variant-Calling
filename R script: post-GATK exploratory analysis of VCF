###############################################################################
# STEP 0: Setup
###############################################################################


setwd("C:/GATK/results")


# Load necessary libraries
library(VariantAnnotation)
library(dplyr)
library(tidyr)

###############################################################################
# STEP 1: Load Funcotator-annotated VCF and extract annotations
###############################################################################

# Load Funcotator-annotated VCF file
vcf_file <- "analysis-ready-snps-filteredGT-functotated.vcf"
vcf <- readVcf(vcf_file, genome="hg38")

# Extract FUNCOTATION field from INFO
funcotator_raw <- info(vcf)$FUNCOTATION

# Convert CharacterList to a standard vector
funcotator_raw <- unlist(funcotator_raw)

# Ensure FUNCOTATION has the same length as VCF
if (length(funcotator_raw) != length(rownames(vcf))) {
  warning("Mismatch detected: VCF has ", length(rownames(vcf)), " variants, FUNCOTATION has ", length(funcotator_raw), " entries.")
  
# Trim or expand FUNCOTATION list to match VCF row count
  funcotator_raw <- funcotator_raw[seq_len(min(length(funcotator_raw), length(rownames(vcf))))]
}

# Split FUNCOTATION (assuming "|" is the delimiter)
funcotator_parsed <- strsplit(funcotator_raw, "\\|")

# Convert to DataFrame
funcotator_df <- do.call(rbind, funcotator_parsed) |> as.data.frame()

# Ensure FUNCOTATION rows match VCF rows
if (nrow(funcotator_df) != length(rownames(vcf))) {
  funcotator_df <- funcotator_df[seq_len(min(nrow(funcotator_df), length(rownames(vcf)))), ]
}

# Assign column names (adjust if needed)
colnames(funcotator_df) <- c("Hugo_Symbol", "Genome_Build", "Chromosome", "Start", "End", 
                             "Region", "Transcript", "Variant_Type", "Reference_Allele", 
                             "Observed_Allele1", "Observed_Allele2", "Protein_Change", 
                             "dbSNP_ID", "Additional_Info")

# Create final dataframe ensuring row count consistency
funcotator_final <- data.frame(
  Variant = rownames(vcf)[seq_len(nrow(funcotator_df))],  
  Gene = funcotator_df$Hugo_Symbol,  
  Variant_Type = funcotator_df$Variant_Type,  
  Protein_Change = funcotator_df$Protein_Change,  
  dbSNP_ID = funcotator_df$dbSNP_ID  
)

# Remove rows with missing Protein Change- Keep only variants with protein changes
funcotator_final <- funcotator_final %>% filter(!is.na(Protein_Change))

# View first few rows
head(funcotator_final)

# Save to CSV
write.csv(funcotator_final, "Funcotator_Extracted_Annotations.csv", row.names = FALSE)




###############################################################################
# STEP 2: Inspect raw VCF structure
###############################################################################

setwd("C:/GATK/results")

#BiocManager::install("VariantAnnotation")


#Load VCF file

library(VariantAnnotation)
vcf <- readVcf("analysis-ready-snps-filteredGT.vcf", "hg38")
vcf

#Overview
header(vcf)

#How many samples are in this VCF?
sampleID <- samples(header(vcf))            # Sample names
sampleID

#What is in the META field?
meta(header(vcf))


#Overview of metadata: Retrieve information in the META field

meta(header(vcf))$fileformat       #File format

meta(header(vcf))$source  #Source used for variant calling


#Retrieve information in the META field
meta(header(vcf))$contig

#Variants information (VRange format)
rd <- rowRanges(vcf)
rd[1:2]

#Retrieving variation information

## Position of the variations

as.vector(seqnames(rd)[1:5])  # Chromosome

start(rd)[1:2]  # Start position

end(rd)[1:2]  # End position

#How to get the reference alleles?

refBase <- ref(vcf)
refBase[1:2]

refBase <- as.character(refBase)
refBase[1:2]

#How to get alternative alleles?
altBase <- alt(vcf)
alt(vcf)[1:2]

altBase <- lapply(altBase, `[[`, 1)    # get the 1st vector of the list
altBase[1:2]

altBase <- unlist(lapply(altBase, as.character))
altBase[1:2]

###############################################################################
# STEP 3: Explore genotype (GT), depth (DP), quality (GQ)
###############################################################################
#INFO section from FIXED field

info(header(vcf))[1:2, ]

info(vcf)[1:2, ]

#GENOTYPE field
geno(header(vcf))[1:2, ]

#Genotype information for each sample
paste0("GT: ", geno(header(vcf))[1, 3])

matGT <- geno(vcf)$GT
matGT[1:2, ]

sampleID <- trimws(sampleID)


#GT Types

tbl <- table(geno(vcf)$GT)
tbl_dat <- as.data.frame(tbl)
tbl


#0/1: heterozygous mutations, one allele is the same as reference sequence
#1/1: homozygous mutations, both alleles are di erent from reference sequence
#1/2: heterozygous mutations, both alleles are di erent from reference sequence


# Plot genotype counts
library(ggplot2)
ggplot(tbl_dat,aes(x=Var1,y=Freq,fill=Var1))+
  geom_bar(stat='Identity')+
  labs(x="",y="Counts",fill="")+
  theme_classic()

#Depth for each sample (DP): Depth distribution

paste0("DP: ", geno(header(vcf))[3, 3])

matDP <- geno(vcf)$DP
matDP[1:2, ]

colnames(matDP) <- trimws(colnames(matDP))


#DP Distribution

summary(as.vector(matDP))

ggplot(as.data.frame(matDP),aes(x=SRR062634))+geom_histogram()+
  labs(x="",y="Counts")+
  scale_x_log10()+
  theme_classic()

#Genotype calling quality (GQ)

paste0("GQ: ", geno(header(vcf))[4, 3])
matGQ <- geno(vcf)$GQ
matGQ[1:2, ]

colnames(matGQ) <- trimws(colnames(matGQ))

#GQ distribution

summary(as.vector(matGQ))

ggplot(as.data.frame(matGQ),aes(x=SRR062634))+geom_histogram()+
  labs(x="",y="Counts")+
  scale_x_log10()+
  theme_classic()


###############################################################################
# STEP 4: Build variant tables for different genotypes
###############################################################################

#Gathering variant information
# Extract heterozygous (0/1, 1/1) and multi-allelic (1/2) variants
var_het_hom <- rownames(matGT)[matGT %in% c("0/1", "1/1")]
var_multi   <- rownames(matGT)[matGT == "1/2"]


##Gathering information ~ GT 0/1 and 1/1
# select variants with GT 0/1 or 1/1

var_1 <- rownames(geno(vcf)$GT)[
  geno(vcf)$GT=="0/1" | 
    geno(vcf)$GT=="1/1"]

#Extract variant information

varTab1 <- data.frame(variant=names(rd)[names(rd) %in% var_1],
                      chr=as.vector(seqnames(rd)[names(rd) %in% var_1]),
                      start=start(rd)[names(rd) %in% var_1],
                      end=end(rd)[names(rd) %in% var_1],
                      stringsAsFactors = FALSE)

#Gathering information ~ GT 0/1 and 1/1
## Ref alleles are retrieved from ref(vcf)

ref_base <- ref(vcf)[rownames(vcf) %in% var_1]
ref_base[1:2]

varTab1$refBase <- as.character(ref_base)


#Gathering information ~ GT 0/1 and 1/1
##Alt alleles are retrieved from alt(vcf)

alt_base <- lapply(alt(vcf)[rownames(vcf) %in% var_1],`[[`,1)
alt_base[1]
varTab1$altBase <- unlist(alt_base)

#Extract counts from AD

adCount <- geno(vcf)$AD[rownames(geno(vcf)$AD) %in% var_1]
adCount[1]

varTab1$refCount <- unlist(lapply(adCount,`[[`,1))
varTab1$altCount <- unlist(lapply(adCount,`[[`,2))


#genoType: genotype (GT)
#gtQuality: genotyping quality (GQ)


varTab1$genoType <- geno(vcf)$GT[rownames(geno(vcf)$GT) %in% var_1]
varTab1$gtQuality <- geno(vcf)$GQ[rownames(geno(vcf)$GQ) %in% var_1]


#Gathering information ~ genotype 1/2

var_2 <- rownames(geno(vcf)$GT)[geno(vcf)$GT=="1/2"]



varTab2 <- data.frame(variant=names(rd)[names(rd) %in% var_2],
                      chr=as.vector(seqnames(rd)[names(rd) %in% var_2]),
                      start=start(rd)[names(rd) %in% var_2],
                      end=end(rd)[names(rd) %in% var_2],
                      refBase=unlist(lapply(lapply(
                        alt(vcf)[rownames(vcf) %in% var_2],`[[`,1),as.character)),
                      altBase=unlist(lapply(lapply(
                        alt(vcf)[rownames(vcf) %in% var_2],`[[`,2),as.character)),
                      refCount=unlist(lapply(
                        geno(vcf)$AD[rownames(geno(vcf)$AD) %in% var_2],`[[`,2)),
                      altCount=unlist(
                        lapply(geno(vcf)$AD[rownames(geno(vcf)$AD) %in% var_2],`[[`,3)),
                      genoType=geno(vcf)$GT[rownames(geno(vcf)$GT) %in% var_2],
                      gtQuality=geno(vcf)$GQ[rownames(geno(vcf)$GQ) %in% var_2],
                      stringsAsFactors = FALSE)


#Merged datatable

varTab <- rbind(varTab1, varTab2)
varTab[1:2, ]

#Differentiate variant types
#SNP: single-nucleotide substitutions
#DEL: deletions
#INS: insertions
#Others: complicated variations, such as Ins/Del or Inversion

library(dplyr)

varTab <- varTab %>%
  mutate(mutType = case_when(
    nchar(refBase) < nchar(altBase) ~ "INS",
    nchar(refBase) > nchar(altBase) ~ "DEL",
    nchar(refBase) == 1 & nchar(altBase) == 1 ~ "SNP",
    TRUE ~ "Others"
  ))

# 
tbl <- table(varTab$mutType)
tbl_dat <- as.data.frame(tbl)
tbl


ggplot(tbl_dat,aes(x=Var1,y=Freq,fill=Var1))+
  geom_bar(stat = 'identity')+
  labs(x="",y="Mutations",fill="")+
  theme_classic()


#______Nucleotide substitution pattern_______#

#only SNPs
#Transition (Ti): purine-to-purine, pyrimidine-to-pyrimidine
#Transversion (Tv): purine-to-pyrimidine, pyrimidine-to-purine

# Transition (Ti)
ti <- c("A>G","G>A","C>T","T>C")
# Transveersion (Tv)
tv <- c("A>T","A>C","G>T","G>C","C>A","C>G","T>A","T>G")
varTab$nuSub <- paste0(varTab$refBase,">",varTab$altBase)
varTab$TiTv[varTab$nuSub %in% ti] <- "Ti"
varTab$TiTv[varTab$nuSub %in% tv] <- "Tv"
varTab[1:2,]

varX <- varTab[varTab$mutType == "SNP", ]
tbl <- table(varX$nuSub)
tbl_dat <- as.data.frame(tbl)
tbl

ggplot(tbl_dat,aes(x=Var1,y=Freq,fill=Var1))+
  geom_bar(stat = 'identity')+
  labs(x="",y="Mutations",fill="")+
  theme(legend.position = "none")


#Ti/Tv
tbl <- table(varX$TiTv)
tbl_dat <- as.data.frame(tbl)
tbl

ggplot(as.data.frame(table(varX$TiTv)),aes(x=Var1,y=Freq,fill=Var1))+
  geom_bar(stat = 'identity')+labs(x="",y="Mutations",fill="")+
  theme(legend.position = "none")


#______Trinucleotide motif analysis______#

#Extract C>T substituion

library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicFeatures)
library(stringr)

rd_idx <- str_split(names(rd), "_", simplify = T)
rd_idx[1:2, ]

rd_sub <- rd[rd_idx[, 2] == "C/T"]
rd_sub[1:2, ]

# Extract sequences beneath the mutation from -1 to +1

rd_sub$triNu <- getSeq(Hsapiens,
                       seqnames(rd_sub),
                       start=start(rd_sub)-1,
                       end=end(rd_sub)+1)
rd_sub[1:2]

#Trinucleotide pattern
tbl <- table(rd_sub$triNu)
tbl_dat <- as.data.frame(tbl)
tbl

ggplot(tbl_dat,aes(x=Var1,y=Freq,fill=Var1))+
  geom_bar(stat='identity')+
  labs(x="",y="Variants",fill="")+
  theme(legend.position = "none")


#________________ APOBEC activation in this sample?_______________#

#APOBEC target: TCW (TCA/TCT)

tbl_dat$APOBEC_target <- tbl_dat$Var1 %in% c("TCA", "TCT")
tbl_dat[1:2]

# collapse by APOBEC_target
apobec_dat <- aggregate(Freq ~ APOBEC_target, tbl_dat, FUN = sum, na.rm = TRUE)
apobec_dat

ggplot(apobec_dat,aes(x=APOBEC_target,y=Freq,fill=APOBEC_target))+
  geom_bar(stat='identity')+
  labs(x="",y="Variants",fill="")+
  theme(legend.position = "none")


###############################################################################
# STEP 5: dbSNP annotation
###############################################################################


library(BSgenome.Hsapiens.UCSC.hg38)  # Genome sequence
library(TxDb.Hsapiens.UCSC.hg38.knownGene)  # Gene annotations

#options(timeout = 5000000)
#BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh38")
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)  # dbSNP locations


#Subset variants in chromosome 1

names(vcf)[1:2]
grepl(names(vcf),pattern = "chr1:")[1:2]

vcf_chr1 <- vcf[grepl(names(vcf),pattern = "chr1:")]
rd_chr1 <- rowRanges(vcf_chr1)

#For all chr
all_chr <- unique(seqnames(vcf))        # Extract chromosome names dynamically
vcf_all <- vcf[seqnames(vcf) %in% all_chr]            # Filter VCF data for all chromosomes
rd_all <- rowRanges(vcf_all)                # Extract rowRanges for all filtered variants


length(vcf_all)  # Number of variants retained
head(rd_all)     # Preview of extracted ranges

#____Annotate rsID from dbSNP____#

# Retrive dbSNP data

all_snps <- SNPlocs.Hsapiens.dbSNP155.GRCh38
all_snps


#Retrieve SNPs by chromosome
tar_chr <- as.vector(seqnames(rd_chr1)@values)
my_snps <- snpsBySeqname(all_snps, c(tar_chr))


#ERROR

## check seqlevels

# seqlevels ~ rd

seqlevels(rd_all)

# seqlevels ~ dbSNP
seqlevels(all_snps)

table(seqnames(rowRanges(vcf)))


#_____seqname in rd_chr1 is chr1 but in all_snps is 1

#check seqlevelStyl
#seqlevelsStyle(rd_chr1)
seqlevelsStyle(rd_all)

seqlevelsStyle(all_snps)

#___________ seqlevelsSytle in rd_chr1 is UCSC, but in all_snps is 


#check genome

#genome(rd_chr1)

genome(rd_all)

genome(all_snps)

#Unify seqnames and Process
tar_chr <- as.vector(seqnames(rd_chr1)@values)
tar_chr <- gsub("chr", "", tar_chr)
tar_chr[grepl(tar_chr, pattern = "M")] <- "MT"
my_snps <- snpsBySeqname(all_snps, c(tar_chr))
my_snps[1:2]



tar_chr <- unique(as.vector(seqnames(rd_all)))
tar_chr <- gsub("chr", "", tar_chr)
tar_chr[grepl(tar_chr, pattern = "M")] <- "MT"
# Keep only standard chromosomes: 1-22, X, Y, MT
tar_chr <- tar_chr[tar_chr %in% c(as.character(1:22), "X", "Y", "MT")]
# Verify the cleaned chromosome list
print(unique(tar_chr))  # Should print only 1-22, X, Y, MT
# Now run snpsBySeqname
my_snps <- snpsBySeqname(all_snps, tar_chr)
my_snps[1:2]



#Convert seqInfo to UCSC style

# change seqlevelsStyle
seqlevelsStyle(my_snps) <- "UCSC"

# change genome
genome(my_snps) <- "hg38"


length(my_snps)
seqlevelsStyle(my_snps)
seqlevelsStyle(rd_chr1)

genome(my_snps)
genome(rd_chr1)

colnames(mcols(my_snps))



#Make rsID table
head(my_snps$RefSNP_id)
head(pos(my_snps))
length(pos(my_snps))

head(seqnames(my_snps))
length(seqnames(my_snps))

snp_ID <- data.frame(
  posIDX = paste0(seqnames(my_snps), ":", pos(my_snps)), 
  rsID = my_snps$RefSNP_id  # Correct column name
)

head(snp_ID)

#Generate Variant table
matV1 <- data.frame(Variant = names(rd_chr1), stringsAsFactors = FALSE)
matV1[1:2, ]


matV1$chromosome <- gsub("(.*):(.*)_(.*)/(.*)", "\\1", matV1$Variant)
matV1$start <- gsub("(.*):(.*)_(.*)/(.*)", "\\2", matV1$Variant)
matV1$end <- gsub("(.*):(.*)_(.*)/(.*)", "\\2", matV1$Variant)
matV1$ref_allele <- gsub("(.*):(.*)_(.*)/(.*)", "\\3", matV1$Variant)
matV1$alt_allele <- gsub("(.*):(.*)_(.*)/(.*)", "\\4", matV1$Variant)
matV1$posIDX <- gsub("(.*)_(.*)", "\\1", matV1$Variant)
matV1[1:2, ]

#Annotation table ~ SNP_ID
matS <- merge(matV1,snp_ID,all.x=TRUE,by="posIDX")
matS <- dplyr::select(matS,-posIDX)
matS[1:2,]

# How many variations in dbSNP
taC2 <- table(!is.na(matS$rsID))
taC2_dat <- as.data.frame(taC2)
taC2

#Variations in dbSNP ~ Plotting
ggplot(taC2_dat,aes(x=Var1,y=Freq,fill=Var1))+
  geom_bar(stat='Identity')+
  labs(x="",y="Counts",fill="in_dbSNP")+
  theme(legend.position = "none")


###############################################################################
# STEP 6: Predict amino acid changes
###############################################################################

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
txdb
#Predict amino acid changes
coding <- predictCoding(vcf_all, txdb, seqSource = Hsapiens)
#Variants and predicted consequence
coding[1]

#Transform into data.frame
matA <- data.frame(Variant=names(coding),
                   chromosome=seqnames(coding),
                   start=start(coding),end=end(coding),
                   ref_allele=as.character(coding$REF),
                   alt_allele=unlist(lapply(lapply(
                     coding$ALT,`[[`,1),as.character)),
                   GeneID=coding$GENEID,
                   TxID=coding$TXID,
                   Protein_posi=unlist(lapply(lapply(
                     coding$PROTEINLOC,`[[`,1),as.integer)),
                   ref_AA=as.character(coding$REFAA),
                   alt_AA=as.character(coding$VARAA),
                   Type=coding$CONSEQUENCE)
matA$aaChange <- paste0("p.",matA$ref_AA,matA$Protein_posi,matA$alt_AA)
matA <- dplyr::select(matA,-Protein_posi,-ref_AA,-alt_AA)

#Annotation table ~ Amino Acid Changes

matA[1:2, ]


# How many variations in coding region

var_in_coding <- data.frame(varName = names(vcf_all), in_coding = names(vcf_all) %in% 
                              matA$Variant, stringsAsFactors = FALSE)
table(var_in_coding$in_coding)


# How many types of mutations in coding region
            
               #nonsense: mutations causing the appearance of stop codon
               #nonsynonymous: mutations causing amino acid changes
               #synonymous: mutations not causing amino acid changes

taC <- table(matA$Type)
taC_dat <- as.data.frame(taC)
taC

# Mutation types in coding region
ggplot(taC_dat,aes(x=Var1,y=Freq,fill=Var1))+
  geom_bar(stat='Identity')+
  labs(x="",y="Counts",fill="")+
  theme(legend.position = "none")


#Integrate SNP and amino acid change into single table
matS$GeneID <- matA$GeneID[match(matS$Variant, matA$Variant)]
matS$AAChange <- matA$GeneID[match(matS$Variant, matA$Variant)]
matS[1:2, ]

