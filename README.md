WGS-Variant-Calling

This repository contains workflows for whole-genome and exome variant calling using GATK4, BWA, Samtools, and R-based post-processing.
Both germline (HaplotypeCaller) and somatic (Mutect2) pipelines are included, following GATK best practices.


🚀 Workflows
   
   🧬 Germline Variant Calling (HaplotypeCaller)
Pipeline steps:
1.	Set up directories
2.	Download reference genome + known sites
3.	FASTQ quality control (FastQC)
4.	Read alignment (BWA-MEM)
5.	Duplicate marking (MarkDuplicatesSpark)
6.	Base Quality Score Recalibration (BQSR)
7.	Collect metrics (alignment, insert size)
8.	Variant calling (HaplotypeCaller)
9.	Hard-filter SNPs/INDELs
10.	Functional annotation (Funcotator)
11.	Extract tabular annotations (VariantsToTable, R scripts)
Output:
•	analysis-ready-snps.vcf
•	analysis-ready-indels.vcf
•	Annotated VCF (*-functotated.vcf)
•	Parsed CSV of gene-level annotations


      🧬 Somatic Variant Calling (Mutect2)
Pipeline steps:
1.	Setup directories
2.	Download reference genome, gnomAD, Panel of Normals (PoN), and exome intervals
3.	FASTQ quality control (FastQC)
4.	Read alignment (BWA-MEM)
5.	Duplicate marking (MarkDuplicatesSpark)
6.	Base Quality Score Recalibration (BQSR)
7.	Somatic variant calling (Mutect2, tumor vs normal)
8.	Contamination estimation (GetPileupSummaries + CalculateContamination)
9.	Orientation bias modeling (LearnReadOrientationModel)
10.	Filter somatic variants (FilterMutectCalls)
11.	Functional annotation (Funcotator)
Output:
•	somatic_raw.vcf.gz
•	somatic_filtered.vcf
•	Annotated VCF (somatic_annotated.vcf)

________________________________________

📊 Downstream R Analysis
The R scripts in postprocessing_analysis.R provide:
•  Extraction of Funcotator annotations into clean CSV tables
•	Exploration of genotype (GT), depth (DP), quality (GQ) distributions
•	Mutation spectrum analysis (SNP, INS, DEL, Ti/Tv ratios)
•	Trinucleotide mutational context & APOBEC signature checks
•	dbSNP annotation and amino acid consequence prediction
Output:
•	Funcotator_Extracted_Annotations.csv
•	Mutation type barplots, Ti/Tv plots, trinucleotide context plots

________________________________________

⚙️ Requirements

**Core tools**
o	GATK 4.2+
o	BWA
o	Samtools
o	FastQC
o	Seqtk (for subsetting reads)

**R packages**
o	VariantAnnotation, dplyr, tidyr, ggplot2
o	BSgenome.Hsapiens.UCSC.hg38, TxDb.Hsapiens.UCSC.hg38.knownGene
o	SNPlocs.Hsapiens.dbSNP155.GRCh38



