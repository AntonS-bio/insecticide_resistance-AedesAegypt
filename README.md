# insecticide_resistance-AedesAegypti
The analysis is based on mapping to Aedes aegypti GCF_002204515.2 reference assembly and uses corresponding GFF file. Neither is included in this repository as both are available via NCBI/ENA. The list of samples used is in ./InputData/SampleRegions.txt

## Dependencies
To download WGS libraries we used [enaBrowserTools ](https://github.com/enasequence/enaBrowserTools) and [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) was used for mapping. Additional processing relied on [Samtools](https://www.htslib.org/), [bedtools](https://bedtools.readthedocs.io/en/latest/index.html), [Biopython](https://biopython.org/).

Main python dependecies are pandas, numpy, sklearn, scipy, matplotlib, umap-learn, and hdbscan.

## The order of scripts:

### Generating vcfs and converting them to pandas DataFrame
align_reads.sh -> maps whole genome sequencing library (short reads) to selection of genes. Each gene is a separte entry in FASTA file.
    This file calls:
+ callSNPs.py -> calls the SNPs using GATK
+ relabelVCFs.py -> relabelles coordinates of VCF file (mapped to file were each gene is separate entry) to match the coordinates of the reference genome which allows snpEff to work 
+ mergeVCFs.sh  -> the process so far has created a VCF file per sample. For analysis these need to be merged into multisample VCF using BCFtools
+ calculateCoverage.py -> calculates coverage per position for each sample
+ generateCoverageDF.py -> convert coverage files into single pandas DataFrame

all_vcf_to_cds_only.ipynb -> removes intronic regions from VCF file. This is necessary because majority of the sequences are intronic with extremely high rate of duplicated sequences
    This generates a AnnotatedResistanceV3_exon_only.vcf file

vcfToPandas.py -> convert the vcf from all_vcf_to_cds_only.ipynb to pandas DataFrame

### Analysing results
resistanceAnalysis.ipynb -> contains the SNP and coverage analysis and generates figures from the paper

## Useful files:
AnnotatedResistanceV3_exon_only.vcf.gz - contains all SNPs identified in non-intronic regions. 