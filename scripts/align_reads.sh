#!/bin/bash

wd="MyDir" #set to directory in which the alingments will be stored

cd ${wd}
mkdir tempDir
mkdir bams
#Fasta file consists of the mRNA regions plus 300nt upstream and downstream 
bedtools getfasta -fi GCF_002204515.2_AaegL5.0_genomic.fna -bed ${wd}/InputData/mappedRegions.bed > mappedRegions.fasta

bowtie2-build ${wd}/InputData/mappedRegions.fasta mappedRegions

#the location of ENA tools (if using) needs to be set.
enaToolsDir="ENAtoolsDir"
aspera_ini="location of aspera.ini file. Otherwise the downloads are likely to take a very long time"


while read -a value; do 

        echo ${value}
        sample=${value[0]}
        python ${enaToolsDir}/enaDataGet.py -a -f fastq -as ${aspera_ini}/aspera.ini $value -d tempDir
        bowtie2  --fast-local --no-unal --no-discordant -p 30 -x mappedRegions  -1 ${wd}/tempDir/${sample}/${sample}_1.fastq.gz -2 ${wd}/tempDir/${sample}/${sample}_2.fastq.gz  | samtools view -@ 10 -F 4 -bS | samtools sort -@ 10 -o ./bams/${sample}.bam
        rm -rf tempDir #to avoid running out of storage

        samtools index ./bams/${value[0]}.bam

done < ${wd}"/InputData/SampleRegions.txt"


#call the SNPs using GATK
python callSNPs.py

#for snpEff to work the coordinates of the fasta file use for alignment need to be changed to the coordinates in the original fasta file. 
#alternatively, the snpEff can be run against custom gff (bedtools intersect -a GCF_002204515.2_AaegL5.0_genomic.gff -b ${wd}/InputData/mappedRegions.bed)
#but the second option would make it difficult to determine the nucleotide coordinates of the substituions.
python relabelVCFs.py

#the process so far has created a VCF file per sample. For analysis these need to be merged into multisample VCF
bash mergeVCFs.sh wd

# calculate coverage based on BAM files
python calculateCoverage.py

# Create dataframe for coverage data
python generateCoverageDF.py
