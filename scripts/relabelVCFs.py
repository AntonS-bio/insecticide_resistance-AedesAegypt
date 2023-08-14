from os import  listdir
from os.path import  isfile, join, splitext

wd="/mnt/storage5/anton/Mosquitoes/ResistanceGenes/"
sourceVcfDir=wd+"/VCFs/"
outputVcfDir=wd+"/relabelled_VCFs/"

vcfFiles = [f for f in listdir(sourceVcfDir) if isfile(join(sourceVcfDir, f)) and (splitext(f)[-1]==".vcf")]

def splitContig(currentLabel):
    chromosome, region = currentLabel.split(":")
    regionStart, regionEnd=region.split("-")
    return {"chr": chromosome, "start": int(regionStart), "end": int(regionEnd)}

#The VCF chromosomal positions needs to be continous as sorted
for file in vcfFiles:
    vcfLines={}
    with open(outputVcfDir+file,"w") as output: 
        with open(sourceVcfDir+file) as source:
            for line in source:
                if line[0]=="#":
                    if line[0:8]=="##contig":
                        pass #the contigs in VCF header are the regions in the bedfile, these need to change to chromosomes.
                    elif line[0:6]=="#CHROM":
                        output.write("##contig=<ID=NC_035109.1,length=409777670>\n")
                        output.write("##contig=<ID=NC_035107.1,length=310827022>\n")
                        output.write("##contig=<ID=NC_035108.1,length=474425716>\n")
                        output.write("##contig=<ID=NC_035109.1,length=409777670>\n")
                        output.write("##contig=<ID=NC_035159.1,length=16790>\n")
                    else:
                        output.write(line)
                else:
                    values=line.strip().split("\t")
                    new_values=splitContig(values[0])
                    values[0]=new_values["chr"]
                    values[1]=new_values["start"]-1+int(values[1]) #There may be an issue here linked to different indexing in bed format and vcf format. 
                    #vcf starts at 1, bed starts at 0. If there is an error, this will be obvious from snpEff resuts because it'll report Ref base mismatches.
                    vcfLines[values[1]]='\t'.join([str(f) for f in values])+"\n"
        for chr in sorted(vcfLines.keys()): #the lines must be sorted by position, which isn't always the case with raw data
            output.write(vcfLines[chr])


# for record in SeqIO.parse("/mnt/storage5/anton/Mosquitoes/GCF_002204515.2_AaegL5.0_genomic.fna", "fasta"):
#     if record.id=="NC_035109.1":
#         for position in referenceNucleotides:
#             print("Reference:"+record.seq[position-1:position]+" vs VCF: "+referenceNucleotides[position])