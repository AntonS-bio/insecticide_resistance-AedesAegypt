import numpy as np
import pandas as pd
import sys


def parseAnnotation(rawAttributes, AltAlleleValue):
    #parses snpEff field in VCF file
    for field in rawAttributes.split(";"):
        if field[0:3]=="ANN":
            annotations=field.replace("ANN=","").split(",")
            for annotation in annotations:
                annotationFields=annotation.split("|")
                if annotationFields[0]==AltAlleleValue:
                    for relevantGene in relevantGeneIDs:
                        if annotationFields[6].find(relevantGene)>-1:
                            if annotationFields[1].find("missense_variant")>-1:
                                return [annotationFields[1], annotationFields[2],annotationFields[10]] #snpEffectOntology, snpEffectImpact
                            else:
                                return [annotationFields[1], annotationFields[2],"NA"] #snpEffectOntology, snpEffectImpact, Amino acid variant
    return ["intergenic_region", "MODIFIER","NA"]

wd="/mnt/storage5/anton/Mosquitoes/ResistanceGenes/"

vcfFileName='/mnt/storage5/anton/Mosquitoes/ResistanceGenes/OutputData/AnnotatedResistanceV3_exon_only.vcf'

#Process bed file used to get subregions for mapping to label each base in the VCF with Gene, Exon/Intron, CDS, Gene Name
positionMetaData={} #key=chromosome, value={key=position, value={"Intron/Exon":, "isCDS":, "GeneID":} }
relevantGeneIDs=set() #these are gene IDs to which reads have been aligned.
with open(f'{wd}/InputData/mappedRegions.bed') as bedFile:
    for line in bedFile:
        values=line.strip().split("\t")
        relevantGeneIDs.add(values[3])
        if values[0] not in positionMetaData:
            positionMetaData[values[0]]={}
        for i in range(int(values[1]), int(values[2])+1):
            positionMetaData[values[0]][i]={"mRNA":False,"exon":False, "CDS":False, "GeneID":values[3]}
##Populate the isExon, isCDS for the positionMetaData
with open("/mnt/storage5/anton/Mosquitoes/GCF_002204515.2_AaegL5.0_genomic.gff") as gffFile:
    for line in gffFile:
        if line[0]!="#":
            values=line.strip().split("\t")
            if values[2] in ["mRNA","exon","CDS"]:
                for relevantGeneID in relevantGeneIDs:
                    if values[8].find(relevantGeneID)>-1: #check if id (ID=rna-XM_021846286.1) is among the target genes:
                        for i in range( int(values[3]),  int(values[4])+1  ):
                            positionMetaData[values[0]][i][values[2]]=True
                        break

##Special case, for COX1, mRNA and exon and CDS are the same
if "NC_035159.1" in positionMetaData:
    for position in positionMetaData["NC_035159.1"]:
        if positionMetaData["NC_035159.1"][position]["CDS"]:
            positionMetaData["NC_035159.1"][position]["mRNA"]=True
            positionMetaData["NC_035159.1"][position]["exon"]=True

# for speed preallocate matrix size
vcfMatrixSize=0
vcf_lines=0
uniquePositions=set()
uniqueAlleles=set()
with open(vcfFileName, "r", newline="\n") as file:
    for line in file:
        vcf_lines+=1
        if line[0:1]!="#":
            values=line.split("\t")
            chr_pos=f'{values[0]}_{values[1]}'
            for allele in values[4].split(","):
                if chr_pos not in uniquePositions:
                    vcfMatrixSize+=3 #Ref, Nocall and ALT
                    uniquePositions.add( chr_pos )
                    uniqueAlleles.add( chr_pos+"_"+allele)
                elif chr_pos+"_"+allele not in uniqueAlleles:
                    vcfMatrixSize+=+1 #ALT
                    uniqueAlleles.add( chr_pos+"_"+allele)

isExon=[]
isCDS=[]
ismRNA=[]
snpType=[]
geneNames=[]
indexVector=[]
alleleSequence=[]
snpPos=[]
aaVariant=[]
chr_pos_with_processed_metadata=set()
chr_pos_alleles_with_processed_metadata={}
snpEffectImpact=[]#High/Moderate/Low/Modifier - output from snpEff
snpEffectOntology=[]#frameshift_variant, synonymous_variant, stop_gained, etc  - output from snpEff
file = open(vcfFileName, "r", newline="\n")
inData=False
aa=[]
processedPositions={} #{key=position number, value=#of alleles at that position}
allele_offsets={}
alt_alleles_at_chr_pos={}
for line in file:
    line_results=[]
    if (line[0:6]=="#CHROM" and 'matrix' not in locals()): ##the second condition only if all regions should form one matrix
        Samples=[] # list of sample names
        line=line.strip().split("\t")
        for i in range(9 , len(line)):
            Samples.append(line[i])
        matrix=np.zeros( (vcfMatrixSize, len(Samples)), dtype=np.byte)
    elif line[0:2]=="##":
        pass
    else:
        inData=True
        line=line.split("\t")
        region=line[0]
        if line[4]!=".": #if site is invariant - skip it.
            alternative_alleles_seqs=line[4].split(",")
            chr_pos=f'{line[0]}_{line[1]}'
            if chr_pos not in alt_alleles_at_chr_pos:
                alt_alleles_at_chr_pos[chr_pos]=set(alternative_alleles_seqs)
                allele_offsets={"0":0,".":1}
            else:
                alt_alleles_at_chr_pos[chr_pos].update(alternative_alleles_seqs)

            for allele_seq in alternative_alleles_seqs:
                if allele_seq not in allele_offsets:
                    allele_offsets[allele_seq]=len(allele_offsets)

            num_all_alleles=len([value for f in alt_alleles_at_chr_pos.values() for value in f])+(len(alt_alleles_at_chr_pos)-1)*2
            #*2 due to REF and NoCall for each chr_pos combination
            #-1 is because the number of elements in alt_alleles_at_chr_pos has already increased by the positions now being added
            num_chr_pos_alleles=len(alt_alleles_at_chr_pos[chr_pos])
            chr_pos_in_matrix=num_all_alleles-num_chr_pos_alleles

            # This fills presence/absence per allele for each sample
            for i in range(9 , len(line)): #loop over all samples for given position line
                line_alleles_nums=line[i].split(":")[0].split("/") #the specific format will vary between files, but is generally 0/0:289,0:289:.:. where / separates alleles
                # ./. means call not made
                for allele_num in set(line_alleles_nums): #the homozyous individuals would have same allele twice
                    if allele_num=="0":
                        #REFERENCE
                        matrix[chr_pos_in_matrix+allele_offsets["0"],i-9]=1  # -alleleOffset because there is only on reference alleleset the 
                    elif allele_num==".": 
                        #MISSING VALUE
                        matrix[chr_pos_in_matrix+allele_offsets["."],i-9]=1  # +1 is due to MISSING being extra allele, but -alleleOffset because there is only on reference alleleset the 
                    else:
                        #ALT
                        allele_seq=alternative_alleles_seqs[int(allele_num)-1] # the value is a number indicating the index of ALT alleles on that line, -1 because 0 is REF
                        matrix[chr_pos_in_matrix+allele_offsets[allele_seq],i-9]=1 # +alleleOffset because ALTs from other line with same position need to be stepped over.


            chr_pos_metadata=positionMetaData[line[0]][int(line[1])]
            if chr_pos not in chr_pos_with_processed_metadata:
                snpEffectOntology+=["REFERENCE", "NoCall"]
                snpEffectImpact+=["REFERENCE", "NoCall"]
                aaVariant+=["REFERENCE", "NoCall"]
                snpType+=["REFERENCE", "NoCall"]

                geneNames+=[chr_pos_metadata["GeneID"],chr_pos_metadata["GeneID"]]
                isExon+=[chr_pos_metadata["exon"],chr_pos_metadata["exon"]]
                isCDS+=[chr_pos_metadata["CDS"],chr_pos_metadata["CDS"]]
                ismRNA+=[chr_pos_metadata["mRNA"],chr_pos_metadata["mRNA"]]

                indexVector+=[region, region]
                snpPos+=[line[1],line[1]]
                alleleSequence+=[line[3],"NoCall"]
                chr_pos_with_processed_metadata.add(chr_pos)
                chr_pos_alleles_with_processed_metadata[chr_pos]=set(["0","."])
            for allele_seq in sorted(allele_offsets.items(), key= lambda x: x[1]):
                if allele_seq[0] not in chr_pos_alleles_with_processed_metadata[chr_pos]:
                    geneNames+=[chr_pos_metadata["GeneID"]]
                    isExon+=[chr_pos_metadata["exon"]]
                    isCDS+=[chr_pos_metadata["CDS"]]
                    ismRNA+=[chr_pos_metadata["mRNA"]]

                    indexVector+=[region]
                    snpPos+=[line[1]]
                    alleleSequence+=[allele_seq[0]]

                    annotationValues=parseAnnotation(line[7], allele_seq[0])
                    snpEffectOntology.append(annotationValues[0])
                    snpEffectImpact.append(annotationValues[1])
                    aaVariant.append(annotationValues[2])
                    snpType.append("ALT")

                    chr_pos_alleles_with_processed_metadata[chr_pos].add(allele_seq[0])


        if len(alt_alleles_at_chr_pos) % 500 == 0:
            print(len(alt_alleles_at_chr_pos)/vcf_lines)
            print(len(indexVector))

print(matrix.shape)
print(len(geneNames))


df=pd.DataFrame(matrix, columns=Samples)
df["GeneticRegion"]=indexVector
df["Pos"]=snpPos    
df['isExon']=isExon 
df['isCDS']=isCDS 
df['ismRNA']=ismRNA
df['Gene']=geneNames
df['SNPType']=snpType
df['AlleleSequence']=alleleSequence
df['snpEffectImpact']=snpEffectImpact
df['snpEffectOntology']=snpEffectOntology
df['AAchange']=aaVariant
df.to_pickle(wd+"OutputData/exonsMatrixV3.pkl")
#np.save(wd+"/exon_matrix", matrix) #will produce API_matrix.npy which can be read via np.load()
file.close()
