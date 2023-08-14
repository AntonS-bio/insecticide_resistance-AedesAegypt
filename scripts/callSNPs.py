import multiprocessing as mp
import subprocess
from os import chdir, remove, listdir, mkdir, isdir
from os.path import exists, isfile, splitext, join

wd="/mnt/storage5/anton/Mosquitoes/ResistanceGenes/"
vcfDir=wd+"/VCFs/"
bamsDir=wd+"/bams/"
tempDir=wd+"/temp/"
referenceFasta="mappedRegions.fasta"


chdir(wd)
if not isdir(vcfDir):
    mkdir(vcfDir)
if not isdir(tempDir):
    mkdir(tempDir)

subprocess.call(f"samtools faidx {referenceFasta}", shell=True, executable="/bin/bash")
if exists(referenceFasta.replace(".fasta",".dict")):
    remove(referenceFasta.replace(".fasta",".dict"))
subprocess.call(f'gatk CreateSequenceDictionary -R {referenceFasta}', shell=True, executable="/bin/bash")

def runCaller(sampleID):
    bamfile=sampleID+".bam"
    subprocess.call(f"samtools addreplacerg  -r \"ID:{sampleID}\\tSM:{sampleID}\\tPL:Illumina\" -O BAM {bamsDir}/{bamfile} > {tempDir}/{bamfile}", shell=True, executable="/bin/bash")
    subprocess.call(f"samtools index {tempDir}/{bamfile}", shell=True, executable="/bin/bash")
    subprocess.call(f"gatk HaplotypeCaller -I {tempDir}/{bamfile} --dont-use-soft-clipped-bases true  -R {referenceFasta} -O {vcfDir}{sampleID}.vcf", shell=True)
    remove(f"{tempDir}/{bamfile}")

bamFiles = [f.replace(".bam","") for f in listdir(bamsDir) if isfile(join(bamsDir, f)) and (splitext(f)[-1]==".bam")]

existingVCFs=[f.replace(".vcf","") for f in listdir(vcfDir) if isfile(join(vcfDir, f)) and (splitext(f)[-1]==".vcf")]

if __name__ == '__main__':

    pool=mp.Pool(10)
    pool.map(runCaller, bamFiles)
    pool.close()
    pool.join()
