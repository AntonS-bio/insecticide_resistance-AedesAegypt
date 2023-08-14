
import multiprocessing as mp
import subprocess
from os import chdir, chdir, listdir, mkdir, isdir
from os.path import isfile, splitext, join

wd="/mnt/storage5/anton/Mosquitoes/ResistanceGenes/"
outputDir=wd+"/coverage/"

chdir(wd)

if not isdir(outputDir):
    mkdir(outputDir)

def getresistancecoverage(sampleID):
    subprocess.call(f'samtools depth -a {wd}/bams/{sampleID}.bam > {outputDir}{sampleID}.tsv', shell=True, executable="/bin/bash")

bam_files=[f.replace(".bam","") for f in listdir(wd+"/bams/") if isfile(join(wd+"/bams/", f)) and (splitext(f)[-1]==".bam")]

print(bam_files)

if __name__ == '__main__':

    pool=mp.Pool(5)
    pool.map(getresistancecoverage, bam_files)
    pool.close()
    pool.join()