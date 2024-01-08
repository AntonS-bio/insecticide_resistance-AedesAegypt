
import multiprocessing as mp
import subprocess
from os import chdir, listdir, mkdir
from os.path import isdir, isfile, join, splitext
from sys import argv

wd=argv[1]
bam_dir=wd
outputDir=wd+"/coverage/"

chdir(wd)

if not isdir(outputDir):
    mkdir(outputDir)

def getresistancecoverage(sampleID):
    subprocess.call(f'samtools depth -a {bam_dir}/{sampleID}.bam > {outputDir}{sampleID}.tsv', shell=True, executable="/bin/bash")

bam_files=[f.replace(".bam","") for f in listdir(bam_dir) if isfile(join(bam_dir, f)) and (splitext(f)[-1]==".bam")]

#print(bam_files)

if __name__ == '__main__':

    pool=mp.Pool(5)
    pool.map(getresistancecoverage, bam_files)
    pool.close()
    pool.join()