import pandas as pd
from os import listdir, mkdir
from os.path import isfile, isdir
from sys import argv

#start with coverage of features by reads
#run CalculateCoverage.py to get the tsv file below
wd=argv[1]
inputData=wd+"/InputData/"
outputData=wd+"/OutputData/"
coverageData=wd+"/coverage/"

if not isdir(coverageData):
    mkdir(coverageData)

coverageFiles=[f for f in listdir(coverageData) if isfile(coverageData+f) and f.split(".")[-1]=="tsv"]
libraryIDs=list(set([f.split("_")[0] for f in coverageFiles]))
validIDs={}
for id in libraryIDs:
    with open(coverageData+id+".tsv") as sampleFile:
        for line in sampleFile:
            if id not in validIDs:
                validIDs[id]=1
                break
            else:
                validIDs[id]+=1
                break
validIDs=[key for key in validIDs if validIDs[key]==2]

#check which samples with coverage are not in SNPs matrix
pdvcf=pd.read_pickle(inputData+"/exonsMatrixV2.pkl")
missingSamples=[f for f in validIDs if f not in pdvcf]


#getLength of resistance and non-resistance files and pre-generate the dataframe
rowIDs=[]
with open(coverageData+validIDs[0]+".tsv") as coverageFile:
    for line in coverageFile:
        values=line.split("\t")
        rowIDs.append(values[0]+"_"+values[1])
coverageDF=pd.DataFrame(index=rowIDs, columns=validIDs)
for sampleID in validIDs:
    print(sampleID)
    rowCounter=0
    row=0
    with open(coverageData+sampleID+".tsv") as coverageFile:
        for line in coverageFile:
            values=line.split("\t")
            rowID=values[0]+"_"+values[1]
            coverageDF.at[rowID,sampleID]=int(values[2])
            row+=1
    if row!=len(rowIDs):
        print("Rows mismatch("+sampleID+"): "+str(row)+" vs "+str(len(rowIDs)))
    
coverageDF.to_pickle(inputData+"coverageDF.pkl")

geneIDs={}
for line in open(f'{inputData}/mappedRegions.bed'):
    values=line.strip().split("\t")
    geneIDs[f'{values[0]}:{values[1]}-{values[2]}'] = values[3]

geneIDsToIndex={}
for key in geneIDs:
    geneIDsToIndex[geneIDs[key]]=key

coverageDF["Gene"]=[  geneIDs['_'.join(f.split("_")[0:-1])] for f in coverageDF.index ]
coverageDF["Pos"]=[int(f.split("_")[-1]) for f in coverageDF.index]

pdvcf=pd.read_pickle(inputData+"/exonsMatrixV2.pkl")

exonIndices=[geneIDsToIndex[pdvcf.at[f,"Gene"]] +"_" + str(int(pdvcf.at[f,"Pos"])+1-int(geneIDsToIndex[pdvcf.at[f,"Gene"]].split(":")[-1].split("-")[0]))    for f in pdvcf.index if pdvcf.at[f,"isExon"]  ]
cdsIndices=[geneIDsToIndex[pdvcf.at[f,"Gene"]] +"_" + str(int(pdvcf.at[f,"Pos"])+1-int(geneIDsToIndex[pdvcf.at[f,"Gene"]].split(":")[-1].split("-")[0]))    for f in pdvcf.index if pdvcf.at[f,"isCDS"]  ]
for indices, column in zip([exonIndices,cdsIndices], ["isExon","isCDS"]):
    indices=set(indices)
    coverageDF[column]=False
    coverageDF.at[indices,column]=True
coverageDF=coverageDF.loc[coverageDF["isExon"]]
coverageDF.to_pickle(inputData+"coverageDFwithExons.pkl")