wd=$1

vcfDir=${wd}"/relabelled_VCFs/"

cd ${vcfDir}

ls -l  | grep "vcf$" | awk '{print $9 }' > ${wd}/VCFs.txt

cd ${wd}
mkdir temp #this probably already exists and can be removed after script has finished.
tempDir=${wd}"/temp/"

rm vcf_files.txt #this files is created as samples are processed. If script is ran multiple times (ex. for failed downloads), the file need to be redone.

while read -r file; do 
    prefix=${file/.vcf/}
    echo $prefix

    #index and gzip with bcftools for merge into multifile vcf
    bcftools view -Oz -o ${tempDir}/${prefix}.vcf.gz ${vcfDir}/${prefix}.vcf
	bcftools index ${tempDir}/${prefix}.vcf.gz

    echo ${tempDir}/${prefix}.vcf.gz >> vcf_files.txt

done < ${wd}/VCFs.txt

bcftools merge -l vcf_files.txt  --missing-to-ref -o temp.vcf

# This step requires some fiddling to set-up the GCF_002204515.2 reference genome as the database for snpEff. Refer to snpEff website to see how to do this.
# At one point there was a bug in snpEff (v5 I think) due to which custom databases weren't working.
#java -Xmx8g -jar ~/snpEff/snpEff/snpEff.jar GCF_002204515.2 CrypticAce.vcf > AnnotatedCrypticAce.vcf