#!/bin/bash

############################################################ 
#        Merge VCFs and annotate them with VEP             #
#                                                          #
#         cankutcubuk [at] {gmail} [dot] {com}             #
#     cankut.cubuk [at] {icr} [dot] {ac} [dot] {uk}        #
#                     2019-2021                            #
#                    @ ICR London                          #
############################################################ 

outputDir=$1 # the cohort/cancer folder: /../tidy/TCGA_BRCA/
tool=$2 # picard or vcftools
cancer=$(basename $outputDir)

module load java/sun8/1.8.0u66
module load gatk/4.1.0.0
module load picard-tools/2.8.1
module load vcftools/0.1.14
module load vep/96

## 1. Merge the sample level VCFs into one

# use picard if all the vcf files have the same entries/variants
# use vcftools if the vcf files have different entries across the samples

inputArgs=''

if [ "$tool" = "picard" ]; then

for sampleID in $outputDir/somatic/variation/*.raw_somatic_mutation_finalOutput.vcf.gz ;do inputArgs=$inputArgs' I='$sampleID; done

java -Xms16G -Xmx32G -Djava.io.tmpdir='/home/ccubuk/gdc/logfiles/picard_tmp' -XX:+UseSerialGC -jar \
/apps/picard-tools/2.8.1/picard.jar MergeVcfs \
$inputArgs \
O=$outputDir'/somatic/'$cancer'_tumor_only_somatic_mutations_finalOutput.vcf.gz'
fi

if [ "$tool" = "vcftools" ]; then

for sampleID in $outputDir/somatic/variation/*.raw_somatic_mutation_finalOutput.vcf.gz ;do inputArgs=$inputArgs' '$sampleID; done

vcf-merge $inputArgs | bgzip -c > $outputDir'/somatic/'$cancer'_tumor_only_somatic_mutations_finalOutput.vcf.gz'
fi

echo >&2 '
************
*** DONE single vcf files of this cohort are merged ***
************
'

## 2. Annotate the variants

# the reference file is different than the reference used in the variant calling step.
# But references are GRCh38
# This is detailed in the TCGA_analysis_folder_tree.docx
 
vep --offline -i  $outputDir'/somatic/'$cancer'_tumor_only_somatic_mutations_finalOutput.vcf.gz' \
-o $outputDir'/somatic/'$cancer'_tumor_only_somatic_mutations_finalOutput.annotated.vcf.gz'  \
--fasta /files/Homo_sapiens_assembly38.fasta --assembly GRCh38 \
--compress_output bgzip --everything --format vcf --no_stats --hgvs --vcf --pick --species homo_sapiens \
--dir_cache /home/ccubuk/.vep

echo >&2 '
************
*** DONE VEP annotation of the somatic variants ***
************
'
