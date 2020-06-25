#!/bin/bash

############################################################ 
#       Index bams and call germline SNPs and indels       #
#                                                          #
#         cankutcubuk [at] {gmail} [dot] {com}             #
#     cankut.cubuk [at] {icr} [dot] {ac} [dot] {uk}        #
#                     2019-2021                            #
#                    @ ICR London                          #
############################################################

sampleID=$1 
sampletype=$2
outputDir=$3
referenceGenome=$4
target=$5

bamDir=$outputDir'/'$sampletype/
varDir=$outputDir'/germline/'
mkdir $outputDir'/germline/'

suffix_vcf='.g.vcf.gz'
fname=$(basename $sampleID)
fbname=${fname%.*}

target_argument='-L '$target

module load samtools/1.5
module load java/sun8/1.8.0u66
module load gatk/4.1.0.0

## 1. Index bam file 

samtools index $bamDir/$sampleID

echo >&2 '
************
*** DONE Bam Indexing ***
************
'

## 2. Call germline SNPs and indels

gatk HaplotypeCaller \
-R $referenceGenome \
-ERC GVCF \
$target_argument -I $bamDir/$sampleID \
--output $varDir/$fbname$suffix_vcf

echo >&2 '
************
*** DONE Variant Calling ***
************
'

