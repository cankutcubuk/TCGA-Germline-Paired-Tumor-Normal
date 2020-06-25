#!/bin/bash

############################################################ 
#        Tumor only variant calling in 8 steps             #
#                                                          #
#     1. Collect bait-bias and pre-adapter artifacts       #
#     2. Generate pileup summaries on tumor sample         #
#     3. Calculate contamination on tumor sample           #
#     4. Find tumor sample name from BAM                   #
#     5. Run MuTect2                                       #
#     6. Sort VCF with Picard                              #
#     7. Filter variant calls from MuTect                  #
#	  8. Filter variants by orientation bias			   #
# 	  													   #
#	  inspired from:									   #
#	  https://docs.gdc.cancer.gov/Data/					   #
#     Bioinformatics_Pipelines/							   #
#	  DNA_Seq_Variant_Calling_Pipeline/					   #
#	  #tumor-only-variant-calling-workflow				   #                                                  
#														   #
#														   #
#         cankutcubuk [at] {gmail} [dot] {com}             #
#     cankut.cubuk [at] {icr} [dot] {ac} [dot] {uk}        #
#                     2019-2021                            #
#                    @ ICR London                          #
###########################################################

sampleID=$1
outputDir=$2 # the cohort/cancer folder: /../../TCGA_BRCA/
referenceGenome=$3 # Do we need to clean the reference genome by removing unknown contigs and keeping only chr1-22, X, Y?
interval=$4
bamDir=$outputDir'/Tumors/'
outputDir=$outputDir'/somatic/'
sampleID=${sampleID%.*}

mkdir $outputDir
mkdir $outputDir'/artifacts/'
mkdir $outputDir'/piles/'
mkdir $outputDir'/filter/'
mkdir $outputDir'/variation/'

module load java/sun8/1.8.0u66
module load gatk/4.1.0.0
module load samtools/1.5
module load picard-tools/2.8.1

samtools index $bamDir/$sampleID'.bam'

## 1. Collect bait-bias and pre-adapter artifacts

# outputs can also be used to generate OXOG metrics

gatk --java-options "-d64 -XX:+UseSerialGC -Xmx8G" CollectSequencingArtifactMetrics \
-I $bamDir'/'$sampleID'.bam' \
-O $outputDir'/artifacts/'$sampleID'.tumor_artifact' \
--FILE_EXTENSION .txt \
-R $referenceGenome  

echo >&2 '
************
*** DONE CollectSequencingArtifactMetrics ***
************
'

## 2. Generate pileup summaries on tumor sample

# the output of this step used to generate contamination.table
# ../files/af-only-gnomad.hg38.vcf.gz is Germline reference from gnomad
# You can use clean reference genome which only contain  chr1-22 + XYM and not contigs

gatk --java-options "-d64 -XX:+UseSerialGC -Xmx8G" GetPileupSummaries \
-I $bamDir'/'$sampleID'.bam' \
-O $outputDir'/piles/'$sampleID'.table' \
-V /files/af-only-gnomad.hg38.vcf.gz \
-L $interval \
-R $referenceGenome

echo >&2 '
************
*** DONE GetPileupSummaries ***
************
'

## 3. Calculate contamination on tumor sample

# '.table' comes from step 2

gatk --java-options "-d64 -XX:+UseSerialGC -Xmx8G" CalculateContamination \
-I $outputDir'/piles/'$sampleID'.table' \
-O $outputDir'/piles/'$sampleID'.contamination.table'

echo >&2 '
************
*** DONE CalculateContamination ***
************
'

## 4. Find tumor sample name from BAM

gatk --java-options "-d64 -XX:+UseSerialGC -Xmx3G" GetSampleName \
-I $bamDir'/'$sampleID'.bam' \
-O $outputDir'/piles/'$sampleID'.sample_name'

echo >&2 '
************
*** DONE GetSampleName ***
************
'

## 5. Run MuTect2

# if you run tumor sample on chromosome level (25 commands with different intervals), merge all chromosome level VCFs into one after this step.
# '.sample_name' comes from step 4
# The filtering will be done in the later stages using germline variants of matched samples. For that reason, here we do not use the arguments of
# --af-of-alleles-not-in-resource, --germline-resource, -pon

gatk --java-options "-Djava.io.tmpdir=/home/ccubuk/gdc/logfiles/tmp -d64 -jar -Xmx8G -XX:+UseSerialGC" Mutect2 \
-R $referenceGenome \
-L $interval \
-I $bamDir'/'$sampleID'.bam' \
-O $outputDir'/filter/'$sampleID'.vcf' \
-tumor $outputDir'/piles/'$sampleID'.sample_name' 
# --af-of-alleles-not-in-resource 2.5e-06 \
# --germline-resource af-only-gnomad.hg38.vcf.gz iz # Germline reference from gnomad
# -pon gatk4_mutect2_4136_pon.vcf.gz # New panel of normal created by 4136 TCGA curated normal samples, using GATK4

echo >&2 '
************
*** DONE Mutect2 Var Calls ***
************
'

## 6. Sort VCF with Picard

# '$sampleID'.vcf' comes from step 5

java -Xms16G -Xmx16G  -Djava.io.tmpdir='/home/ccubuk/gdc/logfiles/picard_tmp' -XX:+UseSerialGC -jar \
/apps/picard-tools/2.8.1/picard.jar SortVcf \
OUTPUT=$outputDir'/filter/'$sampleID'.mutect2.tumor_only.sorted.vcf.gz' \
I=$outputDir'/filter/'$sampleID'.vcf' \
CREATE_INDEX=true

echo >&2 '
************
*** DONE SortVcf ***
************
'

## 7. Filter variant calls from MuTect

# additionally filter on contamination fractions
# '.mutect2.tumor_only.sorted.vcf.gz' comes from step 6
# '.contamination.table' comes from step 3

gatk --java-options "-d64 -XX:+UseSerialGC -Xmx8G" FilterMutectCalls \
-O $outputDir'/filter/'$sampleID'.mutect2.tumor_only.contFiltered.vcf.gz' \
-V $outputDir'/filter/'$sampleID'.mutect2.tumor_only.sorted.vcf.gz' \
--contamination-table $outputDir'/piles/'$sampleID'.contamination.table' \
-L $interval --stats $outputDir'/filter/'$sampleID'.Mutect2FilteringStats.tsv'

echo >&2 '
************
*** DONE FilterMutectCalls ***
************
'
## 8. Filter variants by orientation bias

#'.tumor_only.gatk4_mutect2.raw_somatic_mutation_finalOutput.vcf.gz' is the final output to be annotated with VEP.
# '/artifacts/tumor_artifact.pre_adapter_detail_metrics.txt' comes from step 1
# .mutect2.tumor_only.contFiltered.vcf.gz' is From step 7

gatk --java-options "-d64 -XX:+UseSerialGC -Xmx8G" FilterByOrientationBias \
-O $outputDir'/variation/'$sampleID'.tumor_only.gatk4_mutect2.raw_somatic_mutation_finalOutput.vcf.gz' \
-P $outputDir'/artifacts/'$sampleID'.tumor_artifact.pre_adapter_detail_metrics.txt' \
-V $outputDir'/filter/'$sampleID'.mutect2.tumor_only.contFiltered.vcf.gz' \
-L $interval -R $referenceGenome -AM 'G/T' -AM 'C/T'


echo >&2 '
************
*** DONE FilterByOrientationBias ***
************
'

echo >&2 '
************
*** FINISHED Tumor Only Somatic Call ***
************
'


