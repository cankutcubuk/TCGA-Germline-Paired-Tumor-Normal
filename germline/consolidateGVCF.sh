#!/bin/bash

############################################################ 
#      				 Consolidate GVCFs                     #
#                                                          #
#         cankutcubuk [at] {gmail} [dot] {com}             #
#     cankut.cubuk [at] {icr} [dot] {ac} [dot] {uk}        #
#                     2019-2021                            #
#                    @ ICR London                          #
############################################################

# https://gatkforums.broadinstitute.org/gatk/discussion/3893/calling-variants-on-cohorts-of-samples-using-the-haplotypecaller-in-gvcf-mode
# All case and control samples

outputDir=$1
interval=$2
numThreads=24

module load java/sun8/1.8.0u66
module load gatk/4.1.0.0

varDir=$outputDir'/germline/'
dbDir=$varDir'/GenomicsDB/'

intervalArguments=''
# intervalArguments="$intervalArguments --intervals chr2 --intervals chr3" 
# intervalArguments="$intervalArguments -L chr1 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18- L chr19 -L chr2 -L chr21 -L chr22 -L chr3 -L chr4 -L chr5 -L chr7 -L chr8 -L chr9"
intervalArguments=$intervalArguments' -L '$interval

gVCFarguments=''

for sampleID in $varDir/*.g.vcf.gz ;do gVCFarguments=$gVCFarguments' -V '$sampleID; done
#for sampleID in $varDir/"Tumors"/*".g.vcf.gz" ;do gVCFarguments=$gVCFarguments' -V '$varDir/$sampleID; done

gatk GenomicsDBImport --batch-size $numThreads $gVCFarguments --genomicsdb-workspace-path $dbDir $intervalArguments

echo >&2 '
************
*** DONE ConsolidateGVCFs ***
************
'
