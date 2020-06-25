#!/bin/bash

############################################################ 
#    			     Genotype GVCFs    					   #
#                                                          #
#         cankutcubuk [at] {gmail} [dot] {com}             #
#     cankut.cubuk [at] {icr} [dot] {ac} [dot] {uk}        #
#                     2019-2021                            #
#                    @ ICR London                          #
############################################################

cohortID=$1
outputDir=$2
referenceGenome=$3

module load java/sun8/1.8.0u66
module load gatk/4.1.0.0
module load picard-tools/2.8.1


dbDir=$outputDir'/germline/GenomicsDB/'
outputDir=$outputDir'/germline/'

gatk GenotypeGVCFs \
-R $referenceGenome \
-V gendb://$dbDir \
-O $outputDir/$cohortID'.unsorted.vcf.gz'

# GenotypeGVCF generates unsorted VCFs. Sort them with picard

java -Xms10G -Xmx10G  -Djava.io.tmpdir='/home/ccubuk/gdc/logfiles/picard_tmp' -XX:+UseSerialGC -jar \
/apps/picard-tools/2.8.1/picard.jar SortVcf \
I=$outputDir/$cohortID'.unsorted.vcf.gz' \
O=$outputDir/$cohortID'.vcf.gz'

echo >&2 '
************
*** DONE GenotypeGVCFs ***
************
'
