#!/bin/bash

############################################################ 
#         Variant hard-filtering and annotation            #
#                                                          #
#         cankutcubuk [at] {gmail} [dot] {com}             #
#     cankut.cubuk [at] {icr} [dot] {ac} [dot] {uk}        #
#                     2019-2021                            #
#                    @ ICR London                          #
############################################################

set -x

echo >&2 '
************
*** RUNNING Hard filters Variant filtering ***
************
'
# https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
# https://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set
# https://pmbio.org/module-04-germline/0004/02/02/Germline_SnvIndel_FilteringAnnotationReview/
# SOR can be set greater than 3. Check the link above

hostname >&2
echo $0 >&2
printf START >&2; uptime >&2
date >&2
start=$(date +'%s')

sampleID=$1 # Multi-sample GVCF /../../TCGA_BRCA/normals_TCGA_BRCA...
outputDir=$2 # cohort folder  /../../TCGA_BRCA/
referenceGenome=$3
referenceGenomeVep=$4
DP_threshold=20.0 # ??? I set something very small
SOR_threshold_snvs=3.0
SOR_threshold_indels=10.0

module load java/sun8/1.8.0u66
module load gatk/4.1.0.0
module load vep/96

varDir=$outputDir'/germline/'

## 1. Apply hard filters

# Extract SNPs from the call set
gatk SelectVariants \
-R $referenceGenome \
-V $varDir/$sampleID'.vcf.gz' \
-select-type SNP \
-O $varDir/$sampleID'.snvs.vcf.gz'

# Apply hard filters for SNPs
gatk VariantFiltration \
-R $referenceGenome \
-V $varDir/$sampleID'.snvs.vcf.gz' \
-filter "QD < 2.0" \
--filter-name "QualByDepth" \
-filter "FS > 60.0" \
--filter-name "FisherStrand" \
-filter "MQ < 40.0" \
--filter-name "RMSMappingQuality" \
-filter "MQRankSum < -12.5" \
--filter-name "MappingQualityRankSumTest" \
-filter "ReadPosRankSum < -8.0" \
--filter-name "ReadPosRankSumTest" \
-filter "SOR > $SOR_threshold_snvs" \
--filter-name "StrandOddsRatio" \
-filter "DP < $DP_threshold" \
--filter-name "LowDepth" \
-O $varDir/$sampleID'.snvs.labeled.vcf.gz'


# Extract Indels from the call set
gatk SelectVariants \
-R $referenceGenome \
-V $varDir/$sampleID'.vcf.gz' \
-select-type INDEL \
-O $varDir/$sampleID'.indels.vcf.gz'

# Apply hard filters for Indels
gatk VariantFiltration \
-R $referenceGenome \
-V $varDir/$sampleID'.indels.vcf.gz' \
-filter "QD < 2.0" \
--filter-name "QualByDepth" \
-filter "FS > 200.0" \
--filter-name "FisherStrand" \
-filter "ReadPosRankSum < -20.0" \
--filter-name "ReadPosRankSumTest" \
-filter "SOR > $SOR_threshold_indels" \
--filter-name "StrandOddsRatio" \
-filter "DP < $DP_threshold" \
--filter-name "LowDepth" \
-O $varDir/$sampleID'.indels.labeled.vcf.gz'

# Combine VCFs
gatk MergeVcfs \
-I $varDir/$sampleID'.snvs.labeled.vcf.gz' -I $varDir/$sampleID'.indels.labeled.vcf.gz' \
-O $varDir/$sampleID.labeled.vcf.gz

echo >&2 '
************
*** DONE Hard filters Variant filtering ***
************
'

## 2. Annotate the variants

# download cache files into /home/ccubuk/.vep
# for the cache files and their sources check /mnt/scratch/DGE/MOPOPGEN/ccubuk/gdc/files/NOTE.txt 
vep --offline -i $varDir/$sampleID.labeled.vcf.gz -o $varDir/$sampleID.labeled.annoted.vcf.gz \
--fasta $referenceGenomeVep --assembly GRCh38 --compress_output bgzip --everything --format vcf \
--no_stats --hgvs --vcf --pick --species homo_sapiens --dir_cache /home/ccubuk/.vep

echo >&2 '
************
*** DONE VEP annotation of the variants ***
************
'

printf END >&2; uptime >&2
elapsedTime=$(($(date +'%s') - $start))
echo "Elapsed time $(($elapsedTime / 60)) minutes"


######################################################################

# https://cnfl.extge.co.uk/gere/research-environment-user-guide/7-analysis-scripts-and-workflows/annotate-variants-with-vep-variant-effect-predictor
# This script functionally annotates variants using Ensembl VEP
# Please see https://www.ensembl.org/info/docs/tools/vep/index.html for a full list of parameters
 
# INPUT:
 
        # 1: A VCF file (in place of 'input.vcf.gz'). The test VCF ('input.vcf.gz') is included in this directory
 
# OUTPUT:
 
        # 1: An annotated VCF file (output.vcf.gz) - the output format can be changed (e.g for text output specift '--tab')
 
# STEPS:
        # 1: Annotate variants using Ensembl VEP with the '--everything flag'
        # 2: Annotate variants in the CADD plugin
        # 3: Annotate variants with the ClinVar custom VCF
 
# IMPORTANT!
 
        # 1: If you are using GRCh37, specify '--assembly GRCh37' and change the following parameters:
                # --fasta /public_data_resources/reference/GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa
                # --plugin CADD,/public_data_resources/CADD/v1.4/GRCh37/whole_genome_SNVs.tsv.gz
                # --custom /public_data_resources/clinvar/20190219/clinvar/vcf_GRCh37/clinvar_20190219.vcf.gz
        # 2: Change the file paths to suit your working environment
        # 3: Do not run this script directly on the headnode, submit it as a job to the HPC
 
#!/bin/bash
 
#module load vep/96
 
#vep --input_file input.vcf.gz \
#--output_file output.vcf.gz \
#--vcf \
#--compress_output bgzip \
#--species homo_sapiens \
#--assembly GRCh38 \
#--offline \
#--cache \
#--dir_cache /tools/apps/vep/96/ensembl-vep/.vep \
#--cache_version 96 \
#--force \
#--no_stats \
#--everything \
#--fasta /public_data_resources/reference/GRCh38/GRCh38Decoy_no_alt.fa \
#--plugin CADD,/public_data_resources/CADD/v1.4/GRCh38/whole_genome_SNVs.tsv.gz \
#--custom /public_data_resources/clinvar/20190219/clinvar/vcf_GRCh38/clinvar_20190219.vcf.gz,ClinVar,vcf,exact,0,CLNDN,CLNDNINCL,CLNDISDB,CLNDISDBINCL,CLNHGVS,CLNREVSTAT,CLNSIG,CLNSIGCONF,CLNSIGINCL,CLNVC,CLNVCSO,CLNVI
