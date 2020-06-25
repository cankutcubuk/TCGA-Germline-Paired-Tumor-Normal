#!/bin/bash

genes="gencode=ALK&gencode=APC&gencode=ATM&gencode=BAP1&gencode=BARD1&gencode=BMPR1A&gencode=BRCA1&gencode=BRCA2&gencode=BRIP1&gencode=CDH1&gencode=CDK4&gencode=CDKN2A&gencode=CHEK2&gencode=DICER1&gencode=ERCC3&gencode=FH&gencode=FLCN&gencode=HOXB13&gencode=HRAS&gencode=KIT&gencode=KRAS&gencode=MEN1&gencode=MET&gencode=MITF&gencode=MLH1&gencode=MSH2&gencode=MSH6&gencode=MUTYH&gencode=NBN&gencode=NF1&gencode=NF2&gencode=NRAS&gencode=PALB2&gencode=PDGFRA&gencode=PMS2&gencode=POLE&gencode=PTCH1&gencode=PTEN&gencode=RAD50&gencode=RAD51B&gencode=RAD51C&gencode=RAD51D&gencode=RB1&gencode=RET&gencode=RUNX1&gencode=SDHA&gencode=SDHAF2&gencode=SDHB&gencode=SDHC&gencode=SDHD&gencode=SMAD3&gencode=SMAD4&gencode=SMARCA4&gencode=SMARCB1&gencode=STK11&gencode=SUFU&gencode=TERT&gencode=TGFBR1&gencode=TGFBR2&gencode=TMEM127&gencode=TP53&gencode=TSC1&gencode=TSC2&gencode=VHL&gencode=WT1"

token=$(</home/.../gdc-user-token........txt)
filename=$1

counter=0

while read -r line;
do 
	echo "curl --header \"X-Auth-Token: \$token\" https://api.gdc.cancer.gov/slicing/view/"$line"?"$genes "--output /mnt/scratch/DGE/MOPOPGEN/ccubuk/gdc/bam_splice/allbams/"$line".bam"
	
	if ! (($counter % 5000)); then 
		echo $counter 
	fi

	counter=$((counter+1))


done < $filename

