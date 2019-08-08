#!/bin/bash
#==============================================================================
# HEADER
#==============================================================================
#
#  Description: 
#

#PBS -N sabe_jobs 
#PBS -l nodes=1:ppn=1
#PBS -l mem=24gb
#PBS -l walltime=24:00:00
#PBS -q short
#PBS -j oe
#PBS -o $PBS_JOBNAME.log
#PBS -t 21-22

ROOT_P=/scratch/genevol/sabe/

VCF_P=vcfs_SABE1172/SABE1172.1a1172.CONCAT.hg38.VQSR.PASS.vSR.Biallelic.chr
VCF_S=.recodecleanFBFD.anno.vcf.gz
OUT_P=/scratch/jonatas/admix_diversity_data/

mkdir -p $OUT_P/log

CHR=$PBS_ARRAYID
#CHR=21
POP=BRL
NS=77
POP_SHUFLE=$OUT_P$POP$NS/$POP"_samples.txt"

for CHR in {1..20}
do
$HOME/bin/bcftools view  \
--min-ac 1 \
-m2 -M2  \
--samples-file $POP_SHUFLE \
--output-file $OUT_P$POP$NS/$POP"_chr"$CHR.vcf.gz \
--output-type z \
$ROOT_P$VCF_P$CHR$VCF_S

tabix $OUT_P$POP$NS/$POP"_chr"$CHR.vcf.gz

plink --freq \
  --out $OUT_P$POP$NS/$POP"_chr"$CHR \
  --vcf $OUT_P$POP$NS/$POP"_chr"$CHR.vcf.gz  

plink --hardy \
  --out $OUT_P$POP$NS/$POP"_chr"$CHR \
  --vcf $OUT_P$POP$NS/$POP"_chr"$CHR.vcf.gz  

plink --make-bed \
  --out $OUT_P$POP$NS/$POP"_chr"$CHR \
  --vcf $OUT_P$POP$NS/$POP"_chr"$CHR.vcf.gz  
done
