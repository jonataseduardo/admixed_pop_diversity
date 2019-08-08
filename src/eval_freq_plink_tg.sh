#!/bin/sh
#==============================================================================
# HEADER
#==============================================================================
#
#  Description: 
#

#PBS -N continental_TG
#PBS -l nodes=1:ppn=1
#PBS -l mem=24gb
#PBS -l walltime=24:00:00
#PBS -q short
#PBS -j oe
#PBS -o pbslog/$PBS_JOBNAME.log
#PBS -t 1-22

ROOT_P=/raid/genevol/1kg/phase3/data/

VCF_P=phase3_chr/ALL.chr
VCF_S=.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

OUT_P=/scratch/jonatas/admix_diversity_data/
POP_P=/raid/genevol/1kg/phase3/data/phase3_pops/pop

CHR=$PBS_ARRAYID
NS=61

for POP in ASW YRI CEU
do
  mkdir -p $OUT_P$POP$NS/
  POP_SHUFLE=$OUT_P$POP$NS/$POP"_samples.txt"
  if [ ! -e "$POP_SHUFLE" ]; then 
    shuf -n $NS $POP_P$POP".txt" > $POP_SHUFLE;
  fi;
done


for POP in ASW YRI CEU
do
  POP_SHUFLE=$OUT_P$POP$NS/$POP"_samples.txt"
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
    --keep-allele-order

  plink --hardy \
    --out $OUT_P$POP$NS/$POP"_chr"$CHR \
    --vcf $OUT_P$POP$NS/$POP"_chr"$CHR.vcf.gz  
    --keep-allele-order
    

  plink --make-bed \
    --out $OUT_P$POP$NS/$POP"_chr"$CHR \
    --vcf $OUT_P$POP$NS/$POP"_chr"$CHR.vcf.gz  
    --keep-allele-order
done
