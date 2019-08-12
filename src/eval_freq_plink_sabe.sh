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

mkdir -p $OUT_P$POP$NS/

POP_SHUFLE=$OUT_P$POP$NS/$POP"_samples.txt"
if [ ! -e "$POP_SHUFLE" ]; then 
  shuf -n $NS $POP_P$POP".txt" > $POP_SHUFLE;
fi;
  
VCF_OUT=$OUT_P$POP$NS/$POP"_chr"$CHR.vcf.gz
if [ ! -e "$VCF_OUT" ]; then 
  $HOME/bin/bcftools view  \
  --min-ac 1 \
  -m2 -M2  \
  --samples-file $POP_SHUFLE \
  --output-file $VCF_OUT \
  --output-type z \
  $ROOT_P$VCF_P$CHR$VCF_S;
fi;

TBI_OUT=$VCF_OUT.tbi
if [ ! -e "$TBI_OUT" ]; then 
  tabix VCF_OUT;
fi;

FREQ_OUT=$OUT_P$POP$NS/$POP"_chr"$CHR.afreq
if [ ! -e "$FREQ_OUT" ]; then 
/raid/genevol/users/jonatas/bin/plink2 --freq \
  --out $OUT_P$POP$NS/$POP"_chr"$CHR \
  --vcf $OUT_P$POP$NS/$POP"_chr"$CHR.vcf.gz  \
  --keep-autoconv;
fi;
