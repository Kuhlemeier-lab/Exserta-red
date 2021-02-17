#!/bin/bash
#SBATCH --mail-user=youremail@whatever.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --workdir=/yourdirectory/ASE
#SBATCH --time=8:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --array=4-6
#SBATCH --mem-per-cpu=24G
#SBATCH -J "ase read counting on F1"

module add vital-it/7
module add UHTS/Analysis/picard-tools/2.18.11
module add UHTS/Analysis/GenomeAnalysisTK/3.5.0
module add UHTS/Analysis/samtools/1.8

### ---- Now that we have the variants called (exserta variants against an axillaris reference) we use GATK to count allele-specific expression in the F1 hybrid samples ---- ###
mkdir ASEcounts

#F1 samples separately
#counting reads surrounding SNPs with prev generated filters for each sample
GenomeAnalysisTK -T ASEReadCounter -R Petunia_axillaris.v3.0.3.fasta -I ${SLURM_ARRAY_TASK_ID}.re.bam -o ./ASEcounts/${SLURM_ARRAY_TASK_ID}.strict.counts -U ALLOW_N_CIGAR_READS -sites Pexs_st4_biallelic_snps_selected_strict.vcf -mmq 50 -minDepth 20

#Now do one run with all F1 samples combined for a "total" count file
GenomeAnalysisTK -T ASEReadCounter -R Petunia_axillaris.v3.0.3.fasta -I 4.re.bam -I 5.re.bam -I 6.re.bam -o ./ASEcounts/F1s.strict.counts -U ALLOW_N_CIGAR_READS -sites Pexs_st4_biallelic_snps_selected_strict.vcf -mmq 50 -minDepth 20


#rosetta stone if needed
#1-3 axillaris 1-3
#4-6 f1 1-3
#6-9 exserta 1-3
