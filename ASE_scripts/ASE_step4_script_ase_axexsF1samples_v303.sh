#!/bin/bash
#SBATCH --mail-user=andrea.berardi@ips.unibe.ch
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --workdir=/home/ubelix/ips/berardi/RNAseq_axex/mapped_to_ax303/ASE
#SBATCH --time=8:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --array=4-6
#SBATCH --mem-per-cpu=24G
#SBATCH -J "ase read counting on F1"

module add vital-it/7
module add UHTS/Analysis/picard-tools/2.18.11
module add UHTS/Analysis/GenomeAnalysisTK/3.5.0
module add UHTS/Analysis/samtools/1.8


#F1 samples
#counting reads surrounding SNPs with prev generated filters for each sample
GenomeAnalysisTK -T ASEReadCounter -R Petunia_axillaris.v3.0.3.fasta -I ${SLURM_ARRAY_TASK_ID}.re.bam -o ./ASEcounts/${SLURM_ARRAY_TASK_ID}.strict.counts -U ALLOW_N_CIGAR_READS -sites Pexs_st4_biallelic_snps_selected_strict.vcf -mmq 50 -minDepth 20

#NOW do one run with all combined for a sort of "total" count file
GenomeAnalysisTK -T ASEReadCounter -R Petunia_axillaris.v3.0.3.fasta -I 4.re.bam -I 5.re.bam -I 6.re.bam -o ./ASEcounts/F1s.strict.counts -U ALLOW_N_CIGAR_READS -sites Pexs_st4_biallelic_snps_selected_strict.vcf -mmq 50 -minDepth 20


#rosetta stone if needed
#1-3 ax 1-3
#4-6 f1 1-3
#6-9 exs 1-3