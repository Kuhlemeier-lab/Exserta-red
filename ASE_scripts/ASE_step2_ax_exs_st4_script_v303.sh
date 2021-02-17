#!/bin/bash
#SBATCH --mail-user=youremail@whatever.edu
#SBATCH --time=72:00:00
#SBATCH --mem-per-cpu=60G
#SBATCH --workdir=/yourdirectory/ASE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --array=1-9
#SBATCH -J "realign around snps"


module add vital-it/7
module add UHTS/Analysis/picard-tools/2.18.11

#this script marks duplicates, splits N Cigar reads, and realigns reads around snps and other features
#doing this for both sets of parental reps and the F1s

#genome file = Petunia_axillaris.v3.0.3.fasta

picard-tools MarkDuplicates I=${SLURM_ARRAY_TASK_ID}.sortedByCoord.rg.out.bam O=${SLURM_ARRAY_TASK_ID}.md.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=${SLURM_ARRAY_TASK_ID}.output.metrics TMP_DIR=/yourdirectory/ASE/temp

java -Djava.io.tmpdir=temp -jar ~/gatk/GenomeAnalysisTK.jar -T SplitNCigarReads -R Petunia_axillaris.v3.0.3.fasta -I ${SLURM_ARRAY_TASK_ID}.md.bam -o ${SLURM_ARRAY_TASK_ID}.sp.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

java -Djava.io.tmpdir=temp -jar ~/gatk/GenomeAnalysisTK.jar -T RealignerTargetCreator -R Petunia_axillaris.v3.0.3.fasta -I ${SLURM_ARRAY_TASK_ID}.sp.bam -o ${SLURM_ARRAY_TASK_ID}.list

java -Djava.io.tmpdir=temp -jar ~/gatk/GenomeAnalysisTK.jar -T IndelRealigner -R Petunia_axillaris.v3.0.3.fasta -I ${SLURM_ARRAY_TASK_ID}.sp.bam -targetIntervals ${SLURM_ARRAY_TASK_ID}.list -o ${SLURM_ARRAY_TASK_ID}.re.bam 

#I usually set 72 hrs choose your own time for size of bamfiles etc.
