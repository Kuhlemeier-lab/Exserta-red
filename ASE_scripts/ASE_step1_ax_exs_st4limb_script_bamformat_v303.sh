#!/bin/bash
#SBATCH --mail-user=youremail@wherever.edu
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=60G
#SBATCH --workdir=/your/working/directory
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -J "ASE step 1"

module add vital-it/7
module add UHTS/Analysis/picard-tools/2.18.11



#Step 1: organize bam files (make sure they are soorted by coord) and give symbolic links with numbers so you can use SLURM arrays
mkdir ASE
cd ASE

### ---- Making symbolic links for this analysis so that I can use an array command ---- ###
#usage#
#ln -s source_file symboliclinknameforfile

#list of samples (3 parents species 1, 3 parents species 2, 3 F1 hybrids)
#ax_exF1Limb1_Aligned.sortedByCoord.out.bam  ax_Limb3_Aligned.sortedByCoord.out.bam
#ax_exF1Limb2_Aligned.sortedByCoord.out.bam  exs_Limb1_Aligned.sortedByCoord.out.bam
#ax_exF1Limb3_Aligned.sortedByCoord.out.bam  exs_Limb2_Aligned.sortedByCoord.out.bam
#ax_Limb1_Aligned.sortedByCoord.out.bam      exs_Limb3_Aligned.sortedByCoord.out.bam
#ax_Limb2_Aligned.sortedByCoord.out.bam

#F1 group, name them 4-6
#    ln -s ../ax_exF1Limb1_Aligned.sortedByCoord.out.bam 4.sortedByCoord.out.bam
#    ln -s ../ax_exF1Limb2_Aligned.sortedByCoord.out.bam 5.sortedByCoord.out.bam
#    ln -s ../ax_exF1Limb3_Aligned.sortedByCoord.out.bam 6.sortedByCoord.out.bam
        
#axillaris group, name them 1-3 (parent species 1)
#    ln -s ../ax_Limb1_Aligned.sortedByCoord.out.bam 1.sortedByCoord.out.bam
#    ln -s ../ax_Limb2_Aligned.sortedByCoord.out.bam 2.sortedByCoord.out.bam
#    ln -s ../ax_Limb3_Aligned.sortedByCoord.out.bam 3.sortedByCoord.out.bam

#exserta group, name them 7-9 (parent species 2)
#    ln -s ../exs_Limb1_Aligned.sortedByCoord.out.bam 7.sortedByCoord.out.bam
#    ln -s ../exs_Limb2_Aligned.sortedByCoord.out.bam 8.sortedByCoord.out.bam
#    ln -s ../exs_Limb3_Aligned.sortedByCoord.out.bam 9.sortedByCoord.out.bam
    
    

#I also made a symbolic link (or copied if space) to the reference genome fasta file and gff annotation in this directory
#ln -s /directory_w_genomefiles/peaxi162AQ_PeaxHIC303.cds.swapphase.gff peaxi162AQ_PeaxHIC303.cds.swapphase.gff
#ln -s /directory_w_genomefiles/Petunia_axillaris.v3.0.3.fasta Petunia_axillaris.v3.0.3.fasta


##### things to do before executing #####

#make a temp folder in your working directory if you don't have one already; https://gatkforums.broadinstitute.org/gatk/discussion/6083/picard-tools-markduplicates-spilling-to-disk
#mkdir ./temp

#make sure there is enough space in your directories
    #these wouldn't work for quite a few tries, but then I cleared out some space and then this worked. lesson learned - if getting a confusing error message, clean home directory for more space.

#make sure there is an index AND dictionary for your genome fasta - the samtools faidx command is not sufficient! must use Picardtools
    picard-tools CreateSequenceDictionary R=Petunia_axillaris.v3.0.3.fasta O=Petunia_axillaris.v3.0.3.dict
    samtools faidx Petunia_axillaris.v3.0.3.fasta

#make sure that your bamfiles have read group IDs: https://software.broadinstitute.org/gatk/documentation/article?id=1317
#samtools view -H /path/to/my.bam | grep '^@RG'
#if not, can use picard-tools AddOrReplaceReadGroups. As long as the RGSM is different for each of your samples, things should work.
picard-tools AddOrReplaceReadGroups RGLB=Lane1 RGPL=Illumina RGPU=TTAGGC RGSM=exs1 I=7.sortedByCoord.out.bam O=7.sortedByCoord.rg.out.bam
picard-tools AddOrReplaceReadGroups RGLB=Lane1 RGPL=Illumina RGPU=TTAGGC RGSM=exs2 I=8.sortedByCoord.out.bam O=8.sortedByCoord.rg.out.bam
picard-tools AddOrReplaceReadGroups RGLB=Lane1 RGPL=Illumina RGPU=TTAGGC RGSM=exs3 I=9.sortedByCoord.out.bam O=9.sortedByCoord.rg.out.bam

picard-tools AddOrReplaceReadGroups RGLB=Lane1 RGPL=Illumina RGPU=TTAGGC RGSM=axexs1 I=4.sortedByCoord.out.bam O=4.sortedByCoord.rg.out.bam
picard-tools AddOrReplaceReadGroups RGLB=Lane1 RGPL=Illumina RGPU=TTAGGC RGSM=axexs2 I=5.sortedByCoord.out.bam O=5.sortedByCoord.rg.out.bam
picard-tools AddOrReplaceReadGroups RGLB=Lane1 RGPL=Illumina RGPU=TTAGGC RGSM=axexs3 I=6.sortedByCoord.out.bam O=6.sortedByCoord.rg.out.bam

#doing axillaris just in case we need them
picard-tools AddOrReplaceReadGroups RGLB=Lane1 RGPL=Illumina RGPU=TTAGGC RGSM=ax1 I=1.sortedByCoord.out.bam O=1.sortedByCoord.rg.out.bam
picard-tools AddOrReplaceReadGroups RGLB=Lane1 RGPL=Illumina RGPU=TTAGGC RGSM=ax2 I=2.sortedByCoord.out.bam O=2.sortedByCoord.rg.out.bam
picard-tools AddOrReplaceReadGroups RGLB=Lane1 RGPL=Illumina RGPU=TTAGGC RGSM=ax3 I=3.sortedByCoord.out.bam O=3.sortedByCoord.rg.out.bam
