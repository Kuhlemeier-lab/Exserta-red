#!/bin/bash
#SBATCH --mail-user=andrea.berardi@ips.unibe.ch
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=60G
#SBATCH --workdir=/home/ubelix/ips/berardi/RNAseq_axex/mapped_to_ax303/ASE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -J "ASE step 1"

module add vital-it/7
#module add UHTS/Analysis/picard-tools/1.130;
module add UHTS/Analysis/picard-tools/2.18.11
#module add Development/java_jdk/1.8.0_112;

#after https://wiki.systemsx.ch/display/~tural.yarahmadov@ips.unibe.ch/RNA-seq+data+analysis+pipeline


#Tural helped me to understand this code. Process each individual bamfile through the first four steps, then process them per-group through the last two of this pipeline.

#Making symbolic links for this analysis so that I can use an array command
#ln -s source_file symboliclinknameforfile

#ax_exF1Limb1_Aligned.sortedByCoord.out.bam  ax_Limb3_Aligned.sortedByCoord.out.bam
#ax_exF1Limb2_Aligned.sortedByCoord.out.bam  exs_Limb1_Aligned.sortedByCoord.out.bam
#ax_exF1Limb3_Aligned.sortedByCoord.out.bam  exs_Limb2_Aligned.sortedByCoord.out.bam
#ax_Limb1_Aligned.sortedByCoord.out.bam      exs_Limb3_Aligned.sortedByCoord.out.bam
#ax_Limb2_Aligned.sortedByCoord.out.bam

#F1 group
    #/home/ubelix/ips/berardi/RNAseq_axex/ASE
#    ln -s ../ax_exF1Limb1_Aligned.sortedByCoord.out.bam 4.sortedByCoord.out.bam
#    ln -s ../ax_exF1Limb2_Aligned.sortedByCoord.out.bam 5.sortedByCoord.out.bam
#    ln -s ../ax_exF1Limb3_Aligned.sortedByCoord.out.bam 6.sortedByCoord.out.bam
        
#axillaris group
#    ln -s ../ax_Limb1_Aligned.sortedByCoord.out.bam 1.sortedByCoord.out.bam
#    ln -s ../ax_Limb2_Aligned.sortedByCoord.out.bam 2.sortedByCoord.out.bam
#    ln -s ../ax_Limb3_Aligned.sortedByCoord.out.bam 3.sortedByCoord.out.bam

#exserta group
#    ln -s ../exs_Limb1_Aligned.sortedByCoord.out.bam 7.sortedByCoord.out.bam
#    ln -s ../exs_Limb2_Aligned.sortedByCoord.out.bam 8.sortedByCoord.out.bam
#    ln -s ../exs_Limb3_Aligned.sortedByCoord.out.bam 9.sortedByCoord.out.bam
    
    

#I also made a symbolic link (or copied if space) to the reference genome fasta file in this directory, Petunia_axillaris_v1.6.2_genome.fasta
#ln -s /home/ubelix/ips/berardi/axillaris/Paxillaris_v303/peaxi162AQ_PeaxHIC303.cds.swapphase.gff peaxi162AQ_PeaxHIC303.cds.swapphase.gff
#ln -s /home/ubelix/ips/berardi/axillaris/Paxillaris_v303/Petunia_axillaris.v3.0.3.fasta Petunia_axillaris.v3.0.3.fasta


##### things to do beforehand #####

#make a temp folder in your working directory if you don't have one already; https://gatkforums.broadinstitute.org/gatk/discussion/6083/picard-tools-markduplicates-spilling-to-disk
#mkdir ./temp

#make sure there is enough space in your directories
    #these wouldn't work for quite a few tries, but then I cleared out some space and then this worked. lesson learned - if getting a confusing error message, clean home directory for more space.

#make sure there is an index AND dictionary for your genome fasta - the samtools faidx command is not sufficient! must use Picardtools
    picard-tools CreateSequenceDictionary R=Petunia_axillaris.v3.0.3.fasta O=Petunia_axillaris.v3.0.3.dict
    samtools faidx Petunia_axillaris.v3.0.3.fasta

#make sure that your bamfiles have read group IDs: https://software.broadinstitute.org/gatk/documentation/article?id=1317
#samtools view -H /path/to/my.bam | grep '^@RG'
#if not, can use picard-tools AddOrReplaceReadGroups, where As long as the RGSM is different for each of your samples, things should work
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
