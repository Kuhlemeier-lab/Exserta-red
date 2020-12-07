#!/bin/bash
#SBATCH --mail-user=andrea.berardi@ips.unibe.ch
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=50G
#SBATCH --workdir=/home/ubelix/ips/berardi/RNAseq_axex/mapped_to_ax303/ASE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -J "snpcaller"

module add vital-it/7
module add UHTS/Analysis/GenomeAnalysisTK/3.5.0
module add UHTS/Analysis/vcftools/0.1.15


#note, started this after the step2 command: sbatch --dependency=afterok:34348384 

#use the re.bam files, which are re-aligned around the indels and snps of interest

#variant calling and filtering for bialleles for Pexserta samples (not the F1s!!!), part 1 (without ax)
GenomeAnalysisTK -T HaplotypeCaller -R Petunia_axillaris.v3.0.3.fasta -I 7.re.bam -I 8.re.bam -I 9.re.bam -dontUseSoftClippedBases -stand_emit_conf 10 -stand_call_conf 30 -o Pexs_st4_snps_raw.vcf
GenomeAnalysisTK -T VariantFiltration -R Petunia_axillaris.v3.0.3.fasta -V Pexs_st4_snps_raw.vcf -window 5 -cluster 3 --filterExpression "DP<10" --filterName DP --filterExpression "AF<0.75" --filterName AF --filterExpression "QD<2.0" --filterName QD --filterExpression "MQ<40.0" --filterName MQ --filterExpression "FS>60.0" --filterName FS -o Pexs_st4_snps_filtered.vcf
#first filter, pre-axillaris removal
GenomeAnalysisTK -T SelectVariants -R Petunia_axillaris.v3.0.3.fasta --variant Pexs_st4_snps_filtered.vcf -restrictAllelesTo BIALLELIC -select 'vc.isNotFiltered()' -o Pexs_st4_biallelic_snps_selected_pre-ax.vcf

#variant calling and filtering for bialleles for Pax samples, to make sure we can control for ax-specific SNPs
GenomeAnalysisTK -T HaplotypeCaller -R Petunia_axillaris.v3.0.3.fasta -I 1.re.bam -I 2.re.bam -I 3.re.bam -dontUseSoftClippedBases -stand_emit_conf 10 -stand_call_conf 30 -o Pax_st4_snps_raw.vcf
GenomeAnalysisTK -T VariantFiltration -R Petunia_axillaris.v3.0.3.fasta -V Pax_st4_snps_raw.vcf -window 5 -cluster 3 --filterExpression "DP<10" --filterName DP --filterExpression "AF<0.75" --filterName AF --filterExpression "QD<2.0" --filterName QD --filterExpression "MQ<40.0" --filterName MQ --filterExpression "FS>60.0" --filterName FS -o Pax_st4_snps_filtered.vcf
GenomeAnalysisTK -T SelectVariants -R Petunia_axillaris.v3.0.3.fasta --variant Pax_st4_snps_filtered.vcf -restrictAllelesTo BIALLELIC -select 'vc.isNotFiltered()' -o Pax_st4_biallelic_snps_selected.vcf


#making the exserta VCF even more strict before selecting the variants

#remove any snps between the axillaris RNAseq reads and the reference genome - these will only get in the way and are not informative.
#use vcftools to create a new vcf file (--recode option) for exserta to exclude axillaris SNPs that are shared between ax samples, exserta, and ax ref genome.
vcftools --vcf Pexs_st4_biallelic_snps_selected_pre-ax.vcf --exclude-positions Pax_st4_biallelic_snps_selected.vcf --recode --out Pexs_st4_biallelic_snps_selected_noax.vcf
mv Pexs_st4_biallelic_snps_selected_noax.vcf.recode.vcf Pexs_st4_biallelic_snps_selected_noax.vcf

#Now selecting variants for exserta that have all axillaris-specific ones removed - STRINGENT
#exclude any remaining nonbiallelic snps, exclude indels, exclude nonvariant loci
#(exclude non-variant loci -env and filtered loci -ef (trim remaining alleles by default), excluding indels --selectTypeToExclude INDEL
GenomeAnalysisTK -T SelectVariants -R Petunia_axillaris.v3.0.3.fasta --variant Pexs_st4_biallelic_snps_selected_noax.vcf -restrictAllelesTo BIALLELIC -select 'vc.isNotFiltered()' -env -ef --selectTypeToExclude INDEL -o Pexs_st4_biallelic_snps_selected_strict.vcf

#once more to keep indels, for later and not part of ASE
GenomeAnalysisTK -T SelectVariants -R Petunia_axillaris.v3.0.3.fasta --variant Pexs_st4_biallelic_snps_selected_noax.vcf -restrictAllelesTo BIALLELIC -select 'vc.isNotFiltered()' -env -ef -o Pexs_st4_biallelic_snps_selected_strict_windels.vcf






#rosetta stone if needed
#1-3 ax 1-3
#4-6 f1 1-3
#6-9 exs 1-3