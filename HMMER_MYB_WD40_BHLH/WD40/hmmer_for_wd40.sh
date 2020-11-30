#!/bin/bash
# Slurm options
#SBATCH --mail-user=andrea.berardi@ips.unibe.ch
#SBATCH --mail-type=begin,fail,end
#SBATCH --job-name="hmmer wd40"
#SBATCH --workdir=/home/ubelix/ips/berardi/wd40
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=16GB

module load vital-it/7
module load Blast/ncbi-blast/latest
#module load SequenceAnalysis/GenePrediction/augustus/3.2.3
module load SequenceAnalysis/HMM-Profile/hmmer/3.1b2


#downloaded fasta for solanales/solanaceae (same set of entries) from pfam https://pfam.xfam.org/family/WD40#tabview=tab7

#generate the hmm profile
#hmmbuild WD40.hmm WD40_solanaceae.sto

# hmmbuild :: profile HMM construction from multiple sequence alignments
# HMMER 3.1b2 (February 2015); http://hmmer.org/
# Copyright (C) 2015 Howard Hughes Medical Institute.
# Freely distributed under the GNU General Public License (GPLv3).
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# input alignment file:             WD40_solanaceae.sto
# output HMM file:                  WD40.hmm
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# idx name                  nseq  alen  mlen eff_nseq re/pos description
#---- -------------------- ----- ----- ----- -------- ------ -----------
#1     WD40_solanaceae       5547 12681   501    37.68  0.590

# CPU time: 4.40u 0.28s 00:00:04.68 Elapsed: 00:00:04.70


#makeblastdb first
#makeblastdb -in PeaxHIC303.pep.fasta -dbtype prot -title Pax_303_prot_db -out Pax303_prot -input_type fasta
#makeblastdb -in PeexHIC303.pep.fasta -dbtype prot -title Pexs_303_prot_db -out Pexs303_prot -input_type fasta



#search axillaris, multithread (n+1, here 11+1=12)
hmmsearch -o Petunia_axillaris_303_WD40_search.out -A P_ax_WD40_303_alignment.fasta --cpu 12 WD40.hmm PeaxHIC303.pep.fasta

#search exserta
hmmsearch -o Petunia_exserta_303_WD40_search.out -A P_exs_WD40_303_alignment.fasta --cpu 12 WD40.hmm PeexHIC303.pep.fasta
