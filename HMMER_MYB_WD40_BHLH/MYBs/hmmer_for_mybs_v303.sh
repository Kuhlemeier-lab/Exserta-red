#!/bin/bash
# Slurm options
#SBATCH --mail-user=andrea.berardi@ips.unibe.ch
#SBATCH --mail-type=begin,fail,end
#SBATCH --job-name="hmmer mybs"
#SBATCH --workdir=/home/ubelix/ips/berardi/MYBs/new_annots
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=16GB

module load vital-it/7
module load Blast/ncbi-blast/latest
#module load SequenceAnalysis/GenePrediction/augustus/3.2.3
module load SequenceAnalysis/HMM-Profile/hmmer/3.1b2


#generate the hmm profile
hmmbuild MYBs.hmm mybs_aa.sto

# bash-4.2$ hmmbuild MYBs.hmm mybs_aa.sto
# hmmbuild :: profile HMM construction from multiple sequence alignments
# HMMER 3.1b2 (February 2015); http://hmmer.org/
# Copyright (C) 2015 Howard Hughes Medical Institute.
# Freely distributed under the GNU General Public License (GPLv3).
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# input alignment file:             mybs_aa.sto
# output HMM file:                  MYBs.hmm
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# idx name                  nseq  alen  mlen eff_nseq re/pos description
#---- -------------------- ----- ----- ----- -------- ------ -----------
#1     mybs_aa                179   181   133     1.83  0.590

# CPU time: 0.22u 0.00s 00:00:00.22 Elapsed: 00:00:00.22


#makeblastdb first
makeblastdb -in PeaxHIC303.pep.fasta -dbtype prot -title Pax_303_prot_db -out Pax303_prot -input_type fasta
makeblastdb -in PeexHIC303.pep.fasta -dbtype prot -title Pexs_303_prot_db -out Pexs303_prot -input_type fasta


#had to remove the periods in this file
#sed -e 's/\.//g' P.EXSERTA.contigs.v1.1.3.annotation.v1.proteins.fasta > P.EXSERTA.contigs.v1.1.3.annotation.v1.proteins.fasta
#makeblastdb -in P.EXSERTA.contigs.v1.1.3.annotation.v1.proteins.fasta -dbtype prot -title Pexs_prot_db -out Pexs_prot -input_type fasta

#search axillaris, multithread (n+1, here 11+1=12)
hmmsearch -o Petunia_axillaris_303_MYB_search.out -A P_ax_MYBs_303_alignment.fasta --cpu 12 MYBs.hmm PeaxHIC303.pep.fasta

#search exserta
hmmsearch -o Petunia_exserta_303_MYB_search.out -A P_exs_MYBs_303_alignment.fasta --cpu 12 MYBs.hmm PeexHIC303.pep.fasta
