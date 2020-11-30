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

makeblastdb -in Psecreta_pep.fasta -dbtype prot -title Psec_k119_prot_db -out Psec_k119_prot -input_type fasta

#search secreta
hmmsearch -o Petunia_secreta_k119_MYB_search.out -A P_sec_MYBs_k119_alignment.fasta --cpu 12 MYBs.hmm Psecreta_pep.fasta

