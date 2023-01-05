#!/bin/bash

#SBATCH --partition=shortterm
#SBATCH --mem=90GB
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=6
#SBATCH --nodes=1




#source /etc/profile.d/modules.sh
#module load fastp/v0.20.0
#module load bwa/v0.7.17
#module load samtools/v1.13
#module load picard-tools/2.18.25

# download scRNAseq data from geo data sets 

fastq-dump  --split-3   SRR9036398    #### https://www.ncbi.nlm.nih.gov/sra?term=SRX5813546   
