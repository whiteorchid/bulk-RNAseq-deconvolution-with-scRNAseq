#!/bin/bash
#SBATCH --partition=longterm
#SBATCH --mem=160GB
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=20
#SBATCH --nodes=1

source /etc/profile.d/modules.sh
module load kallisto/0.45.0

### whitelist download from 10x website, output file read by Seurat package and downstream analysis

STAR --runThreadN 20 \
--genomeDir /data/sb_service_01/guo/ref/ref_star \
--readFilesIn  /data/sb_service_01/guo/sg_20221118_singlecell_test/SRR9036398_2.fastq /data/sb_service_01/guo/sg_20221118_singlecell_test/SRR9036398_1.fastq \
--soloType  CB_UMI_Simple   --soloCBwhitelist /data/sb_service_01/guo/sg_20221118_singlecell_test/whitelist.txt

