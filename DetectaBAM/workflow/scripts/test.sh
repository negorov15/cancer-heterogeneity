#!/bin/bash
#SBATCH --job-name=test_cellsnp
#SBATCH --output=/home/nikoe00/DetectaBAM/logs/test_cellsnp_%j.log
#SBATCH --error=/home/nikoe00/DetectaBAM/logs/test_cellsnp_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=6G
#SBATCH --partition=main
#SBATCH --qos=standard

# Testing cellsnp-lite on a single sample to ensure it works before running the full workflow

module load Miniconda3
source /trinity/shared/easybuild/software/Miniconda3/24.7.1-0/bin/activate
conda activate scrna_tools

cd /home/nikoe00/DetectaBAM

cellsnp-lite -s /scratch/nikoe00/bam_subsets/A-253_pseudobulk.bam \
  -b resources/barcodes/a-253_barcodes_suffixed.tsv \
  -R resources/SomaticMutations/A-253/A-253_SomaticMutations_all.vcf \
  -O /tmp/cellsnp_raw_check --cellTAG CB --UMItag UB \
  --minCOUNT 1 --minMAF 0 -p 4 --gzip

zcat /tmp/cellsnp_raw_check/cellSNP.base.vcf.gz | grep -vc '^#'