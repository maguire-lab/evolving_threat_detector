#!/bin/bash
#SBATCH --time=3-00:00:00
#SBATCH --partition=cpubase_bycore_b4
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --account=rrg-rbeiko
#SBATCH --job-name=etd_ngono_query
#SBATCH --output=etd_ngono_query_%j.out
#SBATCH --error=etd_ngono_query_%j.err

module load apptainer
module load nextflow

export NXF_SINGULARITY_CACHEDIR=/project/rrg-rbeiko/$USER/NXF_SINGULARITY_CACHEDIR
export SLURM_ACCOUNT=rrg-rbeiko
export NXF_OPTS='-XX:ActiveProcessorCount=1'
export NXF_OFFLINE='true'

cd /project/rrg-rbeiko/$USER/etd

nextflow run main.nf \
    --mode query_db \
    --input_query /project/rrg-rbeiko/$USER/etd_validation/data/processed/week2/ngono_validation_samplesheets/query_samplesheet.csv \
    --outdir /project/rrg-rbeiko/$USER/ngono_combined_results \
    --amrfinder_db /project/rrg-rbeiko/$USER/amrfinderdb \
    -profile singularity \
    -resume \
    -work-dir /scratch/presh/ngono_query_work_dir
~                                                        
