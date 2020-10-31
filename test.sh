sbatch <<- EOF
#!/bin/bash
#SBATCH -N 1
#SBATCH --partition=batch
#SBATCH -J TCGA
#SBATCH -o slurm/TCGA.%J.out
#SBATCH -e slurm/TCGA.%J.err
#SBATCH --time=1-00:00:00
#SBATCH --mem=100M
#SBATCH --cpus-per-task=32
#run the application:
sleep 30
EOF