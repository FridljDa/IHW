#!/bin/bash

# Example of running R script with a job array

#SBATCH -A huber                # group to which you belong
#SBATCH -N 1                        # number of nodes
#SBATCH -n 1                        # number of cores
#SBATCH --mem 10                    # memory pool for all cores
#SBATCH -t 0-1:00                   # runtime limit (D-HH:MM:SS)
#SBATCH -o measure_time_scripts/out_measure_time/measure_time_out-%j.out
#SBATCH -e measure_time_scripts/err_measure_time/measure_time_er-%j.err          # STDERR
#SBATCH --mail-type=END,FAIL        # notifications for job done & fail
#SBATCH --mail-user=daniel.fridljand@embl.de # send-to address
# Load software

module load R
# Run R script
Rscript measure_time_scripts/measure_time.R 