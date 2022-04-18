#!/bin/bash

#SBATCH -N 1                                                                    
#SBATCH -p htc                                                                  
#SBATCH -t 00:45:00                                                             
#SBATCH -n 1                                                                    
#SBATCH --mem=8G                                                               
#SBATCH --mail-user=daniel.fridljand@embl.de                                                
#SBATCH --mail-type=FAIL,END 

module load R

Rscript Referencemanual.R



