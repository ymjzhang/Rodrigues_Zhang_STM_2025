#!/bin/sh
# request Bourne shell as shell for job
#SBATCH -N 1
#SBATCH -n 1
###################
echo $1
gunzip $1
