#!/bin/bash
#SBATCH --job-name= chase_test_matlab
#SBATCH --output=/home/abram/Discrete_HA/output/run%a.out
#SBATCH --error=/home/abram/Discrete_HA/output/run%a.err
#SBATCH --partition=broadwl
#SBATCH --time=12:00:00
#SBATCH --array=1-44
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=10000

module load matlab/2019b

matlab -nodisplay < /home/abram/AnnualReviewsHA/master.m
