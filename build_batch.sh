#!/bin/sh
#SBATCH -J MC_sim 
#SBATCH -o  1000BL.o%j.txt
#SBATCH -e  1000BL.e%j.txt

#SBATCH -N 1
#SBATCH -n 1 # used for MPI codes, otherwise leave at '1'
##SBATCH --ntasks-per-node=1  # don't trust SLURM to divide the cores evenly
##SBATCH --cpus-per-task=1  # cores per task; set to one if using MPI
##SBATCH --exclusive  # using MPI with 90+% of the cores you should go exclusive
#SBATCH --mem-per-cpu=2G  # memory per core; default is 1GB/core


#SBATCH -t 0-48:00:00  # days-hours:minutes

## mail alert at start, end and abortion of execution
#SBATCH --mail-type=ALL

## send mail to this address
#SBATCH --mail-user=mwsrgf@mail.missouri.edu

## mpirun nrniv -mpi MC_main_small_forBeta_shortburstensamble.hoc #srun
python build_network.py #srun