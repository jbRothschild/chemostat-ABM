#!/bin/bash

# Turn off implicit threading in Python, R
export OMP_NUM_THREADS=1

# EXECUTION COMMAND; ampersand off 40 jobs and wait

# mkdir -p ${RESULTS_DIR}/${SIM_DIR}

parallel --joblog slurm-$SLURM_JOBID.log --sshdelay 0.1 --wd $PWD "./des1015.o 1015 {}" ::: {0..4999}