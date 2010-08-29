#!/bin/bash
# --------------------------
# name
#$ -N ccsfa_traj
#
# Use current working directory
#$ -cwd
#
# pe request 4 slots
#$ -pe mvapich2 901
#
#$ -v MPI_HOME=/usr/local/Packages/mvapich2
# -v MPI_HOME=/usr/local/Packages/openmpi-1.2.8-intel
#
echo "Got $NSLOTS slots."
PATH=$PATH:$MPI_HOME/bin
export MPD_CON_EXT=$JOB_ID
mpiexec -machinefile $TMPDIR/machines  -n $NSLOTS $PWD/test
