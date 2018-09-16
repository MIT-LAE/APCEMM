#!/bin/bash

## SLURM options

# Options to sbatch start with '#SBATCH'. To disable an option, change
# the prefix, e.g. '#xSBATCH'

# Time limit after which the job will be killed. The default time limit
# is 60 minutes. Specified as HH:MM:SS or D-HH:MM
#SBATCH --time=1-00:00

# Prevent multithreading of the cores
#SBATCH --hint=nomultithread

# number of CPUs to utilize
#SBATCH --ntasks=4

# Memory per core. Job will crash if this limit is exceeded.  Default
# is 1000M per allocated core. Use values that will permit multiple
# jobs to run simultaneously when possible, e.g. a memory limit of
# 2000M will allow allow all 12 processors on a node with 24000M of
# RAM to be utilized, while specifying 2G (=2048M) would only allow 11
# of the processors to be used.
#SBATCH --mem-per-cpu=6000M

# Number of nodes to utilize
#SBATCH -N 1

# Partition (queue) for the job
# - "normal" has a time limit of 30 days. Per-user resource limits apply.
# - "debug" has a time limit of 1 hour, and low resource limits,
#   but will pre-empt other jobs (memory permitting).
# - "idle" has no time limit or per-user limit, but jobs can be stopped
# (requeued) by jobs in the normal or debug queues.
#SBATCH --partition=normal
# or, for short, '-p normal'

# Send email when the job begins and ends
#xSBATCH --mail-user=your.address@example.com
#xSBATCH --mail-type=BEGIN,END

# Set a name for the job. This will appear in the output of 'squeue'. The
# default is the name of the job script. This option is probably most useful
# as a command-line argument to sbatch.
#xSBATCH -J some_name

# To submit the job, run:
#
#     sbatch <scriptname.sh>
#
# Additional options or overrides can be specified on the command line:
#
#     sbatch --partition=debug -J case4a scriptname.sh
#
# Command line arguments to your script can be passed in as well:
#
#     sbatch scriptname.sh foo bar
#
# in which case the variables $1 and $2 will be set to 'foo' and 'bar'
# within the script.

# Change this to point to the APCHEM binary file

export exepath=${HOME}/CAPCEMM/build/apps/APCEMM

currDir=${PWD##*/}

echo "Running APCEMM"
echo "Directory: " $currDir
#echo "Running on " $OMP_NUM_THREADS " cores"
echo "Host computer: " `hostname`
echo "Initiation date and time: " `date +%c`

## Job commands

time ${exepath}
#srun --mpi=pmi2 -n 4 ./mpi_hello

# Report additional information if job was killed because of memory limit
oom_check $?
