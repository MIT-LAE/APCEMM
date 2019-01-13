#!/bin/bash

## SLURM options

# Options to sbatch start with '#SBATCH'. To disable an option, change
# the prefix, e.g. '#xSBATCH'

# Time limit after which the job will be killed. The default time limit
# is 60 minutes. Specified as HH:MM:SS or D-HH:MM
#SBATCH --time=30-00:00

# Prevent multithreading of the cores
#SBATCH --hint=nomultithread

# number of CPUs to utilize
#SBATCH -c24

# Memory per core. Job will crash if this limit is exceeded.  Default
# is 1000M per allocated core. Use values that will permit multiple
# jobs to run simultaneously when possible, e.g. a memory limit of
# 2000M will allow allow all 12 processors on a node with 24000M of
# RAM to be utilized, while specifying 2G (=2048M) would only allow 11
# of the processors to be used.
#SBATCH --mem-per-cpu=2000M

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
#SBATCH -J APCEMM

#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err


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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                                                  #
#     Aircraft Plume Chemistry, Emission and Microphysics Model    #
#                             (APCEMM)                             #
#                                                                  #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
#
# DESCRIPTION: This bash script submits an APCEMM simulation
#
#
# REVISION HISTORY:
#  15 Sep 2018 - T. Fritz - Initial version
#  20 Oct 2018 - T. Fritz - Included OMP_NUM_THREADS in slurm
#                           environments
#  21 Oct 2018 - T. Fritz - Pipe output to log
#  11 Dec 2018 - T. Fritz - Created run directories
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Set the proper # of threads for OpenMP
# SLURM_CPUS_PER_TASK ensures this matches the number you set with -c above
if [ -n "$SLURM_CPUS_PER_TASK" ]; then
    omp_threads=$SLURM_CPUS_PER_TASK
else
    omp_threads=1
fi

export OMP_NUM_THREADS=$omp_threads

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Initialize
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Define current directory
APCEMM_runDir=${PWD##}
export APCEMM_runDir

if [[ -e ${APCEMM_runDir}/APCEMM ]]; then
    rm ${APCEMM_runDir}/APCEMM
fi

log=$PWD/log.build
if [[ -e ${APCEMM_runDir}/$log ]]; then
    rm ${APCEMM_runDir}/$log
fi

make realclean
if [ "$OMP_NUM_THREADS" -eq "1" ]; then
    make -j4 spbuild 2>&1 | tee $log
else
    make -j4 mpbuild 2>&1 | tee $log
fi

# Change this to point to the APCEMM binary file
export exepath=${APCEMM_runDir}/APCEMM

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Start the simulation
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

log=$PWD/log.run

# Run APCEMM and pipe output to log file
echo ''
echo 'Running on' $OMP_NUM_THREADS 'core(s)'
echo 'Host computer: ' `hostname`
echo 'Initiation date and time: ' `date +%c`
srun -c $OMP_NUM_THREADS time -p $exepath 2>&1 | tee $log

# Echo end
echo ''
echo 'Run ended at ' `date +%c`

# Report additional information if job was killed because of memory limit
oom_check $?

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Clean up
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Clear variable
unset id
unset log
unset APCEMM_runDir
unset exepath
