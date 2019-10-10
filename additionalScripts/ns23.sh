#!/bin/bash
#PBS -S /bin/bash
#PBS -l nodes=5:ppn=24:cpu24a
#PBS -l walltime=02:00:00
#PBS -j oe
#PBS -e js_err.$PBS_JOBID
#PBS -o js_msgs.$PBS_JOBID
#PBS -N ns23
#PBS -m abe
#PBS -M parthvipulshah@pesu.pes.edu

export NPROCS=`wc -l $PBS_NODEFILE | gawk '//{print $1}'`
echo "NPROCS=$NPROCS"
echo "PBS_NODEFILE name is $PBS_NODEFILE"
echo PBS: Current home directory is $PBS_O_HOME
echo PBS: Working directory is $PBS_O_WORKDIR

cd $PBS_O_WORKDIR
echo pbs nodefile
echo `cat nodefile`

cat $PBS_NODEFILE > nodefile
cat $PBS_NODEFILE | uniq > ./mpd_nodefile_$USER
export NO_OF_COMPUTE_NODES=`wc -l mpd_nodefile_$USER | gawk '//{print $1}'`
export OMP_NUM_THREADS=24
echo mpd nodefile
echo `cat mpd_nodefile_$USER`

PPN=$(( $NPROCS / $NO_OF_COMPUTE_NODES ))
ulimit -c unlimited

source /apps/Intel/bin/compilervars.sh intel64
source /apps/Intel/impi/5.0.3.049/intel64/bin/mpivars.sh 

EXE="/sscu_gpfs/archive/suvo92/VACF/vacf.par2.1"
INPUT="-p 50000,550000,10 -a 500 -i 1000"
OUTPUT="stdout.$PBS_JOBID"

export I_MPI_PIN_PROCESSOR_LIST=0-23
export I_MPI_FABRICS=shm:tmi
export I_MPI_FALLBACK=0

time mpiexec.hydra -np $NPROCS -genvall -ppn $PPN $EXE $INPUT > $OUTPUT 2>&1

