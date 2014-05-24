#!/bin/bash
#SBATCH -J DIRAC
#SBATCH -p sandybridge
#SBATCH --nodes=32
#SBATCH --ntasks=512
#SBATCH --qos=support
#SBATCH --account=spiga-dirac
#SBATCH --time=4:00:00

workdir="$SLURM_SUBMIT_DIR"
cd $workdir

. /etc/profile.d/modules.sh                # Leave this line
module purge
module load default-wilkes
module unload intel/impi cuda intel/mkl intel/cce intel/fce
module load intel/fce/14.0.2.144
module load intel/cce/14.0.2.144
module load intel/mkl/11.1.2.144
module load intel/impi/4.1.3.048
module load fs395/elpa/2013.11-v8/intel-impi-avx

numnodes=$SLURM_JOB_NUM_NODES
numprocs=$SLURM_NTASKS
ppn=$(echo "$SLURM_TASKS_PER_NODE" | sed -e  's/^\([0-9][0-9]*\).*$/\1/')

JOBID=$SLURM_JOB_ID

export IPM_KEYFILE=${IPM_HOME}/etc/ipm_key_mpi
export IPM_LOG=full
ipm_preload="env LD_PRELOAD=${IPM_HOME}/lib/libipm.so"

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=$OMP_NUM_THREADS

cat  << 'EOF' > gen.inp
(large jobs for hot ch4)
GEN-MAT 100000
NROOTS 100000
DIAGONALIZER ELPA-1STAGE (needs to be stated)
GENERATOR SYM-POSITIVE-O2
EOF

application="./diag_IMPI-INTEL.x < gen.inp"
options="| tee out.TESTING.IMPI.ELPA-1STAGE.SYM-POSITIVE-O2.N-100K.$JOBID"

export I_MPI_PIN_DOMAIN=omp:compact # Domains are $OMP_NUM_THREADS cores in size
export I_MPI_PIN_ORDER=scatter # Adjacent domains have minimal sharing of caches/sockets

CMD="mpirun -ppn $mpi_tasks_per_node -np $np $application $options"

###############################################################
### You should not have to change anything below this line ####
###############################################################

echo -e "Changed directory to `pwd`.\n"

echo -e "JobID: $JOBID\n======"
echo "Time: `date`"
echo "Running on master node: `hostname`"
echo "Current directory: `pwd`"

if [ "$SLURM_JOB_NODELIST" ]; then
        #! Create a machine file:
	export NODEFILE=`generate_pbs_nodefile`
	cat $NODEFILE | uniq > machine.file.$JOBID
        #echo $SLURM_JOB_NODELIST | sed -e 's/,/\n/' > machine.file.$JOBID
        
        echo -e "\nNodes allocated:\n================"
        echo `cat machine.file.$JOBID | sed -e 's/\..*$//g'`
fi

echo -e "\nnumprocs=$numprocs, numnodes=$numnodes, ppn=$ppn (OMP_NUM_THREADS=$OMP_NUM_THREADS)"

echo -e "\nExecuting command:\n==================\n$CMD\n"

eval $CMD 

echo "Time: `date`"

