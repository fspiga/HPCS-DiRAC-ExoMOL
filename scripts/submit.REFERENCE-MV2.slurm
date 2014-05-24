#!/bin/bash
#SBATCH -J DIRAC
#SBATCH -p sandybridge
#SBATCH --nodes=16
#SBATCH --ntasks=256
#SBATCH --qos=support
#SBATCH --account=spiga-dirac
#SBATCH --time=6:00:00

workdir="$SLURM_SUBMIT_DIR"
cd $workdir

. /etc/profile.d/modules.sh                # Leave this line
module purge
module load default-wilkes
module unload intel/impi cuda intel/mkl intel/cce intel/fce
module load intel/fce/14.0.2.144
module load intel/cce/14.0.2.144
module load intel/mkl/11.1.2.144
module load fs395/scalapack/2.0.2/intel-mv2
module load fs395/limic2/0.5.6
module load fs395/mvapich2-2.0rc1/intel
module load fs395/elpa/2013.11-v8/intel-mv2-avx

numnodes=$SLURM_JOB_NUM_NODES
numprocs=$SLURM_NTASKS
ppn=$(echo "$SLURM_TASKS_PER_NODE" | sed -e  's/^\([0-9][0-9]*\).*$/\1/')

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=$OMP_NUM_THREADS
export DIRAC_SIZE=60000
export DIRAC_NROOTS=60000
export DIRAC_SOLVER=PDSYEVD

cat  << 'EOF' > gen.inp
(large jobs for hot ch4)
GEN-MAT 60000
NROOTS 60000
DIAGONALIZER PDSYEVD (needs to be stated)
GENERATOR RANDOM-LOCAL
EOF

application="./diag_MV2-INTEL_ELPA.x < gen.inp"
options="| tee out.ELPA-MV2.TESTING-${DIRAC_SOLVER}.N-${DIRAC_SIZE}.ROOTS-${DIRAC_NROOTS}.NODES-$SLURM_JOB_NUM_NODES.$SLURM_JOB_ID"
mv2_affinity="-genv MV2_ENABLE_AFFINITY 1 -genv MV2_CPU_BINDING_LEVEL CORE -genv MV2_CPU_BINDING_POLICY BUNCH"
mv2_shmem="-genv MV2_USE_SHARED_MEM 1 -genv MV2_SMP_USE_LIMIC2 1"
mv2_cuda="-genv MV2_USE_CUDA 0"
CMD="mpirun -np $numprocs -ppn $ppn $mv2_affinity $mv2_shmem $mv2_cuda $application $options"

###############################################################
#### You should not have to change anything below this line ####
###############################################################

echo -e "Changed directory to `pwd`.\n"

echo -e "JobID: $SLURM_JOB_ID\n======"
echo "Time: `date`"
echo "Running on master node: `hostname`"
echo "Current directory: `pwd`"

if [ "$SLURM_JOB_NODELIST" ]; then
        #! Create a machine file:
	export NODEFILE=`generate_pbs_nodefile`
	cat $NODEFILE | uniq > machine.file.$SLURM_JOB_ID
        echo -e "\nNodes allocated:\n================"
        echo `cat machine.file.$SLURM_JOB_ID | sed -e 's/\..*$//g'`
fi

echo -e "\nnumprocs=$numprocs, numnodes=$numnodes, ppn=$ppn (OMP_NUM_THREADS=$OMP_NUM_THREADS)"

echo -e "\nExecuting command:\n==================\n$CMD\n"

eval $CMD 

echo "Time: `date`"
