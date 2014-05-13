#!/bin/bash
#SBATCH -J DIRAC
#SBATCH -p sandybridge
#SBATCH --nodes=4
#SBATCH --ntasks=64
#SBATCH --qos=support
#SBATCH --account=spiga-dirac
#SBATCH --time=2:00:00

workdir="$SLURM_SUBMIT_DIR"
cd $workdir

. /etc/profile.d/modules.sh                # Leave this line
module purge
module load default-wilkes
module unload intel/impi cuda intel/mkl
module load intel/mkl/11.1.2.144
module load fs395/openmpi-1.8.0/intel
module load fs395/scalapack/2.0.2/intel-ompi
module load fs395/elpa/2013.11-v8/intel-ompi
module load fs395/IPM/2.0.3/openmpi-1.8.0_intel
module load fs395/papi/5.2.0

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
GEN-MAT 50000
NROOTS 50000
DIAGONALIZER ELPA-1STAGE (needs to be stated)
GENERATOR SYM-POSITIVE-O2
EOF

application="./diag_OMPI-INTEL.x < gen.inp"
options="| tee out.TESTING.OMPI.ELPA-1STAGE.SYM-POSITIVE-O2.N-100K.$JOBID"

fca_on="-mca coll_fca_enable 1 -mca coll_fca_np 0"
fca_off="-mca coll_fca_enable 0"
hcoll_on="-mca coll_hcoll_enable 1 -mca coll_hcoll_np 0" 
hcoll_off="-mca coll_hcoll_enable 0"
mxm_off="-mca mtl ^mxm -mca btl openib,self,sm"
mxm_on="-mca pml cm -mca mtl mxm -mca mtl_mxm_np 0 -x MXM_TLS=shm,rc,ud,self" 
mpirun_options="--bind-to core $common_params $fca_off $hcoll_off $mxm_off"
CMD="$ipm_preload mpirun --map-by ppr:$ppn:node -np $numprocs $mpirun_options $application $options"

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

