#!/bin/bash -l
#
# Check if our working directory is on the central file server
#
#export verbose

#cd  /shared/scratch/ucapsyu/production/

export exec="diag_pdsyev_darwin.x"

export pwd=`pwd`
echo $pwd

export name=`echo $1 | sed -e 's/\.inp//'`
echo $name

export JOB=$3
echo $JOB

export dir1=`pwd | awk -F/ '{print $2}'`

if [ -e "$name.o" ]; then
   /bin/rm $name.o
fi

if [ -e "$name.e" ]; then
   /bin/rm $name.e
fi

if [ -e "$name.out" ]; then
  if [ -e "$name.tmp" ]; then
    /bin/rm $name.tmp
  fi
  /bin/mv $name.out $name.tmp
fi



export PARNODES=$2
export dmem=`expr $2 \* 63900`
export ntasks=`expr $PARNODES \* 16` 

export nprocs=$ntasks

echo "Nnodes=" $PARNODES, "Nprocs=" $nprocs

export MEM=`echo $nprocs $dmem | awk  '{printf( "%8.0f\n", $2*$1 )}'`


export jobtype="small2"
export wclim=12

if [ $nprocs -lt "8" ]; then
   export jobtype=""
fi

export wclim=$3


#wclim=1

echo "Nnodes=" $PARNODES, "Nprocs=" $nprocs, " Memory = " $dmem, "jobtype = sandybridge", "wclimit = " $wclim
echo "Working dir is " $pwd

#msub -N $name -j oe -e $name.e -q $jobtype -l "walltime=$wclim:00:00,pmem=${dmem}mb,nodes=$PARNODES:ppn=8" \
#     -v "setenv name $name,setenv pwd $pwd" \
#     -v "setenv nprocs $nprocs" \
#     $pwd/run_pdiag.sh


sbatch -A DIRAC-dp020 --nodes=$PARNODES --ntasks=$nprocs --time=$wclim:00:00  -J $name -o $name.o -e $name.e   \
     --workdir=$pwd --hint=compute_bound --no-requeue -p sandybridge \
     $pwd/sub_script.csh $nprocs $name $exec $pwd
     

