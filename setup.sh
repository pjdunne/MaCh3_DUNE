#!/bin/bash

# e.g. this:
#export PATH=/home/cwret/procmail-3.22/new:${PATH}
#source /home/cwret/CMT/setup.sh
#source /home/cwret/root/bin/thisroot.sh
#module load 
#export CUDAPATH=${CUDA_HOME}
#export MACH3_DATA=/home/cwret/P6Data
#export MACH3_MC=/home/cwret/P6MC


#################################################################

if [[ $HOSTNAME == *gpu.rl.ac.uk ]]; then
    # Add lockfile to path (required on emerald!) (commented out for those not on EMERALD)
    export PATH=$PATH:/home/oxford/eisox159/procmail/bin

    # setup root on EMERALD (commented out for those not on EMERALD)
    source /home/stfc/eisext13/root/bin/thisroot.sh
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/stfc/eisext13/root/lib
    export PATH=$PATH:/home/stfc/eisext13/root/bin
    export ROOTLIBDIR=/home/stfc/eisext13/root/lib
    export ROOTINCLUDEDIR=/home/stfc/eisext13/root/include
    module load cuda/8.0.44
    export CUDAPATH=${CUDA_HOME}
fi

# Cuda directory
#   Unset if you don't have a CUDA GPU
#   If you want GPU point CUDAPATH to your CUDA install and makefile will pick this up
if [[ $HOSTNAME == *hep.ph.ic.ac.uk ]]; then
    if [ -z $CUDAPATH ]; then
        source /vols/software/cuda/setup.sh 11.2.0 # (Latest = 11.2.0) Can give this an argument to get a particular CUDA version, e.g. 10.2.2
	export CUDAPATH=$CUDA_PATH
    fi
fi

	# Multithreading? 
	# MULTITHREAD=1 means MP is on, if environment variable doesn't exist it's off
	unset MULTITHREAD
	export MULTITHREAD=1

	# Normal DEBUG flag (internal for MaCh3; mostly timing info)
	# for samplePDFND2014: 1 = print event rates for each iteration broken down by systematic contrib.
	#                      2 = print event by event weights broken down by systematic contrib.
	# Also turns on things like checking TSpline3->Eval() calls are good, that we have no negative weights, etc
	unset DEBUG_ND280
	#export DEBUG_ND280=1
	# Debugging is most verbose for multithreading enabled and OMP_NUM_THREADS=1 is then required. Could be fixed but hasn't happened yet
	export OMP_NUM_THREADS=20

	# For GPU code we can additionally dump CPU and GPU cross-section weight output
	# This checks for biases in the weights
	unset DEBUG_ND280_DUMP
	#export DEBUG_ND280_DUMP=1

	# Debug flags for MaCh3 GPU (internal for MaCh3)
	unset DEBUG_CUDA_ND280
	#export DEBUG_CUDA_ND280=1


	# Automatically set CUDA for Emerald
	if [[ $HOSTNAME == *gpu.rl.ac.uk ]]; then
	  module load cuda/8.0.44
	  export CUDAPATH=${CUDA_HOME}
	# Automatically set up for ComputeCanada
	elif [[ $HOSTNAME == lg-1r[47]* ]]; then
	  module load CUDA/7.5.18
	  export CUDAPATH=${CUDA_HOME}
	  # Also need GSL for some external dependencies
	  module load GSL
        fi


	# Multithreading? (COMMENT OUT TO TURN OFF MULTIPLE CPU)
	#unset MULTITHREAD
	#export MULTITHREAD=1


	# Set the EXPERIMENT MaCh3 directory
	if [ -z $MACH3 ]; then
	  export MACH3=$(pwd)
fi

	# Set the CORE MaCh3 directory
	if [ -z $MACH3CORE ]; then
	  export MACH3CORE=$(pwd)/../MaCh3_refactor
fi

export CMTPATH=$(pwd):CMTPATH
# Add the libraries to lib path
# Add the executables to path
export PATH=${MACH3}/bin:${PATH}

# Add the CUDA libraries and executables if we have CUDA enabled
if [ -d "$CUDAPATH" ]; then
  echo "BLAH"
  # Add the CUDA binaries (e.g. nvcc!)
  export PATH=${CUDAPATH}/bin:${PATH}
  # Add the CUDA libraries (e.g. cudarand)
  export LD_LIBRARY_PATH=${CUDAPATH}/lib64:/usr/lib64/nvidia:${LD_LIBRARY_PATH}
  # Report on what CUDA settings and GPUs we're running
  echo "******************************************"
  echo Running GPU on: $(nvidia-smi --query-gpu=name --format=csv,noheader)
  echo with CUDAPATH = $CUDAPATH
  echo "******************************************"
else
  # Report that we aren't running with GPU
  echo "******************************************"
  echo DID NOT SET CUDAPATH!
  echo This is OK if you do not want to run with a GPU
  echo "******************************************"
fi

mkdir -pv ${MACH3}/bin
mkdir -pv ${MACH3}/lib
#mkdir -pv ${MACH3}/utils/agf/include
mkdir -pv ${MACH3}/plots

echo ========================================================
echo "Are you resetting to change GPU/CPU options?"
echo "If so, don't forget to fully remake!"
echo "'make clean && make'"
echo ========================================================
echo "Finished setup"
