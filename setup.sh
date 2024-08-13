#!/bin/bash
#################################################################

# Cuda directory
#   Unset if you don't have a CUDA GPU
#   If you want GPU point CUDAPATH to your CUDA install and makefile will pick this up
if [[ $HOSTNAME == *hep.ph.ic.ac.uk ]]; then
  if [ -z $CUDAPATH ]; then
	source /vols/software/cuda/setup.sh 11.2.0 # (Latest = 11.2.0) Can give this an argument to get a particular CUDA version, e.g. 10.2.2
	export CUDAPATH=$CUDA_PATH
  fi
fi

# Add the CUDA libraries and executables if we have CUDA enabled
if [ -d "$CUDAPATH" ]; then
  export PATH=${CUDAPATH}/bin:${PATH}
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

mkdir -pv ${MACH3}/plots
echo "Finished setup"
