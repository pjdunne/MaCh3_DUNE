#!/bin/bash
set -x

ONAME=AtmChain.root

FINAL_ODIR=/pnfs/dune/scratch/users/pgranger/mach3
YAML=AtmChain.yaml
ID=$PROCESS

WORKDIR=${_CONDOR_SCRATCH_DIR}/work/
mkdir -p $WORKDIR && cd $WORKDIR
CODEDIR=${INPUT_TAR_DIR_LOCAL}/MaCh3_DUNE/

cd $CODEDIR
source setup_dune_env.sh
LOCAL_YAML=$WORKDIR/$YAML
cp configs/$YAML $LOCAL_YAML

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${CODEDIR}/build/lib/
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${CODEDIR}/build/_deps/yaml-cpp-build/
export LD_LIBRARY_PATH

echo $LD_LIBRARY_PATH
ldd ./build/src/DUNE_atm_chain

sed -i 's#ODIR.*#ODIR: '$WORKDIR'#g' $LOCAL_YAML

./build/src/DUNE_atm_chain $LOCAL_YAML
test -f $WORKDIR/$ONAME && ifdh cp $WORKDIR/$ONAME $FINAL_ODIR/AtmChain_${ID}_real.root
