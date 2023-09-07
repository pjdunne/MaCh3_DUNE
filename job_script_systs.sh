#!/bin/bash
set -x

MAP_FILE=jobs_systs.txt
FINAL_ODIR=/pnfs/dune/scratch/users/pgranger/mach3
YAML=LLH.yaml

cp ${CONDOR_DIR_INPUT}/$MAP_FILE .
LINE_NO=$((PROCESS+2))
LINE=$(sed "${LINE_NO}q;d" $MAP_FILE | tr -s ' ')
ID1=$(echo "$LINE" | cut -d' ' -f1)
NBINS=$(echo "$LINE" | cut -d' ' -f2)
ONAME=$(echo "$LINE" | cut -d' ' -f3)

WORKDIR=${_CONDOR_SCRATCH_DIR}/work/
mkdir -p $WORKDIR && cd $WORKDIR
CODEDIR=${INPUT_TAR_DIR_LOCAL}/MaCh3_DUNE/

cd $CODEDIR
source setup_dune_env.sh
LOCAL_YAML=$WORKDIR/$YAML
cp configs/$YAML  $LOCAL_YAML

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${CODEDIR}/build/lib/
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${CODEDIR}/build/_deps/yaml-cpp-build/
export LD_LIBRARY_PATH

echo $LD_LIBRARY_PATH
ldd ./build/src/DUNE_LLHscan_systs

sed -i 's#ODIR.*#ODIR: '$WORKDIR'#g' $LOCAL_YAML
sed -i 's#nbins.*#nbins: '$NBINS'#g' $LOCAL_YAML
sed -i 's#parId1.*#parId1: '$ID1'#g' $LOCAL_YAML


./build/src/DUNE_LLHscan_systs $LOCAL_YAML
test -f $WORKDIR/TestLLH.root && ifdh cp $WORKDIR/TestLLH.root $FINAL_ODIR/$ONAME