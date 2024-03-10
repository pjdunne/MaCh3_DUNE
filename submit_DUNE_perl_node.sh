#!/bin/bash
#
# Submit job script, bash style
# Because I prefer sed and shizz
#
# This script will run multiple MCMC chains for a SINGLE fit
# on the Emerald cluster. Typically when the jobs finish you can combine
# or reduce the .root files into a single file for analysis.
#
# Current good clusters:  STFC Emerald, 
#                         Compute Canada (Cedar, Graham, Beluga, Niagara, Guillimin, Hades, Helios)
#                         Imperial HPC
#                         Imperial HEP (heppc105, heppc205, lt2gpu00)

# A function to calculate the time from seconds to HH:MM:SS
function SecToH
{
  if [[ "$#" -ne 1 ]]; then
    echo "Second to hour convertor needs one argument"
    exit -1
  fi

  num=$1
  ss=00
  mm=00
  hh=00
  if ((num>59)); then
    ((ss=num%60))
    ((num=num/60))
    if ((num>59)); then
      ((mm=num%60))
      ((num=num/60))
      if ((num>23)); then
        ((hh=num))
      else
        ((hh=num))
      fi
    else
      ((mm=num))
    fi
  else
    ((ss=num))
  fi

  # Now set the variables to 00 form
  ss=$(printf "%.2i" $ss)
  mm=$(printf "%.2i" $mm)
  hh=$(printf "%.2i" $hh)

  echo "$hh:$mm:$ss"
}

if [[ "$#" -ne 6 ]]; then
  echo "I need 6 parameters, you gave $#"
  echo "./submit_DUNE_niagara.sh  EXECUTABLE  CONFIG_FILE  FIT_NAME  N_CHAINS  N_STEPS JOB"
  exit -1
fi

# First argument is executable
EXE=$1

# Second argument is base config file
CFG=$2

# Third argument is output name
FitName=$3

# Fourth argument is number of jobs
NCHAINS=$4
case $NCHAINS in
  ''|[!0-9]*) echo "NCHAINS = $NCHAINS NOT A NUMBER"; exit -1 ;;
esac

# Fifth argument is number of steps
NSTEPS=$5
case $NSTEPS in
  ''|[!0-9]*) echo "N_STEPS = $NSTEPS NOT A NUMBER"; exit -1 ;;
esac

# Fifth argument is number of steps
JOB=$6
case $JOB in
  ''|[!0-9]*) echo "JOB = $JOB NOT A NUMBER"; exit -1 ;;
esac

# Number of threads to use for multi-threading jobs
# Change these if you want!
NTHREADS=32
RAM=14
NGPU=0
SEC_PER_STEP=1.0
CHAINS_PER_NODE=8
NNODES=$(echo "${NCHAINS} / ${CHAINS_PER_NODE}" | bc)
#NNODES=1

CPU_PER_TASK=32


###########################################################
# DO NOT EDIT BELOW HERE
###########################################################

# Different settings for the different clusters
ScratchDir="$SCRATCH/"
##################

echo "---------------"
echo "Job Submission for Perlmutter"

# Check the walltime and calculate a new one
##################
echo "You've given me $NTHREADS CPU threads and ${NGPU} GPUs with $NSTEPS"
echo "---------------"
echo "Setting walltime automatically, assuming $SEC_PER_STEP s/step and 20 min start-up..."

# Calculate recommended walltime
WALLTIME=$(echo "${SEC_PER_STEP} * ${NSTEPS}" | bc)
# Takes about 3hr to set up samples too, just to be sure
WALLTIME=$(echo "${WALLTIME} + 6000" | bc)
# Convert both to ints
WALLTIME_INT=$(printf "%.0f" "$WALLTIME")

# Convert to desired format for ComputeCanada (CC) and Emerald (EM)
# ComputeCanada wants in format HH:MM:SS
# Emerald wants in format HH:MM
WALLTIME_CC=$(SecToH ${WALLTIME_INT})
WALLTIME_EM="${WALLTIME_CC%*:[0-9]*}"

#Print to user
echo "Walltime in seconds: $WALLTIME"
echo "Walltime in hh:mm:s: $WALLTIME_CC"

# Convert the walltime to standard date, needed for qsub/bsub
# Emerald wants hh:mm, Compute Canada wants hh:mm:ss, Imperial same as Compute Canada
# Now use the SecToH converter

# Check the hour of the WALLTIME
WALLTIME_CHECK=${WALLTIME_CC:0:1}

# Find MaCh3 folder and setup script
# Assuming we're running in MaCh3 folder
##################
MaCh3Dir=$(pwd -P)
echo $MaCh3Dir | grep -q 'MaCh3'
greprc=$?
if [ $greprc -ne 0 ]; then
  echo "Coulnd't find MaCh3 in current pwd"
  echo "pwd = $(pwd -P)"
  echo "Have to quit, sorry"
  exit -1
fi
MaCh3Dir=${MaCh3Dir%%/clusters}

##################
# Check that the config file exists
if [ ! -e ${CFG} ]; then
  echo "Did not find config file ${CFG}"
  exit -1
fi
##################


##################
# Output directory
OutputDir=${ScratchDir}/"Chains/"${FitName}


# Make a temp folder that will store all the intermediate scripts, configs and logs
scriptDir=${OutputDir}/scripts
cfgDir=${OutputDir}/cfg
logDir=${OutputDir}/log

# if this is the first job then make the Output Directories
if [[ $JOB -eq 0 ]]
then
  mkdir -p $OutputDir
  mkdir -p $scriptDir
  mkdir -p $cfgDir
  mkdir -p $logDir
fi
##################

declare -i ChainCounter=0
for ((b = 0 ; b < $NNODES ; b++)); do

  SubmitScriptFileName=${scriptDir}/${FitName}_batch_${b}_job_${JOB}.sh
  echo "#!/bin/bash" > $SubmitScriptFileName
  
  # Write the OMP_NUM_THREADS variable
  echo "export OMP_PROC_BIND=spread" >> $SubmitScriptFileName
  echo "export OMP_PLACES=threads" >> $SubmitScriptFileName
  echo "export OMP_NUM_THREADS=${NTHREADS}" >> $SubmitScriptFileName
  # Save the MaCh3 directory into the file
  echo "export MACH3=${MaCh3Dir}" >> $SubmitScriptFileName
  
  # Have the script cat the config so we have a log of the exact config file that was run
  echo "cat ${cfgFileName}" >> $SubmitScriptFileName
  echo "cd \${MACH3}" >> $SubmitScriptFileName
  # Load up cluster defaults
  echo "source setup.sh" >> $SubmitScriptFileName
  echo "source setup_dune_perl_micromamba.sh" >> $SubmitScriptFileName
  echo "PATH=../build/src:$PATH" >> $SubmitScriptFileName
  
  ##################
  # Loop over each job and make the config file and submission files
  for ((i = 0 ; i < $CHAINS_PER_NODE ; i++)); do
    JobName=${FitName}_chain_${ChainCounter}_job_${JOB}
    cfgFileName=${cfgDir}/${JobName}.yaml
    # Copy the "template" config file to new location
    cp $CFG $cfgFileName
  
    # Now over-ride the settings in the card file with user input
    # Replace the output name
    sed -i "s|^    FileName.*|    FileName: \"${OutputDir}/${JobName}.root\"|g" $cfgFileName
    # Replace the number of steps
    sed -i "s|^  NSTEPS.*|  NSTEPS: ${NSTEPS}|g" $cfgFileName
  
    #############################################
    #if you're continuing from the previous chain
    #############################################
    if [[ $JOB -gt 0 ]]
    then
      OLD_JOB=$(( $JOB-1 ))
      sed -i "s|^  StartFromPos.*|  StartFromPos: true|g" $cfgFileName
      OldJobName=${FitName}_chain_${ChainCounter}_job_${OLD_JOB}.root
      sed -i "s|^  PosFileName.*|  PosFileName: \"${OutputDir}/${OldJobName}\"|g" $cfgFileName
    fi
    
    # Write the temporary submission file
    executable="${EXE} ${cfgFileName}"
  
    stdout_chain="${logDir}/${JobName}.log"
    echo "srun --exclusive -n 1 ${executable} &> ${stdout_chain} &" >> $SubmitScriptFileName
    ((ChainCounter++))
  
  done

  chmod 744 $SubmitScriptFileName
  echo "wait" >> $SubmitScriptFileName

  qsubmit="sbatch"
  if [[ $JOB -gt 0 ]]
  then
    OLD_JOB=$(( $JOB-1 ))
    OldJobName=${FitName}_batch_${b}_job_${OLD_JOB}
    PrevJobPID=$(squeue -u lwarsame --format="%.18i %.9P %.60j %.8u %.8T %.10M %.9l %.6D %R" | grep "${OldJobName}" | awk '{print $1}')
    qsubopt="--account=dune --time=${WALLTIME_CC} --qos=regular --nodes=1 --cpus-per-task=${CPU_PER_TASK} --constraint=cpu --dependency=afterany:${PrevJobPID}"
    echo "This job will wait for job with name: ${OldJobName} with PID: "${PrevJobPID}""
  else
    qsubopt="--account=dune --time=${WALLTIME_CC} --qos=regular --nodes=1 --cpus-per-task=${CPU_PER_TASK} --constraint=cpu"
  fi
  
  stdout="--output ${logDir}/${FitName}_batch_${b}_job_${JOB}_%j.log"
  stderr="--error ${logDir}/${FitName}_batch_${b}_job_${JOB}_%j.err"
  qsubmit="${qsubmit} ${qsubopt} ${stdout} ${stderr} ${SubmitScriptFileName}"
  eval $qsubmit
  
done

##################

echo "All ${NCHAINS} chains submitted (This is Job number $JOB for these chains!)"
echo "Used ${CFG} template config, doing ${FitName} with ${NSTEPS} steps each"
echo "Submitted on ${NNODES} Nodes and ${CPU_PER_TASK} threads per Chain "
