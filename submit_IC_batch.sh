#!/bin/bash
#
# Submit job script, bash style

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
  ''|[!0-9]*) echo "N_CHAINS = $NCHAINS NOT A NUMBER"; exit -1 ;;
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
NTHREADS=8
RAM=5
NGPU=0
SEC_PER_STEP=0.4

IsICHEPQ=false
IsICHPC=true
###########################################################
# DO NOT EDIT BELOW HERE
###########################################################

# Different settings for the different clusters
ScratchDir="/vols/dune/ljw20/fit_results"
##################


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

# Check that the executable exists
#if [ ! -e ${EXE} ]; then
#  echo "Did not find executable ${EXE} in ${MaCh3Dir}..."
#  echo "Have you built MaCh3 yet?"
#  exit -1
#fi
##################


##################
# Check that the config file exists
if [ ! -e ${CFG} ]; then
  echo "Did not find config file ${CFG}"
  exit -1
fi
##################


##################
# Output directory
OutputDir=${ScratchDir}/${FitName}

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


##################
# Special treatment of picky task spooler at Imperial
# Essentially set the path, logging directory, socket and simultaneous slots
#if [[ $IsICHEP == true ]]; then
#  # Export the task spooler
#  export PATH=/vols/t2k/users/$USER/software/ts-1.0/bin:${PATH}
#  # Change the TMPDIR env variable which we write to
#  export TMPDIR=${logDir}
#  export TS_SOCKET=/tmp/ts.socket
#  export TS_SLOTS=1
#fi
##################


##################
# Loop over each job and make the config file and submission files
for ((i = 0 ; i < $NCHAINS ; i++)); do
  JobName=${FitName}_chain_${i}_job_${JOB}
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
    OldJobName=${FitName}_chain_${i}_job_${OLD_JOB}.root
    sed -i "s|^  PosFileName.*|  PosFileName: \"${OutputDir}/${OldJobName}\"|g" $cfgFileName
    #if [[ -e "${OutputDir}/${OldJobName}" ]]
    #then
    #  echo "Previous chain exists!"
    #  echo "${OutputDir}/${OldJobName}" 
    #else
    #  echo "Cannot find previous chain ${OutputDir}/${OldJobName}"
    #  exit -1
    #fi
  fi

  # Now make the .sh script that we submit to the cluster
  ScriptFileName=${scriptDir}/${JobName}.sh

  # Write the temporary submission file
  executable="${EXE} ${cfgFileName}"
  echo "#!/bin/bash" > $ScriptFileName

  # Write the OMP_NUM_THREADS variable
  echo "export OMP_NUM_THREADS=${NTHREADS}" >> $ScriptFileName
  # Save the MaCh3 directory into the file
  echo "export MACH3=${MaCh3Dir}" >> $ScriptFileName

  # Have the script cat the config so we have a log of the exact config file that was run
  echo "cat ${cfgFileName}" >> $ScriptFileName
  echo "cd \${MACH3}" >> $ScriptFileName
  # Load up cluster defaults
  echo "source setup.sh" >> $ScriptFileName
  echo "source setup_dune_env.sh" >> $ScriptFileName
  echo "PATH=../build/src:$PATH" >> $ScriptFileName
  # Run exec
  echo "${executable}" >> $ScriptFileName
  # Make executable
  chmod 744 $ScriptFileName

  # Now temp bash file is written and we can submit the contents
  stdout="-o ${logDir}/${JobName}.log"
  stderr="-e ${logDir}/${JobName}.err"

    # For Imperial lt2gpu00 we use qsub now
    # Needs to be submitted from GPU machine
  if [[ ${IsICHEPQ} == true ]]; then
    qsubmit="qsub"
	ICRAM=$(echo "${RAM} * ${NTHREADS} * 3" | bc)
    # Send to GPU queue, choose one GPU, set walltime to HH:MM:SS, request 4GB, request multi-thread
    qsubopt="-q gpu.q -l gpu=1 -l h_rt=${WALLTIME_CC} -l h_vmem=${RAM}G -pe hep.pe ${NTHREADS}"
    # Put together the submit command
    qsubmit="${qsubmit} ${qsubopt} ${stdout} ${stderr} ${ScriptFileName}"
    echo $qsubmit

  elif [[ ${IsICHPC} == true ]]; then
    qsubmit="qsub"
	ICRAM=$(echo "${RAM} * ${NTHREADS} * 3" | bc)
	if [[ $JOB -gt 0 ]]
	then
      OLD_JOB=$(( $JOB-1 ))
      OldJobName=${FitName}_chain_${i}_job_${OLD_JOB}
	  PrevPID=$(qstat -xml | tr '\n' ' ' | sed 's#<job_list[^>]*>#\n#g' | sed 's#<[^>]*>##g' | grep " " | column -t | grep "${OldJobName}" | awk '{print $1}')
      qsubopt="-q hep.q -l h_rt=2:50:00 -l h_vmem=20G -pe hep.pe ${NTHREADS} -hold_jid ${PrevPID}"
	else
      qsubopt="-q hep.q -l h_rt=2:50:00 -l h_vmem=20G -pe hep.pe ${NTHREADS} "
    fi
    # Can also give gpu_type=P100 for spanking new P100 cards
    qsubmit="${qsubmit} ${qsubopt} ${stdout} ${stderr} ${ScriptFileName}"
    eval $qsubmit

  fi

done
##################

echo "All ${NCHAINS} chains submitted (This is Job number $JOB for these chains!)"
echo "Used ${CFG} template config, doing ${FitName} with ${NSTEPS} steps each"
echo "Submitted with GPU and ${NTHREADS} CPUs"
