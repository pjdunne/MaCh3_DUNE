#!/bin/bash

if [[ $HOSTNAME == *heppc105.hep.ph.ic.ac.uk ]];then
  MACHINE="HEPPC105"
  NUM_THREADS=3
elif [[ $HOSTNAME == *lxgpu04.hep.ph.ic.ac.uk ]]; then
  MACHINE="lxgpu04"
  NUM_THREADS=1
elif [[ $HOSTNAME == *lxgpu05.hep.ph.ic.ac.uk ]]; then
  MACHINE="lxgpu05"
  NUM_THREADS=1
elif [[ $HOSTNAME ==  *heppc205.hep.ph.ic.ac.uk ]];then
  MACHINE="HEPPC205"
  NUM_THREADS=1
fi

#Give your resource a name to put in the file
RESOURCE="ICLGPU"

#the day you're submitting the job (nice to just keep track of things)
DATE=$(date "+%x")
DATE=${DATE//\/}
MACH3_DIRECTORY="/vols/t2k/users/ljw20/software/MaCh3_DUNE/shiny_new/MaCh3_DUNE"
#what is the name of the executable?
exe="jointFitDUNE"
#the number of steps you want, this will get used in sed
NSTEPS=100
#OLD DATE 
#OLD_DATE=260922 #use this when continuing from a previous chain so you can find the chain name

#Are you continueing from an old chain?
#If first chain then JOB=0 otherwise increment this!
JOB=0
#Are you running with or without reactor constraints?

FOLDER="DUNE_FD_fits"
SETUPSCRIPTNAME="setup.sh"
original_config="mach3dune_config.yaml"
#What type of fit are you running? Give it a nice clear name!
TYPE="DUNE_FD_thread"

#Name of the directory you want to write things out to
#Again, having a date here is quite useful
INITIALDATE=220223

OUTDIRNAME=MaCh3_${TYPE}_${INITIALDATE}

#########################
#source setup script here
#########################

if [[ ! -e /vols/t2k/users/ljw20/fit_results/${FOLDER}/${OUTDIRNAME} ]]
then
  mkdir -p /vols/t2k/users/ljw20/fit_results/${FOLDER}/${OUTDIRNAME}
fi


#Specify the number of jobs you want to submit in the for loop
#Have to change this manually
for chain in {1..20}
do

  #Name of the script that you will then submit
  #script="MaCh3-t2knovaAsimov1_NOvAlike_FDS_${RESOURCE}_${MACHINE}_${DATE}_chain${chain}_job${JOB}.sh"
  script="MaCh3-${exe}_${RESOURCE}_${MACHINE}_${DATE}_chain${chain}_job${JOB}.sh"

  #sleep for a random time so jobs don't pile up on each other
  #echo "sleep $[ ( $RANDOM % 480 ) +1 ]s" > $script 

  echo "Found Date to be $DATE" 
  DATE=${DATE///}
  echo "Date is now $DATE"
  #Do you need to do some setup?
  echo "cd $MACH3_DIRECTORY" >> $script
  #Turn off core dumps...
  #(I never use them)...
  echo "ulimit -c 0" >> $script
  echo "source ${SETUPSCRIPTNAME}" >> $script
  
  #Make a new config file to run using, good to always have a config file associated with each chain so you
  # can double check that all the correct options were used if you find something weird
  new_config="MaCh3_DUNE_${TYPE}_${RESOURCE}_${MACHINE}_${DATE}_chain${chain}_job${JOB}.yaml"

  echo "cp configs/${original_config} configs/${new_config}" >> $script
  echo "export OMP_NUM_THREADS=${NUM_THREADS}" >> $script

  #Replace the outputname
  echo "sed -i \"s|^    FileName.*|    FileName: \"\\\"/vols/t2k/users/ljw20/fit_results/${FOLDER}/${OUTDIRNAME}/MaCh3_DUNEfit_${TYPE}_${DATE}_${RESOURCE}_${MACHINE}_chain${chain}_job${JOB}.root\\\"\"|g\" configs/$new_config" >> $script

  # Replace the number of steps
  echo "sed -i \"s|^  NSTEPS.*|  NSTEPS: ${NSTEPS}|g\" configs/$new_config" >> $script

  #############################################
  #if you're continuing from the previous chain
  #############################################
  if [[ $JOB -gt 0 ]]
  then
	OLD_JOB=$(( $JOB-1 ))
	echo "sed -i 's/#STARTFROMPOS/STARTFROMPOS/' configs/$new_config" >> $script
	echo "sed -i 's/MaCh3-t2knovajointAsimov2020_startfrompos.root/\/vols\/t2k\/users\/ea2817\/T2K_NOvA_chains\/${FOLDER}\/${OUTDIRNAME}\/MaCh3-t2knova_${TYPE}_${RC}_${OLD_DATE}_${RESOURCE}_${MACHINE}_chain${chain}_job${OLD_JOB}.root/' configs/$new_config" >> $script
	if [[ -e "/vols/t2k/users/ea2817/T2K_NOvA_chains/${FOLDER}/${OUTDIRNAME}/MaCh3-t2knova_${TYPE}_${RC}_${OLD_DATE}_${RESOURCE}_${MACHINE}_chain${chain}_job${OLD_JOB}.root" ]]
	then
	  echo "Previous chain exists!"
	  echo "/vols/t2k/users/ea2817/T2K_NOvA_chains/${FOLDER}/${OUTDIRNAME}/MaCh3-t2knova_${TYPE}_${RC}_${OLD_DATE}_${RESOURCE}_${MACHINE}_chain${chain}_job${OLD_JOB}.root" 
	else
	  echo "Cannot find previous chain /vols/t2k/users/ea2817/T2K_NOvA_chains/${FOLDER}/${OUTDIRNAME}/MaCh3-t2knova_${TYPE}_${RC}_${OLD_DATE}_${RESOURCE}_${MACHINE}_chain${chain}_job${OLD_JOB}.root"
	  return 1
	fi
  fi
 
  echo "$exe configs/$new_config &> /vols/t2k/users/ljw20/fit_results/${FOLDER}/${OUTDIRNAME}/MaCh3_DUNEfit_${TYPE}_${DATE}_${RESOURCE}_${MACHINE}_chain${chain}_job${JOB}.txt" >> $script

  #now submit the chain!
  sub=". $script &> /dev/null"
  echo "submiited job locally to $RESOURCE"
  echo "piping output to /vols/t2k/users/ljw20/${FOLDER}/${OUTDIRNAME}/MaCh3_DUNEfit_${TYPE}_${DATE}_${RESOURCE}_${MACHINE}_chain${chain}_job${JOB}.txt"
  eval $sub
  ((NUM_THREADS++))

done
