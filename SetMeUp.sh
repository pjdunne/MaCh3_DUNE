#!/bin/bash
#
# Script to set up MaCh3+dependencies on external clusters
#
# Tested on: UK STFC Emerald
#            UK STFC SCARF
#            ComputeCanada Guillimin, Hades, Helios, Graham, Cedar, Niagara, Beluga
#            Imperial High Energy Physics cluster
#            Imperial High Performance Computing
#       
#            !!! PLEASE ADD YOURS !!!
#
# Written by cdawg

home=$(pwd -P)
# Needed to get CMT, can use your own if you want or just stick with anoncvs
username=anoncvs

clear

if [[ $home == *MaCh3* ]]; then
  echo -e "Looks like you're running this in a MaCh3 folder"
  echo -e "I'll cd .. and install programs there\n"
  cd ..
  oldhome=${home}
  home=$(pwd -P)
else
  echo -e "Please run this in the MaCh3 folder"
  echo -e "You're currently running in $home"
  exit -1
fi

echo -e "\e[4;1mInstalling to $home\n\e[0m"
sleep 1

echo -e "\e[1m"
# The user specified installs

# Install ROOT?
userROOT=false
while true; do
  read -p "Build ROOT? (y/n)" yn
  case $yn in
    [Yy]* ) userROOT=true; break;;
    [Nn]* ) userROOT=false; break;;
    * ) echo "Please answer...";;
  esac
done

# Install CMT?
userCMT=true
while true; do
  read -p "Build CMT (needed for psyche)? (y/n)" yn
  case $yn in 
    [Yy]* ) userCMT=true; break;;
    [Nn]* ) userCMT=false; break;;
    * ) echo "Please answer...";;
  esac
done

# Install iRODS?
userIRODS=true
while true; do
  read -p "Build iRODS and get files? (y/n)" yn
  case $yn in 
    [Yy]* ) userIRODS=true; break;;
    [Nn]* ) userIRODS=false; break;;
    * ) echo "Please answer...";;
  esac
done

# Setup NIWGReWeight?
userNIWG=true
while true; do
  read -p "Setup NIWGReWeight? (y/n)" yn
  case $yn in
    [Yy]* ) userNIWG=true; break;;
    [Nn]* ) userNIWG=false; break;;
    * ) echo "Please answer...";;
  esac
done
echo -e "\e[0m"

echo -e "\n\e[4;1mThis will build\e[0m: ROOT  \e[1m${userROOT}\e[0m"
echo -e "                 CMT   \e[1m${userCMT}\e[0m"
echo -e "                 iRODS \e[1m${userIRODS}\e[0m"
echo -e "                 NIWG \e[1m${userNIWG}\e[0m\n"
echo -e "I will also need to find GSL, so if you don't have that please modify ${oldhome}/setup.sh"
echo -e "Will try to append the resulting variables to your setup.sh and setup_psyhe.sh\n"
sleep 1

# Find which cluster we're building on
cluster="NULL"
# Emerald
if [[ $HOSTNAME == *gpu.rl.ac.uk ]]; then
  echo -e "\e[1mFound Emerald cluster\e[0m"
  cluster="Emerald"
  cmake="cmake/2.8.12"
  cuda="cuda/8.0.44"

# Guilliumin on Compute Canada
elif [[ $HOSTNAME == lg-1r1[47]-n0[1-4] ]]; then
  echo -e "\e[1mFound Compute Canada cluster\e[0m"
  cluster="ComputeCanada"
  cmake="cmake/2.8.12.2"
  cuda="CUDA/7.5.18"
  extras="module load GSL"

# Imperial HPC
elif [[ $HOSTNAME == *[0-9]-internal ]]; then
  echo -e "\e[1mFound Imperial High Performace Computing cluster\e[0m"
  cluster="ICHPC"
  cmake="cmake/2.8.12.2"
  cuda="cuda/8.0.44"
  extras="module load gsl/2.1 \nexport PATH=/apps/mira/3rc2/build/flex-2.5.35:\${PATH}"

# Imperial HEP
elif [[ $HOSTNAME == @(lt2gpu00|heppc105|heppc205).hep.ph.ic.ac.uk ]]; then
  echo -e "\e[1mFound Imperial HEP cluster\e[0m"
  cluster="ICHEP"
  extras="export PATH=/vols/t2k/users/cvw09/software/ts-1.0/bin:\${PATH} \n\
          export TS_SOCKET=/tmp/ts.socket \n\
          export TS_SLOTS=1 \n\
          export CUDA_HOME=/usr/local/cuda-8.0"

# Helios on Compute Canada
elif [[ $HOSTNAME == helios[0-9] ]]; then
  echo -e "\e[1mFound Helios Computing cluster\e[0m"
  cluster="Helios"
  cmake="apps/cmake/2.8.12.1"
  cuda="cuda/8.0.44"
  extras="export PATH=/software6/libs/gsl/1.16_gcc/bin:\${PATH} \nexport LD_LIBRARY_PATH=/software6/libs/gsl/1.16_gcc/lib:\${LD_LIBRARY_PATH}"

# Hades on Compute Canada
elif [[ $HOSTNAME == hades* ]]; then
  echo -e "\e[1mFound Hades Computing cluster\e[0m"
  cluster="Hades"
  cuda="CUDA/7.5"
  extras="export CUDA_HOME=/usr/local/cuda-7.5 \nmodule load GSL/2.1 \nexport CPATH=/usr/local/cuda-7.5/include:\${CPATH} \nexport PATH=/usr/local/cuda-7.5/open64/bin:/usr/local/cuda-7.5/bin:\${PATH}"

# Graham and Cedar on Compute Canada
elif [[ $HOSTNAME == cedar* || $HOSTNAME == gra-login* || $HOSTNAME == beluga[0-9]* ]]; then
  echo -e "\e[1mFound Cedar, Graham or Beluga Computing cluster\e[0m"
  cluster="Cedar"
  # CMake is already loaded by default on Cedar
  cmake="cmake"
  cuda="cuda"
# Same thing for gcc because this determines the gfortran version
  extras="module load gcc/4.8.5 \nmodule load gsl/1.16"
  module load gcc/4.8.5
  module load gsl/1.16
# Also silly things for ROOT to build!
  rootextras="-DCMAKE_CXX_COMPILER=$(which g++) -DCMAKE_C_COMPILER=$(which gcc) -DCMAKE_Fortran_COMPILER=$(which gfortran)"

# Niagara on Compute Canada (no GPU)
# This cluster needs a lot of hand holding
elif [[ $HOSTNAME == nia-login[0-9]* ]]; then
  echo -e "\e[1mFound Niagara Computing cluster\e[0m"
  cluster="Niagara"
  cmake="cmake"
  cuda=""
  extras="module load gcc/4.9.4"
  module load gcc/4.9.4
  rootextras="-DCMAKE_CXX_COMPILER=$(which g++) -DCMAKE_C_COMPILER=$(which gcc) -DCMAKE_Fortran_COMPILER=$(which gfortran)"

else
  echo -e "\e[1mDidn't find pre-set cluster\e[0m"
  echo -e "\e[1mThis might still work though\e[0m"
  cuda=""
  extras=""
  cmake=""
fi

echo -e "\n"
sleep 1

if [ "${userROOT}" == true ]; then
echo -e "\e[4;1mROOT:\e[0m"
# Check if ROOT is built
buildroot=false
if [ ! -e ${home}/root/bin/thisroot.sh ]; then
  buildroot=true
else
  buildroot=false
  source ${home}/root/bin/thisroot.sh
fi

if [ "${buildroot}" == true ]; then
  echo "Didn't find ROOTSYS, will build ROOT..."
  sleep 1
  # Get ROOT 5.34.34 directly from ROOT rather than CMT
  ROOTVER=root_v5.34.34.source.tar.gz
  if [ ! -e ${ROOTVER} ]; then
    echo -e "\e[1mGetting ROOT ${ROOTVER}\e[0m"
    wget https://root.cern.ch/download/${ROOTVER}
  else 
    echo -e "\e[1mAlready found a ROOT tarball in $ROOTVER, using that instead\e[0m"
  fi

  # Move ROOT source into directory
  if [ ! -d ${home}/root-source ]; then
    echo -e "\e[1mUntarring ROOT install to ${home}/root-source\e[0m"
    mkdir -pv ${home}/root-source
    echo -e "\e[1mPlease be patient, untarring ROOT takes a while...\e[0m"
    tar -xf ${ROOTVER} -C ${home}/root-source
    echo -e "\e[1mUntarred ROOT\e[0m"
  else
    echo -e "\e[1mFound ${home}/root-source directory already\e[0m"
  fi

  echo -e "\e[1mMaking root in ${home}/root with Minuit2, ROOfit, pyROOT with 8 threads\e[0m"
  sleep 2
  mkdir -pv ${home}/root
  echo "Loading cmake module..."

  # Load cmake in different ways
  if ! command -v cmake > /dev/null; then
    echo "Trying to load saved and generic cmake version..."
    module load ${cmake}
    module load cmake
  fi

  # Check if we now have cmake
  CMAKEbuild=true
  if ! command -v cmake > /dev/null; then
    echo "Did not find cmake, using configure instead"
    CMAKEbuild=false
  fi

  # Do the fast ROOT cmake build
  # Documentation for options at https://root.cern.ch/building-root
  # These have been tested, so please don't randomly remove them expecting awesomeness
  if [[ "$CMAKEbuild" == true ]]; then
    echo -e "\e[1mBuilding ROOT with cmake...\e[0m"
    sleep 1
    cd ${home}/root
    cmake ../root-source/root -Dminuit2:BOOL=ON \
                              -Droofit:BOOL=ON \
                              -Dfortran:BOOL=ON \
                              -Dpython:BOOL=OFF \
                              -Dx11:BOOL=OFF \
                              -Dmysql:BOOL=OFF \
                              -Dalien:BOOL=OFF \
                              -Dxml:BOOL=OFF \
                              -Dalien:BOOL=OFF \
                              -Dastiff:BOOL=OFF \
                              -Dbonjour:BOOL=OFF \
                              -Dbuiltin_gsl:BOOL=OFF \
                              -Dopengl:BOOL=OFF \
                              -Dglite:BOOL=OFF \
                              -Dldap:BOOL=OFF \
                              -Dmonalisa:BOOL=OFF \
                              -Dmysql:BOOL=OFF \
                              -Doracle:BOOL=OFF \
                              -Dpgsql:BOOL=OFF \
                              -Dsapdb:BOOL=OFF \
                              -Dsqlite:BOOL=OFF \
                              -Dxft:BOOL=OFF \
                              ${rootextras}
    cmake --build . -- -j8

  # Do the convetional ./configure build
  # Still want to build everything in root directory though
  else
    echo -e "\e[1mBuilding ROOT with configure...\e[0m"
    # Rename root-source to root and build in there
    mv ${home}/root-source/root ${home}
    cd ${home}/root
    ./configure --enable-minuit2 --enable-roofit --enable-python --disable-x11 --disable-mysql
    make
  fi

  # Check if build was successful
  if [[ $? -ne 0 ]]; then
    echo -e "\e[1mFailed ROOT install!\e[0m"
    exit -1
  fi

  echo -e "\e[1mFinished ROOT install!\e[0m"
  sleep 2

  # Source the ROOT environment
  echo "Sourcing ROOT environment..."
  source ${home}/root/bin/thisroot.sh
else 
  echo -e "Found existing ROOT install, skipping getting and installing"
fi
fi
echo -e "Finished with ROOT...\n"

# Check if we wanted CMT (otherwise no point in getting procmail)
if [ "${userCMT}" == true ]; then
echo -e "\e[4;1mprocmail:\e[0m"
sleep 1

# Set up procmail
getprocmail=false
if command -v procmail > /dev/null; then
   echo -e "Found procmail"
   getprocmail=false
 else
   echo -e "Did not find procmail, will check directories"
   getprocmail=true
 fi

if [ "${getprocmail}" == true ]; then
  PROCVER=procmail-3.22.tar.gz
  if [ ! -e ${PROCVER} ]; then
    echo -e "Getting procmail, required by CMT (lockfile)..."
    wget http://www.ring.gr.jp/archives/net/mail/procmail/${PROCVER}
  else
    echo -e "Found procmail tarball already, not downloading..."
  fi

  PROCMAIL="${PROCVER%%.tar.gz}"
  echo $PROCMAIL
  if [ ! -d ${PROCMAIL} ];then
    echo -e "Untarring procmail..."
    mkdir -pv ${home}/${PROCMAIL}
    tar -xf ${PROCVER}
    echo -e "Untarred procmail"
  else
    echo -e "Already had untarred procmail..."
  fi

  echo "Installing procmail..."
  if [ ! -e ${PROCMAIL}/new/procmail ]; then
    # Replace the bloody getline in formisc, fields and formail
    cd ${PROCMAIL}
    sed -i 's/getline/getline1/g' src/formisc.h src/formisc.c src/fields.c src/formail.c
    make
    # Add procmail to path
    cd ${home}
  fi
fi

if [[ $? -ne 0 ]]; then
  echo -e "\e[1mFailed PROCMAIL install!\e[0m"
  exit -1
fi

echo -e "Procmail installed, sourcing...\n"
export PATH=${PATH}:${PROCMAIL}/new

# Now handhold some more, check cvs and GSL (needed for e.g. Niagara)
# GSL (don't use ROOT for this because it points to broken mirrors)
cd ${home}
if command -v gsl-config > /dev/null; then
  echo "Found gsl-config in $(which gsl-config)"
  build_GSL=false
else
  echo "Did not find gsl-config, trying to build GSL"; 
  build_GSL=true;
fi

if [[ ${build_GSL} == true ]]; then
  if [ ! -e gsl-1.16.tar.gz ]; then
    echo -e "Getting gsl tarball..."
    wget https://ftp.gnu.org/gnu/gsl/gsl-1.16.tar.gz
  else
    echo -e "Found gsl tarball, not downloading"
  fi

  if [ ! -d gsl-1.16 ]; then
    echo -e "Untarring gsl..."
    tar -xf gsl-1.16.tar.gz
  else 
    echo -e "Already had untarred gsl"
  fi

  echo "Installing gsl..."
  cd gsl-1.16
  if [ ! -e Linux/bin/gsl-config ]; then
    ./configure --prefix=$(readlink -f Linux)
    make && make install
  else
    echo -e "Already had gsl-config"
  fi

  if [[ $? -ne 0 ]]; then
    echo -e "\e[1mFailed GSL install!\e[0m"
    exit -1
  fi

  gslpath=$(pwd -P)
  export PATH=${gslpath}/Linux/bin:${PATH}
  export LD_LIBRARY_PATH=${gslpath}/Linux/lib:${LD_LIBRARY_PATH}
  # Append to extras
  extras="${extras} \nexport PATH=${gslpath}/Linux/bin:\${PATH}\nexport LD_LIBRARY_PATH=${gslpath}/Linux/lib:\${LD_LIBRARY_PATH}"
  cd ${home}
fi

# CVS
cd ${home}
if command -v cvs > /dev/null; then
  echo "Found cvs in $(which cvs)"
  build_CVS=false
else
  echo "Did not find cvs, trying to build CVS"; 
  build_CVS=true;
fi

if [[ ${build_CVS} == true ]]; then
  if [ ! -e cvs-1.11.23.tar.gz ]; then
    echo -e "Getting cvs tarball..."
    wget https://ftp.gnu.org/non-gnu/cvs/source/stable/1.11.23/cvs-1.11.23.tar.gz
  else 
    echo -e "Found cvs tarball"
  fi

  if [ ! -d cvs-1.11.23 ]; then
    echo -e "Untarring cvs..."
    tar -xf cvs-1.11.23.tar.gz
  else
    echo -e "Already untarred cvs"
  fi

  echo "Installing cvs..."
  cd cvs-1.11.23
  if [ ! -e Linux/bin/cvs ]; then
    ./configure --prefix=$(readlink -f Linux)
    # This sed is needed to use the right getline
    sed -i 's/getline /get_line /' lib/getline.{c,h}
    make && make install
  else
    echo -e "Already had cvs"
  fi

  if [[ $? -ne 0 ]]; then
    echo -e "\e[1mFailed CVS install!\e[0m"
    exit -1
  fi
  cvspath=$(pwd -P)
  export PATH=${cvspath}/Linux/bin:${PATH}
  extras="${extras} \nexport PATH=${cvspath}/Linux/bin:\${PATH}"
  cd ${home}
fi

cd ${home}

## Install CMT
echo -e "\e[4;1mCMT:\e[0m"
sleep 1


# Set CMT version name
CMTNAME=CMTv1r20p20081118.tar.gz
CMTVER="${CMTNAME%%.tar.gz}"
CMTVER="${CMTVER##CMT}"

# First check if CMT has already been installed
if [ -e ${home}/CMT/${CMTVER}/Linux-x86_64/cmt.exe ]; then
    echo "CMT found! Sourcing setup script and moving on"
    source ${home}/CMT/${CMTVER}/mgr/setup.sh
else
    # If a ${home}/CMT folder doesn't exist, check if tarball does
    if [ ! -d ${home}/CMT ]; then

	echo "Checking if tarball exists"

	# If no tarball, download and untarrit
	if [ ! -e ${home}/${CMTNAME} ]; then
	    echo "Getting CMT from T2K software repo"
	    export CVSROOT=:ext:${username}@repo.nd280.org:/home/trt2kmgr/ND280Repository
	    export CVS_RSH=ssh
	    unset CVS_SERVER
	    cvs co CMT/${CMTNAME}
	    mv CMT/${CMTNAME} ${home}
	    rm -rf CMT
	else
	    echo "Found existing CMT tarball"
	fi

	## Untar it
	cd ${home}
	tar -xf ${CMTNAME}
    fi

    ## Now compile CMT
    cd ${home}/CMT/${CMTVER}/mgr
    ./INSTALL
    source setup.sh
    make
    if [[ $? -ne 0 ]]; then
      echo -e "\e[1mFailed CMT install!\e[0m"
      exit -1
    fi
    cp -v setup.sh ${home}/CMT/setup.sh
fi

cd ${home}
echo -e "Built CMT and sourced environment\n"
fi

# Now we can also append stuff in our setup.sh and setup_psyche.sh file
# Like procmail, ROOT, CMT, GSL, hurray!
echo -e "\e[4;1mAppending environment to ${oldhome}/setup.sh and ${oldhome}/setup_psyche.sh\n\e[0m"
# Chuck in the words
# procmail, CMT, ROOT, CUDA
sed -i "/\#\!\/bin\/bash/!b;n;c${extras} \n\
export PATH=${home}/procmail-3.22/new:\${PATH} \n\
source ${home}\/CMT\/setup.sh \n\
source ${home}/root/bin/thisroot.sh \n\
module load ${cuda} \n\
export CUDAPATH=\${CUDA_HOME} \n\
export MACH3_DATA=${home}\/P6Data \n\
export MACH3_MC=${home}\/P6MC" ${oldhome}/setup.sh ${oldhome}/setup_psyche.sh

if [ "${userIRODS}" == true ]; then
echo -e "\e[4;1miRODS:\e[0m\n"
sleep 1

# Write the iRODS file for T2K_ASG_Reader, see http://www.t2k.org/asg/oagroup/gadatastorage/.irodsEnv
function WriteIRODSFile() {

echo "# iRODS server host name:" > .irodsEnv
echo "irodsHost 'hepirods2.ph.qmul.ac.uk'" >> .irodsEnv
echo "# iRODS server port number:" >> .irodsEnv
echo "irodsPort 6835" >> .irodsEnv
echo "# Default storage resource name:" >> .irodsEnv
echo "irodsDefResource 't2kdata1'" >> .irodsEnv
echo "# Home directory in iRODS:" >> .irodsEnv
echo "irodsHome '/QMULZone1/home/asg'" >> .irodsEnv
echo "# Current directory in iRODS:" >> .irodsEnv
echo "irodsCwd '/QMULZone1/home/asg'" >> .irodsEnv
echo "# Account name:" >> .irodsEnv
echo "irodsUserName 'T2K_ASG_Reader'" >> .irodsEnv
echo "# Zone:" >> .irodsEnv
echo "irodsZone 'QMULZone1'" >> .irodsEnv

}

# Also need iRODS to get the files (requires user input so put last)
# http://www.t2k.org/asg/oagroup/gadatastorage/index_html
getirods=false
if [ ! -e ${home}/irods-legacy/iRODS/clients/icommands/bin/ils ]; then
  getirods=true
fi

if [ "$getirods" = true ]; then
  if [ ! -d ${home}/irods-legacy ]; then
    git clone https://github.com/irods/irods-legacy.git
    cd irods-legacy
    git checkout tags/3.2
  else
    cd irods-legacy
  fi

  cd iRODS
  IRODSDIR=$(pwd -P)
  ./irodssetup
  cd ${HOME}
  mkdir -pv .irods
  cd .irods
  # Get the T2K iRODS environment
  # Yuck, this isn't hosted anywhere so literally have to write the file here!
  # This isn't my fault that it isn't encrypted!
  WriteIRODSFile
  cd ${IRODSDIR}
  # Add clients
  echo -e "\e[1;4;92;101mSee http://www.t2k.org/asg/oagroup/gadatastorage/index_html for password\e[0m"
  source add-clients.sh
  export PATH=${PATH}:${home}/irods-legacy/iRODS/clients/icommands/bin
  iinit
fi

# Get the data from iRODS
if [ ! -d ${home}/P6Data ]; then
  echo "iRODS set up, now adding to PATH"
  cd ${home}
  mkdir -pv P6Data
  cd P6Data
  echo "Getting April 24th ND280 data from iRODS..."
  igetwild.sh /QMULZone1/home/asg/asg2016oa_winter/BANFF/prefit/Splines/Data Apr24th m
fi

if [[ ! -d ${home}/P6MC ]]; then
  cd ${home}
  mkdir -pv P6MC
  cd P6MC
  echo "Getting ND280 May15 MC from iRODS..."
  igetwild.sh /QMULZone1/home/asg/asg2016oa_winter/BANFF/prefit/Splines/MC May15th m
fi

fi


echo -e "\e[4;1mDone!\e[0m"
