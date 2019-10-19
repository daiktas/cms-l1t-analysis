#!/usr/bin/env bash
function find_os_version()
{
  try_release=$(lsb_release -rs 2> /dev/null || echo 'not found')
  if [ "${try_release}" != 'not found' ]
  then
    echo $(echo "${try_release}" | cut -f1 -d.)
    return 0
  fi

  try_release=$(rpm -q --queryformat '%{RELEASE}' rpm 2> /dev/null || echo 'not found')
  if [ "${try_release}" != 'not found' ]
  then
    echo $(echo "${try_release}" | grep -o [[:digit:]]*\$)
    return 0
  fi

  try_release=$(cat /etc/redhat-release 2> /dev/null || echo 'not found')
  if [ "${try_release}" != 'not found' ]
  then
    echo $(echo "${try_release}" | grep -o '[[:digit:]]' | head -1)
    return 0
  else
    echo "Unable to find current OS version, aborting"
    return -1
  fi
}

OS_VERSION=centos7
os_ver=$(find_os_version)
if [ "${os_ver}" == "6" ]
then
  OS_VERSION=slc6
fi
export OS_VERSION

#LCG_VERSION=LCG_95
#LCG_ARCH=x86_64-${OS_VERSION}-gcc8-opt
LCG_VERSION=LCG_87
LCG_ARCH=x86_64-${OS_VERSION}-gcc62-opt
export LCG_VERSION
export LCG_ARCH

PROJECT_NAME="cms-l1t-analysis"

if [ -n "${PROJECT_ROOT}" ] ; then
   old_projectbase=${PROJECT_ROOT}
fi

if [ "x${BASH_ARGV[0]}" = "x" ]; then
    if [ ! -f setup.sh ]; then
        echo ERROR: must "cd where/${PROJECT_NAME}/is" before calling ". setup.sh" for this version of bash!
        PROJECT_ROOT=; export PROJECT_ROOT
        return 1
    fi
    PROJECT_ROOT="$PWD"; export PROJECT_ROOT
else
    # get param to "."
    envscript=$(dirname ${BASH_ARGV[0]})
    PROJECT_ROOT=$(cd ${envscript};pwd); export PROJECT_ROOT
fi

# clean PATH and PYTHONPATH
if [ -n "${old_projectbase}" ] ; then
  PATH=`python ${PROJECT_ROOT}/bin/remove_from_path.py "$PATH" "${old_projectbase}"`
  PYTHONPATH=`python ${PROJECT_ROOT}/bin/remove_from_path.py "$PYTHONPATH" "${old_projectbase}"`
fi

# add project to PYTHONPATH
if [ -z "${PYTHONPATH}" ]; then
   PYTHONPATH=${PROJECT_ROOT}; export PYTHONPATH
else
   PYTHONPATH=${PROJECT_ROOT}:$PYTHONPATH; export PYTHONPATH
fi

# add project to PATH
PATH=${PROJECT_ROOT}/bin:$PATH; export PATH
# add local bin to PATH (for local flake8 installations, etc)
export PATH=~/.local/bin:$PATH

unset old_projectbase
unset envscript

if [[ -z "${NO_CVMFS}" ]]
then
  echo "Getting dependencies from CVMFS"
  # this gives you voms-proxy-*, xrdcp and other grid tools
  source /cvmfs/grid.cern.ch/etc/profile.d/setup-cvmfs-ui.sh
  # ROOT 6, voms-proxy-init and other things
  source /cvmfs/sft.cern.ch/lcg/views/${LCG_VERSION}/${LCG_ARCH}/setup.sh
  # to fix java for the hadoop commands:
  unset JAVA_HOME
  pip install --user -r requirements.txt
else
  echo "No CVMFS available, setting up Anaconda Python"
  if [ ! -d "${CMSL1T_CONDA_PATH}" ]
  then
    wget -nv https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh -O /tmp/miniconda.sh
    bash /tmp/miniconda.sh -b -p ${CMSL1T_CONDA_PATH}
    PATH=${CMSL1T_CONDA_PATH}/bin:$PATH; export PATH
    rm -f miniconda.sh
    echo "Finished conda installation, creating new conda environment"
    conda config --add channels http://conda.anaconda.org/NLeSC
    conda config --set show_channel_urls yes
    conda create -n cms python=2.7 -yq
    source activate cms
    conda install -y -q \
      matplotlib \
      numpy \
      root>=6.04 \
      rootpy
    echo "Created conda environment, installing basic dependencies"
    pip install -U pip
    pip install -r requirements.txt
    conda clean -t -y
  fi
  source activate cms
fi

# Capture the user's site-packages directory:
USER_SITE_PACKAGES="$(python -c "import site; print(site.USER_SITE)")"
# add project to PYTHONPATH
PYTHONPATH="${USER_SITE_PACKAGES}:$PYTHONPATH"

git submodule init
git submodule update

ulimit -c 0

echo "Environment for ${PROJECT_NAME} is ready"
