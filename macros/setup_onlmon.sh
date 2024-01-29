#! /bin/bash

if [[ ! "$0" != "$BASH_SOURCE" ]]
then
  echo "setup_onlmon.sh needs to be sourced"
  exit
fi

opt_n=0

for arg in "$@"
do
    case "$arg" in
    -n)
        opt_n=1
        ;;
    esac
done

if [[ $opt_n != 0 ]]
then
  tmponlmonmain=$ONLMON_MAIN
  unset ${!ONLMON_*}
  export ONLMON_MAIN=$tmponlmonmain
fi

if [[ -z "$ONLMON_MAIN" ]]
then
  echo "ONLMON_MAIN environment variable not set, doing nothing"
  return
fi

# set the macros directory to the current dir (they
# reside in the same directory as this setup script)
if [[ -z "$ONLMON_MACROS" ]]
then
  export ONLMON_MACROS=`pwd`
fi

# create save directory if not exist
if [ ! -n "$ONLMON_SAVEDIR" ]
then
#  [ -d $ONLMON_SAVEDIR ] || mkdir -p $ONLMON_SAVEDIR
#else
  export ONLMON_SAVEDIR=`pwd`
fi

# create save directory if not exist
if [ ! -n "$ONLMON_LOGDIR" ]
then
#  if [ ! -d $ONLMON_LOGDIR ]
#  then
#    mkdir -p $ONLMON_LOGDIR
#  fi
#else
  export ONLMON_LOGDIR=`pwd`
fi

if [[ -z "$ONLMON_HTMLDIR" ]]
then
  export ONLMON_HTMLDIR=/sphenix/WWW/subsystem/OnlMonHtmlTest
fi


# create history directory if not exist
#if [ -n "$ONLMON_HISTORYDIR" ]
#then
#  if [ ! -d $ONLMON_HISTORYDIR ]
#  then
#    mkdir -p $ONLMON_HISTORYDIR
#  fi
#fi

#if [[ -z "$ONLINE_LOG" ]] then
#  export ONLINE_LOG=/export/data1/log
#fi

#if [[ -z "$ONLINE_CONFIGURATION" ]] then
#  export ONLINE_CONFIGURATION=/export/software/oncs/online_configuration
#fi

if [[ -z "$ONLMON_BIN" ]]
then
  export ONLMON_BIN=$ONLMON_MAIN/bin
fi

if [[ -z "$ONLMON_CALIB" ]]
then
  export ONLMON_CALIB=$ONLMON_MAIN/share/onlmon
fi

if [[ -z "$ONLMON_RUNDIR" ]]
then
  export ONLMON_RUNDIR=$ONLMON_MAIN/share
fi

if [ -f ${OPT_SPHENIX}/bin/setup_local.sh ]
then
  source ${OPT_SPHENIX}/bin/setup_local.sh ${ONLMON_MAIN}
else
  ROOT_INCLUDE_PATH=./:$ONLINE_MAIN/include
  for local_incdir in $ONLINE_MAIN/include/*
    do
    if [ -d $local_incdir ]
    then
	ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:$local_incdir
    fi
  done
  export ROOT_INCLUDE_PATH=${ONLMON_MAIN}/include:${ROOT_INCLUDE_PATH}
  export LD_LIBRARY_PATH=${ONLMON_MAIN}/lib:$LD_LIBRARY_PATH
  export PATH=${ONLMON_MAIN}/bin:${PATH}
fi
# all subsystems scripts end in Setup.csh
for script in $ONLMON_BIN/*Setup.sh
do
  source $script
done
