#! /bin/bash

if [[ ! "$0" != "$BASH_SOURCE" ]]
then
  echo "setup_onlmon.sh needs to be sourced"
  exit
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
if [[ $ONLMON_SAVEDIR ]]
then
  if  ! -d $ONLMON_SAVEDIR ]
  then
    mkdir -p $ONLMON_SAVEDIR
  fi
fi

# create save directory if not exist
if [[ ! -z "$ONLMON_LOGDIR" ]]
then
  if [ ! -d $ONLMON_LOGDIR ]
  then
    mkdir -p $ONLMON_LOGDIR
  fi
else
  export ONLMON_LOGDIR=`pwd`
fi

if [[ -z "$ONLMON_HTMLDIR" ]]
then
  export ONLMON_HTMLDIR=/sphenix/WWW/subsystem/OnlMonHtmlTest
fi


# create history directory if not exist
if [[ ! -z "$ONLMON_HISTORYDIR" ]]
then
  if [ ! -d $ONLMON_HISTORYDIR ]
  then
    mkdir -p $ONLMON_HISTORYDIR
  fi
fi

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

source /opt/sphenix/core/bin/setup_local.sh $ONLMON_MAIN
# all subsystems scripts end in Setup.csh
for script in $ONLMON_BIN/*Setup.sh
do
  source $script
done
