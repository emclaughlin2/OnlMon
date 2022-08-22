#! /bin/csh

if (! $?ONLMON_MAIN ) then
  echo "ONLMON_MAIN environment variable not set, exiting"
  exit
endif

# set the macros directory to the current dir (they
# reside in the same directory as this setup script)
if (! $?ONLMON_MACROS) then
  setenv ONLMON_MACROS `pwd`
endif

# create save directory if not exist
if ($?ONLMON_SAVEDIR) then
  if (! -d $ONLMON_SAVEDIR) then
    mkdir -p $ONLMON_SAVEDIR
  endif
endif

# create save directory if not exist
if ($?ONLMON_LOGDIR) then
  if (! -d $ONLMON_LOGDIR) then
    mkdir -p $ONLMON_LOGDIR
  endif
else
  setenv ONLMON_LOGDIR `pwd`
endif

if (! $?ONLMON_HTMLDIR) then
  setenv ONLMON_HTMLDIR /sphenix/WWW/subsystem/OnlMonHtmlTest
endif


# create history directory if not exist
if ($?ONLMON_HISTORYDIR) then
  if (! -d $ONLMON_HISTORYDIR) then
    mkdir -p $ONLMON_HISTORYDIR
  endif
endif

#if (! $?ONLINE_LOG) then
#  setenv ONLINE_LOG /export/data1/log
#endif

#if (! $?ONLINE_CONFIGURATION) then
#  setenv ONLINE_CONFIGURATION /export/software/oncs/online_configuration
#endif

if (! $?ONLMON_BIN) then
  setenv ONLMON_BIN $ONLMON_MAIN/bin
endif

if (! $?ONLMON_CALIB) then
  setenv ONLMON_CALIB $ONLMON_MAIN/share/onlmon
endif

if (! $?ONLMON_RUNDIR) then
  setenv ONLMON_RUNDIR $ONLMON_MAIN/share
endif

source /opt/sphenix/core/bin/setup_local.csh $ONLMON_MAIN
# all subsystems scripts end in Setup.csh
foreach script ($ONLMON_BIN/*Setup.csh)
  source $script
end
