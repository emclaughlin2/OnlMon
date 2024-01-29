#! /bin/csh
# no echoss it interferes with the monitor watcher
if (! $?TPOTCALIBREF ) then
  setenv TPOTCALIBREF /home/phnxrc/operations/TPOT/onlmon_ref
#  echo "TPOTCALIBREF environment variable set to $TPOTCALIBREF"
  exit
endif

setenv TPOTCALIB $ONLMON_CALIB/tpot
