# OnlMon

[![Build Status](https://web.sdcc.bnl.gov/jenkins-sphenix/buildStatus/icon?job=sPHENIX%2FsPHENIX_OnlMon_MasterBranch)](https://web.sdcc.bnl.gov/jenkins-sphenix/job/sPHENIX/job/sPHENIX_OnlMon_MasterBranch/)


To build it, create a build area and run

\<srcdir\>/autogen.sh --prefix=\<installdir\>

make install

then set ONLMON_MAIN to installdir, cd to \<srcdir\>/macros where the macros and scripts are located

source setup_onlmon.sh

(or for csh source setup_onlmon.csh)

## You need three sessions to run the show, two for the servers, one for the client. In the server windows start root and run the example server2:

root.exe

.x run_example_server0.C

and

root.exe

.x run_example_server1.C

run 100 events to get ths show going (otherwise your histos will be empty)

prun(100) 

## in the other start the client

root.exe

.L run_example_client.C

exampleDrawInit()

exampleDraw()

