# OnlMon

[![Build Status](https://web.sdcc.bnl.gov/jenkins-sphenix/buildStatus/icon?job=sPHENIX%2FsPHENIX_OnlMon_MasterBranch)](https://web.sdcc.bnl.gov/jenkins-sphenix/job/sPHENIX/job/sPHENIX_OnlMon_MasterBranch/)


To build it, create a build area and run

\<srcdir\>/autogen.sh --prefix=\<installdir\>

make install

then set ONLMON_MAIN to installdir, cd to \<srcdir\>/macros where the macros and scripts are located

source setup_onlmon.sh

(or for csh source setup_onlmon.csh)

## You need two sessions to run the show, one for the server, one for the client. In the server window start root and run the example server:

root.exe

.x run_example_server.C

prun(200)

## in the other start the client

root.exe

.L run_example_client.C

exampleDrawInit()

exampleDraw()

