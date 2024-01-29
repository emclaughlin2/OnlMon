#include "ServerFuncs.C"

#include <onlmon/intt/InttMonConstants.h>
#include <onlmon/intt/InttMon.h>
#include <onlmon/OnlMonServer.h>

// cppcheck-suppress unknownMacro
R__LOAD_LIBRARY(libonlinttmon_server.so)

void run_intt_server(const std::string &name = "INTTMON", unsigned int serverid = 0, const std::string &prdffile = "/sphenix/u/jbertaux/evt_files/intt_intt0-00000025-0000.evt")
{
  OnlMon *m = new InttMon(name);      // create subsystem Monitor object
  m->SetMonitorServerId(serverid);
 
  OnlMonServer *se = OnlMonServer::instance(); // get pointer to Server Framework
  se->registerMonitor(m);       // register subsystem Monitor with Framework
  start_server(prdffile);
  return;
}
