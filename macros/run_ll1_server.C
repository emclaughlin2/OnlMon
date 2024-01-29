#include "ServerFuncs.C"

#include <onlmon/ll1/LL1Mon.h>

#include <onlmon/OnlMonServer.h>

// cppcheck-suppress unknownMacro
R__LOAD_LIBRARY(libonlll1mon_server.so)

void run_ll1_server(const std::string &name = "LL1MON", unsigned int serverid = 0, const std::string &prdffile = "/sphenix/lustre01/sphnxpro/rawdata/commissioning/LL1/calib/calib_LL1-00013020-0000.prdf")
{
  OnlMon *m = new LL1Mon(name);                  // create subsystem Monitor object
  m->SetMonitorServerId(serverid);
  
  OnlMonServer *se = OnlMonServer::instance();  // get pointer to Server Framework
  se->registerMonitor(m);                       // register subsystem Monitor with Framework
  start_server(prdffile);
  return;
}
