#include "ServerFuncs.C"

#include <onlmon/tpc/TpcMon.h>

#include <onlmon/OnlMonServer.h>

// cppcheck-suppress unknownMacro
R__LOAD_LIBRARY(libonltpcmon_server.so)

void run_tpc_server(const std::string &name = "TPCMON", unsigned int serverid = 0, const std::string &prdffile = "/sphenix/data/data02/sphnxpro/tpc/chughes/prdf/00010169/TPC_ebdc00_pedestal-00010169-0000.prdf")
{
  OnlMon *m = new TpcMon(name);                     // create subsystem Monitor object
  m->SetMonitorServerId(serverid);
                                                //  m->AddTrigger("PPG(Laser)");  // high efficiency triggers selection at et pool
                                                //  m->AddTrigger("ONLMONBBCLL1"); // generic bbcll1 minbias trigger (defined in ServerFuncs.C)
  OnlMonServer *se = OnlMonServer::instance();  // get pointer to Server Framework
  se->registerMonitor(m);                       // register subsystem Monitor with Framework
  start_server(prdffile);
  return;
}
