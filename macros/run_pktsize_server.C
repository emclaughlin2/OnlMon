#include "ServerFuncs.C"

#include <onlmon/pktsize/PktSizeMon.h>

#include <onlmon/OnlMonServer.h>

// cppcheck-suppress unknownMacro
R__LOAD_LIBRARY(libonlpktsizemon_server.so)

void run_pktsize_server(const std::string &name = "PKTSIZEMON", unsigned int serverid = 0, const std::string &prdffile = "/sphenix/lustre01/sphnxpro/commissioning/aligned/beam-00022982-0000.prdf")
{
  OnlMon *m = new PktSizeMon("PKTSIZEMON");                     // create subsystem Monitor object
  m->SetMonitorServerId(serverid);
                                                //  m->AddTrigger("PPG(Laser)");  // high efficiency triggers selection at et pool
                                                //  m->AddTrigger("ONLMONBBCLL1"); // generic bbcll1 minbias trigger (defined in ServerFuncs.C)
  OnlMonServer *se = OnlMonServer::instance();  // get pointer to Server Framework
  se->registerMonitor(m);                       // register subsystem Monitor with Framework
  start_server(prdffile);
  return;
}
