#include "ServerFuncs.C"

#include <onlmon/daq/DaqMon.h>

#include <onlmon/OnlMonServer.h>

// cppcheck-suppress unknownMacro
R__LOAD_LIBRARY(libonldaqmon_server.so)

void run_daq_server(const std::string &name = "DAQMON", unsigned int serverid = 0, const std::string &prdffile = "/sphenix/data/data02/sphenix/t1044/rcdaq-00000221-0000.prdf")
{
  OnlMon *m = new DaqMon(name);                     // create subsystem Monitor object
  m->SetMonitorServerId(serverid);
                                                //  m->AddTrigger("PPG(Laser)");  // high efficiency triggers selection at et pool
                                                //  m->AddTrigger("ONLMONBBCLL1"); // generic bbcll1 minbias trigger (defined in ServerFuncs.C)
  OnlMonServer *se = OnlMonServer::instance();  // get pointer to Server Framework
  se->registerMonitor(m);                       // register subsystem Monitor with Framework
  start_server(prdffile);
  return;
}
