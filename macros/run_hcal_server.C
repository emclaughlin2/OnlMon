#include "ServerFuncs.C"

#include <onlmon/hcal/HcalMon.h>

#include <onlmon/OnlMonServer.h>

// cppcheck-suppress unknownMacro
R__LOAD_LIBRARY(libonlhcalmon_server.so)

void run_hcal_server(const std::string &name = "HCALMON", unsigned int serverid = 0, const std::string &prdffile = "/bbox/commissioning/HCal/calib/calib_West-00011989-0000.prdf")
{
  char OHCALMON_STR[100];
  sprintf(OHCALMON_STR, "O%s", name.c_str());
  char IHCALMON_STR[100];
  sprintf(IHCALMON_STR, "I%s", name.c_str());
  OnlMon *m = new HcalMon(OHCALMON_STR);                    // create subsystem Monitor object
  m->SetMonitorServerId(serverid);
  
                                                //  m->AddTrigger("PPG(Laser)");  // high efficiency triggers selection at et pool
                                                //  m->AddTrigger("ONLMONBBCLL1"); // generic bbcll1 minbias trigger (defined in ServerFuncs.C)
  OnlMonServer *se = OnlMonServer::instance();  // get pointer to Server Framework
  se->registerMonitor(m);                       // register subsystem Monitor with Framework
  m = new HcalMon(IHCALMON_STR);                  
  m->SetMonitorServerId(serverid);
  se->registerMonitor(m);                       
  start_server(prdffile);
  return;
}
