#include "ServerFuncs.C"

#include <onlmon/tpot/TpotMon.h>

#include <onlmon/OnlMonServer.h>

// cppcheck-suppress unknownMacro
R__LOAD_LIBRARY(libonltpotmon_server.so)

void run_tpot_server(
  const std::string &name = "TPOTMON",
  unsigned int serverid = 0,
  const std::string &prdffile = "/sphenix/lustre01/sphnxpro/commissioning/TPOT/junk/TPOT_ebdc39_junk-00020121-0000.prdf"
  )
{
  // create subsystem Monitor object
  auto m = new TpotMon(name);
  m->SetMonitorServerId(serverid);

  // get pointer to Server Framework
  auto se = OnlMonServer::instance();

  // register subsystem Monitor with Framework
  se->registerMonitor(m);

  // and start
  start_server(prdffile);
  return;
}
