// This starts the example server with id=0

#include <ServerFuncs.C>

#include <onlmon/example/MyMon.h>

#include <onlmon/OnlMonServer.h>

// cppcheck-suppress unknownMacro
R__LOAD_LIBRARY(libonlmymon_server.so)

void run_example_server0(const std::string &name = "MYMON", unsigned int serverid = 0, const std::string &prdffile = "/sphenix/data/data02/sphenix/t1044/rcdaq-00000221-0000.prdf")
{
// create subsystem Monitor object
  OnlMon *m = new MyMon(name);
// set server id needed for running multiple servers
  m->SetMonitorServerId(serverid);
// get pointer to Server Framework
  OnlMonServer *se = OnlMonServer::instance();
//  se->Verbosity(3);
  se->registerMonitor(m);
  start_server(prdffile);
  prun(100);
  return;
}
