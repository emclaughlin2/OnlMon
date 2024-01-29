#include <CommonFuncs.C>

#include <onlmon/cemc/CemcMonDraw.h>

#include <onlmon/OnlMonClient.h>

// cppcheck-suppress unknownMacro
R__LOAD_LIBRARY(libonlcemcmon_client.so)

std::string monitorname("CEMCMON");

void cemcDrawInit(const int online = 0)
{
  OnlMonClient *cl = OnlMonClient::instance();
  cl->AddServerHost("localhost");  // stand-in for SEB00
  cl->AddServerHost("rcas2068");   //  stand-in for SEB01
  // CreateHostList(online);
  // register histos we want with monitor name
  for (unsigned int i = 0; i < 2; i++)
  {
  std::string monitor = monitorname + string("_") + to_string(i);
  cl->registerHisto("cemc_occupancy", monitor);
  cl->registerHisto("cemc_runningmean", monitor);
  cl->requestHistoBySubSystem(monitor, 1);
  }
  //  get my histos from server, the second parameter = 1
  //  says I know they are all on the same node
  OnlMonDraw *cemcmon = new CemcMonDraw("CEMCMONDRAW");  // create Drawing Object
  cl->registerDrawer(cemcmon);              // register with client framework
}

void cemcDraw(const char *what = "Standard")
{
  OnlMonClient *cl = OnlMonClient::instance();  // get pointer to framewrk
  for (unsigned int i = 0; i < 2; i++)
  {
  std::string monitor = monitorname + string("_") + to_string(i);
  cl->requestHistoBySubSystem(monitor);       // update histos
  }
  cl->Draw("CEMCMONDRAW", what);                    // Draw Histos of registered Drawers
}

void cemcPS()
{
  OnlMonClient *cl = OnlMonClient::instance();  // get pointer to framewrk
  cl->MakePS("CEMCMONDRAW");                        // Create PS files
  return;
}

void cemcHtml()
{
  OnlMonClient *cl = OnlMonClient::instance();  // get pointer to framewrk
  cl->MakeHtml("CEMCMONDRAW");                      // Create html output
  return;
}
