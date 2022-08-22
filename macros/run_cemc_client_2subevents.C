#include <CommonFuncs.C>

#include <onlmon/cemc/CemcMonDraw.h>

#include <onlmon/OnlMonClient.h>

R__LOAD_LIBRARY(libonlcemcmon_client.so)

void cemcDrawInit(const int online = 0)
{
  OnlMonClient *cl = OnlMonClient::instance();
  // register histos we want with monitor name
  cl->registerHisto("cemc_occupancy_SEB00", "CEMCMON");
  cl->registerHisto("cemc_occupancy_SEB01", "CEMCMON");
  cl->registerHisto("cemc_runningmean_SEB00", "CEMCMON");
  cl->registerHisto("cemc_runningmean_SEB01", "CEMCMON");
  cl->AddServerHost("localhost");   // stand-in for SEB00
  cl->AddServerHost("rcas2068");   //  stand-in for SEB01
  //CreateHostList(online);
  // get my histos from server, the second parameter = 1
  // says I know they are all on the same node
  cl->requestHistoBySubSystem("CEMCMON", 1);
  OnlMonDraw *cemcmon = new CemcMonDraw();    // create Drawing Object
  cl->registerDrawer(cemcmon);              // register with client framework
}

void cemcDraw(const char *what="Standard")
{
  OnlMonClient *cl = OnlMonClient::instance();  // get pointer to framewrk
  cl->requestHistoBySubSystem("CEMCMON");         // update histos
  cl->Draw("CEMCMON",what);                       // Draw Histos of registered Drawers

}

void cemcPS()
{
  OnlMonClient *cl = OnlMonClient::instance();  // get pointer to framewrk
  cl->MakePS("CEMCMON");                          // Create PS files
  return;
}

void cemcHtml()
{
  OnlMonClient *cl = OnlMonClient::instance();  // get pointer to framewrk
  cl->MakeHtml("CEMCMON");                        // Create html output
  return;
}
