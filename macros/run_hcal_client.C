#include <CommonFuncs.C>

#include <onlmon/hcal/HcalMonDraw.h>

#include <onlmon/OnlMonClient.h>

R__LOAD_LIBRARY(libonlhcalmon_client.so)

void hcalDrawInit(const int online = 0)
{
  OnlMonClient *cl = OnlMonClient::instance();
  // register histos we want with monitor name
  cl->registerHisto("h2_hcal_hits", "HCALMON");
  cl->registerHisto("h2_hcal_mean", "HCALMON");
  cl->registerHisto("h2_hcal_history", "HCALMON");
  cl->AddServerHost("localhost");   // check local host first
  CreateHostList(online);
  // get my histos from server, the second parameter = 1
  // says I know they are all on the same node
  cl->requestHistoBySubSystem("HCALMON", 1);
  OnlMonDraw *hcalmon = new HcalMonDraw();    // create Drawing Object
  cl->registerDrawer(hcalmon);              // register with client framework
}

void hcalDraw(const char *what="ALL")
{
  OnlMonClient *cl = OnlMonClient::instance();  // get pointer to framewrk
  cl->requestHistoBySubSystem("HCALMON");         // update histos
  cl->Draw("HCALMON",what);                       // Draw Histos of registered Drawers
}

void hcalPS()
{
  OnlMonClient *cl = OnlMonClient::instance();  // get pointer to framewrk
  cl->MakePS("HCALMON");                          // Create PS files
  return;
}

void hcalHtml()
{
  OnlMonClient *cl = OnlMonClient::instance();  // get pointer to framewrk
  cl->MakeHtml("HCALMON");                        // Create html output
  return;
}
