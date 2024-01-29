#include "CommonFuncs.C"

#include <onlmon/daq/DaqMonDraw.h>

#include <onlmon/OnlMonClient.h>

// cppcheck-suppress unknownMacro
R__LOAD_LIBRARY(libonldaqmon_client.so)

void daqDrawInit(const int online = 0)
{
  OnlMonClient *cl = OnlMonClient::instance();
  // register histos we want with monitor name
  cl->registerHisto("daqmon_hist1", "DAQMON_0");
  cl->registerHisto("daqmon_hist2", "DAQMON_0");
  cl->AddServerHost("localhost");  // check local host first
  CreateHostList(online);
  // get my histos from server, the second parameter = 1
  // says I know they are all on the same node
  cl->requestHistoBySubSystem("DAQMON_0", 1);
  OnlMonDraw *daqmon = new DaqMonDraw("DAQMONDRAW");  // create Drawing Object
  cl->registerDrawer(daqmon);             // register with client framework
}

void daqDraw(const char *what = "ALL")
{
  OnlMonClient *cl = OnlMonClient::instance();  // get pointer to framewrk
  cl->requestHistoBySubSystem("DAQMON_0",1);        // update histos
  cl->Draw("DAQMONDRAW", what);                     // Draw Histos of registered Drawers
}

void daqSavePlot()
{
  OnlMonClient *cl = OnlMonClient::instance();  // get pointer to framewrk
  cl->SavePlot("DAQMONDRAW");                         // Save Plots
  return;
}

void daqHtml()
{
  OnlMonClient *cl = OnlMonClient::instance();  // get pointer to framewrk
  cl->MakeHtml("DAQMONDRAW");                       // Create html output
  return;
}
