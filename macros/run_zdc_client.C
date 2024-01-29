#include "CommonFuncs.C"

#include <onlmon/zdc/ZdcMonDraw.h>

#include <onlmon/OnlMonClient.h>

// cppcheck-suppress unknownMacro
R__LOAD_LIBRARY(libonlzdcmon_client.so)

void zdcDrawInit(const int online = 0)
{

  OnlMonClient *cl = OnlMonClient::instance();
  // register histos we want with monitor name
  cl->registerHisto("zdc_adc_north","ZDCMON_0");
  cl->registerHisto("zdc_adc_south","ZDCMON_0");

  cl->AddServerHost("localhost");  // check local host first
  CreateHostList(online);

 // get my histos from server, the second parameter = 1 
// says I know they are all on the same node
  cl->requestHistoBySubSystem("ZDCMON_0", 1);
  OnlMonDraw *zdcmon = new ZdcMonDraw("ZDCMONDRAW");    // create Drawing Object
  cl->registerDrawer(zdcmon);              // register with client framework
}

void zdcDraw(const char *what="ALL")
{
  OnlMonClient *cl = OnlMonClient::instance();  // get pointer to framewrk
  cl->requestHistoBySubSystem("ZDCMON_0",1);      // update histos
  cl->Draw("ZDCMONDRAW",what);                      // Draw Histos of registered Drawers
}

void zdcSavePlot()
{
    OnlMonClient *cl = OnlMonClient::instance();  // get pointer to framewrk
    cl->SavePlot("ZDCMONDRAW");                  // Save Plots
    return;
}

void zdcHtml()
{
  OnlMonClient *cl = OnlMonClient::instance();  // get pointer to framewrk
  cl->MakeHtml("ZDCMONDRAW");                       // create html output
  return;
}
