#include "CommonFuncs.C"

#include <onlmon/bbc/BbcMonDraw.h>

#include <onlmon/OnlMonClient.h>

// cppcheck-suppress unknownMacro
R__LOAD_LIBRARY(libonlbbcmon_client.so)

void bbcDrawInit(const int online = 0)
{
  OnlMonClient *cl = OnlMonClient::instance();
  // register histos we want with monitor name
  cl->registerHisto("bbc_adc", "BBCMON_0");
  cl->registerHisto("bbc_tdc", "BBCMON_0");
  cl->registerHisto("bbc_tdc_overflow", "BBCMON_0");
  cl->registerHisto("bbc_tdc_armhittime", "BBCMON_0");
  cl->registerHisto("bbc_zvertex", "BBCMON_0");
  cl->registerHisto("bbc_zvertex_bbll1", "BBCMON_0");
  cl->registerHisto("bbc_nevent_counter", "BBCMON_0");
  cl->registerHisto("bbc_tzero_zvtx", "BBCMON_0");
  cl->registerHisto("bbc_prescale_hist", "BBCMON_0");
  cl->registerHisto("bbc_avr_hittime", "BBCMON_0");
  cl->registerHisto("bbc_north_hittime", "BBCMON_0");
  cl->registerHisto("bbc_south_hittime", "BBCMON_0");
  cl->registerHisto("bbc_north_chargesum", "BBCMON_0");
  cl->registerHisto("bbc_south_chargesum", "BBCMON_0");
  cl->registerHisto("bbc_time_wave", "BBCMON_0");
  cl->registerHisto("bbc_charge_wave", "BBCMON_0");
  cl->registerHisto("bbc_south_hitmap", "BBCMON_0");
  cl->registerHisto("bbc_north_hitmap", "BBCMON_0");

  cl->AddServerHost("localhost");  // check local host first
  CreateHostList(online);

  // get my histos from server, the second parameter = 1
  // says I know they are all on the same node
  cl->requestHistoBySubSystem("BBCMON_0", 1);
  OnlMonDraw *bbcmon = new BbcMonDraw("BBCMONDRAW");  // create Drawing Object
  cl->registerDrawer(bbcmon);             // register with client framework
}

void bbcDraw(const char *what = "ALL")
{
  OnlMonClient *cl = OnlMonClient::instance();  // get pointer to framewrk
  cl->requestHistoBySubSystem("BBCMON_0",1);        // update histos
  cl->Draw("BBCMONDRAW", what);                     // Draw Histos of registered Drawers
}

void bbcSavePlot()
{
  OnlMonClient *cl = OnlMonClient::instance();  // get pointer to framewrk
  cl->SavePlot("BBCMONDRAW");                         // Save Plots
  return;
}

void bbcHtml()
{
  OnlMonClient *cl = OnlMonClient::instance();  // get pointer to framewrk
  cl->MakeHtml("BBCMONDRAW");                       // Create html output
  return;
}
