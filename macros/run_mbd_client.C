#include <CommonFuncs.C>

#include <onlmon/mbd/MbdMonDraw.h>

#include <onlmon/OnlMonClient.h>

R__LOAD_LIBRARY(libonlmbdmon_client.so)

void mbdDrawInit(const int online = 0)
{
  OnlMonClient *cl = OnlMonClient::instance();
  // register histos we want with monitor name
  cl->registerHisto("bbc_adc", "MBDMON");
  cl->registerHisto("bbc_tdc", "MBDMON");
  cl->registerHisto("bbc_tdc_overflow", "MBDMON");
  cl->registerHisto("bbc_tdc_armhittime", "MBDMON");
  cl->registerHisto("bbc_zvertex", "MBDMON");
  cl->registerHisto("bbc_zvertex_bbll1", "MBDMON");
  //cl->registerHisto("bbc_zvertex_zdc", "MBDMON");
  //cl->registerHisto("bbc_zvertex_zdc_scale3", "MBDMON");
  cl->registerHisto("bbc_zvertex_bbll1_novtx", "MBDMON");
  cl->registerHisto("bbc_zvertex_bbll1_narrowvtx", "MBDMON");
  cl->registerHisto("bbc_nevent_counter", "MBDMON");
  cl->registerHisto("bbc_tzero_zvtx", "MBDMON");
  cl->registerHisto("bbc_prescale_hist", "MBDMON");
  cl->registerHisto("bbc_avr_hittime", "MBDMON");
  cl->registerHisto("bbc_north_hittime", "MBDMON");
  cl->registerHisto("bbc_south_hittime", "MBDMON");
  cl->registerHisto("bbc_north_chargesum", "MBDMON");
  cl->registerHisto("bbc_south_chargesum", "MBDMON");
  //cl->registerHisto("bbc_zvertex_bbll1_zdc", "MBDMON");

  cl->AddServerHost("localhost");   // check local host first
  CreateHostList(online);

  // get my histos from server, the second parameter = 1
  // says I know they are all on the same node
  cl->requestHistoBySubSystem("MBDMON", 1);
  OnlMonDraw *mbdmon = new MbdMonDraw();    // create Drawing Object
  cl->registerDrawer(mbdmon);              // register with client framework
}

void mbdDraw(const char *what="ALL")
{
  OnlMonClient *cl = OnlMonClient::instance();  // get pointer to framewrk
  cl->requestHistoBySubSystem("MBDMON");         // update histos
  cl->Draw("MBDMON",what);                       // Draw Histos of registered Drawers
}

void mbdPS()
{
  OnlMonClient *cl = OnlMonClient::instance();  // get pointer to framewrk
  cl->MakePS("MBDMON");                          // Create PS files
  return;
}

void mbdHtml()
{
  OnlMonClient *cl = OnlMonClient::instance();  // get pointer to framewrk
  cl->MakeHtml("MBDMON");                        // Create html output
  return;
}
