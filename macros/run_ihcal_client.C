#include "CommonFuncs.C"

#include <onlmon/hcal/HcalMonDraw.h>

#include <onlmon/OnlMonClient.h>

// cppcheck-suppress unknownMacro
R__LOAD_LIBRARY(libonlhcalmon_client.so)

void ihcalDrawInit(const int online = 0)
{
  OnlMonClient *cl = OnlMonClient::instance();
  // register histos we want with monitor name
  char HCALMON_STR[100];
  for (int i = 0; i < 2; i++)
  {
    sprintf(HCALMON_STR, "IHCALMON_%i", i);
    cl->registerHisto("h2_hcal_hits", HCALMON_STR);
    cl->registerHisto("h2_hcal_rm", HCALMON_STR);
    cl->registerHisto("h2_hcal_mean", HCALMON_STR);
    cl->registerHisto("h2_hcal_waveform", HCALMON_STR);
    cl->registerHisto("h2_hcal_correlation", HCALMON_STR);
    cl->registerHisto("h_event", HCALMON_STR);
    cl->registerHisto("h_sectorAvg_total", HCALMON_STR);
    cl->registerHisto("h_waveform_twrAvg", HCALMON_STR);
    cl->registerHisto("h_waveform_time", HCALMON_STR);
    cl->registerHisto("h_waveform_pedestal", HCALMON_STR);
    cl->registerHisto("h_ntower", HCALMON_STR);
    cl->registerHisto("h1_packet_chans", HCALMON_STR);
    cl->registerHisto("h1_packet_length", HCALMON_STR);
    cl->registerHisto("h1_packet_number", HCALMON_STR);
    cl->registerHisto("h1_packet_event", HCALMON_STR);

    for (int ih = 0; ih < 32; ih++)
    {
      cl->registerHisto(Form("h_rm_sectorAvg_s%d", ih), HCALMON_STR);
    }
    for (int ieta = 0; ieta < 24; ieta++)
    {
      for (int iphi = 0; iphi < 64; iphi++)
      {
        cl->registerHisto(Form("h_rm_tower_%d_%d", ieta, iphi), HCALMON_STR);
      }
    }
  }
  cl->AddServerHost("localhost");  // check local host first
  CreateHostList(online);
  // get my histos from server, the second parameter = 1
  // says I know they are all on the same node
  for (int i = 0; i < 2; i++)
  {
    sprintf(HCALMON_STR, "IHCALMON_%i", i);
    cl->requestHistoBySubSystem(HCALMON_STR, 1);
  }
  OnlMonDraw *hcalmon = new HcalMonDraw("IHCALMONDRAW");  // create Drawing Object
  
  cl->registerDrawer(hcalmon);                           // register with client framework
}

void ihcalDraw(const char *what = "ALL")
{
  OnlMonClient *cl = OnlMonClient::instance();  // get pointer to framewrk
  char HCALMON_STR[100];
  for (int i = 0; i < 2; i++)
  {
    sprintf(HCALMON_STR, "IHCALMON_%i", i);
    cl->requestHistoBySubSystem(HCALMON_STR, 1);
  }
  cl->Draw("IHCALMONDRAW", what);  // Draw Histos of registered Drawers
}

void ihcalSavePlot()
{
  OnlMonClient *cl = OnlMonClient::instance();  // get pointer to framewrk
  cl->SavePlot("IHCALMONDRAW");                    // Save Plots
  return;
}

void ihcalHtml()
{
  OnlMonClient *cl = OnlMonClient::instance();  // get pointer to framewrk
  cl->MakeHtml("IHCALMONDRAW");                  // Create html output
  return;
}
