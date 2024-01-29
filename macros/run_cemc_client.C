#include "CommonFuncs.C"

#include <onlmon/cemc/CemcMonDraw.h>

#include <onlmon/OnlMonClient.h>

// cppcheck-suppress unknownMacro
R__LOAD_LIBRARY(libonlcemcmon_client.so)
//int nserver = 0;
void cemcDrawInit(const int online = 0)
{
  OnlMonClient *cl = OnlMonClient::instance();
  cl -> Verbosity(0);
  // register histos we want with monitor name
  for (int serverid = 0; serverid <= 7; serverid++)
    {
      cl->registerHisto("h2_cemc_rm", Form("CEMCMON_%d",serverid));
      cl->registerHisto("h2_cemc_mean", Form("CEMCMON_%d",serverid));
      cl->registerHisto("h1_event", Form("CEMCMON_%d",serverid));
      cl->registerHisto("h2_waveform_twrAvg", Form("CEMCMON_%d",serverid));
      cl->registerHisto("h1_waveform_time", Form("CEMCMON_%d",serverid));
      cl->registerHisto("h1_waveform_pedestal", Form("CEMCMON_%d",serverid));
      cl->registerHisto("h1_fitting_sigDiff",Form("CEMCMON_%d",serverid));
      cl->registerHisto("h1_fitting_pedDiff",Form("CEMCMON_%d",serverid));
      cl->registerHisto("h1_fitting_timeDiff",Form("CEMCMON_%d",serverid));
      cl->registerHisto("h1_packet_number",Form("CEMCMON_%d",serverid));
      cl->registerHisto("h1_packet_length",Form("CEMCMON_%d",serverid));
      cl->registerHisto("h1_packet_chans",Form("CEMCMON_%d",serverid));
      cl->registerHisto("h1_cemc_adc",Form("CEMCMON_%d",serverid));
		    
      for (int ih=0; ih<32; ih++){
	cl->registerHisto(Form("h1_rm_sectorAvg_s%d",ih), Form("CEMCMON_%d",serverid));
      }
    }
  CreateHostList(online);
  //  get my histos from server, the second parameter = 1
  //  says I know they are all on the same node
  //cl->requestHistoBySubSystem(Form("CEMCMON_%d",serverid), 1);
  for(int i = 0; i <= 7; i++)
     {
       cl->requestHistoBySubSystem(Form("CEMCMON_%d",i), 1);
     }
  CemcMonDraw *cemcmon = new CemcMonDraw("CEMCMONDRAW");  // create Drawing Object
  cemcmon -> setSave(0);
  cl->registerDrawer(cemcmon);              // register with client framework
}

void cemcDraw(const char *what = "ALL")
{
  OnlMonClient *cl = OnlMonClient::instance();  // get pointer to framewrk
  for(int i = 0; i <= 7; i++)
    {
      cl->requestHistoBySubSystem(Form("CEMCMON_%d",i), 1);
    }
  cl->Draw("CEMCMONDRAW", what);  // Draw Histos of registered Drawers
}

void cemcSavePlot()
{
  OnlMonClient *cl = OnlMonClient::instance();  // get pointer to framewrk
  cl->SavePlot("CEMCMONDRAW");                         // Save Plots
  return;
}

void cemcHtml()
{
  OnlMonClient *cl = OnlMonClient::instance();  // get pointer to framewrk
  cl->MakeHtml("CEMCMONDRAW");                      // Create html output
  return;
}
