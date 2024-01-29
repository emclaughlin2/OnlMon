#include "CommonFuncs.C"

#include <onlmon/tpc/TpcMonDraw.h>

#include <onlmon/OnlMonClient.h>

// cppcheck-suppress unknownMacro
R__LOAD_LIBRARY(libonltpcmon_client.so)

void tpcDrawInit(const int online = 0)
{
  OnlMonClient *cl = OnlMonClient::instance();
  // register histos we want with monitor name
  cl->registerHisto("tpcmon_hist1", "TPCMON_0");
  cl->registerHisto("tpcmon_hist2", "TPCMON_0");

  char TPCMON_STR[100];
  // TPC ADC pie chart
  for( int i=0; i<24; i++ )
  //for( int i: {1,20} )
  {
    sprintf(TPCMON_STR,"TPCMON_%i",i);
    //std::cout<<"You registered the NSIDEADC/SSIDEADC "<<i<<" histo"<<std::endl;

    cl->registerHisto("NorthSideADC", TPCMON_STR);

    cl->registerHisto("NorthSideADC_clusterXY_R1", TPCMON_STR);
    cl->registerHisto("NorthSideADC_clusterXY_R2", TPCMON_STR);
    cl->registerHisto("NorthSideADC_clusterXY_R3", TPCMON_STR);

    cl->registerHisto("NorthSideADC_clusterXY_R1_unw", TPCMON_STR);
    cl->registerHisto("NorthSideADC_clusterXY_R2_unw", TPCMON_STR);
    cl->registerHisto("NorthSideADC_clusterXY_R3_unw", TPCMON_STR);

    cl->registerHisto("SouthSideADC", TPCMON_STR);

    cl->registerHisto("SouthSideADC_clusterXY_R1", TPCMON_STR);
    cl->registerHisto("SouthSideADC_clusterXY_R2", TPCMON_STR);
    cl->registerHisto("SouthSideADC_clusterXY_R3", TPCMON_STR);

    cl->registerHisto("SouthSideADC_clusterXY_R1_unw", TPCMON_STR);
    cl->registerHisto("SouthSideADC_clusterXY_R2_unw", TPCMON_STR);
    cl->registerHisto("SouthSideADC_clusterXY_R3_unw", TPCMON_STR);

    cl->registerHisto("sample_size_hist",TPCMON_STR);
    cl->registerHisto("Check_Sum_Error",TPCMON_STR);
    cl->registerHisto("Check_Sums",TPCMON_STR);
    cl->registerHisto("ADC_vs_SAMPLE",TPCMON_STR); 
    cl->registerHisto("ADC_vs_SAMPLE_large",TPCMON_STR); 
    cl->registerHisto("MAXADC",TPCMON_STR);

    cl->registerHisto("RAWADC_1D_R1",TPCMON_STR);
    cl->registerHisto("MAXADC_1D_R1",TPCMON_STR);
    cl->registerHisto("RAWADC_1D_R2",TPCMON_STR);
    cl->registerHisto("MAXADC_1D_R2",TPCMON_STR);
    cl->registerHisto("RAWADC_1D_R3",TPCMON_STR);
    cl->registerHisto("MAXADC_1D_R3",TPCMON_STR);

    cl->registerHisto("NorthSideADC_clusterZY", TPCMON_STR);
    cl->registerHisto("SouthSideADC_clusterZY", TPCMON_STR);

    cl->registerHisto("NorthSideADC_clusterZY_unw",TPCMON_STR);
    cl->registerHisto("SouthSideADC_clusterZY_unw",TPCMON_STR);
  } //



  cl->AddServerHost("localhost");  // check local host first
  CreateHostList(online);
  // get my histos from server, the second parameter = 1
  // says I know they are all on the same node

  for( int i=0; i<24; i++ )
  //for( int i: {1,20} )
  {
    sprintf(TPCMON_STR,"TPCMON_%i",i);
    cl->requestHistoBySubSystem(TPCMON_STR, 1);
  }

  OnlMonDraw *tpcmon = new TpcMonDraw("TPCMONDRAW");  // create Drawing Object
  cl->registerDrawer(tpcmon);             // register with client framework
}

void tpcDraw(const char *what = "ALL")
{
  OnlMonClient *cl = OnlMonClient::instance();  // get pointer to framewrk

  char TPCMON_STR[100];

  for( int i=0; i<24; i++ )
  //for( int i: {1,20} )
  {
    sprintf(TPCMON_STR,"TPCMON_%i",i);
    cl->requestHistoBySubSystem(TPCMON_STR, 1);
  }


  cl->Draw("TPCMONDRAW", what);                     // Draw Histos of registered Drawers
}

void tpcSavePlot()
{
  OnlMonClient *cl = OnlMonClient::instance();  // get pointer to framewrk
  cl->SavePlot("TPCMONDRAW");                         // Save Plots
  return;
}

void tpcHtml()
{
  OnlMonClient *cl = OnlMonClient::instance();  // get pointer to framewrk
  cl->MakeHtml("TPCMONDRAW");                       // Create html output
  return;
}
