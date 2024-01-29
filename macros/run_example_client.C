//
// You need to run 2 servers by executing (in 2 separate root sessions)
// .x run_example_server0.C
// .x run_example_server1.C
// and then start the client
// if you only start one server only one set of histograms will be displayed

#include "CommonFuncs.C"

#include <onlmon/example/MyMonDraw.h>

#include <onlmon/OnlMonClient.h>

// cppcheck-suppress unknownMacro
R__LOAD_LIBRARY(libonlmymon_client.so)

void exampleDrawInit(const int online = 0)
{
  OnlMonClient *cl = OnlMonClient::instance();
  // register histos we want with monitor name
  cl->registerHisto("mymon_hist1", "MYMON_0");
  cl->registerHisto("mymon_hist2", "MYMON_0");
  cl->registerHisto("mymon_hist1", "MYMON_1");
  cl->registerHisto("mymon_hist2", "MYMON_1");
  cl->AddServerHost("localhost");  // check local host first
  CreateHostList(online);
  // get my histos from server, the second parameter = 1
  // says I know they are all on the same node
  cl->requestHistoBySubSystem("MYMON_0", 1);
  cl->requestHistoBySubSystem("MYMON_1", 1);
  OnlMonDraw *mymon = new MyMonDraw("MYMONDRAW");  // create Drawing Object
  cl->registerDrawer(mymon);            // register with client framework
}

void exampleDraw(const char *what = "ALL")
{
  OnlMonClient *cl = OnlMonClient::instance();  // get pointer to framewrk
  cl->requestHistoBySubSystem("MYMON_1",1);         // update histos
  cl->requestHistoBySubSystem("MYMON_0",1);         // update histos
  cl->Draw("MYMONDRAW", what);                      // Draw Histos of registered Drawers
}

void exampleSavePlot()
{
  OnlMonClient *cl = OnlMonClient::instance();  // get pointer to framewrk
  cl->SavePlot("MYMONDRAW");                          // Save Plots
  return;
}

void exampleHtml()
{
  OnlMonClient *cl = OnlMonClient::instance();  // get pointer to framewrk
  cl->MakeHtml("MYMONDRAW");                        // Create html output
  return;
}
