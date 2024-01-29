#include "CommonFuncs.C"

#include <onlmon/pktsize/PktSizeMonDraw.h>

#include <onlmon/OnlMonClient.h>

// cppcheck-suppress unknownMacro
R__LOAD_LIBRARY(libonlpktsizemon_client.so)

void pktsizeDrawInit(const int online = 0)
{
  OnlMonClient *cl = OnlMonClient::instance();
  // register histos we want with monitor name
  cl->registerHisto("pktsize_hist", "PKTSIZEMON_0");
  cl->AddServerHost("localhost");  // check local host first
//  CreateHostList(online);
  // get my histos from server, the second parameter = 1
  // says I know they are all on the same node
  cl->requestHistoBySubSystem("PKTSIZEMON_0", 1);
  OnlMonDraw *pktsizemon = new PktSizeMonDraw("PKTSIZEMONDRAW");  // create Drawing Object
  cl->registerDrawer(pktsizemon);             // register with client framework
}

void pktsizeDraw(const char *what = "ALL")
{
  OnlMonClient *cl = OnlMonClient::instance();  // get pointer to framewrk
  cl->requestHistoBySubSystem("PKTSIZEMON_0",1);        // update histos
  cl->Draw("PKTSIZEMONDRAW", what);                     // Draw Histos of registered Drawers
}

void pktsizeSavePlot()
{
  OnlMonClient *cl = OnlMonClient::instance();  // get pointer to framewrk
  cl->SavePlot("PKTSIZEMONDRAW");                         // Save Plots
  return;
}

void pktsizeHtml()
{
  OnlMonClient *cl = OnlMonClient::instance();  // get pointer to framewrk
  cl->MakeHtml("PKTSIZEMONDRAW");                       // Create html output
  return;
}
