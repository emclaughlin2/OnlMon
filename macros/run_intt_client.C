#include "CommonFuncs.C"

#include <onlmon/OnlMonClient.h>
#include <onlmon/intt/InttMonDraw.h>

// cppcheck-suppress unknownMacro
R__LOAD_LIBRARY(libonlinttmon_client.so)

void inttDrawInit(const int online = 0)
{
  OnlMonClient *cl = OnlMonClient::instance();
  // register histos we want with monitor name

  cl->AddServerHost("localhost");	// check local host first
  CreateHostList(online);

  for(int felix = 0; felix < INTT::FELIX; ++felix)
  {
    cl->registerHisto("InttNumEvents",	Form("INTTMON_%d", felix));
    cl->registerHisto("InttMap",	Form("INTTMON_%d", felix));
  }
  //cl->registerHisto("InttHitMapRef",	"INTTMON_0");

  // get my histos from server, the second parameter = 1
  // says I know they are all on the same node
  for(int felix = 0; felix < INTT::FELIX; ++felix)
  {
    cl->requestHistoBySubSystem(Form("INTTMON_%d", felix), 1);
  }
  OnlMonDraw *inttmon = new InttMonDraw("INTTMONDRAW");    // create Drawing Object
  inttmon->Init(); //registers the hists it will need to the OnlMonClient
  cl->registerDrawer(inttmon);              // register with client framework
}

void inttDraw(const char *what="ALL")
{
  OnlMonClient *cl = OnlMonClient::instance();		// get pointer to framewrk
  for(int felix = 0; felix < INTT::FELIX; ++felix)	// update histos
  {
    cl->requestHistoBySubSystem(Form("INTTMON_%d", felix),1);
  }
  cl->Draw("INTTMONDRAW",what);				// Draw Histos of registered Drawers
}

void inttSavePlot()
{
  OnlMonClient *cl = OnlMonClient::instance();  // get pointer to framewrk
  cl->SavePlot("INTTMONDRAW");                          // Save Plots
  return;
}

void inttHtml()
{
  OnlMonClient *cl = OnlMonClient::instance();  // get pointer to framewrk
  cl->MakeHtml("INTTMONDRAW");                        // Create html output
  return;
}
