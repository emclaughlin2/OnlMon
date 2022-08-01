// Phenix Online Monitoring System
//
// Macro loader for starting Online Monitoring GUI Framework
//
// Author:  Cody McCain (mcm99c@acu.edu)
// Date:    January 2, 2003
//
// See run_Poms.C.README

#include <onlmon/Poms.h>

R__LOAD_LIBRARY(libpoms.so)

void StartPoms()
{
  PomsMainFrame* pmf;
  pmf = PomsMainFrame::Instance();
  
  // Register SubSystems
  pmf->RegisterSubSystem("EXAMPLE", "example");
  
  SubSystem *cemc = new SubSystem("CEMC", "cemc");
  cemc->AddAction("cemcDraw(\"FIRST\")", "CEMC 1st Monitor");
  cemc->AddAction("cemcDraw(\"SECOND\")", "CEMC 2nd Monitor");
  cemc->AddAction("cemcDraw(\"HISTORY\")", "CEMC History Monitor");
  cemc->AddAction(new SubSystemActionDrawPS(cemc));
  pmf->RegisterSubSystem(cemc);  
  
  
  SubSystem *daq = new SubSystem("DAQ", "daq");
  daq->AddAction("daqDraw(\"FIRST\")", "DAQ 1st Monitor");
  daq->AddAction("daqDraw(\"SECOND\")", "DAQ 2nd Monitor");
  daq->AddAction("daqDraw(\"HISTORY\")", "DAQ History Monitor");
  daq->AddAction(new SubSystemActionDrawPS(daq));
  pmf->RegisterSubSystem(daq);  

  SubSystem *epd = new SubSystem("EPD", "epd");
  epd->AddAction("epdDraw(\"FIRST\")", "EPD 1st Monitor");
  epd->AddAction("epdDraw(\"SECOND\")", "EPD 2nd Monitor");
  epd->AddAction("epdDraw(\"HISTORY\")", "EPD History Monitor");
  epd->AddAction(new SubSystemActionDrawPS(epd));
  pmf->RegisterSubSystem(epd);  

  SubSystem *hcal = new SubSystem("HCAL", "hcal");
  hcal->AddAction("hcalDraw(\"FIRST\")", "HCAL 1st Monitor");
  hcal->AddAction("hcalDraw(\"SECOND\")", "HCAL 2nd Monitor");
  hcal->AddAction("hcalDraw(\"HISTORY\")", "HCAL History Monitor");
  hcal->AddAction(new SubSystemActionDrawPS(hcal));
  pmf->RegisterSubSystem(hcal);  

  SubSystem *intt = new SubSystem("INTT", "intt");
  intt->AddAction("inttDraw(\"FIRST\")", "INTT 1st Monitor");
  intt->AddAction("inttDraw(\"SECOND\")", "INTT 2nd Monitor");
  intt->AddAction("inttDraw(\"HISTORY\")", "INTT History Monitor");
  intt->AddAction(new SubSystemActionDrawPS(intt));
  pmf->RegisterSubSystem(intt);  

  SubSystem *mbd = new SubSystem("MBD", "mbd");
  mbd->AddAction("mbdDraw(\"FIRST\")", "MBD 1st Monitor");
  mbd->AddAction("mbdDraw(\"SECOND\")", "MBD 2nd Monitor");
  mbd->AddAction("mbdDraw(\"HISTORY\")", "MBD History Monitor");
  mbd->AddAction(new SubSystemActionDrawPS(mbd));
  pmf->RegisterSubSystem(mbd);  

  SubSystem *mvtx = new SubSystem("MVTX", "mvtx");
  mvtx->AddAction("mvtxDraw(\"FIRST\")", "MVTX 1st Monitor");
  mvtx->AddAction("mvtxDraw(\"SECOND\")", "MVTX 2nd Monitor");
  mvtx->AddAction("mvtxDraw(\"HISTORY\")", "MVTX History Monitor");
  mvtx->AddAction(new SubSystemActionDrawPS(mvtx));
  pmf->RegisterSubSystem(mvtx);  
  
  SubSystem *mbdll1 = new SubSystem("MBDLL1", "mbdll1");
  mbdll1->AddAction("mbdll1Draw(\"FIRST\")", "MBDLL1 1st Monitor");
  mbdll1->AddAction("mbdll1Draw(\"SECOND\")", "MBDLL1 2nd Monitor");
  mbdll1->AddAction("mbdll1Draw(\"HISTORY\")", "MBDLL1 History Monitor");
  mbdll1->AddAction(new SubSystemActionDrawPS(mbdll1));
  pmf->RegisterSubSystem(mbdll1);  

  SubSystem *tpc = new SubSystem("TPC", "tpc");
  tpc->AddAction("tpcDraw(\"FIRST\")", "TPC 1st Monitor");
  tpc->AddAction("tpcDraw(\"SECOND\")", "TPC 2nd Monitor");
  tpc->AddAction("tpcDraw(\"HISTORY\")", "TPC History Monitor");
  tpc->AddAction(new SubSystemActionDrawPS(tpc));
  pmf->RegisterSubSystem(tpc);  

  SubSystem *tpot = new SubSystem("TPOT", "tpot");
  tpot->AddAction("tpotDraw(\"FIRST\")", "TPOT 1st Monitor");
  tpot->AddAction("tpotDraw(\"SECOND\")", "TPOT 2nd Monitor");
  tpot->AddAction("tpotDraw(\"HISTORY\")", "TPOT History Monitor");
  tpot->AddAction(new SubSystemActionDrawPS(tpot));
  pmf->RegisterSubSystem(tpot);  

  pmf->Draw();
  
}
