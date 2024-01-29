// Phenix Online Monitoring System
//
// Macro loader for starting Online Monitoring GUI Framework
//
// Author:  Cody McCain (mcm99c@acu.edu)
// Date:    January 2, 2003
//
// See run_Poms.C.README

#include <onlmon/Poms.h>

// cppcheck-suppress unknownMacro
R__LOAD_LIBRARY(libpoms.so)

void StartPoms()
{
  PomsMainFrame *pmf;
  pmf = PomsMainFrame::Instance();

  // Register SubSystems
  //  pmf->RegisterSubSystem("EXAMPLE", "example");
  SubSystem *subsys = nullptr;

  subsys = new SubSystem("MBD", "bbc");
  subsys->AddAction("bbcDraw(\"FIRST\")", "MBD Vertex Monitor");
  subsys->AddAction("bbcDraw(\"SECOND\")", "MBD Timing Monitor");
  subsys->AddAction(new SubSystemActionSavePlot(subsys));
  pmf->RegisterSubSystem(subsys);

  subsys = new SubSystem("CEMC", "cemc");
  subsys->AddAction("cemcDraw(\"FIRST\")", "Towers");
  subsys->AddAction("cemcDraw(\"SECOND\")", "Packet Health");
  subsys->AddAction("cemcDraw(\"THIRD\")", "Wave Forms");
  subsys->AddAction("cemcDraw(\"FOURTH\")", "Wave vs Fast Fitting");
  subsys->AddAction(new SubSystemActionSavePlot(subsys));
  pmf->RegisterSubSystem(subsys);

  // subsys = new SubSystem("DAQ", "daq");
  // subsys->AddAction("daqDraw(\"FIRST\")", "DAQ 1st Monitor");
  // subsys->AddAction("daqDraw(\"SECOND\")", "DAQ 2nd Monitor");
  // subsys->AddAction("daqDraw(\"HISTORY\")", "DAQ History Monitor");
  // subsys->AddAction(new SubSystemActionSavePlot(subsys));
  // pmf->RegisterSubSystem(subsys);

  // subsys = new SubSystem("EPD", "epd");
  // subsys->AddAction("epdDraw(\"FIRST\")", "EPD 1st Monitor");
  // subsys->AddAction("epdDraw(\"SECOND\")", "EPD 2nd Monitor");
  // subsys->AddAction("epdDraw(\"THIRD\")", "EPD 3rd Monitor");
  // subsys->AddAction(new SubSystemActionSavePlot(subsys));
  // pmf->RegisterSubSystem(subsys);

  subsys = new SubSystem("Inner HCAL", "ihcal");
  subsys->AddAction("ihcalDraw(\"FIRST\")", "Towers");
  subsys->AddAction("ihcalDraw(\"SECOND\")", "Sector Average");
  subsys->AddAction("ihcalDraw(\"THIRD\")", "Wave Form");
  subsys->AddAction("ihcalDraw(\"FOURTH\")", "Packet Health");
  subsys->AddAction(new SubSystemActionSavePlot(subsys));
  pmf->RegisterSubSystem(subsys);

  subsys = new SubSystem("Outer HCAL", "ohcal");
  subsys->AddAction("ohcalDraw(\"FIRST\")", "Towers");
  subsys->AddAction("ohcalDraw(\"SECOND\")", "Sector Average");
  subsys->AddAction("ohcalDraw(\"THIRD\")", "Wave Form");
  subsys->AddAction("ohcalDraw(\"FOURTH\")", "Packet Health");
  subsys->AddAction(new SubSystemActionSavePlot(subsys));
  pmf->RegisterSubSystem(subsys);

  subsys = new SubSystem("INTT", "intt");
  subsys->AddAction("inttDraw(\"chip_hitmap\")", "Chip Hitmap");
  subsys->AddAction("inttDraw(\"ladder_hitmap\")", "Ladder Hitmap");
  subsys->AddAction("inttDraw(\"chip_nll\")", "Chip NLL");
  subsys->AddAction(new SubSystemActionSavePlot(subsys));
  pmf->RegisterSubSystem(subsys);

  subsys = new SubSystem("Local Level 1", "ll1");
  subsys->AddAction("ll1Draw(\"FIRST\")", "Hits");
  subsys->AddAction("ll1Draw(\"SECOND\")", "Hits Correlation");
  subsys->AddAction("ll1Draw(\"THIRD\")", "Hits Lineup");
  subsys->AddAction(new SubSystemActionSavePlot(subsys));
  pmf->RegisterSubSystem(subsys);
  
  // subsys = new SubSystem("MVTX", "mvtx");
  // subsys->AddAction("mvtxDraw(\"GENERAL\")", "General Monitor");
  // subsys->AddAction("mvtxDraw(\"FEE\")", "FEE");
  // subsys->AddAction("mvtxDraw(\"FHR\")", "FHR");
  // subsys->AddAction("mvtxDraw(\"OCC\")", "OCC");
  // subsys->AddAction(new SubSystemActionSavePlot(subsys));
  // pmf->RegisterSubSystem(subsys);

  // subsys = new SubSystem("SEPD", "sepd");
  // subsys->AddAction("sepdDraw(\"FIRST\")", "Towers");
  // subsys->AddAction("sepdDraw(\"SECOND\")", "Packet Health");
  // subsys->AddAction("sepdDraw(\"FOURTH\")", "Wave Forms");
  // subsys->AddAction("sepdDraw(\"FIFTH\")", "Wave vs Fast Fitting");
  // subsys->AddAction(new SubSystemActionSavePlot(subsys));
  // pmf->RegisterSubSystem(subsys);

  subsys = new SubSystem("TPC", "tpc");
  subsys->AddAction("tpcDraw(\"TPCMODULE\")", "TPC SUM[ADC]");
  subsys->AddAction("tpcDraw(\"TPCSAMPLESIZE\")", "TPC Sample Size");
  subsys->AddAction("tpcDraw(\"TPCCHECKSUMERROR\")", "TPC Checksum Error Prob.");
  subsys->AddAction("tpcDraw(\"TPCADCVSSAMPLE\")", "TPC ADC vs. Sample");
  subsys->AddAction("tpcDraw(\"TPCADCVSSAMPLELARGE\")", "TPC ADC vs. Large Sample");
  subsys->AddAction("tpcDraw(\"TPCMAXADCMODULE\")", "TPC MAX ADC");
  subsys->AddAction("tpcDraw(\"TPCRAWADC1D\")", "TPC RAW ADC 1D");
  subsys->AddAction("tpcDraw(\"TPCMAXADC1D\")", "TPC (WindowMAX-Pedestal) ADC 1D");
  subsys->AddAction("tpcDraw(\"TPCCLUSTERSXYWEIGTHED\")","TPC Cluster XY MaxADC-Pedestal,w");
  subsys->AddAction("tpcDraw(\"TPCCLUSTERSXYUNWEIGTHED\")","TPC Cluster XY MaxADC-Pedestal,u");
  subsys->AddAction(new SubSystemActionSavePlot(subsys));
  pmf->RegisterSubSystem(subsys);

  subsys = new SubSystem("TPOT", "tpot");
  subsys->AddAction("tpotDraw(\"TPOT_counters\")", "Event counters");
  subsys->AddAction("tpotDraw(\"TPOT_detector_occupancy\")", "Detector Occupancy");
  subsys->AddAction("tpotDraw(\"TPOT_resist_occupancy\")", "Resist Occupnacy");
  subsys->AddAction("tpotDraw(\"TPOT_adc_vs_sample\")", "ADC Sample");
  subsys->AddAction("tpotDraw(\"TPOT_counts_vs_sample\")", "Counts Sample");
  subsys->AddAction("tpotDraw(\"TPOT_hit_charge\")", "Hit Charge");
  subsys->AddAction("tpotDraw(\"TPOT_hit_multiplicity\")", "Hit Multiplicity");
  subsys->AddAction("tpotDraw(\"TPOT_hit_vs_channel\")", "Hit Channel");
  subsys->AddAction(new SubSystemActionSavePlot(subsys));
  pmf->RegisterSubSystem(subsys);

  // subsys = new SubSystem("ZDC", "zdc");
  // subsys->AddAction("zdcDraw(\"FIRST\")", "ZDC 1st Monitor");
  // subsys->AddAction("zdcDraw(\"SECOND\")", "ZDC 2nd Monitor");
  // subsys->AddAction("zdcDraw(\"THIRD\")", "ZDC 3rd Monitor");
  // subsys->AddAction(new SubSystemActionSavePlot(subsys));
  // pmf->RegisterSubSystem(subsys);


  pmf->Draw();
}
