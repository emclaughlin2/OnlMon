// use #include "" only for your local include and put
// those in the first line(s) before any #include <>
// otherwise you are asking for weird behavior
// (more info - check the difference in include path search when using "" versus <>)

#include "TpcMon.h"

#include <onlmon/OnlMon.h>  // for OnlMon
#include <onlmon/OnlMonDB.h>
#include <onlmon/OnlMonServer.h>

#include <Event/Event.h>
#include <Event/msg_profile.h>
#include <Event/oncsSubConstants.h>

#include <tpc/TpcMap.h> 

#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TTree.h>
#include <TLatex.h>

#include <cmath>
#include <vector>
#include <cstdio>  // for printf
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>  // for allocator, string, char_traits

enum
{
  TRGMESSAGE = 1,
  FILLMESSAGE = 2
};

TpcMon::TpcMon(const std::string &name)
  : OnlMon(name)
{
  // leave ctor fairly empty, its hard to debug if code crashes already
  // during a new TpcMon()

  serverid = 0; //default case - initializing in constructor
  //BCO initialization in TPCMon
  starting_BCO = -1;
  rollover_value = 0;
  current_BCOBIN = 0;
  M.setMapNames("AutoPad-R1-RevA.sch.ChannelMapping.csv", "AutoPad-R2-RevA-Pads.sch.ChannelMapping.csv", "AutoPad-R3-RevA.sch.ChannelMapping.csv");

  return;
}

TpcMon::~TpcMon()
{
  // you can delete NULL pointers it results in a NOOP (No Operation)
  return;
}

int TpcMon::Init()
{
  // read our calibrations from TpcMonData.dat
  const char *tpccalib = getenv("TPCCALIB");
  if (!tpccalib)
  {
    std::cout << "TPCCALIB environment variable not set" << std::endl;
    exit(1);
  }
  std::string fullfile = std::string(tpccalib) + "/" + "TpcMonData.dat";
  std::ifstream calib(fullfile);
  calib.close();
  // use printf for stuff which should go the screen but not into the message
  // system (all couts are redirected)
  printf("doing the Init\n");
  tpchist1 = new TH1F("tpcmon_hist1", "test 1d histo", 101, 0., 100.);
  tpchist2 = new TH2F("tpcmon_hist2", "test 2d histo", 101, 0., 100., 101, 0., 100.);

  //TPC GEM Module Displays
  NorthSideADC = new TH2F("NorthSideADC" , "ADC Counts North Side", N_thBins, -TMath::Pi()/12. , 23.*TMath::Pi()/12. , N_rBins , rBin_edges );
  SouthSideADC = new TH2F("SouthSideADC" , "ADC Counts South Side", N_thBins, -TMath::Pi()/12. , 23.*TMath::Pi()/12. , N_rBins , rBin_edges );
  //

  //TPC "cluster" XY heat maps WEIGHTED
  //NorthSideADC_clusterXY_R1 = new TH2F("NorthSideADC_clusterXY_R1" , "ADC Peaks North Side", N_phi_binx_XY_R1, -TMath::Pi(), TMath::Pi(), N_rBins_XY, r_bins);
  NorthSideADC_clusterXY_R1 = new TH2F("NorthSideADC_clusterXY_R1" , "(ADC-Pedestal) > 20 North Side", 400, -800, 800, 400, -800, 800);
  NorthSideADC_clusterXY_R1->SetXTitle("X [mm]");
  NorthSideADC_clusterXY_R1->SetYTitle("Y [mm]");

  //NorthSideADC_clusterXY_R2 = new TH2F("NorthSideADC_clusterXY_R2" , "ADC Peaks North Side", N_phi_binx_XY_R2, -TMath::Pi(), TMath::Pi(), N_rBins_XY, r_bins);
  NorthSideADC_clusterXY_R2 = new TH2F("NorthSideADC_clusterXY_R2" , "(ADC-Pedestal) > 20 North Side", 400, -800, 800, 400, -800, 800);
  NorthSideADC_clusterXY_R2->SetXTitle("X [mm]");
  NorthSideADC_clusterXY_R2->SetYTitle("Y [mm]");

  //NorthSideADC_clusterXY_R3 = new TH2F("NorthSideADC_clusterXY_R3" , "ADC Peaks North Side", N_phi_binx_XY_R3, -TMath::Pi(), TMath::Pi(), N_rBins_XY, r_bins);
  NorthSideADC_clusterXY_R3 = new TH2F("NorthSideADC_clusterXY_R3" , "(ADC-Pedestal) > 20 North Side", 400, -800, 800, 400, -800, 800);
  NorthSideADC_clusterXY_R3->SetXTitle("X [mm]");
  NorthSideADC_clusterXY_R3->SetYTitle("Y [mm]");


  //SouthSideADC_clusterXY_R1 = new TH2F("SouthSideADC_clusterXY_R1" , "ADC Peaks South Side", N_phi_binx_XY_R1, -TMath::Pi(), TMath::Pi(), N_rBins_XY, r_bins);
  SouthSideADC_clusterXY_R1 = new TH2F("SouthSideADC_clusterXY_R1" , "(ADC-Pedestal) >20 South Side", 400, -800, 800, 400, -800, 800);
  SouthSideADC_clusterXY_R1->SetXTitle("X [mm]");
  SouthSideADC_clusterXY_R1->SetYTitle("Y [mm]");

  //SouthSideADC_clusterXY_R2 = new TH2F("SouthSideADC_clusterXY_R2" , "ADC Peaks South Side", N_phi_binx_XY_R2, -TMath::Pi(), TMath::Pi(), N_rBins_XY, r_bins);
  SouthSideADC_clusterXY_R2 = new TH2F("SouthSideADC_clusterXY_R2" , "(ADC-Pedestal) > 20 South Side", 400, -800, 800, 400, -800, 800);
  SouthSideADC_clusterXY_R2->SetXTitle("X [mm]");
  SouthSideADC_clusterXY_R2->SetYTitle("Y [mm]");

  //SouthSideADC_clusterXY_R3 = new TH2F("SouthSideADC_clusterXY_R3" , "ADC Peaks South Side", N_phi_binx_XY_R3, -TMath::Pi(), TMath::Pi(), N_rBins_XY, r_bins);  
  SouthSideADC_clusterXY_R3 = new TH2F("SouthSideADC_clusterXY_R3" , "(ADC-Pedestal) > 20 South Side", 400, -800, 800, 400, -800, 800);
  SouthSideADC_clusterXY_R3->SetXTitle("X [mm]");
  SouthSideADC_clusterXY_R3->SetYTitle("Y [mm]");

  //____________________________________________________________________//

  //TPC "cluster" XY heat maps UNWEIGHTED
  //NorthSideADC_clusterXY_R1 = new TH2F("NorthSideADC_clusterXY_R1" , "ADC Peaks North Side", N_phi_binx_XY_R1, -TMath::Pi(), TMath::Pi(), N_rBins_XY, r_bins);
  NorthSideADC_clusterXY_R1_unw = new TH2F("NorthSideADC_clusterXY_R1_unw" , "(ADC-Pedestal) > 20 North Side", 400, -800, 800, 400, -800, 800);
  NorthSideADC_clusterXY_R1_unw->SetXTitle("X [mm]");
  NorthSideADC_clusterXY_R1_unw->SetYTitle("Y [mm]");

  //NorthSideADC_clusterXY_R2 = new TH2F("NorthSideADC_clusterXY_R2" , "ADC Peaks North Side", N_phi_binx_XY_R2, -TMath::Pi(), TMath::Pi(), N_rBins_XY, r_bins);
  NorthSideADC_clusterXY_R2_unw = new TH2F("NorthSideADC_clusterXY_R2_unw" , "(ADC-Pedestal) > 20 North Side", 400, -800, 800, 400, -800, 800);
  NorthSideADC_clusterXY_R2_unw->SetXTitle("X [mm]");
  NorthSideADC_clusterXY_R2_unw->SetYTitle("Y [mm]");

  //NorthSideADC_clusterXY_R3 = new TH2F("NorthSideADC_clusterXY_R3" , "ADC Peaks North Side", N_phi_binx_XY_R3, -TMath::Pi(), TMath::Pi(), N_rBins_XY, r_bins);
  NorthSideADC_clusterXY_R3_unw = new TH2F("NorthSideADC_clusterXY_R3_unw" , "(ADC-Pedestal) > 20 North Side", 400, -800, 800, 400, -800, 800);
  NorthSideADC_clusterXY_R3_unw->SetXTitle("X [mm]");
  NorthSideADC_clusterXY_R3_unw->SetYTitle("Y [mm]");


  //SouthSideADC_clusterXY_R1 = new TH2F("SouthSideADC_clusterXY_R1" , "ADC Peaks South Side", N_phi_binx_XY_R1, -TMath::Pi(), TMath::Pi(), N_rBins_XY, r_bins);
  SouthSideADC_clusterXY_R1_unw = new TH2F("SouthSideADC_clusterXY_R1_unw" , "(ADC-Pedestal) >20 South Side", 400, -800, 800, 400, -800, 800);
  SouthSideADC_clusterXY_R1_unw->SetXTitle("X [mm]");
  SouthSideADC_clusterXY_R1_unw->SetYTitle("Y [mm]");

  //SouthSideADC_clusterXY_R2 = new TH2F("SouthSideADC_clusterXY_R2" , "ADC Peaks South Side", N_phi_binx_XY_R2, -TMath::Pi(), TMath::Pi(), N_rBins_XY, r_bins);
  SouthSideADC_clusterXY_R2_unw = new TH2F("SouthSideADC_clusterXY_R2_unw" , "(ADC-Pedestal) > 20 South Side", 400, -800, 800, 400, -800, 800);
  SouthSideADC_clusterXY_R2_unw->SetXTitle("X [mm]");
  SouthSideADC_clusterXY_R2_unw->SetYTitle("Y [mm]");

  //SouthSideADC_clusterXY_R3 = new TH2F("SouthSideADC_clusterXY_R3" , "ADC Peaks South Side", N_phi_binx_XY_R3, -TMath::Pi(), TMath::Pi(), N_rBins_XY, r_bins);  
  SouthSideADC_clusterXY_R3_unw = new TH2F("SouthSideADC_clusterXY_R3_unw" , "(ADC-Pedestal) > 20 South Side", 400, -800, 800, 400, -800, 800);
  SouthSideADC_clusterXY_R3_unw->SetXTitle("X [mm]");
  SouthSideADC_clusterXY_R3_unw->SetYTitle("Y [mm]");

  //____________________________________________________________________//

  //TPC "cluster" ZY heat maps WEIGHTED
   NorthSideADC_clusterZY = new TH2F("NorthSideADC_clusterZY" , "(ADC-Pedestal) > 20 North Side", 206, -1030, 1030, 400, -800, 800);
   SouthSideADC_clusterZY = new TH2F("SouthSideADC_clusterZY" , "(ADC-Pedestal) > 20 North Side", 206, -1030, 1030, 400, -800, 800);

  //____________________________________________________________________//

  //TPC "cluster" ZY heat maps UNWEIGHTED
   NorthSideADC_clusterZY_unw = new TH2F("NorthSideADC_clusterZY_unw" , "(ADC-Pedestal) > 20 South Side", 206, -1030, 1030, 400, -800, 800);
   SouthSideADC_clusterZY_unw = new TH2F("SouthSideADC_clusterZY_unw" , "(ADC-Pedestal) > 20 South Side", 206, -1030, 1030, 400, -800, 800);

  //____________________________________________________________________//

  //

  // ADC vs Sample (small)
  char ADC_vs_SAMPLE_str[100];
  char ADC_vs_SAMPLE_xaxis_str[100];
  sprintf(ADC_vs_SAMPLE_str,"ADC Counts vs Sample: SECTOR %i",MonitorServerId());
  sprintf(ADC_vs_SAMPLE_xaxis_str,"Sector %i: ADC Time bin [1/20MHz]",MonitorServerId());
  ADC_vs_SAMPLE = new TH2F("ADC_vs_SAMPLE", ADC_vs_SAMPLE_str, 360, 0, 360, 256, 0, 1024);
  ADC_vs_SAMPLE -> SetXTitle(ADC_vs_SAMPLE_xaxis_str);
  ADC_vs_SAMPLE -> SetYTitle("ADC [ADU]");

  // ADC vs Sample (large)
  char ADC_vs_SAMPLE_large_str[100];
  char ADC_vs_SAMPLE_xaxis_large_str[100];
  sprintf(ADC_vs_SAMPLE_large_str,"ADC Counts vs Large Sample: SECTOR %i",MonitorServerId());
  sprintf(ADC_vs_SAMPLE_xaxis_large_str,"Sector %i: ADC Time bin [1/20MHz]",MonitorServerId());
  ADC_vs_SAMPLE_large = new TH2F("ADC_vs_SAMPLE_large", ADC_vs_SAMPLE_large_str, 1080, 0, 1080, 256, 0, 1024);
  ADC_vs_SAMPLE_large -> SetXTitle(ADC_vs_SAMPLE_xaxis_large_str);
  ADC_vs_SAMPLE_large -> SetYTitle("ADC [ADU]");

  // Sample size distribution 1D histogram
  char sample_size_title_str[100];
  sprintf(sample_size_title_str,"Distribution of Sample Sizes in Events: SECTOR %i",MonitorServerId());
  sample_size_hist = new TH1F("sample_size_hist" , sample_size_title_str, 1000, 0.5, 1000.5);
  sample_size_hist->SetXTitle("sample size");
  sample_size_hist->SetYTitle("counts");

  // entries vs FEE*8 + SAMPA Number
  Check_Sums = new TH1F("Check_Sums" , "Entries vs Fee*8 + SAMPA in Events",208,-0.5, 207.5);
  Check_Sums->SetXTitle("FEE_NUM*8 + SAMPA_ADRR");
  Check_Sums->SetYTitle("Entries");
  Check_Sums->Sumw2(kFALSE); //explicity turn off Sumw2 - we do not want it

  // checksum error vs FEE*8 + SAMPA Number
  char checksum_title_str[100];
  sprintf(checksum_title_str,"Check Sum Error Probability vs Fee*8 + SAMPA in Events: SECTOR %i",MonitorServerId());
  Check_Sum_Error = new TH1F("Check_Sum_Error" , checksum_title_str,208,-0.5, 207.5);
  Check_Sum_Error->SetXTitle("FEE_NUM*8 + SAMPA_ADDR");
  Check_Sum_Error->SetYTitle("Prob. Check. Sum. Err.");
  Check_Sum_Error->Sumw2(kFALSE); //explicity turn off Sumw2 - we do not want it

  // Max ADC per waveform dist for each module (R1, R2, R3)
  char MAXADC_str[100];
  char YLabel_str[5];

  sprintf(MAXADC_str,"MAX ADC per Waveform in SLIDING WINDOW: SECTOR %i",MonitorServerId());

  MAXADC = new TH2F("MAXADC" , MAXADC_str,1025,-0.5, 1024.5,3,-0.5,2.5);
  MAXADC->SetXTitle("LocalMAX ADC in Waveform [ADU]");
  MAXADC->SetYTitle("Entries");
  MAXADC->Sumw2(kFALSE); //explicity turn off Sumw2 - we do not want it

  for(int i = 0; i < 3; i++ )
  {
    sprintf(YLabel_str,"R%i",i+1);
    MAXADC->GetYaxis()->SetBinLabel(i+1,YLabel_str);
    MAXADC->GetYaxis()->SetLabelSize(0.12);
    MAXADC->GetXaxis()->SetLabelSize(0.04);
  }

  //________For 1D ADC spectra - Doing this the hard way because I'm not sure if I can have a TH1 array
  char RAWADC_1D_titlestr[100];
  char MAXADC_1D_titlestr[100];

  sprintf(RAWADC_1D_titlestr,"RAW ADC for Sector %i R1",MonitorServerId());
  sprintf(MAXADC_1D_titlestr,"MAX ADC in SLIDING WINDOW for Sector %i R1",MonitorServerId());

  RAWADC_1D_R1 = new TH1F("RAWADC_1D_R1",RAWADC_1D_titlestr,1025,-0.5,1024.5);
  MAXADC_1D_R1 = new TH1F("MAXADC_1D_R1",MAXADC_1D_titlestr,1025,-0.5,1024.5);

  RAWADC_1D_R1->SetYTitle("Entries");
  RAWADC_1D_R1->SetXTitle("ADC [ADU]");

  MAXADC_1D_R1->SetYTitle("Entries");
  MAXADC_1D_R1->SetXTitle("(MAXADC-pedestal) [ADU]");

  MAXADC_1D_R1->Sumw2(kFALSE);
  RAWADC_1D_R1->Sumw2(kFALSE);

  MAXADC_1D_R1->SetLineColor(2);
  RAWADC_1D_R1->SetLineColor(2);

  sprintf(RAWADC_1D_titlestr,"RAW ADC for Sector %i R2",MonitorServerId());
  sprintf(MAXADC_1D_titlestr,"MAX ADC for Sector %i R2",MonitorServerId());

  RAWADC_1D_R2 = new TH1F("RAWADC_1D_R2",RAWADC_1D_titlestr,1025,-0.5,1024.5);
  MAXADC_1D_R2 = new TH1F("MAXADC_1D_R2",MAXADC_1D_titlestr,1025,-0.5,1024.5);

  RAWADC_1D_R2->SetYTitle("Entries");
  RAWADC_1D_R2->SetXTitle("ADC [ADU]");

  MAXADC_1D_R2->SetYTitle("Entries");
  MAXADC_1D_R2->SetXTitle("(MAXADC-pedestal) [ADU]");

  MAXADC_1D_R2->Sumw2(kFALSE);
  RAWADC_1D_R2->Sumw2(kFALSE);

  MAXADC_1D_R2->SetLineColor(3);
  RAWADC_1D_R2->SetLineColor(3);

  sprintf(RAWADC_1D_titlestr,"RAW ADC for Sector %i R3",MonitorServerId());
  sprintf(MAXADC_1D_titlestr,"MAX ADC for Sector %i R3",MonitorServerId());

  RAWADC_1D_R3 = new TH1F("RAWADC_1D_R3",RAWADC_1D_titlestr,1025,-0.5,1024.5);
  MAXADC_1D_R3 = new TH1F("MAXADC_1D_R3",MAXADC_1D_titlestr,1025,-0.5,1024.5);

  RAWADC_1D_R3->SetYTitle("Entries");
  RAWADC_1D_R3->SetXTitle("ADC [ADU]");

  MAXADC_1D_R3->SetYTitle("Entries");
  MAXADC_1D_R3->SetXTitle("(MAXADC-pedestal) [ADU]");

  MAXADC_1D_R3->Sumw2(kFALSE);  // ADC vs Sample (small)
  RAWADC_1D_R3->Sumw2(kFALSE);

  MAXADC_1D_R3->SetLineColor(4);
  RAWADC_1D_R3->SetLineColor(4);

  OnlMonServer *se = OnlMonServer::instance();
  // register histograms with server otherwise client won't get them
  se->registerHisto(this, tpchist1);  // uses the TH1->GetName() as key
  se->registerHisto(this, tpchist2);
  se->registerHisto(this, NorthSideADC);
  se->registerHisto(this, SouthSideADC);
  se->registerHisto(this, sample_size_hist);
  se->registerHisto(this, Check_Sum_Error);
  se->registerHisto(this, Check_Sums);
  se->registerHisto(this, ADC_vs_SAMPLE);
  se->registerHisto(this, ADC_vs_SAMPLE_large);
  se->registerHisto(this, MAXADC);
  se->registerHisto(this, RAWADC_1D_R1);
  se->registerHisto(this, MAXADC_1D_R1);
  se->registerHisto(this, RAWADC_1D_R2);
  se->registerHisto(this, MAXADC_1D_R2);
  se->registerHisto(this, RAWADC_1D_R3);
  se->registerHisto(this, MAXADC_1D_R3);

  se->registerHisto(this, NorthSideADC_clusterXY_R1);
  se->registerHisto(this, NorthSideADC_clusterXY_R2);
  se->registerHisto(this, NorthSideADC_clusterXY_R3);

  se->registerHisto(this, SouthSideADC_clusterXY_R1);
  se->registerHisto(this, SouthSideADC_clusterXY_R2);
  se->registerHisto(this, SouthSideADC_clusterXY_R3);

  se->registerHisto(this, NorthSideADC_clusterXY_R1_unw);
  se->registerHisto(this, NorthSideADC_clusterXY_R2_unw);
  se->registerHisto(this, NorthSideADC_clusterXY_R3_unw);

  se->registerHisto(this, SouthSideADC_clusterXY_R1_unw);
  se->registerHisto(this, SouthSideADC_clusterXY_R2_unw);
  se->registerHisto(this, SouthSideADC_clusterXY_R3_unw);

  se->registerHisto(this, NorthSideADC_clusterZY);
  se->registerHisto(this, SouthSideADC_clusterZY);

  se->registerHisto(this, NorthSideADC_clusterZY_unw);
  se->registerHisto(this, SouthSideADC_clusterZY_unw);


  Reset();
  return 0;
}

int TpcMon::BeginRun(const int /* runno */)
{
  // if you need to read calibrations on a run by run basis
  // this is the place to do it

  // we reset the BCO for the new run
  starting_BCO = -1;
  rollover_value = 0;
  current_BCOBIN = 0;

  return 0;
}

int TpcMon::process_event(Event *evt/* evt */)
{

  //std::cout << "TpcMon::process_event(Event * evt) Processing Event" << std::endl;

  if (evt == nullptr)
  {
    std::cout << "TpcMon::process_event - Event not found" << std::endl;
    return -1;
  }

  if (evt->getEvtType() >= 8)  /// special events
  {
    std::cout << "TpcMon::process_event - Special Event type >= 8, moving on" << std::endl;
    return -1;
  }

  //reset these each event
  float North_Side_Arr[36] = {0};
  float South_Side_Arr[36] = {0};

  std::vector<int> store_ten; 
  std::vector<int> mean_and_stdev_vec;


  // we check if we have legacy data and start with packet 4000
  // the range for the TPC is really 4001...4032
  // we assume we start properly at 4001, but check if not
  
  int firstpacket=4001;
  if (evt->existPacket(4000))
    {
      Packet *p = evt->getPacket(4000);
      if (p->getIdentifier() == IDTPCFEEV3 ) firstpacket = 4000;
      delete p;
    }
  int lastpacket = firstpacket+232;

  
  for( int packet = firstpacket; packet < lastpacket; packet++) //packet 4001 or 4002 = Sec 00, packet 4231 or 4232 = Sec 23
  {
    Packet* p = evt->getPacket(packet);
    if (!p)
    {
      //std::cout << "TpcMon::process_event - No packet numbered " << packet << " in this event!!" << std::endl;
      continue;
    }
    else
    {
      //std::cout << "____________________________________" << std::endl;
      //std::cout << "Packet # " << packet << std::endl;
      int nr_of_waveforms = p->iValue(0, "NR_WF");
      //std::cout << "Hello Waveforms ! - There are " << nr_of_waveforms << " of you !" << std::endl;

      bool is_channel_stuck = 0;

      for( int wf = 0; wf < nr_of_waveforms; wf++)
      {

        int current_BCO = p->iValue(wf, "BCO") + rollover_value;
        if (starting_BCO < 0)
        {
          starting_BCO = current_BCO;
        }

        if (current_BCO < starting_BCO)  // we have a rollover
        {
          rollover_value += 0x100000;
          current_BCO = p->iValue(wf, "BCO") + rollover_value;
          starting_BCO = current_BCO;
          current_BCOBIN++;
        }


        int fee = p->iValue(wf, "FEE");
        int sampaAddress = p->iValue(wf, "SAMPAADDRESS");
        int checksumError = p->iValue(wf, "CHECKSUMERROR");
        int channel = p->iValue(wf, "CHANNEL");

        Check_Sums->Fill(fee*8 + sampaAddress); 
        if( checksumError == 1){Check_Sum_Error->Fill(fee*8 + sampaAddress);}

        int nr_Samples = p->iValue(wf, "SAMPLES");
        sample_size_hist->Fill(nr_Samples);


        // clockwise FEE mapping
        //int FEE_map[26]={5, 6, 1, 3, 2, 12, 10, 11, 9, 8, 7, 1, 2, 4, 8, 7, 6, 5, 4, 3, 1, 3, 2, 4, 6, 5};
        int FEE_R[26]={2, 2, 1, 1, 1, 3, 3, 3, 3, 3, 3, 2, 2, 1, 2, 2, 1, 1, 2, 2, 3, 3, 3, 3, 3, 3};
        // counter clockwise FEE mapping (From Takao - DEPRECATED AS OF 08.29)
        //int FEE_map[26]={3, 2, 5, 3, 4, 0, 2, 1, 3, 4, 5, 7, 6, 2, 0, 1, 0, 1, 4, 5, 11, 9, 10, 8, 6, 7};

        // FEE mapping from Jin
        int FEE_map[26]={4, 5, 0, 2, 1, 11, 9, 10, 8, 7, 6, 0, 1, 3, 7, 6, 5, 4, 3, 2, 0, 2, 1, 3, 5, 4};

        
        //int pads_per_sector[3] = {96, 128, 192};


        serverid = MonitorServerId();

        // setting the mapp of the FEE
        int feeM = FEE_map[fee];
        if(FEE_R[fee]==2) feeM += 6;
        if(FEE_R[fee]==3) feeM += 14;

        // getting R and Phi coordinates
        double R = M.getR(feeM, channel);
        int layer = M.getLayer(feeM, channel);
        double phi = 0;

        //double phi = M.getPhi(feeM, channel) + (serverid - 12*side(serverid)) * M_PI / 6 ;

        if( side(serverid) == 0 ) //NS
        {
          phi = M.getPhi(feeM, channel) + (serverid ) * M_PI / 6 ; 
        }
        else if( side(serverid) == 1 ) //SS
        {
          phi = M.getPhi(feeM, channel) + (18 - serverid ) * M_PI / 6 ; 
        }


        //std::cout<<"Sector = "<< serverid <<" FEE = "<<fee<<" channel = "<<channel<<std::endl;

        int mid = floor(nr_Samples/2); //get median sample

        if( nr_Samples > 9)
        {
          if( (p->iValue(wf,mid) == p->iValue(wf,mid-1)) && (p->iValue(wf,mid) == p->iValue(wf,mid-2)) && (p->iValue(wf,mid) == p->iValue(wf,mid+1)) && (p->iValue(wf,mid) == p->iValue(wf,mid+2)) )     
          {
            is_channel_stuck = 1;
          }
          for( int si=0;si < 10; si++ ) //get pedestal and noise before hand
          {
            mean_and_stdev_vec.push_back(p->iValue(wf,si));
          }
        } //Compare 5 values to determine stuck !!

        std::pair<float, float> result = calculateMeanAndStdDev(mean_and_stdev_vec);
        int pedestal = result.first; //average/pedestal
        //int noise = result.second; //stdev/noise

        int wf_max = 0;
        int t_max = 0;

        for( int s =0; s < nr_Samples ; s++ )
        {
          
          //int t = s + 2 * (current_BCO - starting_BCO);

          int adc = p->iValue(wf,s);       

          if( adc > wf_max){ wf_max = adc; t_max = s; }

          if( s >= 10 && s <= 19) // get first 10-19
          {
            store_ten.push_back(adc);
          }
          else if( s > 19 )  
          {

            //nine_max = Max_Nine(p->iValue(wf,s-9),p->iValue(wf,s-8),p->iValue(wf,s-7),p->iValue(wf,s-6),p->iValue(wf,s-5),p->iValue(wf,s-4),p->iValue(wf,s-3),p->iValue(wf,s-2),p->iValue(wf,s-1)); //take the previous 9 numbers

            //VERY IMPORTANT, erase first entry, push_back current
            store_ten.erase(store_ten.begin());
            store_ten.push_back(adc);
             
            int max_of_previous_10 = *max_element(store_ten.begin(), store_ten.end());

            if(adc == max_of_previous_10 && (checksumError == 0 && is_channel_stuck == 0)) //if the new value is greater than the previous 9
            {
               MAXADC->Fill(adc - pedestal,Module_ID(fee)); 
               if(Module_ID(fee)==0){MAXADC_1D_R1->Fill(adc - pedestal);} //Raw 1D for R1
               else if(Module_ID(fee)==1){MAXADC_1D_R2->Fill(adc - pedestal);} //Raw 1D for R2
               else if(Module_ID(fee)==2){MAXADC_1D_R3->Fill(adc - pedestal);} //Raw 1D for R3
            }

          }

          if( checksumError == 0 && is_channel_stuck == 0)
          {
            ADC_vs_SAMPLE -> Fill(s, adc);
            ADC_vs_SAMPLE_large -> Fill(s, adc);

            if(Module_ID(fee)==0){RAWADC_1D_R1->Fill(adc);} //Raw 1D for R1
            else if(Module_ID(fee)==1){RAWADC_1D_R2->Fill(adc);} //Raw 1D for R2
            else if(Module_ID(fee)==2){RAWADC_1D_R3->Fill(adc);} //Raw 1D for R3
          }

          //increment 
          if(serverid >= 0 && serverid < 12 ){ North_Side_Arr[ Index_from_Module(serverid,fee) ] += adc;}
          else {South_Side_Arr[ Index_from_Module(serverid,fee)-36 ] += adc;}

        } //nr samples

        //for complicated XY stuff ____________________________________________________
        //20 = 3-5 * sigma - hard-coded
        // OR 10*noise = 10 sigma

        float z = 0; //mm

        if( (serverid < 12 && (wf_max - pedestal) > 20) && layer != 0 )
        {
          if(Module_ID(fee)==0){NorthSideADC_clusterXY_R1->Fill(R*cos(phi),R*sin(phi),wf_max - pedestal);NorthSideADC_clusterXY_R1_unw->Fill(R*cos(phi),R*sin(phi));} //Raw 1D for R1
          else if(Module_ID(fee)==1){NorthSideADC_clusterXY_R2->Fill(R*cos(phi),R*sin(phi),wf_max - pedestal);NorthSideADC_clusterXY_R2_unw->Fill(R*cos(phi),R*sin(phi));} //Raw 1D for R2
          else if(Module_ID(fee)==2){NorthSideADC_clusterXY_R3->Fill(R*cos(phi),R*sin(phi),wf_max - pedestal);NorthSideADC_clusterXY_R3_unw->Fill(R*cos(phi),R*sin(phi));} //Raw 1D for R3

          if( t_max >= 10 && t_max <=255 ){z = 1030 - (t_max - 10)*(50 * 0.084);NorthSideADC_clusterZY->Fill(z,R*sin(phi),wf_max - pedestal);NorthSideADC_clusterZY_unw->Fill(z,R*sin(phi));}
        }
        else if( (serverid >=12 && (wf_max - pedestal) > 20) && layer != 0)
        {
          if(Module_ID(fee)==0){SouthSideADC_clusterXY_R1->Fill(R*cos(phi),R*sin(phi),wf_max - pedestal);SouthSideADC_clusterXY_R1_unw->Fill(R*cos(phi),R*sin(phi));} //Raw 1D for R1
          else if(Module_ID(fee)==1){SouthSideADC_clusterXY_R2->Fill(R*cos(phi),R*sin(phi),wf_max - pedestal);SouthSideADC_clusterXY_R2_unw->Fill(R*cos(phi),R*sin(phi));} //Raw 1D for R2
          else if(Module_ID(fee)==2){SouthSideADC_clusterXY_R3->Fill(R*cos(phi),R*sin(phi),wf_max - pedestal);SouthSideADC_clusterXY_R3_unw->Fill(R*cos(phi),R*sin(phi));} //Raw 1D for R3

          if( t_max >= 10 && t_max <=255 ){z = -1030 + (t_max - 10)*(50 * 0.084);SouthSideADC_clusterZY->Fill(z,R*sin(phi),wf_max - pedestal);SouthSideADC_clusterZY_unw->Fill(z,R*sin(phi));}
        }
        //________________________________________________________________________________

        is_channel_stuck = 0; //reset after looping through waveform samples

        store_ten.clear(); //clear this after every waveform
        
        mean_and_stdev_vec.clear(); //clear this after every waveform


      } //nr waveforms
      delete p;
    }
  } //packet loop

  evtcnt++;

  // get temporary pointers to histograms
  // one can do in principle directly se->getHisto("tpchist1")->Fill()
  // but the search in the histogram Map is somewhat expensive and slows
  // things down if you make more than one operation on a histogram
  tpchist1->Fill((float) idummy);
  tpchist2->Fill((float) idummy, (float) idummy, 1.);

  //fill the TPC module displays
  float r, theta;

  //dummy data
  //float North_Side_Arr[36] = { 12, 8, 40, 39, 80, 50, 12, 8, 40, 39, 80, 50, 12, 8, 40, 39, 80, 50, 12, 8, 40, 39, 80, 50, 12, 8, 40, 39, 80, 50, 12, 8, 40, 39, 80, 50 };
  //float South_Side_Arr[36] = { 12, 8, 40, 39, 80, 50, 12, 8, 40, 39, 80, 50, 12, 8, 40, 39, 80, 50, 12, 8, 40, 39, 80, 50, 12, 8, 40, 39, 80, 50, 12, 8, 40, 39, 80, 50 };

  for(int tpciter = 1; tpciter < 73 ; tpciter++){

    Locate(tpciter, &r, &theta);
    //std::cout << "r is: "<< r <<" theta is: "<< theta <<"\n";
    if(tpciter < 37){ //South side
      NorthSideADC->Fill(theta,r, North_Side_Arr[tpciter-1]); //fill South side with the weight = bin content
    }
    else { //North side
      SouthSideADC->Fill(theta,r,South_Side_Arr[tpciter-37]); //fill North side with the weight = bin content
    }
  }
  //

  return 0;
}

int TpcMon::Module_ID(int fee_id) //for simply determining which module you are in (doesn't care about sector)
{
  int mod_id;
 
  if( fee_id == 2 || fee_id == 3 || fee_id == 4 || fee_id == 13 || fee_id == 16 || fee_id == 17 ){mod_id = 0;} // R1 

  else if( fee_id == 0 || fee_id == 1 || fee_id == 11 || fee_id == 12 || fee_id == 14 || fee_id == 15 || fee_id == 18 || fee_id == 19 ){mod_id = 1;} // R2

  else if( fee_id == 5 || fee_id == 6 || fee_id ==7 || fee_id == 8 || fee_id == 9 || fee_id == 10 || fee_id == 20 || fee_id == 21 || fee_id == 22 || fee_id == 23 || fee_id == 24 || fee_id == 25 ){mod_id = 2;} // R3

  else mod_id = 0;

  return mod_id;
}


int TpcMon::Index_from_Module(int sec_id, int fee_id) //for placing in the array (takes into account sector)
{
  int mod_id;
 
  if( fee_id == 2 || fee_id == 3 || fee_id == 4 || fee_id == 13 || fee_id == 16 || fee_id == 17 ){mod_id = 3*sec_id + 0;} // R1 

  else if( fee_id == 0 || fee_id == 1 || fee_id == 11 || fee_id == 12 || fee_id == 14 || fee_id == 15 || fee_id == 18 || fee_id == 19 ){mod_id = 3*sec_id + 1;} // R2

  else if( fee_id == 5 || fee_id == 6 || fee_id ==7 || fee_id == 8 || fee_id == 9 || fee_id == 10 || fee_id == 20 || fee_id == 21 || fee_id == 22 || fee_id == 23 || fee_id == 24 || fee_id == 25 ){mod_id = 3*sec_id + 2;} // R3

  else mod_id = 0;

  return mod_id;
}

void TpcMon::Locate(int id, float *rbin, float *thbin)
{
   float CSIDE_angle_bins[12] = { 0.1*2.*TMath::Pi()/12 , 1.1*2.*TMath::Pi()/12 , 2.1*2.*TMath::Pi()/12 , 3.1*2.*TMath::Pi()/12 , 4.1*2.*TMath::Pi()/12 , 5.1*2.*TMath::Pi()/12 , 6.1*2.*TMath::Pi()/12 , 7.1*2.*TMath::Pi()/12 , 8.1*2.*TMath::Pi()/12 , 9.1*2.*TMath::Pi()/12 , 10.1*2.*TMath::Pi()/12 , 11.1*2.*TMath::Pi()/12 }; //CCW from x = 0 (RHS horizontal)

   float ASIDE_angle_bins[12] = { 6.1*2.*TMath::Pi()/12 , 5.1*2.*TMath::Pi()/12 , 4.1*2.*TMath::Pi()/12 , 3.1*2.*TMath::Pi()/12 , 2.1*2.*TMath::Pi()/12 , 1.1*2.*TMath::Pi()/12 , 0.1*2.*TMath::Pi()/12 , 11.1*2.*TMath::Pi()/12 , 10.1*2.*TMath::Pi()/12 , 9.1*2.*TMath::Pi()/12 , 8.1*2.*TMath::Pi()/12 , 7.1*2.*TMath::Pi()/12  }; //CCW from x = 0 (RHS horizontal)

   int modid3 = id % 3;

   switch(modid3) {
     case 1:
       *rbin = 0.4; //R1
       break;
     case 2:
       *rbin = 0.6; //R2
       break;
     case 0:
       *rbin = 0.8; //R3
       break;
   }

  if( id < 37){
    *thbin = CSIDE_angle_bins[TMath::FloorNint((id-1)/3)];
  }
  else if( id >= 37){
    *thbin = ASIDE_angle_bins[TMath::FloorNint((id-37)/3)];
  }
}

int TpcMon::Max_Nine(int one, int two, int three, int four, int five, int six, int seven, int eight, int nine)
{ 
  int max = 0;
  int nine_array[9] = {one, two, three, four, five, six, seven, eight, nine};

  for( int i = 0; i < 9; i++ )
  {
    if( nine_array[i] > max ){max = nine_array[i];}
  }

  return max;
}

bool TpcMon::side(int server_id)
{
  bool side_id = 0; // side = 0 for NS, side = 1 for SS
  if(server_id >= 12){side_id=1;} // side = 1 when serverid 12 or more
  
  return side_id;
}

std::pair<float, float> TpcMon::calculateMeanAndStdDev(const std::vector<int>& values) {
    // Calculate the mean
    float sum = 0.0;
    for (const auto& value : values) {
        sum += value;
    }
    float mean = sum / values.size();

    // Calculate the standard deviation
    float squaredSum = 0.0;
    for (const auto& value : values) {
        double diff = value - mean;
        squaredSum += diff * diff;
    }
    float variance = squaredSum / values.size();
    float stdDev = std::sqrt(variance);

    return std::make_pair(mean, stdDev);
}



int TpcMon::Reset()
{
  // reset our internal counters
  evtcnt = 0;
  idummy = 0;
  return 0;
}

