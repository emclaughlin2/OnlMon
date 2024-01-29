// use #include "" only for your local include and put
// those in the first line(s) before any #include <>
// otherwise you are asking for weird behavior
// (more info - check the difference in include path search when using "" versus <>)

#include "MvtxMon.h"

//#include <fun4allraw/SingleMvtxInput.h>
//#include <ffarawobjects/MvtxRawHitContainerv1.h>
//#include <ffarawobjects/MvtxRawHitv1.h>
//#include <ffarawobjects/MvtxRawEvtHeaderv1.h>

#include <onlmon/OnlMon.h>  // for OnlMon
#include <onlmon/OnlMonDB.h>
#include <onlmon/OnlMonServer.h>

#include <Event/msg_profile.h>

#include <TH1.h>
#include <TString.h>
#include <TList.h>
#include <TLatex.h>
#include <TH2.h>
#include <TH2Poly.h>
#include <TLine.h>
#include <TH3.h>

#include <Event/Event.h>
#include <Event/packet.h>

#include <cmath>
#include <cstdio>  // for printf
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>  // for allocator, string, char_traits
#include <utility>



enum
{
  TRGMESSAGE = 1,
  FILLMESSAGE = 2
};

MvtxMon::MvtxMon(const std::string &name)
  : OnlMon(name)
{
  // leave ctor fairly empty, its hard to debug if code crashes already
  // during a new MvtxMon()
  return;
}

MvtxMon::~MvtxMon()
{
  // you can delete NULL pointers it results in a NOOP (No Operation)
  return;
}

int MvtxMon::Init()
{
  // read our calibrations from MvtxMonData.dat
  const char *mvtxcalib = getenv("MVTXCALIB");
  if (!mvtxcalib)
  {
    std::cout << "MVTXCALIB environment variable not set" << std::endl;
    exit(1);
  }
  std::string fullfile = std::string(mvtxcalib) + "/" + "MvtxMonData.dat";
  std::ifstream calib(fullfile);
  calib.close();
  // use printf for stuff which should go the screen but not into the message
  // system (all couts are redirected)
  printf("doing the Init\n");
  OnlMonServer *se = OnlMonServer::instance();

  // register histograms with server otherwise client won't get them
  mvtxmon_ChipStaveOcc = new TH2D("MVTXMON_RawData_ChipStaveOcc","Average Occupancy: Stave Number and Chip Number",NCHIP,-0.5,NCHIP-0.5,NSTAVE,-0.5,NSTAVE-0.5);
  mvtxmon_ChipStave1D = new TH1D("MVTXMON_RawData_ChipStave1D","Average Occupancy per Chip Stave",NCHIP*NSTAVE,-0.5,NCHIP * NSTAVE-0.5);
  mvtxmon_ChipFiredHis= new TH1D("MVTXMON_RawData_ChipFiredHis","Number of Chips Fired per Felix per RCDAQ event Distribution",NCHIP*NSTAVE/6,-0.5,NCHIP*NSTAVE/6 - 0.5);
  mvtxmon_EvtHitChip= new TH1D("MVTXMON_RawData_EvtHitChip","Number of Hits Per Strobe Per Chip Distribution",101,-0.5,100.5);
  mvtxmon_EvtHitDis = new TH1D("MVTXMON_RawData_EvtHitDis","Number of Hits Per Strobe Distribution",7201,-0.5,7200.5);

  mvtxmon_ChipStaveOcc->SetStats(0);
  mvtxmon_ChipStave1D->SetStats(0);
  mvtxmon_ChipFiredHis->SetStats(0);
  mvtxmon_EvtHitChip->SetStats(0);
  mvtxmon_EvtHitDis->SetStats(0);

  mvtxmon_ChipStave1D->GetXaxis()->SetTitle("Chip*Stave");
  mvtxmon_ChipStave1D->GetYaxis()->SetTitle("Occupancy");

  mvtxmon_ChipFiredHis->GetXaxis()->SetTitle("Number of Chips");
  mvtxmon_ChipFiredHis->GetYaxis()->SetTitle("Counts");

  mvtxmon_EvtHitChip->GetXaxis()->SetTitle("Number of Hits");
  mvtxmon_EvtHitChip->GetYaxis()->SetTitle("Counts");

  mvtxmon_EvtHitDis->GetXaxis()->SetTitle("Number of Hits");
  mvtxmon_EvtHitDis->GetYaxis()->SetTitle("Counts");

  se->registerHisto(this, mvtxmon_ChipStaveOcc);
  se->registerHisto(this, mvtxmon_ChipStave1D);
  se->registerHisto(this, mvtxmon_ChipFiredHis);
  se->registerHisto(this, mvtxmon_EvtHitChip);
  se->registerHisto(this, mvtxmon_EvtHitDis);


  mRCDAQevt = new TH1I("RCDAQ_evt", "Number of RCDAQ events processed", 6,-0.5,5.5);
  mRCDAQevt->GetXaxis()->SetTitle("FELIX server");
  mRCDAQevt->GetYaxis()->SetTitle("Counts");
  mRCDAQevt->SetStats(0);
  se->registerHisto(this,  mRCDAQevt);

  for (int i = 0; i < NFlags; i++) {
    mvtxmon_LaneStatusOverview[i] = new TH2Poly();
    mvtxmon_LaneStatusOverview[i]->SetName(Form("MVTXMON_LaneStatus_laneStatusOverviewFlag%s", mLaneStatusFlag[i].c_str()));
    TString title = Form("Fraction of lanes into %s", mLaneStatusFlag[i].c_str());
    //title += ";mm (IB 3x);mm (IB 3x)";
     title += ";mm;mm";
    mvtxmon_LaneStatusOverview[i]->SetTitle(title);
    createPoly(mvtxmon_LaneStatusOverview[i]);
    se->registerHisto(this, mvtxmon_LaneStatusOverview[i]); // mLaneStatusOverview
  }

  //fee
  mTriggerVsFeeId = new TH2I("MVTXMON_FEE_TriggerVsFeeid", "Trigger count vs Trigger ID and Fee ID", NFees, 0, NFees, NTrigger, 0.5, NTrigger + 0.5);
  mTriggerVsFeeId->GetXaxis()->SetTitle("FeeID");
  mTriggerVsFeeId->GetYaxis()->SetTitle("Trigger ID");
  mTriggerVsFeeId->SetStats(0);
  for (int i = 0; i < NTrigger; i++) {
    mTriggerVsFeeId->GetYaxis()->SetBinLabel(i + 1, mTriggerType[i]);
  }
  se->registerHisto(this, mTriggerVsFeeId);

  mTrigger = new TH1I("MVTXMON_FEE_TriggerFlag", "Trigger vs counts", NTrigger, 0.5, NTrigger + 0.5);
  mTrigger->GetXaxis()->SetTitle("Trigger ID");
  mTrigger->GetYaxis()->SetTitle("Counts");
  mTrigger->SetStats(0);
  for (int i = 0; i < NTrigger; i++) {
    mTrigger->GetXaxis()->SetBinLabel(i + 1, mTriggerType[i]);
  }
  se->registerHisto(this, mTrigger);

 /* mLaneInfo = new TH2I("MVTXMON/FEE/LaneInfo", "Lane Information", NLanesMax, -.5, NLanesMax - 0.5, NFlags, -.5, NFlags - 0.5);
  mLaneInfo->GetXaxis()->SetTitle("Lane");
  mLaneInfo->GetYaxis()->SetTitle("Flag");
  mLaneInfo->SetStats(0);
  se->registerHisto(this, mLaneInfo);*/

  for (int i = 0; i < NFlags; i++) {
    mLaneStatus[i] = new TH2I(Form("MVTXMON_LaneStatus_laneStatusFlag%s", mLaneStatusFlag[i].c_str()), Form("Lane Status Flag: %s", mLaneStatusFlag[i].c_str()), NFees, 0, NFees, NLanesMax, 0, NLanesMax);
    mLaneStatus[i]->GetXaxis()->SetTitle("FEE ID");
    mLaneStatus[i]->GetYaxis()->SetTitle("Lane");
    mLaneStatus[i]->SetStats(0);
    for (const int& lay : LayerBoundaryFEE) {
      auto l = new TLine(lay, 0, lay, mLaneStatus[i]->GetNbinsY());
      mLaneStatus[i]->GetListOfFunctions()->Add(l);
    }
    se->registerHisto(this,  mLaneStatus[i]);

    mLaneStatusCumulative[i] = new TH2I(Form("MVTXMON_LaneStatus_laneStatusFlagCumulative%s", mLaneStatusFlag[i].c_str()), Form("Lane Status Flags since SOX: %s", mLaneStatusFlag[i].c_str()), NFees, 0, NFees, NLanesMax, 0, NLanesMax);
    mLaneStatusCumulative[i]->GetXaxis()->SetTitle("FEE ID");
    mLaneStatusCumulative[i]->GetYaxis()->SetTitle("Lane");
    mLaneStatusCumulative[i]->SetStats(0);
    for (const int& lay : LayerBoundaryFEE) {
      auto l = new TLine(lay, 0, lay, mLaneStatusCumulative[i]->GetNbinsY());
      mLaneStatusCumulative[i]->GetListOfFunctions()->Add(l);
    }
    se->registerHisto(this,  mLaneStatusCumulative[i]);
  }


  for (int i = 0; i < NLAYERS; i++) {
    mLaneStatusSummary[i] = new TH1I(Form("MVTXMON_LaneStatusSummary_LaneStatusSummaryL%i", i), Form("Lane Status Summary L%i", i), 3, 0, 3);
    mLaneStatusSummary[i]->GetYaxis()->SetTitle("#Lanes");
    for (int j = 0; j < NFlags; j++) {
      mLaneStatusSummary[i]->GetXaxis()->SetBinLabel(j + 1, mLaneStatusFlag[j].c_str());
    }
    mLaneStatusSummary[i]->GetXaxis()->CenterLabels();
    mLaneStatusSummary[i]->SetStats(0);
    se->registerHisto(this,  mLaneStatusSummary[i]);
  }

  mLaneStatusSummaryIB = new TH1I("MVTXMON_LaneStatusSummary_LaneStatusSummary", "Lane Status Summary", 3, 0, 3);
  mLaneStatusSummaryIB->GetYaxis()->SetTitle("#Lanes");
  for (int j = 0; j < NFlags; j++) {
    mLaneStatusSummaryIB->GetXaxis()->SetBinLabel(j + 1, mLaneStatusFlag[j].c_str());
  }
  mLaneStatusSummaryIB->GetXaxis()->CenterLabels();
  mLaneStatusSummaryIB->SetStats(0);
  se->registerHisto(this,  mLaneStatusSummaryIB);

  //raw task
  hErrorPlots = new TH1D("MVTXMON_General_ErrorPlots", "Decoding Errors", NError, 0.5, NError + 0.5);
  hErrorPlots->GetYaxis()->SetTitle("Counts");
  hErrorPlots->GetXaxis()->SetTitle("Error ID");
  hErrorPlots->SetMinimum(0);  //remove
  hErrorPlots->SetFillColor(kRed); //remove
  hErrorPlots->SetStats(0);

  TPaveText* pt[NError]={nullptr};
  for (int i = 0; i < NError; i++) {
    pt[i] = new TPaveText(0.20, 0.80 - i * 0.05, 0.85, 0.85 - i * 0.05, "NDC");
    pt[i]->SetTextSize(0.04);
    pt[i]->SetTextAlign(12);
    pt[i]->SetFillColor(0);
    pt[i]->SetTextColor(2);
    pt[i]->AddText(ErrorType[i].Data());
    hErrorPlots->GetListOfFunctions()->Add(pt[i]);
  }
  se->registerHisto(this,  hErrorPlots);

  hErrorFile = new TH2D("MVTXMON_General_ErrorFile", "Decoding Errors vs Packet ID", 6*8, 0, 6*8+1, NError, 0.5, NError + 0.5);
  hErrorFile->GetYaxis()->SetTitle("Error ID");
  hErrorFile->GetXaxis()->SetTitle("Packet ID");
  hErrorFile->GetZaxis()->SetTitle("Counts");
  hErrorFile->SetMinimum(0);
  hErrorFile->SetStats(0);

  for (int i = 1; i < 6; i++) {
    auto l = new TLine(i*8+0.5+((i-3)*0.15), 0.5, i*8+0.5+((i-3)*0.15), hErrorFile->GetNbinsY()+0.5);
    hErrorFile->GetListOfFunctions()->Add(l);
  }

  se->registerHisto(this,  hErrorFile);

  for (int aLayer = 0; aLayer < 3; aLayer++) {
    hOccupancyPlot[aLayer] = new TH1D(Form("MVTXMON_Occupancy_Layer%dOccupancy", aLayer), Form("MVTX Layer %d, Occupancy Distribution", aLayer), 301, -3, 0.05);
    hOccupancyPlot[aLayer]->GetYaxis()->SetTitle("Counts");
    hOccupancyPlot[aLayer]->GetXaxis()->SetTitle("log10(Pixel Occupancy)");
    hOccupancyPlot[aLayer]->SetStats(0);
    se->registerHisto(this,  hOccupancyPlot[aLayer]);

    hChipStaveOccupancy[aLayer] = new TH2D(Form("MVTXMON_Occupancy_Layer%d_Layer%dChipStave", aLayer, aLayer), Form("MVTX Layer%d, Occupancy vs Chip and Stave", aLayer), 9, -.5, 9 - .5,     NStaves[aLayer], -.5, NStaves[aLayer] - .5);
    hChipStaveOccupancy[aLayer]->GetYaxis()->SetTitle("Stave Number");
    hChipStaveOccupancy[aLayer]->GetXaxis()->SetTitle("Chip Number");
    hChipStaveOccupancy[aLayer]->GetZaxis()->SetTitle("Number of Hits");
    hChipStaveOccupancy[aLayer]->SetStats(0);
    se->registerHisto(this,  hChipStaveOccupancy[aLayer]);

    TString tmp = Form("MVTX Layer %d", aLayer);
    tmp += ", Number of noisy pixels (>20\% stobes in event) per RCDAQ event";

    hChipStaveNoisy[aLayer] = new TH2D(Form("MVTXMON_Noisy_Layer%d_ChipStave", aLayer), tmp, 9, -.5, 9 - .5,     NStaves[aLayer], -.5, NStaves[aLayer] - .5);
    hChipStaveNoisy[aLayer]->GetYaxis()->SetTitle("Stave Number");
    hChipStaveNoisy[aLayer]->GetXaxis()->SetTitle("Chip Number");
    //hChipStaveNoisy[aLayer]->GetZaxis()->SetTitle("Number of Noisy Pixels");
    hChipStaveNoisy[aLayer]->SetStats(0);
    se->registerHisto(this,  hChipStaveNoisy[aLayer]);
  }

  hChipHitmap = new TH3I(Form("MVTXMON_chipHitmapFLX%d", this->MonitorServerId()), Form("MVTXMON_chipHitmapFLX%d", this->MonitorServerId()), 1024, -.5, 1023.5, 512, -.5, 511.5,8*9*6,-.5,8*9*6-0.5);
  hChipHitmap_evt = new TH3I(Form("MVTXMON_chipHitmapFLX%d_evt", this->MonitorServerId()), Form("MVTXMON_chipHitmapFLX%d_evt", this->MonitorServerId()), 1024, -.5, 1023.5, 512, -.5, 511.5,8*9*6,-.5,8*9*6-0.5);
  hChipHitmap->GetXaxis()->SetTitle("Col");
  hChipHitmap->GetYaxis()->SetTitle("Row");
  hChipHitmap->SetStats(0);
  se->registerHisto(this, hChipHitmap);

  hChipStrobes = new TH1I("hChipStrobes", "Chip Strobes vs Chip*Stave", 8*9*6,-.5,8*9*6-0.5);
  hChipStrobes->GetXaxis()->SetTitle("Chip*Stave");
  hChipStrobes->GetYaxis()->SetTitle("Counts");
  hChipStrobes->SetStats(0);
  se->registerHisto(this, hChipStrobes);

  hChipL1 = new TH1I("hChipL1", "L1 triggers vs Chip*Stave", 8*9*6,-.5,8*9*6-0.5);
  hChipL1->GetXaxis()->SetTitle("Chip*Stave");
  hChipL1->GetYaxis()->SetTitle("Counts");
  hChipL1->SetStats(0);
  se->registerHisto(this, hChipL1);

  //fhr
  mErrorVsFeeid = new TH2I("MVTXMON_General_ErrorVsFeeid", "Error count vs Error id and Fee id", 3 * StaveBoundary[3], 0, 3 * StaveBoundary[3], NErrorExtended, 0.5, NErrorExtended + 0.5);
  mErrorVsFeeid->GetXaxis()->SetTitle("FEE ID");
  mErrorVsFeeid->GetYaxis()->SetTitle("Error ID");
  mErrorVsFeeid->SetStats(0);
  se->registerHisto(this,  mErrorVsFeeid);

  mGeneralOccupancy = new TH2Poly();
  mGeneralOccupancy->SetTitle("General Occupancy;mm;mm");
  mGeneralOccupancy->SetName("MVTXMON_General_Occupancy");
  //mGeneralOccupancy->GetXaxis()->SetTitle("");
  //mGeneralOccupancy->GetYaxis()->SetTitle("");
  mGeneralOccupancy->SetStats(0);

  mGeneralNoisyPixel = new TH2Poly();
  mGeneralNoisyPixel->SetTitle("Noisy Pixel Number;mm (IB 3x);mm (IB 3x)");
  mGeneralNoisyPixel->SetName("MVTXMON_General_Noisy_Pixel");
  mGeneralNoisyPixel->GetXaxis()->SetTitle("");
  mGeneralNoisyPixel->GetYaxis()->SetTitle("");
  mGeneralNoisyPixel->SetStats(0);

  createPoly(mGeneralOccupancy);
  createPoly(mGeneralNoisyPixel);
  se->registerHisto(this,  mGeneralOccupancy);
  se->registerHisto(this,  mGeneralNoisyPixel);

  for (int mLayer = 0; mLayer < 3; mLayer++) {
    mDeadChipPos[mLayer] = new TH2D(Form("MVTXMON_Occupancy_Layer%d_Layer%dDeadChipPos", mLayer, mLayer), Form("DeadChipPos on Layer %d", mLayer), 9, -0.5, 9 - 0.5, NStaves[mLayer], -0.5, NStaves[mLayer] - 0.5);
    mDeadChipPos[mLayer] ->GetXaxis()->SetTitle("Chip Number");
    mDeadChipPos[mLayer] ->GetYaxis()->SetTitle("Stave Number");
    mDeadChipPos[mLayer]->SetStats(0);
    mAliveChipPos[mLayer] = new TH2D(Form("MVTXMON_Occupancy_Layer%d_Layer%dAliveChipPos", mLayer, mLayer), Form("Fraction of RCDAQ events a chip received data on Layer %d", mLayer), 9, -0.5, 9 - 0.5, NStaves[mLayer], -0.5, NStaves[mLayer] - 0.5);
    mAliveChipPos[mLayer]->GetXaxis()->SetTitle("Chip Number");
    mAliveChipPos[mLayer]->GetYaxis()->SetTitle("Stave Number");
    mAliveChipPos[mLayer]->SetStats(0);

    mDeadChipPos[mLayer]->SetStats(0);
    mAliveChipPos[mLayer]->SetStats(0);

    for (int binx = 0; binx < mDeadChipPos[mLayer]->GetNbinsX(); binx++) {
      for (int biny = 0; biny < mDeadChipPos[mLayer]->GetNbinsY(); biny++) {
        mDeadChipPos[mLayer]->SetBinContent(binx+1,biny+1,1);
      }
    }


    mOccupancyPlot[mLayer] = new TH1D(Form("MVTXMON_Occupancy_Layer%dOccupancy_LOG", mLayer), Form("MVTX Layer %d Noise pixels occupancy distribution", mLayer), 301, -3, 0.05);
    mOccupancyPlot[mLayer]->GetXaxis()->SetTitle("log(Occupancy)");
    mOccupancyPlot[mLayer]->GetYaxis()->SetTitle("");
    mOccupancyPlot[mLayer]->SetStats(0);
    se->registerHisto(this,  mDeadChipPos[mLayer]);
    se->registerHisto(this,  mAliveChipPos[mLayer]);
    se->registerHisto(this,  mOccupancyPlot[mLayer]);
  }

  mTotalDeadChipPos = new TH2D(Form("MVTXMON_Occupancy_TotalDeadChipPos"), Form("TotalDeadChipPos "), 9, -0.5, 9 - 0.5, NStaves[2], -0.5, NStaves[2] - 0.5);
  mTotalDeadChipPos->GetXaxis()->SetTitle("Chip Number");
  mTotalDeadChipPos->GetYaxis()->SetTitle("Stave Number");
  se->registerHisto(this,  mTotalDeadChipPos);
  
  mTotalAliveChipPos = new TH2D(Form("MVTXMON_Occupancy_TotalAliveChipPos"), Form("TotalAliveChipPos "), 9, -0.5, 9 - 0.5, NStaves[2], -0.5, NStaves[2] - 0.5);
  mTotalAliveChipPos->GetXaxis()->SetTitle("Chip Number");
  mTotalAliveChipPos->GetYaxis()->SetTitle("Stave Number");
  se->registerHisto(this,  mTotalAliveChipPos);  


  Reset();
  return 0;
}

int MvtxMon::BeginRun(const int /* runno */)
{
  // if you need to read calibrations on a run by run basis
  // this is the place to do it
  return 0;
}

int MvtxMon::process_event(Event *evt)
{
  evtcnt++;
  mRCDAQevt->Fill(this->MonitorServerId());
  //std::cout << "Processing Event " << evtcnt << std::endl;
  OnlMonServer *se = OnlMonServer::instance();

  plist = new Packet *[2];
   //SingleMvtxInput *reader = new SingleMvtxInput("onlmonreader");

  for(int l = 0; l < NLAYERS; l++){
    for(int i = 0; i < NStaves[l]; i++){
      for(int j = 0; j < NCHIP; j++){
        mHitPerChip[l][i][j] = 0;
      }
    }
  }
  hChipHitmap_evt->Reset("ICESM");

   int nChipStrobes[8*9*6] = {0};


   int npackets = evt->getPacketList(plist, 2);

    if (npackets > 2)
    {
      exit(1);
    }

    for (int i = 0; i < npackets; i++)
    {
      // Ignoring packet not from MVTX detector
      if ( (plist[i]->getIdentifier() < 2001) || (plist[i]->getIdentifier() > 2052) )
      {
        continue;
      }
      if (Verbosity() > 1)
      {
        plist[i]->identify();
      }
      int num_feeId = plist[i]->iValue(-1, "NR_LINKS");
      if (Verbosity()  > 1)
      {
        std::cout << "Number of feeid in RCDAQ events: " << num_feeId << " for packet "
          << plist[i]->getIdentifier() << std::endl;
      }
      if (num_feeId > 0)
      {
        for (int i_fee{0}; i_fee < num_feeId; ++i_fee)
        {
          auto feeId = plist[i]->iValue(i_fee, "FEEID");
          auto link = DecodeFeeid(feeId);
          auto num_strobes = plist[i]->iValue(feeId, "NR_STROBES");
          ntriggers = num_strobes;
          auto num_L1Trgs = plist[i]->iValue(feeId, "NR_PHYS_TRG");
          for ( int iL1 = 0; iL1 < num_L1Trgs; ++iL1 )
          {
            //auto l1Trg_bco = plist[i]->lValue(feeId, iL1, "L1_IR_BCO");
            hChipL1->Fill((StaveBoundary[link.layer]+link.stave)*9 + 3 * link.gbtid + 0); //same for chip id 0 1 and 2
	    hChipL1->Fill((StaveBoundary[link.layer]+link.stave)*9 + 3 * link.gbtid + 1);
	    hChipL1->Fill((StaveBoundary[link.layer]+link.stave)*9 + 3 * link.gbtid + 2);
          }

          //m_FeeStrobeMap[feeId] += num_strobes;
          for (int i_strb{0}; i_strb < num_strobes; ++i_strb)
          {
            auto strb_bco = plist[i]->lValue(feeId, i_strb, "TRG_IR_BCO");
            //auto strb_bc  = plist[i]->iValue(feeId, i_strb, "TRG_IR_BC");
            auto num_hits = plist[i]->iValue(feeId, i_strb, "TRG_NR_HITS");
            if (Verbosity() > 4)
            {
 	      if(link.layer == 0){
              std::cout << "evtno: " << ", Fee: " << feeId;
              std::cout << " Layer: " << link.layer << " Stave: " << link.stave;
              std::cout << " GBT: " << link.gbtid << ", bco: 0x" << std::hex << strb_bco << std::dec;
              std::cout << ", n_hits: " << num_hits << std::endl;
	      }
            }
            hChipStrobes->Fill((StaveBoundary[link.layer]+link.stave%20)*9 + 3 * link.gbtid + 0); //same for chip id 0 1 and 2
	    hChipStrobes->Fill((StaveBoundary[link.layer]+link.stave%20)*9 + 3 * link.gbtid + 1);
	    hChipStrobes->Fill((StaveBoundary[link.layer]+link.stave%20)*9 + 3 * link.gbtid + 2);
	    nChipStrobes[(StaveBoundary[link.layer]+link.stave%20)*9 + 3 * link.gbtid + 0]++;
	    nChipStrobes[(StaveBoundary[link.layer]+link.stave%20)*9 + 3 * link.gbtid + 1]++;
            nChipStrobes[(StaveBoundary[link.layer]+link.stave%20)*9 + 3 * link.gbtid + 2]++;

           
            

            for (int i_hit{0}; i_hit < num_hits; ++i_hit)
            {
              //auto chip_bc = plist[i]->iValue(feeId, i_strb, i_hit, "HIT_BC");
              auto chip_id = plist[i]->iValue(feeId, i_strb, i_hit, "HIT_CHIP_ID");
              auto chip_row = plist[i]->iValue(feeId, i_strb, i_hit, "HIT_ROW");
              auto chip_col = plist[i]->iValue(feeId, i_strb, i_hit, "HIT_COL");

              mHitPerChip[link.layer][link.stave%20][3 * link.gbtid + chip_id]++;
              hChipHitmap->Fill(chip_col,chip_row,(StaveBoundary[link.layer]+link.stave%20)*9 + 3 * link.gbtid + chip_id);
              hChipStaveOccupancy[link.layer]->Fill(3 * link.gbtid + chip_id, link.stave%20);
	      hChipHitmap_evt->Fill(chip_col,chip_row,(StaveBoundary[link.layer]+link.stave%20)*9 + 3 * link.gbtid + chip_id);

            }

	    
		

            //m_BeamClockFEE[strb_bco].insert(feeId);
            //m_BclkStack.insert(strb_bco);
           // m_FEEBclkMap[feeId] = strb_bco;
          }
        }
      }
      delete plist[i];
    }

	int firedChips = 0;
        int firedPixels = 0;
        int sumstrobes = 0;
        int nstrobes = 0;
            for(int l = 0; l < NLAYERS; l++){
    	      for(int s = 0; s < NStaves[l]; s++){
      		for(int j = 0; j < NCHIP; j++){
        	  if(mHitPerChip[l][s][j] > 0){
                     firedChips++;
	             mvtxmon_EvtHitChip->Fill(mHitPerChip[l][s][j]/nChipStrobes[(StaveBoundary[l]+s)*9+j]);
                     firedPixels+=mHitPerChip[l][s][j];
		     sumstrobes+=nChipStrobes[(StaveBoundary[l]+s)*9+j];
                     nstrobes++;
		  }
      		}
              }
            }
	    mvtxmon_ChipFiredHis->Fill(firedChips);
           mvtxmon_EvtHitDis->Fill((double)firedPixels/((double)sumstrobes/(double)nstrobes));

  double pixelOccupancy, chipOccupancy;
  for (int iLayer = 0; iLayer < 3; iLayer++) {
    for (int iStave = 0; iStave < NStaves[iLayer]; iStave++) {
      for (int iChip = 0; iChip < 9; iChip++) {
        chipOccupancy = hChipHitmap->Integral(0,-1,0,-1,(StaveBoundary[iLayer]+iStave)*9+iChip+1,(StaveBoundary[iLayer]+iStave)*9+iChip+1); //scale at client
        double chipOccupancyNorm = chipOccupancy/hChipStrobes->GetBinContent((StaveBoundary[iLayer]+iStave)*9+iChip+1)/1024/512;
         //if (chipOccupancyNorm > 0) mvtxmon_ChipStave1D->SetBinContent((iLayer==0?iStave:NStaves[iLayer]+iStave)*9+iChip+1,chipOccupancyNorm); //need to remember total number of occ and events and scale here
         if (chipOccupancyNorm > 0)mvtxmon_ChipStave1D->SetBinContent((StaveBoundary[iLayer]+iStave)*9 + iChip +1,chipOccupancyNorm); 
         if (chipOccupancyNorm > 0)mGeneralOccupancy->SetBinContent(mapstave[iLayer][iStave], chipOccupancyNorm);
        int nTrg = nChipStrobes[(StaveBoundary[iLayer]+iStave)*9+iChip];
        for (int iCol = 0; iCol < NCols; iCol++) {
          for (int iRow = 0; iRow < NRows; iRow++) {
            pixelOccupancy = hChipHitmap_evt->GetBinContent(iCol + 1, iRow + 1,(StaveBoundary[iLayer]+iStave)*9+iChip+1);
            if (pixelOccupancy > 0) {
              if (pixelOccupancy/(double)nTrg > mOccupancyCutForNoisyPixel) {
                  mNoisyPixelNumber[iLayer][iStave][iChip]++;
                  mOccupancyPlot[iLayer]->Fill(log10(pixelOccupancy/(double)nTrg));
              }                         
              hOccupancyPlot[iLayer]->Fill(log10(pixelOccupancy/(double)nTrg));
            }
          }
        }
      }      
    }
  }

  for (int iLayer = 0; iLayer < 3; iLayer++) {
    for (int iStave = 0; iStave < NStaves[iLayer]; iStave++) {
      for (int iChip = 0; iChip < 9; iChip++) {
        double noisy = *(std::max_element(mNoisyPixelNumber[iLayer][iStave], mNoisyPixelNumber[iLayer][iStave] + 9));
        mGeneralNoisyPixel->SetBinContent(mapstave[iLayer][iStave], noisy);
        if(mNoisyPixelNumber[iLayer][iStave][iChip]>0) hChipStaveNoisy[iLayer]->SetBinContent(iChip+1,iStave+1,mNoisyPixelNumber[iLayer][iStave][iChip]);
      }
    }
  }

for (int iLayer = 0; iLayer < 3; iLayer++) {
    for (int iStave = 0; iStave < NStaves[iLayer]; iStave++) {
      for (int iChip = 0; iChip < 9; iChip++) {
        if (!mHitPerChip[iLayer][iStave][iChip]) {
          if(mDeadChipPos[iLayer]->GetBinContent(mDeadChipPos[iLayer]->GetXaxis()->FindBin(iChip), mDeadChipPos[iLayer]->GetYaxis()->FindBin(iStave)) > 0.5) mDeadChipPos[iLayer]->SetBinContent(mDeadChipPos[iLayer]->GetXaxis()->FindBin(iChip), mDeadChipPos[iLayer]->GetYaxis()->FindBin(iStave), 1);
          if(mTotalDeadChipPos->GetBinContent(mTotalDeadChipPos->GetXaxis()->FindBin(iChip), mTotalDeadChipPos->GetYaxis()->FindBin(iStave)) > 0.5) mTotalDeadChipPos->SetBinContent(mTotalDeadChipPos->GetXaxis()->FindBin(iChip), mTotalDeadChipPos->GetYaxis()->FindBin(iStave), 1);
        } else {
	  mAliveChipPos[iLayer]->Fill(iChip, iStave);
          mTotalAliveChipPos->Fill(iChip, iStave);
          mDeadChipPos[iLayer]->SetBinContent(mDeadChipPos[iLayer]->GetXaxis()->FindBin(iChip), mDeadChipPos[iLayer]->GetYaxis()->FindBin(iStave), 0);                // not dead
          mTotalDeadChipPos->SetBinContent(mTotalDeadChipPos->GetXaxis()->FindBin(iChip), mTotalDeadChipPos->GetYaxis()->FindBin(iStave), 0); // not dead
        }
      }
    }
  }

  
  int packet_init = 2001;
  //std::cout<<"processing rcdaq event"<<std::endl;
  for ( int iPkt = 0; iPkt < 8; iPkt++ ){
    Packet *p = evt->getPacket(packet_init + iPkt);  
    if (p){
      //std::cout<<"PACKET: "<<iPkt<<std::endl;
     // std::map<mvtx::InteractionRecord, std::vector<mvtx::ChipPixelData>> *data = reinterpret_cast<std::map<mvtx::InteractionRecord, std::vector<mvtx::ChipPixelData>>*>(p->pValue(-1, "ChipData"));
     // int nevents_packet = 0; 
    if(false){
    //std::vector<mvtx::GBTLinkDecodingStat> *linkErrors = reinterpret_cast<std::vector<mvtx::GBTLinkDecodingStat>*>(p->pValue(-1, "linkErrors"));
/*
      if(linkErrors){
        for (auto error : *linkErrors){
          for (int ierror = 0; ierror < mvtx::GBTLinkDecodingStat::NErrorsDefined; ierror++) {
            if (error.errorCounts[ierror] <= 0) {
              continue;
            }
            //mErrorCount[istave][ilink][ierror] = GBTLinkInfo->statistics.errorCounts[ierror];
            std::cout<<error.feeID<<" "<<ierror<<" "<<error.errorCounts[ierror]<<std::endl;
          }

      //    for (int i = 0; i < NError; i++) {
       //     hErrorPlots->SetBinContent(i + 1, mErrors[i]);
          //hErrorFile->SetBinContent((FileID + 1 + (EPID - 4) * 12), i + 1, mErrorPerFile[i]);
       //   }
        }
      }
*/	
  delete p;
}
  }// if(p)
}// for packet loop


  // get temporary pointers to histograms
  // one can do in principle directly se->getHisto("mvtxhist1")->Fill()
  // but the search in the histogram Map is somewhat expensive and slows
  // things down if you make more than one operation on a histogram


    std::ostringstream msg;
    msg << "Filling Histos";
    se->send_message(this, MSG_SOURCE_UNSPECIFIED, MSG_SEV_INFORMATIONAL, msg.str(), FILLMESSAGE);
    idummy = 0;
  return 0;
}

int MvtxMon::Reset()
{
  // reset our internal counters
  evtcnt = 0;
  idummy = 0;
  return 0;
}


void MvtxMon::getStavePoint(int layer, int stave, double* px, double* py)
{
  float stepAngle = M_PI * 2 / NStaves[layer];              // the angle between to stave
  float midAngle = StartAngle[layer] + (stave * stepAngle); // mid point angle
  float staveRotateAngle = M_PI / 2 - (stave * stepAngle);  // how many angle this stave rotate(compare with first stave)
  px[1] = MidPointRad[layer] * std::cos(midAngle);          // there are 4 point to decide this TH2Poly bin
                                                            // 0:left point in this stave;
                                                            // 1:mid point in this stave;
                                                            // 2:right point in this stave;
                                                            // 3:higher point int this stave;
  py[1] = MidPointRad[layer] * std::sin(midAngle);          // 4 point calculated accord the blueprint
                                                            // roughly calculate
  px[0] = 7.7 * std::cos(staveRotateAngle) + px[1];
  py[0] = -7.7 * std::sin(staveRotateAngle) + py[1];
  px[2] = -7.7 * std::cos(staveRotateAngle) + px[1];
  py[2] = 7.7 * std::sin(staveRotateAngle) + py[1];
  px[3] = 5.623 * std::sin(staveRotateAngle) + px[1];
  py[3] = 5.623 * std::cos(staveRotateAngle) + py[1];
}

void MvtxMon::createPoly(TH2Poly *h){
  for (int ilayer = 0; ilayer < NLAYERS; ilayer++) {
    for (int istave = 0; istave < NStaves[ilayer]; istave++) {
      double* px = new double[4];
      double* py = new double[4];
      getStavePoint(ilayer, istave, px, py);
      for (int icoo = 0; icoo < 4; icoo++) {
        px[icoo] *= 3.;
        py[icoo] *= 3.;
      }
      h->AddBin(4, px, py);
    }
  }
}


