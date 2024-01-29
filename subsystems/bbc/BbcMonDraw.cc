#include "BbcMonDraw.h"

#include <onlmon/OnlMonClient.h>

#include <TArc.h>
#include <TArrow.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TDatime.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TLatex.h>
#include <TLine.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TSpectrum.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TText.h>

#include <algorithm>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include <chrono>

#define DEBUG
#ifdef DEBUG
#define PRINT_DEBUG(x) std::cout<<x<<std::endl
#else
#define PRINT_DEBUG(x) {};
#endif

#ifdef DEBUGNEW
#define ifdelete(x) if(x!=nullptr){ std::cout << "Delete " << #x << std::endl; delete x;x=nullptr;}
#define ifnew(t,x) {if(x!=nullptr){ std::cout << "Delete " << #x << std::endl; delete x;}std::cout << "New "<< #x << std::endl;x = new t;}
#else
#define ifdelete(x) if(x!=nullptr){ delete x;x=nullptr;}
#define ifnew(t,x) {if(x!=nullptr){ delete x;}x = new t;}
#endif

// x position of trigger, scale factor and vtx mean
// static float xpos[3] = {0.30, 0.50, 0.75};
static float xpos[4] = {0.20, 0.35, 0.55, 0.75};

BbcMonDraw::BbcMonDraw(const std::string &name)
  : OnlMonDraw(name)
{
  PRINT_DEBUG("In BbcMonDraw::BbcMonDraw()");
  std::cout << "Done with BbcMonDraw::BbcMonDraw" << std::endl;
  return;
}

BbcMonDraw::~BbcMonDraw()
{
  PRINT_DEBUG("In BbcMonDraw::~BbcMonDraw()");
  // ifdelete( bbccalib );

  ifdelete( bbcStyle );

  // ------------------------------------------------------
  // Canvas and Histogram

  ifdelete(PaveTop);
  ifdelete(TextTop);

  for (int i = 0; i < MAX_WARNING; i++)
  {
    ifdelete(PaveWarning[i]);
    ifdelete(ArcWarning[i]);
  }

  // ------------------------------------------------------
  // Graph for visualization

  // for the 1st Page

  for (int side = 0; side < nSIDE; side++)
  {
    ifdelete(TdcOver[side]);
    ifdelete(LineTdcOver[side][0]);
    ifdelete(LineTdcOver[side][1]);
    ifdelete(BoxTdcOver[side]);

    for (int trig = 0; trig < nTRIGGER; trig++)
    {
      ifdelete(nHit[trig][side]);
      ifdelete(LinenHit[trig][side][0]);
      ifdelete(LinenHit[trig][side][1]);
      ifdelete(BoxnHit[trig][side]);
    }
  }
  ifdelete(TextnHitStatus);

  // for 2nd Page

  for (int side = 0; side < nSIDE; side++)
  {
    ifdelete(HitTime[side]);
    ifdelete(LineHitTime[side][0]);
    ifdelete(LineHitTime[side][1]);
    ifdelete(ArrowHitTime[side]);
    ifdelete(TextHitTime[side]);
    ifdelete(FitHitTime[side]);
  }

  ifdelete(AvrHitTime);
  ifdelete(LineAvrHitTime[0]);
  ifdelete(LineAvrHitTime[1]);
  ifdelete(ArrowAvrHitTime);
  ifdelete(TextAvrHitTime);
  ifdelete(FitAvrHitTime);

  ifdelete(NorthHitTime);
  ifdelete(LineNorthHitTime[0]);
  ifdelete(LineNorthHitTime[1]);
  ifdelete(ArrowNorthHitTime);
  ifdelete(TextNorthHitTime);
  ifdelete(FitNorthHitTime);
  ifdelete(SouthHitTime);
  ifdelete(LineSouthHitTime[0]);
  ifdelete(LineSouthHitTime[1]);
  ifdelete(ArrowSouthHitTime);
  ifdelete(TextSouthHitTime);
  ifdelete(FitSouthHitTime);

  ifdelete(NorthChargeSum);
  ifdelete(SouthChargeSum);
  ifdelete(TextSouthChargeSum);
  ifdelete(TextNorthChargeSum);

  ifdelete(TextBbcSummaryHitTime[0]);
  ifdelete(TextBbcSummaryHitTime[1]);
  ifdelete(TextBbcSummaryGlobalOffset[0]);
  ifdelete(TextBbcSummaryGlobalOffset[1]);
  ifdelete(TextBbcSummaryZvertex);
  ifdelete(TextBbcSummaryTrigRate);

  ifdelete(Zvtx);
  ifdelete(FitZvtx);
  ifdelete(LineZvtx[0]);
  ifdelete(LineZvtx[1]);
  ifdelete(ArrowZvtx);
  ifdelete(TextZvtx);
  ifdelete(TextZvtxNorth);
  ifdelete(TextZvtxSouth);
  ifdelete(TextZvtxStatus[0]);
  ifdelete(TextZvtxStatus[1]);

  ifdelete(ArmHit);
  ifdelete(TextArmHit);
  ifdelete(ArcArmHit);

  // for 3rd Page
  ifdelete(Adc);
  ifdelete(MultiView1F);
  ifdelete(MultiView2F);
  ifdelete(PaveWarnings);

  // for 4th Page
  ifdelete(Zvtx_bbll1);
  ifdelete(Zvtx);
  ifdelete(TzeroZvtx);
  ifdelete(TextZVertexNotice);

  ifdelete(LineTzeroZvtx[0]);
  ifdelete(LineTzeroZvtx[1]);
  ifdelete(LineTzeroZvtx[2]);
  ifdelete(LineTzeroZvtx[3]);
  ifdelete(TextTzeroZvtx);

  // for (int i = 0; i < 4; i++)//need a look
  for (int i = 0; i < 5; i++)
  {
    ifdelete(TextZVertex[i]);
    ifdelete(TextZVertex_scale[i]);
    ifdelete(TextZVertex_mean[i]);
  }
  ifdelete( FitZvtx );
  ifdelete( TextZvtxStatus[0] );
  ifdelete( TextZvtxStatus[1] );

  for (int icv=0; icv<nCANVAS; icv++)
  {
    ifdelete( TC[icv] );
  }

  ifdelete(Prescale_hist);
  ifdelete(tspec);

  return;
}

int BbcMonDraw::Init()
{
  PRINT_DEBUG("In BbcMonDraw::Init()");
  TStyle *oldStyle = gStyle;
  bbcStyle = new TStyle("bbcStyle", "BBC/MBD Online Monitor Style");
  bbcStyle->SetOptStat(0);
  bbcStyle->SetTitleH(0.075);
  bbcStyle->SetTitleW(0.98);
  bbcStyle->SetPalette(1);
  bbcStyle->SetFrameBorderMode(0);
  bbcStyle->SetPadBorderMode(0);
  bbcStyle->SetCanvasBorderMode(0);
  oldStyle->cd();

  // ------------------------------------------------------
  // Canvas and Histogram

  memset(TC, 0, sizeof(TC));
  memset(transparent, 0, sizeof(transparent));
  memset(PadTop, 0, sizeof(PadTop));
  memset(PaveWarning, 0, sizeof(PaveWarning));
  memset(ArcWarning, 0, sizeof(ArcWarning));
  memset(PadWarning, 0, sizeof(PadWarning));
  memset(nPadWarning, 0, sizeof(nPadWarning));
  // ------------------------------------------------------
  // Graph for visualization

  // for the 1st Page
  memset(TdcOver, 0, sizeof(TdcOver));
  memset(PadTdcOver, 0, sizeof(PadTdcOver));
  memset(LineTdcOver, 0, sizeof(LineTdcOver));
  memset(BoxTdcOver, 0, sizeof(BoxTdcOver));
  memset(FrameTdcOver, 0, sizeof(FrameTdcOver));
  memset(PadnHit, 0, sizeof(PadnHit));
  memset(nHit, 0, sizeof(nHit));
  memset(LinenHit, 0, sizeof(LinenHit));
  memset(BoxnHit, 0, sizeof(BoxnHit));
  memset(FramenHit, 0, sizeof(FramenHit));

  // for 2nd Page

  memset(HitTime, 0, sizeof(HitTime));
  memset(PadHitTime, 0, sizeof(PadHitTime));
  memset(LineHitTime, 0, sizeof(LineHitTime));
  memset(ArrowHitTime, 0, sizeof(ArrowHitTime));
  memset(TextHitTime, 0, sizeof(TextHitTime));
  memset(FitHitTime, 0, sizeof(FitHitTime));
  memset(LineAvrHitTime, 0, sizeof(LineAvrHitTime));

  memset(LineNorthHitTime, 0, sizeof(LineNorthHitTime));
  memset(LineSouthHitTime, 0, sizeof(LineSouthHitTime));

  memset(TextBbcSummaryHitTime, 0, sizeof(TextBbcSummaryHitTime));
  memset(TextBbcSummaryGlobalOffset, 0, sizeof(TextBbcSummaryGlobalOffset));

  memset(LineZvtx, 0, sizeof(LineZvtx));

  memset(LineTzeroZvtx, 0, sizeof(LineTzeroZvtx));

  memset(TextZVertex, 0, sizeof(TextZVertex));
  memset(TextZVertex_scale, 0, sizeof(TextZVertex_scale));
  memset(TextZVertex_mean, 0, sizeof(TextZVertex_mean));

  FitZvtx = new TF1("FitZvtx", "gaus");
  TextZvtxStatus[0] = new TLatex;
  TextZvtxStatus[1] = new TLatex;

  tspec = new TSpectrum(5);  // 5 peaks is enough - we have 4

  return 0;
}

int BbcMonDraw::ClearWarning()
{
  for (int i = 0; i < nWarning; i++)
  {
    ifdelete(PaveWarning[i]);
    ifdelete(ArcWarning[i]);
  }
  nWarning = 0;

  if (PaveWarnings)
  {
    PaveWarnings->Clear();
  }
  // if( PadWarnings )
  // PaveWarnings->Draw();
  return 0;
}

int BbcMonDraw::Warning(TPad *pad, const float x, const float y, const int r, const std::string &msg)
{
  static int brink = 1;
  float x1, x2, y1, y2;

  x1 = pad->GetX1() * 0.97 + pad->GetX2() * 0.03;
  x2 = pad->GetX1() * 0.03 + pad->GetX2() * 0.97;
  // pave rises on top if error point is lower
  if (y > pad->GetY1() * 0.5 + pad->GetY2() * 0.5)
  {
    y1 = pad->GetY1() + (pad->GetY2() - pad->GetY1()) * 0.10;
    y2 = pad->GetY1() + (pad->GetY2() - pad->GetY1()) * 0.30;
  }
  else
  {
    y1 = pad->GetY1() + (pad->GetY2() - pad->GetY1()) * 0.70;
    y2 = pad->GetY1() + (pad->GetY2() - pad->GetY1()) * 0.90;
  }

  int newWarning = nWarning;
  for (int i = 0; i < nWarning; i++)
  {
    if (PadWarning[i] == pad)
    {
      newWarning = i;
    }
  }

  ifnew(TArc(x, y, r), ArcWarning[newWarning]);
  ArcWarning[newWarning]->SetLineWidth(r);
  ArcWarning[newWarning]->SetLineColor(2);

  ifnew(TPaveText(x1, y1, x2, y2), PaveWarning[newWarning]);
  PaveWarning[newWarning]->AddText(msg.c_str());
  ArcWarning[newWarning]->SetLineColor(2 + 3 * brink);
  PaveWarning[newWarning]->SetFillColor(5 + 5 * brink);
  PaveWarning[newWarning]->SetTextColor(2);
  PaveWarning[newWarning]->SetLineColor(2);
  PaveWarning[newWarning]->SetLineWidth(2);

  if (newWarning == nWarning)
  {
    PadWarning[newWarning] = pad;
    nPadWarning[newWarning] = 1;
    nWarning++;
  }
  else
  {
    std::ostringstream wmsg;
    wmsg << "... and other " << nPadWarning[newWarning] << " warnings";
    nPadWarning[newWarning]++;
    PaveWarning[newWarning]->AddText(wmsg.str().c_str());
    wmsg.str("");
  }
  // PaveWarning[newWarning]->AddText( "Call BBC/MBD Expert. If this run is a PHYSICS run " );
  PaveWarning[newWarning]->AddText("Do Quick-feed BBC/MBD");
  PaveWarning[newWarning]->AddText("If it is NOT fixed by feed at PHYSICS run, call BBC/MBD experts");
  // PaveWarning[newWarning]->AddText( "If it is NOT fixed by quick feed at PHYSICS run, conform to the BBC/MBD standing orders" );

  if (PaveWarnings)
  {
    PaveWarnings->AddText(msg.c_str());
  }

  if (PadWarnings && PaveWarnings)
  {
    PaveWarnings->Draw();
  }
  // ArcWarning[newWarning]->Draw();
  PaveWarning[newWarning]->Draw();

  // brink = 1 - brink;

  if (nWarning == MAX_WARNING - 2)
  {
    std::string bmsg = "Too Many Warnings";
    Warning(pad, x, y, 0, bmsg);
    bmsg.erase();
  }

  return 0;
}

int BbcMonDraw::MakeCanvas(const std::string &name)
{
  PRINT_DEBUG("In BbcMonDraw::MakeCanvas()");
  OnlMonClient *cl = OnlMonClient::instance();
  int xsize = cl->GetDisplaySizeX();
  int ysize = cl->GetDisplaySizeY();

  if (name == "BbcMon1")
  {
    std::cout << "Creating Canvas BbcMon1..." << std::endl;

    TC[0] = new TCanvas("BbcMon1", "BBC/MBD Z-Vertex View for Shift crew", -1, 0, xsize / 2, ysize);

    // root is pathetic, whenever a new TCanvas is created root piles up
    // 6kb worth of X11 events which need to be cleared with
    // gSystem->ProcessEvents(), otherwise your process will grow and
    // grow and grow but will not show a definitely lost memory leak
    gSystem->ProcessEvents();

    TC[0]->cd();
    PadTop[0] = new TPad("PadTop0", "PadTop0", 0.00, 0.90, 1.00, 1.00, 0, 0, 0);
    PadZVertex = new TPad("PadZVertex", "PadZVertex", 0.00, 0.60, 1.00, 0.90, 0, 0, 0);
    PadZVertexSummary = new TPad("PadZVertexSummary", "PadZVertexSummary", 0.00, 0.40, 1.00, 0.60, 0, 0, 0);
    PadSouthHitMap = new TPad("PadSouthHitMap", "PadSouthHitMap", 0.00, 0.00, 0.495, 0.40, 0, 0, 0);
    PadNorthHitMap = new TPad("PadNorthHitMap", "PadNorthHitMap", 0.505, 0.00, 1.0, 0.40, 0, 0, 0);
    //PadTzeroZVertex = new TPad("PadTzeroZVertex", "PadTzeroZVertex", 0.00, 0.00, 1.00, 0.40, 0, 0, 0);

    PadTop[0]->Draw();
    // PadZVertex->SetLogy();
    PadZVertex->Draw();
    //PadTzeroZVertex->Draw();
    PadSouthHitMap->Draw();
    PadNorthHitMap->Draw();
    PadZVertexSummary->Draw();

    PadZVertexSummary->cd();

    // ifnew( TText, TextZVertexNotice );
    // TextZVertexNotice->SetTextSize(0.08);
    // TextZVertexNotice->SetText(0.05, 0.75, "< Z vertex deviation may NOT be due to BBC/MBD, don't page expert easily! >");

    // for (int i = 0; i < 4; i++)
    for (int i = 0; i < 1; i++)  //  pp
    {
      ifnew(TText, TextZVertex[i]);
      TextZVertex[i]->SetTextColor(bbc_onlmon::BBC_COLOR_ZVTX[i]);
      // TextZVertex[i]->SetTextSize(0.15);
      TextZVertex[i]->SetTextSize(0.08);

      ifnew(TText, TextZVertex_scale[i]);
      ifnew(TText, TextZVertex_mean[i]);
      TextZVertex_scale[i]->SetTextColor(bbc_onlmon::BBC_COLOR_ZVTX[i]);
      TextZVertex_mean[i]->SetTextColor(bbc_onlmon::BBC_COLOR_ZVTX[i]);
      TextZVertex_scale[i]->SetTextSize(0.08);
      TextZVertex_mean[i]->SetTextSize(0.08);
    }
    /*chiu
    TextZVertex[0]->SetText(xpos[0], 0.65, "Zbbc [BBLL1]");       // RUN11 pp
    */
    TC[0]->cd();
    transparent[0] = new TPad("transparent3", "this does not show", 0, 0, 1, 1, 0, 0);
    transparent[0]->SetFillStyle(4000);
    transparent[0]->Draw();

    /*
    ifnew(TLine(bbc_onlmon::BBC_MIN_REGULAR_ZVERTEX_MEAN, -6,
                bbc_onlmon::BBC_MIN_REGULAR_ZVERTEX_MEAN, 16),
          LineTzeroZvtx[0]);
    ifnew(TLine(bbc_onlmon::BBC_MAX_REGULAR_ZVERTEX_MEAN, -6,
                bbc_onlmon::BBC_MAX_REGULAR_ZVERTEX_MEAN, 16),
          LineTzeroZvtx[1]);
    ifnew(TLine(-200, bbc_onlmon::BBC_MIN_REGULAR_TDC1_MEAN - 5,
                200, bbc_onlmon::BBC_MIN_REGULAR_TDC1_MEAN - 5),
          LineTzeroZvtx[2]);
    ifnew(TLine(-200, bbc_onlmon::BBC_MAX_REGULAR_TDC1_MEAN - 5,
                200, bbc_onlmon::BBC_MAX_REGULAR_TDC1_MEAN - 5),
          LineTzeroZvtx[3]);
    ifnew(TText, TextTzeroZvtx);
    */
  }

  else if (name == "BbcMon2")
  {
    std::cout << "Creating Canvas BbcMon2..." << std::endl;

    TC[1] = new TCanvas("BbcMon2", "Beam status view for Shift crew", -xsize / 2, 0, xsize / 2, ysize);
    gSystem->ProcessEvents();
    TC[1]->cd();
    PadTop[1] = new TPad("PadTop1", "PadTop1", 0.00, 0.90, 1.00, 1.00, 0, 0, 0);
    PadTimeWave = new TPad("PadTimeWave", "PadTimeWave", 0.00, 0.30, 1.00, 0.90, 0, 0, 0);
    PadSouthHitTime = new TPad("PadHitTimeSouth1", "PadHitTimeSouth1", 0.00, 0.00, 0.50, 0.30, 0, 0, 0);
    PadNorthHitTime = new TPad("PadHitTimeNorth1", "PadHitTimeNorth1", 0.50, 0.00, 1.00, 0.30, 0, 0, 0);

    //PadArmHit = new TPad("PadArmHit", "PadArmHit", 0.50, 0.00, 1.00, 0.40, 0, 0, 0);
    //PadBbcSummary = new TPad("PadBbcSummary", "PadBbcSummary", 0.00, 0.40, 1.00, 0.55, 0, 0, 0);
    //PadAvrHitTime = new TPad("PadAvrHitTime", "PadAvrHitTime", 0.00, 0.00, 0.50, 0.40, 0, 0, 0);

    //PadArmHit->SetLogz();

    PadTop[1]->Draw();
    if (PadTimeWave)   PadTimeWave->Draw();
    if (PadNorthHitTime) PadNorthHitTime->Draw();
    if (PadSouthHitTime) PadSouthHitTime->Draw();
    if (PadArmHit) PadArmHit->Draw();
    if (PadBbcSummary) PadBbcSummary->Draw();
    if (PadAvrHitTime) PadAvrHitTime->Draw();

    ifnew(TText, TextBbcSummaryHitTime[0]);
    ifnew(TText, TextBbcSummaryHitTime[1]);
    ifnew(TText, TextBbcSummaryZvertex);
    ifnew(TText, TextBbcSummaryGlobalOffset[0]);
    ifnew(TText, TextBbcSummaryGlobalOffset[1]);
    // TextBbcSummaryHitTime->SetTextSize(0.2);
    TextBbcSummaryHitTime[0]->SetTextSize(0.2);
    TextBbcSummaryHitTime[1]->SetTextSize(0.2);
    TextBbcSummaryZvertex->SetTextSize(0.2);
    TextBbcSummaryGlobalOffset[0]->SetTextSize(0.2);
    TextBbcSummaryGlobalOffset[1]->SetTextSize(0.2);

    ifnew(TArc(10.5, 10.5, 4, 0, 360), ArcArmHit);
    ArcArmHit->SetFillStyle(4000);
    // ifnew( TArc(10.5, 10.5, 4), ArcArmHit );
    ifnew(TText, TextArmHit);
    TextArmHit->SetTextSize(0.08);
    TextArmHit->SetText(10, 4, "Good region");

    // PadSouthHitTime = new TPad("PadSouthHitTime", "PadSouthHitTime", 0.00, 0.55, 0.50, 0.90, 0);
    // PadNorthHitTime = new TPad("PadNorthHitTime", "PadNorthHitTime", 0.50, 0.55, 1.00, 0.90, 0);

    if (PadAvrHitTime)
    {
      ifnew(TLine(bbc_onlmon::BBC_MIN_REGULAR_TDC1_MEAN, 0,
                  bbc_onlmon::BBC_MIN_REGULAR_TDC1_MEAN, 10000),
            LineAvrHitTime[0]);
      ifnew(TLine(bbc_onlmon::BBC_MAX_REGULAR_TDC1_MEAN, 0,
                  bbc_onlmon::BBC_MAX_REGULAR_TDC1_MEAN, 10000),
            LineAvrHitTime[1]);

      ifnew(TArrow(bbc_onlmon::BBC_MIN_REGULAR_TDC1_MEAN, 10000,
                   bbc_onlmon::BBC_MAX_REGULAR_TDC1_MEAN, 10000),
            ArrowAvrHitTime);
      ArrowAvrHitTime->SetOption("<>");

      ifnew(TText, TextAvrHitTime);
      TextAvrHitTime->SetTextSize(0.08);
      TextAvrHitTime->SetText(bbc_onlmon::BBC_MAX_REGULAR_TDC1_MEAN, 10000, " good mean");

      LineAvrHitTime[0]->SetLineColor(9);
      LineAvrHitTime[1]->SetLineColor(9);
      ArrowAvrHitTime->SetLineColor(9);
      TextAvrHitTime->SetTextColor(9);
    }

    if (PadSouthHitTime)
    {
      ifnew(TLine(bbc_onlmon::BBC_MIN_REGULAR_TDC1_MEAN, 0,
                  bbc_onlmon::BBC_MIN_REGULAR_TDC1_MEAN, 10000),
            LineSouthHitTime[0]);
      ifnew(TLine(bbc_onlmon::BBC_MAX_REGULAR_TDC1_MEAN, 0,
                  bbc_onlmon::BBC_MAX_REGULAR_TDC1_MEAN, 10000),
            LineSouthHitTime[1]);

      ifnew(TArrow(bbc_onlmon::BBC_MIN_REGULAR_TDC1_MEAN, 10000,
                   bbc_onlmon::BBC_MAX_REGULAR_TDC1_MEAN, 10000),
            ArrowSouthHitTime);
      ArrowSouthHitTime->SetOption("<>");

      ifnew(TText, TextSouthHitTime);
      TextSouthHitTime->SetTextSize(0.08);
      TextSouthHitTime->SetText(bbc_onlmon::BBC_MAX_REGULAR_TDC1_MEAN, 10000, " good mean");

      LineSouthHitTime[0]->SetLineColor(9);
      LineSouthHitTime[1]->SetLineColor(9);
      ArrowSouthHitTime->SetLineColor(9);
      TextSouthHitTime->SetTextColor(9);
    }

    if (PadNorthHitTime)
    {
      ifnew(TLine(bbc_onlmon::BBC_MIN_REGULAR_TDC1_MEAN, 0,
                  bbc_onlmon::BBC_MIN_REGULAR_TDC1_MEAN, 10000),
            LineNorthHitTime[0]);
      ifnew(TLine(bbc_onlmon::BBC_MAX_REGULAR_TDC1_MEAN, 0,
                  bbc_onlmon::BBC_MAX_REGULAR_TDC1_MEAN, 10000),
            LineNorthHitTime[1]);

      ifnew(TArrow(bbc_onlmon::BBC_MIN_REGULAR_TDC1_MEAN, 10000,
                   bbc_onlmon::BBC_MAX_REGULAR_TDC1_MEAN, 10000),
            ArrowNorthHitTime);
      ArrowNorthHitTime->SetOption("<>");

      ifnew(TText, TextNorthHitTime);
      TextNorthHitTime->SetTextSize(0.08);
      TextNorthHitTime->SetText(bbc_onlmon::BBC_MAX_REGULAR_TDC1_MEAN, 10000, " good mean");

      LineNorthHitTime[0]->SetLineColor(9);
      LineNorthHitTime[1]->SetLineColor(9);
      ArrowNorthHitTime->SetLineColor(9);
      TextNorthHitTime->SetTextColor(9);
    }

    transparent[1] = new TPad("transparent1", "this does not show", 0, 0, 1, 1, 0, 0);
    transparent[1]->SetFillStyle(4000);
    transparent[1]->Draw();
  }

  else if (name == "BbcMon3")
  {
    std::cout << "Creating Canvas BbcMon3..." << std::endl;

    TC[2] = new TCanvas("BbcMon3", "BBC/MBD status view for Expert", -xsize / 2, 0, xsize / 2, ysize);
    gSystem->ProcessEvents();
    TC[2]->cd();


    PadTop[2] = new TPad("PadTop0", "PadTop0", 0.00, 0.90, 1.00, 1.00, 0, 0, 0);
    PadTdcOver[0] = new TPad("PadTdcOverSouth0", "PadTdcOverSouth0", 0.00, 0.30, 0.50, 0.60, 0, 0, 0);
    PadTdcOver[1] = new TPad("PadTdcOverNorth1", "PadTdcOverNorth1", 0.50, 0.60, 1.00, 0.90, 0, 0, 0);
    PadnHit[1] = new TPad("PadnHitNorth", "PadnHitNorth", 0.50, 0.03, 1.00, 0.30, 0, 0, 0);
    PadnHit[0] = new TPad("PadnHitSouth", "PadnHitSouth", 0.00, 0.03, 0.50, 0.30, 0);
    PadnHitStatus = new TPad("PadnHitStatus", "PadnHitStatus", 0.00, 0.00, 1.00, 0.03, 0, 0, 0);

    PadTop[2]->Draw();
    for (int side = 0; side < nSIDE; side++)
    {
      PadTdcOver[side]->Draw();
      PadnHit[side]->Draw();
      PadnHitStatus->Draw();

      ifnew(TBox(0.5, bbc_onlmon::BBC_nHIT_MB_MIN[side], nPMT_1SIDE_BBC + .5, bbc_onlmon::BBC_nHIT_MB_MAX[side]), BoxnHit[0][side]);
      BoxnHit[0][side]->SetFillColor(5);
      BoxnHit[0][side]->SetLineColor(3);

      ifnew(TBox(0.5, bbc_onlmon::BBC_nHIT_LASER_MIN[side], nPMT_1SIDE_BBC + .5, bbc_onlmon::BBC_nHIT_LASER_MAX[side]), BoxnHit[1][side]);
      BoxnHit[1][side]->SetFillColor(7);
      BoxnHit[1][side]->SetLineColor(4);
    }

    for (auto &side : BoxTdcOver)
    {
      ifnew(TBox(0.5, bbc_onlmon::BBC_TDC_OVERFLOW_REGULAR_MIN, nPMT_1SIDE_BBC + .5, bbc_onlmon::BBC_TDC_OVERFLOW_REGULAR_MAX), side);
      side->SetFillColor(5);
      side->SetLineColor(3);
    }
    ifnew(TText, TextnHitStatus);
    TextnHitStatus->SetTextSize(0.7);
    TextnHitStatus->SetText(0.05, 0.5, "Red Square : Collision Event  / Blue Triangle : Laser Event");
    transparent[2] = new TPad("transparent2", "this does not show", 0, 0, 1, 1, 0, 0);
    transparent[2]->SetFillStyle(4000);
    transparent[2]->Draw();
  }

  //
  // 4th Page
  else if (name == "BbcMon4")
  {
    std::cout << "Creating Canvas BbcMon4..." << std::endl;

    // ifnew( TText, TextBbcSummaryTrigRate );
    ifnew(TLatex, TextBbcSummaryTrigRate);
    // TextBbcSummaryTrigRate->SetTextSize(0.109);
    TextBbcSummaryTrigRate->SetTextSize(0.08);
    // TextBbcSummaryTrigRate->SetTextSize(0.075);

    TC[3] = new TCanvas("BbcMon4", "Beam/Trigger status view for Shift clew", -1, 0, xsize / 2, ysize);
    gSystem->ProcessEvents();
    TC[3]->cd();

    PadTop[3] = new TPad("PadTop3", "PadTop3", 0.00, 0.90, 1.00, 1.00, 0, 0, 0);
    PadTop[3]->Draw();

    PadHitTime[0] = new TPad("PadHitTimeSouth1", "PadHitTimeSouth1", 0.00, 0.60, 0.50, 0.90, 0, 0, 0);
    PadHitTime[1] = new TPad("PadHitTimeNorth1", "PadHitTimeNorth1", 0.50, 0.60, 1.00, 0.90, 0, 0, 0);
    PadHitTime[0]->SetLeftMargin(0.17);
    PadHitTime[1]->SetLeftMargin(0.17);
    PadHitTime[0]->Draw();
    PadHitTime[1]->Draw();

    for (int side = 0; side < nSIDE; side++)
    {
      ifnew(TLine(bbc_onlmon::BBC_MIN_REGULAR_TDC1_MEAN, 0, bbc_onlmon::BBC_MIN_REGULAR_TDC1_MEAN, 10000), LineHitTime[side][0]);
      ifnew(TLine(bbc_onlmon::BBC_MAX_REGULAR_TDC1_MEAN, 0, bbc_onlmon::BBC_MAX_REGULAR_TDC1_MEAN, 10000), LineHitTime[side][1]);

      ifnew(TArrow(bbc_onlmon::BBC_MIN_REGULAR_TDC1_MEAN, 10000,
                   bbc_onlmon::BBC_MAX_REGULAR_TDC1_MEAN, 10000),
            ArrowHitTime[side]);
      ArrowHitTime[side]->SetOption("<>");

      ifnew(TText, TextHitTime[side]);
      TextHitTime[side]->SetTextSize(0.08);
      TextHitTime[side]->SetText(bbc_onlmon::BBC_MAX_REGULAR_TDC1_MEAN, 10000, " good mean");

      LineHitTime[side][0]->SetLineColor(9);
      LineHitTime[side][1]->SetLineColor(9);
      ArrowHitTime[side]->SetLineColor(9);
      TextHitTime[side]->SetTextColor(9);
    }
    PadZvtx = new TPad("PadZvtx", "PadZvtx", 0.00, 0.30, 0.50, 0.60, 0, 0, 0);

    // for Zvtx
    if (PadZvtx)
    {
      PadZvtx->SetLeftMargin(0.17);
      PadZvtx->Draw();
      ifnew(TLine(bbc_onlmon::BBC_MIN_REGULAR_ZVERTEX_MEAN, 0,
                  bbc_onlmon::BBC_MIN_REGULAR_ZVERTEX_MEAN, 10000),
            LineZvtx[0]);
      ifnew(TLine(bbc_onlmon::BBC_MAX_REGULAR_ZVERTEX_MEAN, 0,
                  bbc_onlmon::BBC_MAX_REGULAR_ZVERTEX_MEAN, 10000),
            LineZvtx[1]);

      ifnew(TArrow(bbc_onlmon::BBC_MIN_REGULAR_ZVERTEX_MEAN, 10000,
                   bbc_onlmon::BBC_MAX_REGULAR_ZVERTEX_MEAN, 10000),
            ArrowZvtx);
      ArrowZvtx->SetOption("<>");

      ifnew(TText, TextZvtx);
      TextZvtx->SetTextSize(0.08);
      TextZvtx->SetText(bbc_onlmon::BBC_MAX_REGULAR_ZVERTEX_MEAN, 10000, " good mean");
      ifnew(TText, TextZvtxNorth);
      TextZvtxNorth->SetTextSize(0.05);
      TextZvtxNorth->SetText(130, 0, "North");
      ifnew(TText, TextZvtxSouth);
      TextZvtxSouth->SetTextSize(0.05);
      TextZvtxSouth->SetText(-160, 0, "South");

      LineZvtx[0]->SetLineColor(9);
      LineZvtx[1]->SetLineColor(9);
      ArrowZvtx->SetLineColor(9);
      TextZvtx->SetTextColor(9);
    }

    PadChargeSum = new TPad("PadChargeSum", "PadCHargeSum", 0.50, 0.30, 1.00, 0.60, 0, 0, 0);
    PadChargeSum->SetLogy();
    PadChargeSum->Draw();

    ifnew(TText, TextNorthChargeSum);
    // SetNDC(): this will use normalized coordinates (0-1/0-1) for drawing no matter what the histogram uses
    TextNorthChargeSum->SetNDC();
    TextNorthChargeSum->SetTextSize(0.08);
    TextNorthChargeSum->SetText(0.3, 0.8, "--North");  // for p+p
    TextNorthChargeSum->SetTextColor(4);

    ifnew(TText, TextSouthChargeSum);
    // SetNDC(): this will use normalized coordinates (0-1/0-1) for drawing no matter what the histogram uses
    TextSouthChargeSum->SetNDC();
    TextSouthChargeSum->SetTextSize(0.08);
    TextSouthChargeSum->SetText(0.6, 0.8, "--South");  // for p+p
    TextSouthChargeSum->SetTextColor(2);

    PadAdc = new TPad("PadAdc", "PadAdc", 0.00, 0.00, 1.00, 0.30, 0, 0, 0);
    PadAdc->SetLogz();
    PadAdc->Draw();
    transparent[3] = new TPad("transparent0", "this does not show", 0, 0, 1, 1);
    transparent[3]->SetFillStyle(4000);
    transparent[3]->Draw();

    // PadHitTime[1][0]->SetLogy();
    // PadHitTime[1][1]->SetLogy();

    /*
       ifnew( TPad("PadButton" , "PadButton" , 0.00, 0.70, 1.00, 0.90, 0), PadButton);
       ifnew( TPad("PadMultiView" , "PadMultiView" , 0.10, 0.20, 1.00, 0.70, 0), PadMultiView);
       ifnew( TPad("PadWarnings" , "PadWarnings" , 0.00, 0.00, 1.00, 0.20, 0), PadWarnings);
       PadMultiView->Draw();
       PadWarnings->Draw();
       PadButton->Draw();
       ifnew( TPaveText(0.05, 0.05, 0.95, 0.95) , PaveWarnings );
       */

    // TC[3]->SetEditable(0);
  }

  //

  /*
     if ( ( strlen(name) < 1 || strcmp(name, "BBCMon1") == 0 ) ||
     ( strlen(name) < 1 || strcmp(name, "BBCMon2") == 0 ) ||
     ( strlen(name) < 1 || strcmp(name, "BBCMon3") == 0 ) ||
     ( strlen(name) < 1 || strcmp(name, "BBCMon4") == 0 ) )
     */
  {
    ifnew(TPaveText(0.05, 0.65, 0.70, 0.95), PaveTop);
    PaveTop->AddText("BBC/MBD ONLINE MONITOR");
    ifnew(TText, TextTop);
    TextTop->SetTextSize(0.25);
  }

  std::cout << "Initialize completed" << std::endl;

  return 0;
}

int BbcMonDraw::Draw(const std::string &what)
{
  PRINT_DEBUG("In BbcMonDraw::Draw()");
  std::cout << what << std::endl;

  TStyle *oldStyle = gStyle;
  bbcStyle->cd();

  if (!gROOT->FindObject("BbcMon1"))
  {
    PRINT_DEBUG("In BbcMonDraw::Draw(), BbcMon1");
    TC[0] = nullptr;
    if (what == "ALL" || what == "FIRST" || what == "BbcMon1" || what == "BbcMonitor")
    {
      MakeCanvas("BbcMon1");
    }
  }
  if (!gROOT->FindObject("BbcMon2"))
  {
    PRINT_DEBUG("In BbcMonDraw::Draw(), BbcMon2");
    TC[1] = nullptr;
    if (what == "ALL" || what == "SECOND" || what == "BbcMon2" || what == "VertexMonitor")
    {
      MakeCanvas("BbcMon2");
    }
  }
  if (!gROOT->FindObject("BbcMon3"))
  {
    TC[2] = nullptr;
    if (what == "BbcMon3" || what == "BbcMonitor")
    {
      MakeCanvas("BbcMon3");
    }
  }
  // 4th canvas
  if (!gROOT->FindObject("BbcMon4"))
  {
    TC[3] = nullptr;
    if (what == "BbcMon4" || what == "VertexMonitor")
    {
      MakeCanvas("BbcMon4");
    }
  }
  // Histogram
  TH2 *bbc_adc{nullptr};
  TH2 *bbc_tdc{nullptr};
  // TH2 *bbc_tdc_overflow;
  // TH1 *bbc_tdc_overflow_each[nPMT_BBC];
  //TH1 *bbc_nhit[nTRIGGER];
  TH2 *bbc_tdc_armhittime{nullptr};
  TH1 *bbc_zvertex{nullptr};
  TH1 *bbc_zvertex_bbll1{nullptr};
  TH1 *bbc_nevent_counter{nullptr};
  TH2 *bbc_tzero_zvtx{nullptr};
  TH1 *bbc_avr_hittime{nullptr};
  TH1 *bbc_south_hittime{nullptr};
  TH1 *bbc_north_hittime{nullptr};
  TH1 *bbc_south_chargesum{nullptr};
  TH1 *bbc_north_chargesum{nullptr};
  TH1 *bbc_prescale_hist{nullptr};
  TH2 *bbc_time_wave{nullptr};
  TH2 *bbc_charge_wave{nullptr};
  TH2 *bbc_south_hitmap{nullptr};
  TH2 *bbc_north_hitmap{nullptr};

  ClearWarning();

  std::ostringstream otext;
  std::string text;
  std::string textok;
  // std::ostringstream textok;
  std::ostringstream name;

  // ---------------------------------------------------------------------------------
  // get pointer for each histrams

  OnlMonClient *cl = OnlMonClient::instance();

  PRINT_DEBUG("Start Getting Histogram");

  name << "bbc_tdc";
  bbc_tdc = static_cast<TH2 *>(cl->getHisto("BBCMON_0",name.str().c_str()));
  name.str("");

  /*
  name << "bbc_tdc_overflow" ;
  bbc_tdc_overflow =  static_cast<TH2 *> (cl->getHisto("BBCMON_0",name.str().c_str()));
  name.str("");

  for ( int i = 0 ; i < nPMT_BBC ; i++ )
  {
    name << "bbc_tdc_overflow_" << setw(3) << setfill('0') << i ;
    bbc_tdc_overflow_each[i] = cl->getHisto("BBCMON_0",name.str().c_str());
    name.str("");
  }
  */

  bbc_adc = static_cast<TH2 *>(cl->getHisto("BBCMON_0","bbc_adc"));
  ifdelete(Adc);
  for (int i = 0; i < nCANVAS; i++)
  {
    if (TC[i])
    {
      transparent[i]->Clear();  // clear dead server msg if it was printed before
    }
  }

  if (bbc_adc)
  {
    Adc = static_cast<TH2 *>(bbc_adc->Clone());
  }
  else
  {
    // here I assume server is dead and there must be a message saying so
    std::cout << "SERVER IS DEAD, no bbc_adc" << std::endl;
    for (auto &i : TC)
    {
      if (i)
      {
        i->cd();
        i->Clear("D");
        //        DrawDeadServer(transparent[i]);
        i->Update();
      }
    }
    return -1;
  }

  bbc_tdc_armhittime = static_cast<TH2 *>(cl->getHisto("BBCMON_0","bbc_tdc_armhittime"));
  ifdelete(ArmHit);
  ArmHit = static_cast<TH2 *>(bbc_tdc_armhittime->Clone());

  bbc_zvertex = cl->getHisto("BBCMON_0","bbc_zvertex");
  ifdelete(Zvtx);
  Zvtx = static_cast<TH1 *>(bbc_zvertex->Clone());

  bbc_zvertex_bbll1 = cl->getHisto("BBCMON_0","bbc_zvertex_bbll1");
  ifdelete(Zvtx_bbll1);
  Zvtx_bbll1 = static_cast<TH1 *>(bbc_zvertex_bbll1->Clone());

  bbc_nevent_counter = cl->getHisto("BBCMON_0","bbc_nevent_counter");

  bbc_tzero_zvtx = static_cast<TH2 *>(cl->getHisto("BBCMON_0","bbc_tzero_zvtx"));
  ifdelete(TzeroZvtx);
  TzeroZvtx = static_cast<TH2 *>(bbc_tzero_zvtx->Clone());

  bbc_avr_hittime = cl->getHisto("BBCMON_0","bbc_avr_hittime");
  ifdelete(AvrHitTime);
  AvrHitTime = static_cast<TH1 *>(bbc_avr_hittime->Clone());

  bbc_north_hittime = cl->getHisto("BBCMON_0","bbc_north_hittime");
  ifdelete(NorthHitTime);
  NorthHitTime = static_cast<TH1 *>(bbc_north_hittime->Clone());

  bbc_south_hittime = cl->getHisto("BBCMON_0","bbc_south_hittime");
  ifdelete(SouthHitTime);
  SouthHitTime = static_cast<TH1 *>(bbc_south_hittime->Clone());

  bbc_south_chargesum = cl->getHisto("BBCMON_0","bbc_south_chargesum");
  ifdelete(SouthChargeSum);
  SouthChargeSum = static_cast<TH1 *>(bbc_south_chargesum->Clone());

  bbc_north_chargesum = cl->getHisto("BBCMON_0","bbc_north_chargesum");
  ifdelete(NorthChargeSum);
  NorthChargeSum = static_cast<TH1 *>(bbc_north_chargesum->Clone());

  bbc_prescale_hist = cl->getHisto("BBCMON_0","bbc_prescale_hist");
  ifdelete(Prescale_hist);
  Prescale_hist = static_cast<TH1 *>(bbc_prescale_hist->Clone());

  bbc_time_wave = static_cast<TH2 *>( cl->getHisto("BBCMON_0","bbc_time_wave") );
  ifdelete(TimeWave);
  TimeWave = static_cast<TH2 *>(bbc_time_wave->Clone());

  bbc_charge_wave = static_cast<TH2 *>( cl->getHisto("BBCMON_0","bbc_charge_wave") );
  ifdelete(ChargeWave);
  ChargeWave = static_cast<TH2 *>(bbc_charge_wave->Clone());

  bbc_south_hitmap = static_cast<TH2 *>( cl->getHisto("BBCMON_0","bbc_south_hitmap") );
  ifdelete(SouthHitMap);
  SouthHitMap = static_cast<TH2 *>(bbc_south_hitmap->Clone());

  bbc_north_hitmap = static_cast<TH2 *>( cl->getHisto("BBCMON_0","bbc_north_hitmap") );
  ifdelete(NorthHitMap);
  NorthHitMap = static_cast<TH2 *>(bbc_north_hitmap->Clone());

  PRINT_DEBUG("Start Creating graphs");

  // Create HitTime projection ------------------------------------------

  for (int side = 0; side < nSIDE; side++)
  {
    ifdelete(HitTime[side]);

    name << "FitHitTime" << bbc_onlmon::SIDE_Str[side];
    HitTime[side] = bbc_tdc->ProjectionY(name.str().c_str(), side * 64, side * 64 + 63);
    name.str("");

    name << bbc_onlmon::SIDE_Str[side] << " BBC/MBD TDC Distribution";
    // name << bbc_onlmon::SIDE_Str[side] << " BBC/MBD TDC" << tdc << " Distribution(Trigger:BBLL1)" ; // Run14 AuAu 15GeV 2014.02.23
    HitTime[side]->SetTitle(name.str().c_str());
    ifnew(TF1(name.str().c_str(), "gaus"), FitHitTime[side]);
    name.str("");
    FitHitTime[side]->SetLineWidth(1);  // 0.05 was used fro run3
    FitHitTime[side]->SetLineColor(4);
  }

  // Create ZVertex Fit

  ifnew(TF1("FitZvtx", "gaus"), FitZvtx);
  ifnew(TF1("FitAvrHitTime", "gaus"), FitAvrHitTime);
  ifnew(TF1("FitNorthHitTime", "gaus"), FitNorthHitTime);
  ifnew(TF1("FitSouthHitTime", "gaus"), FitSouthHitTime);
  FitZvtx->SetLineWidth(1);  // 0.05 was used for run3
  FitZvtx->SetLineColor(4);

  // ------------------------------------------------------------------------------
  // making nHit TGraph

  double pmt[nPMT_1SIDE_BBC] = {0.};
  double zero[nPMT_1SIDE_BBC] = {0.};
  double nhitPmt[nTRIGGER][nSIDE][nPMT_1SIDE_BBC];
  double nhit_total = bbc_nevent_counter->GetBinContent(2);
  double nhit[nTRIGGER];
  nhit[0] = bbc_nevent_counter->GetBinContent(2);
  nhit[1] = bbc_nevent_counter->GetBinContent(3);

  for (int i = 0; i < nPMT_1SIDE_BBC; i++)
  {
    pmt[i] = i + 1;
    // zero[i] = 0;
  }

  PRINT_DEBUG("Creating nHit Graph");

  for (int side = 0; side < nSIDE; side++)
  {
    for (int trig = 0; trig < nTRIGGER; trig++)
    {
      for (int i = 0; i < nPMT_1SIDE_BBC; i++)
      {
        if (nhit[trig] != 0)
        {
          // RUN12: to ignore laser rate for ch8, fiber is broken.
          if (side == 0 && trig == 1 && i == 7)
          {
            nhitPmt[trig][side][i] = 999;
          }
          else
          {
            //nhitPmt[trig][side][i] = bbc_nhit[trig]->GetBinContent(i + side * nPMT_1SIDE_BBC + 1) / nhit[trig];
            nhitPmt[trig][side][i] = 999;
          }
        }
        else
        {
          nhitPmt[trig][side][i] = 0;
        }
      }
      ifnew(TGraph(nPMT_1SIDE_BBC, pmt, nhitPmt[trig][side]), nHit[trig][side]);
      // nHit[trig][side]->GetXaxis()->SetTitleSize(  0.05);
      // nHit[trig][side]->GetYaxis()->SetTitleSize(  0.05);
      // nHit[trig][side]->GetXaxis()->SetTitleOffset(0.70);
      // nHit[trig][side]->GetYaxis()->SetTitleOffset(1.75);
    }
  }
  PRINT_DEBUG("Creating OverFlow Grapth");

  // Create TDC Overflow Graph

  double tdcOverMean[nSIDE][nPMT_1SIDE_BBC];
  double tdcOverErrY[nSIDE][nPMT_1SIDE_BBC];
  for (int iside = 0; iside < nSIDE; iside++)
  {
    for (int i = 0; i < nPMT_1SIDE_BBC; i++)
    {
      tdcOverMean[iside][i] = 0.0;
      tdcOverErrY[iside][i] = 0.0;
    }
  }

  for (int side = 0; side < nSIDE; side++)
  {
    for (int i = 0; i < nPMT_1SIDE_BBC; i++)
    {
      /*
      tdcOverMean[0][side][i] = bbc_tdc_overflow_each[0][i + side * nPMT_1SIDE_BBC]->GetMean() /
        bbccalib->getOverflow0()->getCalibPar(i)->getDeviation();
      tdcOverErrY[0][side][i] = bbc_tdc_overflow_each[0][i + side * nPMT_1SIDE_BBC]->GetRMS() /
        bbccalib->getOverflow0()->getCalibPar(i)->getDeviation();

      tdcOverMean[1][side][i] = bbc_tdc_overflow_each[1][i + side * nPMT_1SIDE_BBC]->GetMean() /
        bbccalib->getOverflow1()->getCalibPar(i)->getDeviation();
      tdcOverErrY[1][side][i] = bbc_tdc_overflow_each[1][i + side * nPMT_1SIDE_BBC]->GetRMS() /
        bbccalib->getOverflow1()->getCalibPar(i)->getDeviation();
        */

      ifnew(TGraphErrors(nPMT_1SIDE_BBC, pmt, tdcOverMean[side], zero,
                         tdcOverErrY[side]),
            TdcOver[side]);
    }
  }

  // Redraw each Pad

  // ------------------------------------------------------------------------
  // Draw 1st Page
  // ------------------------------------------------------------------------
  std::ostringstream msg;
  PRINT_DEBUG("Drawing Graphs on Canvas");

  // Make TopPave
  time_t evttime = cl->EventTime("CURRENT");

  otext.str("");
  otext << "Run #" << cl->RunNumber();
  otext << " Events: " << nhit_total;
  otext << " Date:" << ctime(&evttime);
  text = otext.str();
  TextTop->SetText(0.01, 0.25, text.c_str());

  if (TC[0])
  {
    TC[0]->cd();

    PadTop[0]->cd();
    PaveTop->Draw();
    TextTop->Draw();

    PadZVertexSummary->cd();

    Zvtx->SetLineColor(4);
    Zvtx->SetFillColor(4);

    Zvtx_bbll1->SetLineColor(4);
    Zvtx_bbll1->SetFillColor(4);

    // Get Maximum at the inside of BBC which is 130cm from center;
    float maxEntries = 10;

    for (int i = 0; i < Zvtx->GetNbinsX(); i++)
    {
      if (fabs(Zvtx->GetBinCenter(i)) < 130)
      {
        if (maxEntries < Zvtx->GetBinContent(i + 1))
        {
          maxEntries = Zvtx->GetBinContent(i + 1);
        }
      }
    }

    // Fit No-Vertex Distribution
    FitZvtx->SetRange(-75, 75);
    FitZvtx->SetLineColor(7);
    //Zvtx->Fit("FitZvtx", "LRQ");
    Zvtx->Fit("FitZvtx", "R");

    // here we get the relative scaling right to put all on the same plot
    // the binning might be different, so we first find the bins corresponding to
    // +- 20cm and then get the Integral corrected for the binwidth
    // this is then used to get the histograms on the same scale
    /*
       int lowbin = Zvtx_zdc->GetXaxis()->FindBin(-20.);
       int hibin = Zvtx_zdc->GetXaxis()->FindBin(20.);
       double zdc_count = Zvtx_zdc->Integral(lowbin, hibin, "width");
       */

    /*
       int lowbin = Zvtx->GetXaxis()->FindBin(-20.);
       int hibin = Zvtx->GetXaxis()->FindBin(20.);
       double bbc_count = Zvtx->Integral(lowbin, hibin, "width");
       */

    // std::cout << "the ratio of integral (-30cm < ZVertex < 30cm) between BBLL1 without vtx cut and ZDC : " << bbc_count_novtx / zdc_count << std::endl ;
    // std::cout << "the ratio of integral (-30cm < ZVertex < 30cm) between BBLL1 without vtx cut  and BBLL1 with BBCZ < |30cm|  : " << bbc_count_novtx/bbc_count << std::endl ;

    // Draw ZVertex triggerd variable trigger
    Zvtx->SetMaximum(maxEntries * 1.05);
    Zvtx->SetTitle("BBC/MBD ZVertex (south<-->north)");
    // PadZVertex->DrawFrame(-160,0,160,maxEntries*1.05,"Bbc ZVertex (south<-->north)");
    // std::cout << "maxEntries " << maxEntries << std::endl;
    // Zvtx->Draw("hist");
    // Zvtx->Scale(bbc_count_novtx/bbc_count);

    // just in case the bbcll1 with vtx cut histo is empty
    // draw the novertex cut (happens during setup)
    if (Zvtx->GetEntries() > 0)
    {
      Zvtx->Draw("hist");
    }
    // trigger rate between BBCLL1 and Zvertex within +- BBC_ZVERTEX_CUT_FOR_TRIG_RATE
    // bbll1, zdc, bbll1_novtx

    float trig_rate[3], trig_rate_err[3];
    double nevent[3];
    memset(trig_rate, 0, sizeof(trig_rate));
    memset(trig_rate_err, 0, sizeof(trig_rate_err));
    memset(nevent, 0, sizeof(nevent));

    /*
       CalculateTriggerRates(trig_rate[0], trig_rate_err[0], nevent[0],
       trig_rate[1], trig_rate_err[1], nevent[1],
       trig_rate[2], trig_rate_err[2], nevent[2]);
       */

/* chiu
    float sigma_zdc = 0.0;
    // float sigma_zdc_err = 0.0;
    float effic_bbc = 0.0;
    // float effic_bbc_err = 0.0;
    float beamInZdc = 0.0;
    // float beamInZdc_err = 0.0;
    float beamInBbc = 0.0;
    // float beamInBbc_err = 0.0;
*/

    /*
       BeamMonitoring( sigma_zdc, sigma_zdc_err, effic_bbc, effic_bbc_err,
       beamInZdc, beamInZdc_err, beamInBbc, beamInBbc_err);
       */

    /*chiu
    std::ostringstream trig_rate_text;
    trig_rate_text << " #sigma_{ZDC} = " << std::fixed << std::setprecision(3) << sigma_zdc
                   << " #epsilon_{BBC} = " << effic_bbc
                   << " {}^{beam in acceptance}_{       ZDC      } = " << beamInZdc
                   << " {}^{beam in acceptance}_{       BBC      } = " << beamInBbc;
    TextBbcSummaryTrigRate->SetText(0.02, 0.05, trig_rate_text.str().c_str());
    trig_rate_text.str("");
    */

    // float xpos[3] = {0.1, 0.34, 0.58};
    //     float xpos[3] = {0.30, 0.50, 0.75};
    // TextZVertexNotice->Draw();

    TH1 *Zvtx_array[4];  // with narrow
    // TH1 *Zvtx_array[3];
    Zvtx_array[0] = Zvtx;
    //Zvtx_array[1] = Zvtx;
    //Zvtx_array[2] = Zvtx;
    // Zvtx_array[1] = Zvtx_zdc;

    // Show status of ZVertex
    int i = 0;

    //TextZVertex[i]->Draw();

    // scale factor ---------------------------------------------------------------------
    //	std::cout << "  " << i << " " << Prescale_hist->GetBinContent(i + 1) << std::endl;
    otext.str("");
    // otext << "( "<< Prescale_hist->GetBinContent(i+1)<<" )";
    // otext << "( "<< Prescale_hist->GetBinContent(i+1)<<" )"<<" ";
    // otext << nevent[i] ;
    // otext.precision(8);
    otext << " ( " << Prescale_hist->GetBinContent(i + 1) << " ) "
        << " ";
    // otext << nevent[i]/Prescale_hist->GetBinContent(i+1) ;
    otext << Zvtx_array[i]->GetEntries();

    text = otext.str();
    TextZVertex_scale[i]->SetText(xpos[i], 0.50, text.c_str());
    //TextZVertex_scale[i]->Draw();

    // mean and RMS ---------------------------------------------------------------------
    otext.str("");
    otext << ((float) int(Zvtx_array[i]->GetMean() * 10)) / 10.0 << "cm ( "
        << ((float) int(Zvtx_array[i]->GetRMS() * 10)) / 10.0 << " cm) ";
    text = otext.str();

    TextZVertex_mean[i]->SetText(xpos[i], 0.25, text.c_str());
    //TextZVertex_mean[i]->Draw();
    // otext.precision(6);

    //TextZVertex[i]->Draw();

    /*
    TextZVertex_scale[i]->SetText(0.00, 0.50, "(Scale Fac.) #Evt.");
    TextZVertex_scale[i]->Draw();
    // TextZVertex_mean[i]->SetText(0.05, 0.25, "Vertex Mean (RMS)");
    TextZVertex_mean[i]->SetText(0.00, 0.25, "Vertex Mean (RMS)");
    TextZVertex_mean[i]->Draw();
    */

    // Draw Status
    otext.str("");
    otext << "Z_{All Trigs}^{Fit}= " << ((float) int(FitZvtx->GetParameter(1) * 10)) / 10.0 << " cm";
        //	      << " #pm " << ((float)int(FitZvtx->GetParError(1)*10))/10.0
        //	      << " #pm " << ((float)int(FitZvtx->GetParError(2)*10))/10.0

    text = otext.str();
    //TextZvtxStatus[0]->SetText(0.0, 0.85, text.c_str());
    TextZvtxStatus[0]->SetText(-230., maxEntries*0.8, text.c_str());
    TextZvtxStatus[0]->SetTextSize(0.10);
    TextZvtxStatus[0]->Draw();

    otext.str("");
    otext << "#sigma = " << int(FitZvtx->GetParameter(2)) << " cm";
    text = otext.str();
    TextZvtxStatus[1]->SetText(100., maxEntries*0.8, text.c_str());
    TextZvtxStatus[1]->SetTextSize(0.10);
    TextZvtxStatus[1]->Draw();

    //chiu TextBbcSummaryTrigRate->Draw();

    PadZVertex->cd();

    if (Zvtx_bbll1->GetEntries() > 0)
    {
      Zvtx_bbll1->GetXaxis()->SetRangeUser(-30,30);
      Zvtx_bbll1->Draw("hist");
    }

    /*
    // replaced with hitmap
    PadTzeroZVertex->cd();
    TzeroZvtx->Draw("colz");

    // insert line
    LineTzeroZvtx[0]->Draw();
    LineTzeroZvtx[1]->Draw();
    LineTzeroZvtx[2]->Draw();
    LineTzeroZvtx[3]->Draw();

    // insert text
    TextTzeroZvtx->SetText(10, 4, "Good region");
    TextTzeroZvtx->Draw();
    */

    double nevents = bbc_nevent_counter->GetBinContent(2);
    PadSouthHitMap->cd();
    SouthHitMap->Scale(1.0/nevents);
    SouthHitMap->Draw("colz");

    PadNorthHitMap->cd();
    NorthHitMap->Scale(1.0/nevents);
    NorthHitMap->Draw("colz");
  }


  //  bbc_t0_pave->Draw("same");

  // ------------------------------------------------------------------------
  // Draw 2nd Page
  // ------------------------------------------------------------------------
  if (TC[1])
  {
    TC[1]->cd();

    PadTop[1]->cd();
    PaveTop->Draw();
    TextTop->Draw();

    if (PadAvrHitTime)
    {
      PadAvrHitTime->cd();
      AvrHitTime->Draw();

      float rangemin;
      float rangemax;
      int npeak = tspec->Search(AvrHitTime, 2, "goff");  // finds the highest peak, draws marker
      if (npeak < 3)                                     // no center peak
      {
        AvrHitTime->Fit("FitAvrHitTime", "QN0L");
        rangemin = FitAvrHitTime->GetParameter(1) - FitAvrHitTime->GetParameter(2);
        rangemax = FitAvrHitTime->GetParameter(1) + FitAvrHitTime->GetParameter(2);
      }
      else
      {
        double *peakpos = tspec->GetPositionX();
        float centerpeak = peakpos[0];
        float sidepeak[2];
        if (peakpos[2] > peakpos[1])
        {
          sidepeak[0] = peakpos[1];
          sidepeak[1] = peakpos[2];
        }
        else
        {
          sidepeak[1] = peakpos[1];
          sidepeak[0] = peakpos[2];
        }
        rangemin = centerpeak - (centerpeak - sidepeak[0]) / 2.;
        rangemax = centerpeak + (sidepeak[1] - centerpeak) / 2.;
      }

      FitAvrHitTime->SetRange(rangemin, rangemax);
      AvrHitTime->Fit("FitAvrHitTime", "QRL");
      FitAvrHitTime->Draw("same");

      float height = AvrHitTime->GetMaximum();
      FitAvrHitTime->Draw("same");

      LineAvrHitTime[1]->SetY2(height);
      LineAvrHitTime[0]->SetY2(height);
      ArrowAvrHitTime->SetY1(height * 0.90);
      ArrowAvrHitTime->SetY2(height * 0.90);
      TextAvrHitTime->SetY(height * 0.88);
      LineAvrHitTime[0]->Draw();
      LineAvrHitTime[1]->Draw();
      ArrowAvrHitTime->Draw();
      TextAvrHitTime->Draw();
    }

    if (PadTimeWave)
    {
      PadTimeWave->cd();
      TimeWave->Draw("colz");
    }

    if (PadSouthHitTime)
    {
      PadSouthHitTime->cd();
      SouthHitTime->Draw();
      float rangemin;
      float rangemax;
      int npeak = tspec->Search(SouthHitTime, 2, "goff");  // finds the highest peak, draws marker
      if (npeak < 3)                                       // no center peak
      {
        SouthHitTime->Fit("FitSouthHitTime", "QN0L");
        rangemin = FitSouthHitTime->GetParameter(1) - FitSouthHitTime->GetParameter(2);
        rangemax = FitSouthHitTime->GetParameter(1) + FitSouthHitTime->GetParameter(2);
      }
      else
      {
        double *peakpos = tspec->GetPositionX();
        float centerpeak = peakpos[0];
        float sidepeak[2];
        if (peakpos[2] > peakpos[1])
        {
          sidepeak[0] = peakpos[1];
          sidepeak[1] = peakpos[2];
        }
        else
        {
          sidepeak[1] = peakpos[1];
          sidepeak[0] = peakpos[2];
        }
        rangemin = centerpeak - (centerpeak - sidepeak[0]) / 2.;
        rangemax = centerpeak + (sidepeak[1] - centerpeak) / 2.;
      }

      FitSouthHitTime->SetRange(rangemin, rangemax);
      SouthHitTime->Fit("FitSouthHitTime", "QRL");
      FitSouthHitTime->Draw("same");

      FitSouthHitTime->Draw("same");

      /*
      // Lines to indicate good mean
      float height = SouthHitTime->GetMaximum();
      LineSouthHitTime[1]->SetY2(height);
      LineSouthHitTime[0]->SetY2(height);
      ArrowSouthHitTime->SetY1(height * 0.90);
      ArrowSouthHitTime->SetY2(height * 0.90);
      TextSouthHitTime->SetY(height * 0.88);
      LineSouthHitTime[0]->Draw();
      LineSouthHitTime[1]->Draw();
      ArrowSouthHitTime->Draw();
      TextSouthHitTime->Draw();
      */
    }

    if (PadNorthHitTime)
    {
      PadNorthHitTime->cd();
      NorthHitTime->Draw();
      float rangemin;
      float rangemax;
      int npeak = tspec->Search(NorthHitTime, 2, "goff");  // finds the highest peak, draws marker
      if (npeak < 3)                                       // no center peak
      {
        NorthHitTime->Fit("FitNorthHitTime", "QN0L");
        rangemin = FitNorthHitTime->GetParameter(1) - FitNorthHitTime->GetParameter(2);
        rangemax = FitNorthHitTime->GetParameter(1) + FitNorthHitTime->GetParameter(2);
      }
      else
      {
        double *peakpos = tspec->GetPositionX();
        float centerpeak = peakpos[0];
        float sidepeak[2];
        if (peakpos[2] > peakpos[1])
        {
          sidepeak[0] = peakpos[1];
          sidepeak[1] = peakpos[2];
        }
        else
        {
          sidepeak[1] = peakpos[1];
          sidepeak[0] = peakpos[2];
        }
        rangemin = centerpeak - (centerpeak - sidepeak[0]) / 2.;
        rangemax = centerpeak + (sidepeak[1] - centerpeak) / 2.;
      }

      FitNorthHitTime->SetRange(rangemin, rangemax);
      NorthHitTime->Fit("FitNorthHitTime", "QRL");
      FitNorthHitTime->Draw("same");

      FitNorthHitTime->Draw("same");

      /*
      // Lines to indicate good mean
      float height = NorthHitTime->GetMaximum();
      LineNorthHitTime[1]->SetY2(height);
      LineNorthHitTime[0]->SetY2(height);
      ArrowNorthHitTime->SetY1(height * 0.90);
      ArrowNorthHitTime->SetY2(height * 0.90);
      TextNorthHitTime->SetY(height * 0.88);
      LineNorthHitTime[0]->Draw();
      LineNorthHitTime[1]->Draw();
      ArrowNorthHitTime->Draw();
      TextNorthHitTime->Draw();
      */
    }

    /*
    PadArmHit->cd();
    if (ArmHit)
    {
      ArmHit->Draw("colz");
      ArcArmHit->Draw();
      TextArmHit->Draw();
    }
    */

    if ( PadBbcSummary )
    {
      PadBbcSummary->cd();
      otext.str("");
      otext << "South:" << ((float) int(FitSouthHitTime->GetParameter(1) * 10)) / 10 << "[ns]  ";
      otext << "North:" << ((float) int(FitNorthHitTime->GetParameter(1) * 10)) / 10 << "[ns]  ";
      otext << "...  ";
      if (bbc_onlmon::BBC_MIN_REGULAR_TDC0_MEAN < FitNorthHitTime->GetParameter(1) &&
              bbc_onlmon::BBC_MAX_REGULAR_TDC0_MEAN > FitNorthHitTime->GetParameter(1) &&
              bbc_onlmon::BBC_MIN_REGULAR_TDC0_MEAN < FitSouthHitTime->GetParameter(1) &&
              bbc_onlmon::BBC_MAX_REGULAR_TDC0_MEAN > FitSouthHitTime->GetParameter(1))
      {
          // otext << "OK";
          textok = "                                                         OK";
      }
      else
      {
          textok = " ";
          if (bbc_onlmon::BBC_MIN_WORNING_STATISTICS_FOR_ZVERTEX_MEAN > Zvtx->GetEntries())
          {
              otext << "Too low statistics";
          }
          else
          {
              otext << "Change Global-Offset on V124";
          }
      }
      text = otext.str();
      TextBbcSummaryHitTime[0]->SetText(0.01, 0.75, text.c_str());
      TextBbcSummaryHitTime[0]->Draw();
      text = textok;
      // TextBbcSummaryHitTime[1]->SetText(0.65, 0.75, text.c_str() );
      TextBbcSummaryHitTime[1]->SetText(0.01, 0.75, text.c_str());
      TextBbcSummaryHitTime[1]->SetTextColor(3);
      TextBbcSummaryHitTime[1]->Draw();
    }

    // Global offset
    float delay = (bbc_onlmon::BBC_DEFAULT_OFFSET -
                   ((FitNorthHitTime->GetParameter(1) + FitSouthHitTime->GetParameter(1)) * 0.5));  //[ns]
    otext.str("");
    otext << "   Global offset : ";
    if (int(-2 * delay) == 0)
    {
      // otext << "need not to move ... OK";
      // otext << " ... OK";
      //textok = "                                    ... OK";
    }
    else
    {
      //textok = " ";
      if (!(bbc_onlmon::BBC_MIN_REGULAR_TDC0_MEAN < FitNorthHitTime->GetParameter(1) &&
            bbc_onlmon::BBC_MAX_REGULAR_TDC0_MEAN > FitNorthHitTime->GetParameter(1) &&
            bbc_onlmon::BBC_MIN_REGULAR_TDC0_MEAN < FitSouthHitTime->GetParameter(1) &&
            bbc_onlmon::BBC_MAX_REGULAR_TDC0_MEAN > FitSouthHitTime->GetParameter(1)))
      {
        // otext << "Need to move " << int( -2*delay) << "count";
        otext << "Need to call BBC expert: " << int(-2 * delay) << "count shifted";
      }
      else
      {
        // otext << " ... OK";
        //textok = "                            OK";
      }
    }

/*chiu
    // otext << "...     )" ;
    text = otext.str();
    TextBbcSummaryGlobalOffset[0]->SetText(0.01, 0.50, text.c_str());
    TextBbcSummaryGlobalOffset[0]->Draw();
    text = textok;
    // TextBbcSummaryGlobalOffset[1]->SetText(0.35, 0.50, text.c_str() );
    TextBbcSummaryGlobalOffset[1]->SetText(0.01, 0.50, text.c_str());
    TextBbcSummaryGlobalOffset[1]->SetTextColor(3);
    TextBbcSummaryGlobalOffset[1]->Draw();
    // textok = " ";
*/

    // ZVertex
    /*
       otext.str("");
       otext << "Mean ZVertex:" << ((float)int(FitZvtx->GetParameter(1)*10)) / 10 << "[cm] ";
       if(Zvtx->GetEntries() > bbc_onlmon::BBC_MIN_WORNING_STATISTICS_FOR_ZVERTEX_MEAN ) {
       otext << "( sigma " << ((float)int(FitZvtx->GetParameter(2)*10)) / 10 << "cm)";
       }else{
    // otext << "( RMS " << ((float)int(Zvtx->GetRMS()*10)) / 10 << "cm)";
    }
    otext << " ... ";
    if ( bbc_onlmon::BBC_MIN_REGULAR_ZVERTEX_MEAN < FitZvtx->GetParameter(1) &&
    bbc_onlmon::BBC_MAX_REGULAR_ZVERTEX_MEAN > FitZvtx->GetParameter(1) )
    otext << "OK";
    else {
    if( bbc_onlmon::BBC_MIN_WORNING_STATISTICS_FOR_ZVERTEX_MEAN > Zvtx->GetEntries() )
    otext << "Too low statistics";
    else
    otext << "Ask SL to contact MCR";
    }
    text = otext.str();
    TextBbcSummaryZvertex->SetText(0.01, 0.25, text.c_str() );
    */

/*chiu
    TextBbcSummaryZvertex->SetText(0.01, 0.25, "Shown data are triggered by BBLL1 |z|<130cm");
    TextBbcSummaryZvertex->Draw();
*/
  }

  // ------------------------------------------------------------------------
  // Draw 3rd Page
  // ------------------------------------------------------------------------
  if (TC[2])
  {
    TC[2]->cd();

    // -------------------------------------------------------------------------
    PadTop[2]->cd();
    PaveTop->Draw();
    TextTop->Draw();

    // TDC

    for (int side = 0; side < nSIDE; side++)
    {
      PadTdcOver[side]->cd();
      // ifdelete( FrameTdcOver[side] );
      FrameTdcOver[side] = TC[0]->DrawFrame(0.5, -5, 64.5, 5);
      // FrameTdcOver[side] = gPad->DrawFrame( 0.5, -5, 64.5, 5);

      std::cout << "FrameTdcOver[" << side << "] = " << (unsigned long) FrameTdcOver[side] << std::endl;
      BoxTdcOver[side]->Draw();

      name << bbc_onlmon::SIDE_Str[side] << " BBC/MBD TDC Distribution";
      FrameTdcOver[side]->SetTitle(name.str().c_str());
      name.str("");

      FrameTdcOver[side]->GetXaxis()->SetTitle("PMT Number");
      FrameTdcOver[side]->GetXaxis()->SetTitleSize(0.05);
      FrameTdcOver[side]->GetXaxis()->SetLabelSize(0.05);
      FrameTdcOver[side]->GetYaxis()->SetTitle("Deviation [#sigma]");
      FrameTdcOver[side]->GetYaxis()->SetTitleSize(0.07);
      FrameTdcOver[side]->GetYaxis()->SetTitleOffset(0.50);
      FrameTdcOver[side]->GetYaxis()->SetLabelSize(0.05);
      TdcOver[side]->SetMarkerStyle(21);
      TdcOver[side]->SetMarkerSize(0.5);
      TdcOver[side]->SetMarkerColor(2);
      TdcOver[side]->Draw("P");

      // Check Warning
      if (nhit[0] > 100)
      {
        for (int i = 0; i < nPMT_1SIDE_BBC; i++)
        {
          if (tdcOverMean[side][i] == 0 && tdcOverErrY[side][i] == 0)
          {
            msg.str("");
            msg << "Stop the run (ch " << i + 1 << " is out of the range)";
            std::string wmsg = msg.str();
            Warning(PadTdcOver[side], i, tdcOverMean[side][i], 1, wmsg);
            wmsg.erase();
            msg.str("");
          }
          if (tdcOverMean[side][i] < bbc_onlmon::BBC_TDC_OVERFLOW_REGULAR_MIN)
          {
            msg.str("");
            msg << "ch " << i + 1 << " is too low ( " << std::fixed << std::setprecision(1) << tdcOverMean[side][i] << " #sigma)";
            std::string wmsg = msg.str();
            Warning(PadTdcOver[side], i, tdcOverMean[side][i], 1, wmsg);
            wmsg.erase();
            msg.str("");
          }
          if (tdcOverMean[side][i] > bbc_onlmon::BBC_TDC_OVERFLOW_REGULAR_MAX)
          {
            msg.str("");
            msg << "ch " << i + 1 << " is too high ( " << std::fixed << std::setprecision(1) << tdcOverMean[side][i] << " #sigma)";
            std::string wmsg = msg.str();
            Warning(PadTdcOver[side], i, tdcOverMean[side][i], 1, wmsg);
            wmsg.erase();
            msg.str("");
          }
          if (tdcOverMean[side][i] > bbc_onlmon::BBC_TDC_OVERFLOW_REGULAR_RMS_MAX)
          {
            msg.str("");
            msg << "ch " << i + 1 << " is too wide ( " << std::fixed << std::setprecision(1) << tdcOverErrY[side][i] << " #sigma)";
            std::string wmsg = msg.str();
            Warning(PadTdcOver[side], i, tdcOverMean[side][i], 1, wmsg);
            wmsg.erase();
            msg.str("");
          }
        }
      }

      PadnHit[side]->cd();
      // ifdelete( FramenHit[side] );
      // FramenHit[side] = TC[0]->DrawFrame( 0.5, 0, 64.5, 1);
      FramenHit[side] = gPad->DrawFrame(0.5, 0, 64.5, 1);
      BoxnHit[0][side]->Draw();
      BoxnHit[1][side]->Draw();

      name << bbc_onlmon::SIDE_Str[side] << " BBC/MBD number of Hit per Event";
      FramenHit[side]->SetTitle(name.str().c_str());
      name.str("");

      FramenHit[side]->GetXaxis()->SetTitle("PMT Number");
      FramenHit[side]->GetXaxis()->SetTitleSize(0.05);
      FramenHit[side]->GetXaxis()->SetLabelSize(0.05);
      FramenHit[side]->GetYaxis()->SetTitleSize(0.05);
      FramenHit[side]->GetYaxis()->SetLabelSize(0.05);
      nHit[0][side]->SetMarkerStyle(21);
      nHit[0][side]->SetMarkerSize(0.5);
      nHit[0][side]->SetMarkerColor(2);
      nHit[0][side]->Draw("P");
      nHit[1][side]->SetMarkerStyle(22);
      nHit[1][side]->SetMarkerSize(0.5);
      nHit[1][side]->SetMarkerColor(4);
      nHit[1][side]->Draw("P");

      for (int i = 0; i < nPMT_1SIDE_BBC; i++)
      {
        if (nhit[0] > 100)
        {
          if (nhitPmt[0][side][i] < bbc_onlmon::BBC_nHIT_MB_MIN[side])
          {
            // RUN11: to ignore hit rate since ch29 before FEM input is dead.
            // RUN11: to ignore hit rate since the gain for ch40 is unstable.
            // if( (side==0)&&(i==29 || i==40))
            //{
            // if(i == 29){
            //    std::cout << "ch29 : Ignore hit rate into ch29" << std::endl;
            // }
            // else {
            //    std::cout << "ch40 : Ignore hit rate into ch40" << std::endl;
            // }
            //	continue;
            // }

            msg.str("");
            msg << "Too low hit-rate into ch " << i + 1 << " ("
                << std::fixed << std::setprecision(2) << nhitPmt[0][side][i]
                << "/" << std::setprecision(0) << nhit[0] << "evt)";
            std::string wmsg = msg.str();
            Warning(PadnHit[side], i, nhitPmt[0][side][i], 1, wmsg);
            wmsg.erase();
            msg.str("");
          }
          if (nhitPmt[0][side][i] > bbc_onlmon::BBC_nHIT_MB_MAX[side])
          {
            msg.str("");
            msg << "Too high hit-rate into ch " << i + 1 << " ("
                << std::fixed << std::setprecision(2) << nhitPmt[0][side][i]
                << "/" << std::setprecision(0) << nhit[0] << "evt)";
            std::string wmsg = msg.str();
            Warning(PadnHit[side], i, nhitPmt[0][side][i], 1, wmsg);
            wmsg.erase();
            msg.str("");
          }
        }

        if (nhit[1] > 0)
        {
          if (nhitPmt[1][side][i] < bbc_onlmon::BBC_nHIT_LASER_MIN[side])
          {
            if (side == 0)
            {
              // RUN12: to ignore lack of laser rate, since fiber of ch7 is broken.
              if (i == 7)
              {
                std::cout << "ch7(S8) : Ignore hit rate of laser" << std::endl;
                continue;
              }
              msg.str("");
              msg << "Lack of laser's hit into ch " << i + 1 << " ("
                  << std::fixed << std::setprecision(2) << nhitPmt[1][side][i]
                  << "/" << std::setprecision(0) << nhit[1] << "evt)";
              std::string wmsg = msg.str();
              Warning(PadnHit[side], i, nhitPmt[1][side][i], 1, wmsg);
              wmsg.erase();
              msg.str("");
            }
            else  // North
            {
              msg.str("");
              msg << "Lack of laser's hit into ch " << i + 1 << " ("
                  << std::fixed << std::setprecision(2) << nhitPmt[1][side][i]
                  << "/" << std::setprecision(0) << nhit[1] << "evt)";
              std::string wmsg = msg.str();
              Warning(PadnHit[side], i, nhitPmt[1][side][i], 1, wmsg);
              wmsg.erase();
              msg.str("");
            }
          }
        }
      }
      // if ( nhit[0] == 0 )
      // Warning( PadnHit[side], 1, 0, 1, "No Event" );
      // if ( nhit[1] == 0 )
      // Warning( PadnHit[side], 1, 0, 1, "No Laser Event" );
    }

    PadnHitStatus->cd();
    TextnHitStatus->Draw();
  }

  //  std::cout << "Got Histgram Got-Pad 0" << std::endl;

  // ------------------------------------------------------------------------
  // Draw 4th Page
  // ------------------------------------------------------------------------
  bbcStyle->cd();
  if (TC[3])
  {
    PadTop[3]->cd();
    PaveTop->Draw();
    TextTop->Draw();

    for (int side = 0; side < nSIDE; side++)
    {
      PadHitTime[side]->cd();

      HitTime[side]->Draw();
      float rangemin;
      float rangemax;
      int npeak = tspec->Search(HitTime[side], 2, "goff");  // finds the highest peak, draws marker
      if (npeak < 3)                                        // no center peak
      {
        FitHitTime[side]->SetRange(bbc_onlmon::TDC_FIT_MIN, bbc_onlmon::TDC_FIT_MAX);
        HitTime[side]->Fit(FitHitTime[side]->GetName(), "QRNL");
        rangemax = std::min(bbc_onlmon::TDC_FIT_MAX,
                            FitHitTime[side]->GetParameter(1) + FitHitTime[side]->GetParameter(2));
        rangemin = std::max(bbc_onlmon::TDC_FIT_MIN,
                            FitHitTime[side]->GetParameter(1) - FitHitTime[side]->GetParameter(2));
      }
      else
      {
        double *peakpos = tspec->GetPositionX();
        float centerpeak = peakpos[0];
        float sidepeak[2];
        if (peakpos[2] > peakpos[1])
        {
          sidepeak[0] = peakpos[1];
          sidepeak[1] = peakpos[2];
        }
        else
        {
          sidepeak[1] = peakpos[1];
          sidepeak[0] = peakpos[2];
        }
        rangemin = centerpeak - (centerpeak - sidepeak[0]) / 2.;
        rangemax = centerpeak + (sidepeak[1] - centerpeak) / 2.;
      }

      FitHitTime[side]->SetRange(rangemin, rangemax);

      HitTime[side]->Fit(FitHitTime[side]->GetName(), "QRL");
      float height = HitTime[side]->GetMaximum();
      FitHitTime[side]->Draw("same");

      LineHitTime[side][1]->SetY2(height);
      LineHitTime[side][0]->SetY2(height);
      ArrowHitTime[side]->SetY1(height * 0.90);
      ArrowHitTime[side]->SetY2(height * 0.90);
      TextHitTime[side]->SetY(height * 0.88);
      LineHitTime[side][0]->Draw();
      LineHitTime[side][1]->Draw();
      ArrowHitTime[side]->Draw();
      TextHitTime[side]->Draw();
    }
    // PadWarnings->cd();
    // PaveWarnings->Draw();

    if (PadZvtx)
    {
      PadZvtx->cd();
      Zvtx->Draw();
      Zvtx->Fit("FitZvtx", "QN0L");
      FitZvtx->SetRange(FitZvtx->GetParameter(1) - FitZvtx->GetParameter(2) * 2,
              FitZvtx->GetParameter(1) + FitZvtx->GetParameter(2) * 2);
      Zvtx->Fit("FitZvtx", "QRL");
      FitZvtx->Draw("same");

      float height = Zvtx->GetMaximum();
      FitZvtx->Draw("same");

      LineZvtx[1]->SetY2(height);
      LineZvtx[0]->SetY2(height);
      ArrowZvtx->SetY1(height * 0.90);
      ArrowZvtx->SetY2(height * 0.90);
      TextZvtx->SetY(height * 0.88);
      LineZvtx[0]->Draw();
      LineZvtx[1]->Draw();
      ArrowZvtx->Draw();
      TextZvtx->Draw();
      TextZvtxNorth->Draw();
      TextZvtxSouth->Draw();
    }

    if (PadChargeSum)
    {
      PadChargeSum->cd();
      NorthChargeSum->SetLineColor(4);
      SouthChargeSum->SetLineColor(2);
      NorthChargeSum->Draw();
      SouthChargeSum->Draw("same");

      //       float height = NorthChargeSum->GetMaximum();
      //       TextNorthChargeSum->SetY( height*0.88);
      //       TextSouthChargeSum->SetY( height*0.88);
      TextNorthChargeSum->SetTextColor(4);
      TextSouthChargeSum->SetTextColor(2);
      TextNorthChargeSum->Draw();
      TextSouthChargeSum->Draw();
    }

    if (PadAdc)
    {
       PadAdc->cd();
       Adc->Draw("colz");
    }
  }  // TC[3]

  if (TC[0])
  {
    TC[0]->Update();
  }
  if (TC[1])
  {
    TC[1]->Update();
  }
  if (TC[2])
  {
    TC[2]->Update();
  }
  if (TC[3])
  {
    TC[3]->Update();
  }

  // std::cout << "Got Histgram Got-Pad 1" << std::endl;

  oldStyle->cd();

  //******************************************
  int iret = 0;
  //  int idraw = 0;
  /*
     if (what == "ALL" || what == "FIRST")
     {
     iret += DrawFirst(what);
     idraw++;
     }
     if (what == "ALL" || what == "SECOND")
     {
     iret += DrawSecond(what);
     idraw++;
     }
  if (what == "ALL" || what == "HISTORY")
  {
    iret += DrawHistory(what);
    //    idraw++;
  }
     */

  /*
     if (!idraw)
     {
     std::cout << __PRETTY_FUNCTION__ << " Unimplemented Drawing option: " << what << std::endl;
     iret = -1;
     }
     */

  //******************************************
  return iret;
}

int BbcMonDraw::DrawFirst(const std::string & )
{
  PRINT_DEBUG("In BbcMonDraw::DrawFirst()");
  OnlMonClient *cl = OnlMonClient::instance();
  TH1 *bbcmon_hist1 = cl->getHisto("BBCMON_0","bbc_zvertex");
  TH2 *bbcmon_hist2 = (TH2*)cl->getHisto("BBCMON_0","bbc_tzero_zvtx");

  if (!gROOT->FindObject("BbcMon1"))
  {
      MakeCanvas("BbcMon1");
  }
  TC[0]->SetEditable(1);
  TC[0]->Clear("D");
  Pad[0]->cd();
  if (bbcmon_hist1)
  {
      bbcmon_hist1->DrawCopy();
  }
  else
  {
      DrawDeadServer(transparent[0]);
      TC[0]->SetEditable(0);
      return -1;
  }

  Pad[1]->cd();
  if (bbcmon_hist2)
  {
      bbcmon_hist2->DrawCopy();
  }
  TText PrintRun;
  PrintRun.SetTextFont(62);
  PrintRun.SetTextSize(0.04);
  PrintRun.SetNDC();          // set to normalized coordinates
  PrintRun.SetTextAlign(23);  // center/top alignment
  std::ostringstream runnostream;
  std::string runstring;
  time_t evttime = cl->EventTime("CURRENT");
  // fill run number and event time into string
  runnostream << ThisName << "_1 Run " << cl->RunNumber()
      << ", Time: " << ctime(&evttime);
  runstring = runnostream.str();
  transparent[0]->cd();
  PrintRun.DrawText(0.5, 1., runstring.c_str());
  TC[0]->Update();
  TC[0]->Show();
  TC[0]->SetEditable(0);
  return 0;
}

int BbcMonDraw::DrawSecond(const std::string & )
{
  PRINT_DEBUG("In BbcMonDraw::DrawSecond()");

  OnlMonClient *cl = OnlMonClient::instance();
  TH1 *bbcmon_hist1 = cl->getHisto("BBCMON_0","bbc_zvertex");
  TH2 *bbcmon_hist2 = (TH2*)cl->getHisto("BBCMON_0","bbc_tzero_zvtx");
  if (!gROOT->FindObject("BbcMon2"))
  {
      MakeCanvas("BbcMon2");
  }
  TC[1]->SetEditable(1);
  TC[1]->Clear("D");
  Pad[2]->cd();
  if (bbcmon_hist1)
  {
      bbcmon_hist1->DrawCopy();
  }
  else
  {
      DrawDeadServer(transparent[1]);
      TC[1]->SetEditable(0);
      return -1;
  }
  Pad[3]->cd();
  if (bbcmon_hist2)
  {
      bbcmon_hist2->DrawCopy();
  }
  TText PrintRun;
  PrintRun.SetTextFont(62);
  PrintRun.SetTextSize(0.04);
  PrintRun.SetNDC();          // set to normalized coordinates
  PrintRun.SetTextAlign(23);  // center/top alignment
  std::ostringstream runnostream;
  std::string runstring;
  time_t evttime = cl->EventTime("CURRENT");

  // fill run number and event time into string
  runnostream << ThisName << "_2 Run " << cl->RunNumber()
      << ", Time: " << ctime(&evttime);
  runstring = runnostream.str();
  transparent[1]->cd();
  PrintRun.DrawText(0.5, 1., runstring.c_str());
  TC[1]->Update();
  TC[1]->Show();
  TC[1]->SetEditable(0);
  return 0;
}

int BbcMonDraw::SavePlot(const std::string &what, const std::string &type)
{

  OnlMonClient *cl = OnlMonClient::instance();
  int iret = Draw(what);
  if (iret)  // on error no png files please
  {
      return iret;
  }
  int icnt = 0;
  for (TCanvas *canvas : TC)
  {
    if (canvas == nullptr)
    {
      continue;
    }
    icnt++;
    std::string filename = ThisName + "_" + std::to_string(icnt) + "_" +
      std::to_string(cl->RunNumber()) + "." + type;
    cl->CanvasToPng(canvas, filename);
  }
  return 0;
}

int BbcMonDraw::MakeHtml(const std::string &what)
{
  int iret = Draw(what);
  if (iret)  // on error no html output please
  {
      return iret;
  }

  OnlMonClient *cl = OnlMonClient::instance();

  int icnt = 0;
  for (TCanvas *canvas : TC)
  {
    if (canvas == nullptr)
    {
      continue;
    }
    icnt++;
    // Register the canvas png file to the menu and produces the png file.
    std::string pngfile = cl->htmlRegisterPage(*this, canvas->GetTitle(), std::to_string(icnt), "png");
    cl->CanvasToPng(canvas, pngfile);
  }

  // Now register also EXPERTS html pages, under the EXPERTS subfolder.
  // std::string logfile = cl->htmlRegisterPage(*this, "EXPERTS/Log", "log", "html");
  // std::ofstream out(logfile.c_str());
  // out << "<HTML><HEAD><TITLE>Log file for run " << cl->RunNumber()
  //     << "</TITLE></HEAD>" << std::endl;
  // out << "<P>Some log file output would go here." << std::endl;
  // out.close();

  // std::string status = cl->htmlRegisterPage(*this, "EXPERTS/Status", "status", "html");
  // std::ofstream out2(status.c_str());
  // out2 << "<HTML><HEAD><TITLE>Status file for run " << cl->RunNumber()
  //     << "</TITLE></HEAD>" << std::endl;
  // out2 << "<P>Some status output would go here." << std::endl;
  // out2.close();

  // cl->SaveLogFile(*this);

  return 0;
}

int BbcMonDraw::DrawHistory(const std::string & /* what */)
{
  // you need to provide the following vectors
  // which are filled from the db
  std::vector<float> var;
  std::vector<float> varerr;
  std::vector<time_t> timestamp;
  /*chiu
    std::vector<int> runnumber;
  //std::string varname = "bbcmondummy";
  // this sets the time range from whihc values should be returned
  time_t begin = 0;            // begin of time (1.1.1970)
  time_t end = time(nullptr);  // current time (right NOW)
  int iret = dbvars->GetVar(begin, end, varname, timestamp, runnumber, var, varerr);
  if (iret)
  {
  std::cout << __PRETTY_FUNCTION__ << " Error in db access" << std::endl;
  return iret;
  }
  */
  if (!gROOT->FindObject("BbcMon3"))
  {
      MakeCanvas("BbcMon3");
  }
  // timestamps come sorted in ascending order
  float *x = new float[var.size()];
  float *y = new float[var.size()];
  float *ex = new float[var.size()];
  float *ey = new float[var.size()];
  // int n = var.size();
  for (unsigned int i = 0; i < var.size(); i++)
  {
      //       std::cout << "timestamp: " << ctime(&timestamp[i])
      // 	   << ", run: " << runnumber[i]
      // 	   << ", var: " << var[i]
      // 	   << ", varerr: " << varerr[i]
      // 	   << std::endl;
      x[i] = timestamp[i] - TimeOffsetTicks;
      y[i] = var[i];
      ex[i] = 0;
      ey[i] = varerr[i];
  }

  /* need to implement history for BBC
     Pad[4]->cd(); // neeed to fix
     if (gr[0])
     {
     delete gr[0];
     }
     gr[0] = new TGraphErrors(n, x, y, ex, ey);
     gr[0]->SetMarkerColor(4);
     gr[0]->SetMarkerStyle(21);
     gr[0]->Draw("ALP");
     gr[0]->GetXaxis()->SetTimeDisplay(1);
     gr[0]->GetXaxis()->SetLabelSize(0.03);
  // the x axis labeling looks like crap
  // please help me with this, the SetNdivisions
  // don't do the trick
  gr[0]->GetXaxis()->SetNdivisions(-1006);
  gr[0]->GetXaxis()->SetTimeOffset(TimeOffsetTicks);
  gr[0]->GetXaxis()->SetTimeFormat("%Y/%m/%d %H:%M");
  */
  delete[] x;
  delete[] y;
  delete[] ex;
  delete[] ey;

  /*chiu
    varname = "bbcmoncount";
    iret = dbvars->GetVar(begin, end, varname, timestamp, runnumber, var, varerr);
    if (iret)
    {
    std::cout << __PRETTY_FUNCTION__ << " Error in db access" << std::endl;
    return iret;
    }
    */
  x = new float[var.size()];
  y = new float[var.size()];
  ex = new float[var.size()];
  ey = new float[var.size()];
  // n = var.size();
  for (unsigned int i = 0; i < var.size(); i++)
  {
      //       std::cout << "timestamp: " << ctime(&timestamp[i])
      // 	   << ", run: " << runnumber[i]
      // 	   << ", var: " << var[i]
      // 	   << ", varerr: " << varerr[i]
      // 	   << std::endl;
      x[i] = timestamp[i] - TimeOffsetTicks;
      y[i] = var[i];
      ex[i] = 0;
      ey[i] = varerr[i];
  }

  /* Need to implement
     Pad[5]->cd();
     if (gr[1])
     {
     delete gr[1];
     }
     gr[1] = new TGraphErrors(n, x, y, ex, ey);
     gr[1]->SetMarkerColor(4);
     gr[1]->SetMarkerStyle(21);
     gr[1]->Draw("ALP");
     gr[1]->GetXaxis()->SetTimeDisplay(1);
  // TC[2]->Update();
  //    h1->GetXaxis()->SetTimeDisplay(1);
  //    h1->GetXaxis()->SetLabelSize(0.03);
  gr[1]->GetXaxis()->SetLabelSize(0.03);
  gr[1]->GetXaxis()->SetTimeOffset(TimeOffsetTicks);
  gr[1]->GetXaxis()->SetTimeFormat("%Y/%m/%d %H:%M");
  //    h1->Draw();
  */

  delete[] x;
  delete[] y;
  delete[] ex;
  delete[] ey;

  TC[2]->Update();
  return 0;
}

int BbcMonDraw::DrawDeadServer(TPad *transparent_pad)
{
  transparent_pad->cd();
  TText FatalMsg;
  FatalMsg.SetTextFont(62);
  FatalMsg.SetTextSize(0.1);
  FatalMsg.SetTextColor(4);
  FatalMsg.SetNDC();  // set to normalized coordinates
  FatalMsg.SetTextAlign(23); // center/top alignment
  FatalMsg.DrawText(0.5, 0.9, "BBCMONITOR");
  FatalMsg.SetTextAlign(22); // center/center alignment
  FatalMsg.DrawText(0.5, 0.5, "SERVER");
  FatalMsg.SetTextAlign(21); // center/bottom alignment
  FatalMsg.DrawText(0.5, 0.1, "DEAD");
  transparent_pad->Update();
  return 0;
}


