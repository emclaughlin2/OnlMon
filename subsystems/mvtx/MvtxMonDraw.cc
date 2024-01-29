#include "MvtxMonDraw.h"

#include <onlmon/OnlMonClient.h>
#include <onlmon/OnlMonDB.h>

#include <TAxis.h>  // for TAxis
#include <TCanvas.h>
#include <TDatime.h>
#include <TGaxis.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TH2Poly.h>
#include <TPad.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TText.h>
#include <TStyle.h>
#include <TPaveText.h>
#include <TLatex.h>
#include <TColor.h>
#include <TEllipse.h>

#include <cstring>  // for memset
#include <ctime>
#include <fstream>
#include <iostream>  // for operator<<, basic_ostream, basic_os...
#include <sstream>
#include <vector>  // for vector

MvtxMonDraw::MvtxMonDraw(const std::string &name)
  : OnlMonDraw(name)
{
  // this TimeOffsetTicks is neccessary to get the time axis right
  TDatime T0(2003, 01, 01, 00, 00, 00);
  TimeOffsetTicks = T0.Convert();
  //dbvars = new OnlMonDB(ThisName);
  for (int iFelix = 0; iFelix < NFlx; iFelix++){
    dbvars[iFelix] = new OnlMonDB(Form("MVTXMON_%d",iFelix));
  }
  return;
}

int MvtxMonDraw::Init()
{
  return 0;
}

int MvtxMonDraw::MakeCanvas(const std::string &name)
{
  OnlMonClient *cl = OnlMonClient::instance();
  int xsize = cl->GetDisplaySizeX();
  int ysize = cl->GetDisplaySizeY();
 /* if (name == "MvtxMon1")
  {
    // xpos (-1) negative: do not draw menu bar
    TC[0] = new TCanvas(name.c_str(), "MvtxMon Example Monitor", -1, 0, xsize / 2, ysize);
    // root is pathetic, whenever a new TCanvas is created root piles up
    // 6kb worth of X11 events which need to be cleared with
    // gSystem->ProcessEvents(), otherwise your process will grow and
    // grow and grow but will not show a definitely lost memory leak
    gSystem->ProcessEvents();
    Pad[0] = new TPad("mvtxpad1", "who needs this?", 0.1, 0.7, 0.9, 0.9, 0);
    Pad[1] = new TPad("mvtxpad2", "who needs this?", 0.1, 0.3, 0.9, 0.7, 0);
    Pad[2] = new TPad("mvtxpad3", "who needs this?", 0.1, 0.05, 0.9, 0.3, 0);
    Pad[0]->Draw();
    Pad[1]->Draw();
      Pad[2]->Draw();
    // this one is used to plot the run number on the canvas
    transparent[0] = new TPad("transparent0", "this does not show", 0, 0, 1, 1);
    transparent[0]->SetFillStyle(4000);
    transparent[0]->Draw();
    TC[0]->SetEditable(false);
  }
  else if (name == "MvtxMon2")
  {
    // xpos negative: do not draw menu bar
    TC[1] = new TCanvas(name.c_str(), "MvtxMon2 Example Monitor", -xsize / 2, 0, xsize / 2, ysize);
    gSystem->ProcessEvents();
    Pad[3] = new TPad("mvtxpad4", "who needs this?", 0.1, 0.5, 0.9, 0.9, 0);
    Pad[4] = new TPad("mvtxpad5", "who needs this?", 0.1, 0.05, 0.9, 0.45, 0);
    Pad[3]->Draw();
    Pad[4]->Draw();
    // this one is used to plot the run number on the canvas
    transparent[1] = new TPad("transparent1", "this does not show", 0, 0, 1, 1);
    transparent[1]->SetFillStyle(4000);
    transparent[1]->Draw();
    TC[1]->SetEditable(false);
  }
  else if (name == "MvtxMon3")
  {
    TC[2] = new TCanvas(name.c_str(), "MvtxMon3 Example Monitor", xsize / 2, 0, xsize / 2, ysize);
    gSystem->ProcessEvents();
    Pad[5] = new TPad("mvtxpad6", "who needs this?", 0.1, 0.1, 0.9, 0.9, 0);
    //Pad[6] = new TPad("mvtxpad7", "who needs this?", 0.1, 0.05, 0.9, 0.45, 0);
    Pad[5]->Draw();
    //Pad[6]->Draw();
    // this one is used to plot the run number on the canvas
    //        transparent[2] = new TPad("transparent2", "this does not show", 0, 0, 1, 1);
    //        transparent[2]->SetFillStyle(4000);
    //        transparent[2]->Draw();
    //      TC[2]->SetEditable(0);
  }*/
  if (name == "MvtxMon_HitMap")
  {
    TC[0] = new TCanvas(name.c_str(), "MvtxMon Hitmaps", -1, 0, xsize, ysize);
    gSystem->ProcessEvents();
    //TC[3]->Divide(2,1/*NSTAVE*/);
    Pad[0] = new TPad("mvtxpad0", "who needs this?", 0.05, 0.02, 0.98, 0.92, 0);  
    Pad[0]->Draw();
    // this one is used to plot the run number on the canvas
    transparent[0] = new TPad("transparent0", "this does not show", 0, 0, 1, 1);
    transparent[0]->SetFillStyle(4000);
    transparent[0]->Draw();
    TC[0]->SetEditable(false);
    TC[0]->SetTopMargin(0.05);
    TC[0]->SetBottomMargin(0.05);
  }
  else if (name == "MvtxMon_General")
  {
    TC[1] = new TCanvas(name.c_str(), "MVTX Monitoring General", -1, 0, xsize, ysize);
    gSystem->ProcessEvents();
    TC[1]->cd();
    Pad[1] = new TPad("mvtxpad1", "who needs this?", 0., 0.02, 0.98, 0.97, 0);  
    Pad[1]->Draw();
    //Pad[7]->Divide(4,2);  
    // this one is used to plot the run number on the canvas
    transparent[0] = new TPad("transparent0", "this does not show", 0, 0, 1, 1);
    transparent[0]->SetFillStyle(4000);
    transparent[0]->Draw();
    TC[1]->SetEditable(false);
    TC[1]->SetTopMargin(0.05);
    TC[1]->SetBottomMargin(0.05);
//Pad[8]->SetRightMargin(0.13);
  }
  else if (name == "MvtxMon_FEE")
  {
    TC[2] = new TCanvas(name.c_str(), "MVTX Monitoring FEE", -1, 0, xsize, ysize);
    gSystem->ProcessEvents();
    Pad[2] = new TPad("mvtxpad2", "who needs this?", 0., 0.02, 0.98, 0.97, 0);  
    Pad[2]->Draw();
    
    // this one is used to plot the run number on the canvas
    transparent[0] = new TPad("transparent0", "this does not show", 0, 0, 1, 1);
    transparent[0]->SetFillStyle(4000);
    transparent[0]->Draw();
    TC[2]->SetEditable(false);
    TC[2]->SetTopMargin(0.05);
    TC[2]->SetBottomMargin(0.05);
  }
  else if (name == "MvtxMon_OCC")
  {
    TC[3] = new TCanvas(name.c_str(), "MVTX Monitoring Occupancy", -1, 0, xsize, ysize);
    gSystem->ProcessEvents();
    Pad[3] = new TPad("mvtxpad3", "who needs this?", 0., 0.02, 0.98, 0.97, 0);  
    Pad[3]->Draw();
    
    // this one is used to plot the run number on the canvas
    transparent[0] = new TPad("transparent0", "this does not show", 0, 0, 1, 1);
    transparent[0]->SetFillStyle(4000);
    transparent[0]->Draw();
    TC[3]->SetEditable(false);
    TC[3]->SetTopMargin(0.05);
    TC[3]->SetBottomMargin(0.05);
  }
  else if (name == "MvtxMon_FHR")
  {
    TC[4] = new TCanvas(name.c_str(), "MVTX Monitoring FHR", -1, 0, xsize, ysize);
    gSystem->ProcessEvents();
    Pad[4] = new TPad("mvtxpad4", "who needs this?", 0., 0.02, 0.98, 0.97, 0);  
    Pad[4]->Draw();
    
    // this one is used to plot the run number on the canvas
    transparent[0] = new TPad("transparent0", "this does not show", 0, 0, 1, 1);
    transparent[0]->SetFillStyle(4000);
    transparent[0]->Draw();
    TC[4]->SetEditable(false);
    TC[4]->SetTopMargin(0.05);
    TC[4]->SetBottomMargin(0.05);
  }
  else if (name == "MvtxMon3")
  {
    TC[5] = new TCanvas(name.c_str(), "MvtxMon3 Example Monitor", xsize / 2, 0, xsize / 2, ysize);
    gSystem->ProcessEvents();
    Pad[5] = new TPad("mvtxpad5", "who needs this?", 0.1, 0.1, 0.9, 0.9, 0);
    //Pad[6] = new TPad("mvtxpad7", "who needs this?", 0.1, 0.05, 0.9, 0.45, 0);
    Pad[5]->Draw();
    transparent[0] = new TPad("transparent0", "this does not show", 0, 0, 1, 1);
    transparent[0]->SetFillStyle(4000);
    transparent[0]->Draw();
    TC[5]->SetEditable(false);
    TC[5]->SetTopMargin(0.05);
    TC[5]->SetBottomMargin(0.05);
    //Pad[6]->Draw();
    // this one is used to plot the run number on the canvas
    //        transparent[2] = new TPad("transparent2", "this does not show", 0, 0, 1, 1);
    //        transparent[2]->SetFillStyle(4000);
    //        transparent[2]->Draw();
    //      TC[2]->SetEditable(0);
  }
  return 0;
}

int MvtxMonDraw::Draw(const std::string &what)
{
  int iret = 0;
  int idraw = 0;
  if (what == "ALL" || what == "L0L"|| what == "L0R"|| what == "L1L"|| what == "L1R"|| what == "L2L"|| what == "L2R")
  {
    iret += DrawHitMap(what);
    idraw++;
  }
  if (what == "ALL" || what == "GENERAL")
  {
    iret += DrawGeneral(what);
    idraw++;
  }
  if (what == "ALL" || what == "FEE")
  {
    iret += DrawFEE(what);
    idraw++;
  }
  if (what == "ALL" || what == "OCC")
  {
    iret += DrawOCC(what);
    idraw++;
  }
  if (what == "ALL" || what == "FHR")
  {
    iret += DrawFHR(what);
    idraw++;
  }
  if (what == "HISTORY")
  {
    iret += DrawHistory(what);
    idraw++;
  }
  if (!idraw)
  {
    std::cout << __PRETTY_FUNCTION__ << " Unimplemented Drawing option: " << what << std::endl;
    iret = -1;
  }
  return iret;
}

int MvtxMonDraw::DrawHitMap(const std::string &what )
{
  OnlMonClient *cl = OnlMonClient::instance();

  const int canvasID = 0;
  const int padID = 0;
  TH3 *mvtxmon_HitMap[NFlx+1] = {nullptr};
  if (!gROOT->FindObject("MvtxMon_HitMap"))
  {
    MakeCanvas("MvtxMon_HitMap");
  }
  TC[canvasID]->SetEditable(true);
  TC[canvasID]->Clear("D");

  int aLs = 0;
  int aLe = 0;
  int aSs = 0;
  int aSe = 0;
  if (what == "L0R" || what == "ALL"){
    aLs = 0; aLe = 1;
    aSs = 0; aSe = 6;
  }
  else if (what == "L0L"){
    aLs = 0; aLe = 1;
    aSs = 6; aSe = 12;
  }
  else if (what == "L1R"){
    aLs = 1; aLe = 2;
    aSs = 0; aSe = 8;
  }
  else if (what == "L1L"){
    aLs = 1; aLe = 2;
    aSs = 8; aSe = 16;
  }
  else if (what == "L2R"){
    aLs = 2; aLe = 3;
    aSs = 0; aSe = 10;
  }
  else if (what == "L2L"){
    aLs = 2; aLe = 3;
    aSs = 10; aSe = 20;
  }
  Pad[padID]->Divide(NCHIP,aSe-aSs/*NSTAVE*/);
  
  int ipad = 0;
  int returnCode = 0;

  for (int iFelix = 0; iFelix <NFlx; iFelix++){
    mvtxmon_HitMap[iFelix] = dynamic_cast<TH3*>(cl->getHisto(Form("MVTXMON_%d",iFelix),Form("MVTXMON_chipHitmapFLX%d", iFelix)));
  }
  MergeServers<TH3*>(mvtxmon_HitMap);
  if(mvtxmon_HitMap[NFlx]){
    mvtxmon_HitMap[NFlx]->GetXaxis()->CenterTitle();
    mvtxmon_HitMap[NFlx]->GetYaxis()->CenterTitle();
    mvtxmon_HitMap[NFlx]->GetYaxis()->SetTitleOffset(1.4);
    mvtxmon_HitMap[NFlx]->GetXaxis()->SetTitleOffset(0.75);
    mvtxmon_HitMap[NFlx]->GetXaxis()->SetTitleSize(0.06);
    mvtxmon_HitMap[NFlx]->GetYaxis()->SetTitleOffset(0.75);
    mvtxmon_HitMap[NFlx]->GetYaxis()->SetTitleSize(0.06);
  }

  for (int aLayer = aLs; aLayer < aLe; aLayer++) {
    for (int aStave = aSs; aStave < aSe; aStave++) {
      for (int iChip = 0; iChip < 9; iChip++) {
        //int stave = aLayer==0?aStave:NStaves[aLayer]+aStave;
	TString prefix = Form("%d%d%d",aLayer,aStave,iChip);
        mvtxmon_HitMap[NFlx]->GetZaxis()->SetRange(((chipmapoffset[aLayer]+aStave)*9+iChip+1),((chipmapoffset[aLayer]+aStave)*9+iChip+1));
        returnCode += PublishHistogram(Pad[padID],ipad*9+iChip+1,mvtxmon_HitMap[NFlx]->Project3D(prefix+"yx"),"colz"); //publish merged one
	gStyle->SetOptStat(0);
      }
      ipad++;
    }
  }

  TC[canvasID]->cd();
  TPaveText *ptchip[9] = {nullptr};
  for(int i = 0; i < 9; i++){
    ptchip[i] = new TPaveText(0.05+0.1033*i,.92,0.05+0.1033*(i+1),0.97,"blNDC");
    ptchip[i]->SetTextSize(0.03);
    ptchip[i]->SetFillColor(4000);
    ptchip[i]->SetLineColor(0);
    ptchip[i]->SetBorderSize(1);
    ptchip[i]->AddText(Form("Chip %d",i));
    ptchip[i]->Draw();
  }

  TPaveText *ptstave[10] = {nullptr};
  for(int i = 0; i < aSe - aSs; i++){
    ptstave[i] = new TPaveText(0.0,0.93- ((i+1) * 0.92/(aSe - aSs)),0.05,0.93- (i * 0.92/(aSe - aSs)),"blNDC");
    ptstave[i]->SetTextSize(0.02);
    ptstave[i]->SetFillColor(4000);
    ptstave[i]->SetLineColor(0);
    ptstave[i]->SetBorderSize(1);
    ptstave[i]->AddText(Form("Stave %d",aSs + i));
    ptstave[i]->Draw();
  }

  TPaveText *ptlayer = {nullptr};
  ptlayer = new TPaveText(0.01,.92,0.06,0.97,"blNDC");
  ptlayer->SetTextSize(0.035);
  ptlayer->SetFillColor(4000);
  ptlayer->SetLineColor(0);
  ptlayer->SetBorderSize(1);
  if (what == "L0R" || what == "L0R"|| what == "ALL") ptlayer->AddText("Layer 0");
  if (what == "L1R" || what == "L1R") ptlayer->AddText("Layer 1");
  if (what == "L2R" || what == "L2R") ptlayer->AddText("Layer 2");
  ptlayer->Draw();

  PublishStatistics(TC[canvasID],cl);

  //TC[canvasID]->Modified();
  TC[canvasID]->Update();
  TC[canvasID]->SetEditable(false);
  return returnCode < 0 ? -1 : 0;
}

int MvtxMonDraw::DrawGeneral(const std::string & /* what */)
{
  OnlMonClient *cl = OnlMonClient::instance();
  const int canvasID = 1;
  const int padID = 1;

  TH2Poly* mvtxmon_LaneStatusOverview[3][NFlx+1] = {nullptr};  
  TH2Poly* mvtxmon_mGeneralOccupancy[NFlx+1] = {nullptr}; 
  TH1D* mvtxmon_mGeneralErrorPlots[NFlx+1]  = {nullptr}; 
  TH2D* mvtxmon_mGeneralErrorFile[NFlx+1]  = {nullptr}; 

  for (int iFelix = 0; iFelix < NFlx; iFelix++){
    mvtxmon_LaneStatusOverview[0][iFelix] = dynamic_cast<TH2Poly*>(cl->getHisto(Form("MVTXMON_%d",iFelix),"MVTXMON_LaneStatus_laneStatusOverviewFlagWARNING"));
    mvtxmon_LaneStatusOverview[1][iFelix] = dynamic_cast<TH2Poly*>(cl->getHisto(Form("MVTXMON_%d",iFelix),"MVTXMON_LaneStatus_laneStatusOverviewFlagERROR"));
    mvtxmon_LaneStatusOverview[2][iFelix] = dynamic_cast<TH2Poly*>(cl->getHisto(Form("MVTXMON_%d",iFelix),"MVTXMON_LaneStatus_laneStatusOverviewFlagFAULT"));
    mvtxmon_mGeneralOccupancy[iFelix] = dynamic_cast<TH2Poly*>(cl->getHisto(Form("MVTXMON_%d",iFelix),"MVTXMON_General_Occupancy"));
    mvtxmon_mGeneralErrorPlots[iFelix] = dynamic_cast<TH1D*>(cl->getHisto(Form("MVTXMON_%d",iFelix),"MVTXMON_General_ErrorPlots"));
    mvtxmon_mGeneralErrorFile[iFelix] = dynamic_cast<TH2D*>(cl->getHisto(Form("MVTXMON_%d",iFelix),"MVTXMON_General_ErrorFile"));
  }
  
  //injecting errors here
  //mvtxmon_LaneStatusOverview[0][1] = mvtxmon_LaneStatusOverview[0][0];
  //injecting errors done

  int bitset = 0;
  int bitsetOR = 0;
  int bitsetAND = 63;
  bitset = MergeServers<TH2Poly*>(mvtxmon_LaneStatusOverview[0]);
  bitsetOR |= bitset; bitsetAND &= bitset;
  bitset = MergeServers<TH2Poly*>(mvtxmon_LaneStatusOverview[1]);
  bitsetOR |= bitset; bitsetAND &= bitset;
  bitset = MergeServers<TH2Poly*>(mvtxmon_LaneStatusOverview[2]);
  bitsetOR |= bitset; bitsetAND &= bitset;
  bitset = MergeServers<TH2Poly*>(mvtxmon_mGeneralOccupancy);
  bitsetOR |= bitset; bitsetAND &= bitset;
  bitset = MergeServers<TH1D*>(mvtxmon_mGeneralErrorPlots);
  bitsetOR |= bitset; bitsetAND &= bitset;
  bitset = MergeServers<TH2D*>(mvtxmon_mGeneralErrorFile);
  bitsetOR |= bitset; bitsetAND &= bitset;
  
  for(int i = 0; i < 3; i++){
    if(mvtxmon_LaneStatusOverview[i][NFlx]){
      mvtxmon_LaneStatusOverview[i][NFlx]->SetStats(0);
      mvtxmon_LaneStatusOverview[i][NFlx]->GetYaxis()->SetTitleOffset(0.6);
      mvtxmon_LaneStatusOverview[i][NFlx]->SetMinimum(0);
      mvtxmon_LaneStatusOverview[i][NFlx]->SetMaximum(1);
      mvtxmon_LaneStatusOverview[i][NFlx]->SetBit(TH1::kIsAverage);
    }
  }
  if(mvtxmon_mGeneralOccupancy[NFlx]) mvtxmon_mGeneralOccupancy[NFlx]->GetYaxis()->SetTitleOffset(0.6);

  if(mvtxmon_mGeneralErrorFile[NFlx]){
    for(int bin = 1; bin < mvtxmon_mGeneralErrorFile[NFlx]->GetNbinsX()+1;bin++ ){
      mvtxmon_mGeneralErrorFile[NFlx]->GetXaxis()->SetBinLabel(bin,Form("%d",(bin-1)%8));
    }
    mvtxmon_mGeneralErrorFile[NFlx]->GetXaxis()->SetTitleSize(0.045);
    mvtxmon_mGeneralErrorFile[NFlx]->GetYaxis()->SetTitleSize(0.045);
  }

  if (!gROOT->FindObject("MvtxMon_General"))
  {
    MakeCanvas("MvtxMon_General");
  }

  TC[canvasID]->SetEditable(true);
  TC[canvasID]->Clear("D");
  
  Pad[padID]->Divide(4,2);
  int returnCode = 0;

  Pad[padID]->cd(4)->SetRightMargin(0.12);
  Pad[padID]->cd(5)->SetLeftMargin(0.16);
  Pad[padID]->cd(6)->SetTopMargin(0.16);

  //injecting errors here
/*  mvtxmon_LaneStatusOverview[1][NFlx]->SetBinContent(5,0.5);
  mvtxmon_LaneStatusOverview[2][NFlx]->SetBinContent(11,0.5);
  mvtxmon_LaneStatusOverview[2][NFlx]->SetBinContent(12,0.5);
  mvtxmon_LaneStatusOverview[2][NFlx]->SetBinContent(13,0.5);
  mvtxmon_LaneStatusOverview[2][NFlx]->SetBinContent(14,0.5);
  mvtxmon_LaneStatusOverview[2][NFlx]->SetBinContent(15,0.5);
  mvtxmon_LaneStatusOverview[2][NFlx]->SetBinContent(16,0.5);
  mvtxmon_LaneStatusOverview[2][NFlx]->SetBinContent(17,0.5);
  mvtxmon_LaneStatusOverview[2][NFlx]->SetBinContent(18,0.5);
  mvtxmon_LaneStatusOverview[2][NFlx]->SetBinContent(19,0.5);
  mvtxmon_mGeneralErrorFile[NFlx]->SetBinContent(17,5,20);*/
  //injecting errors done

  if(mvtxmon_mGeneralErrorPlots[NFlx]){
    mvtxmon_mGeneralErrorPlots[NFlx]->GetXaxis()->SetTitleSize(0.045);
    mvtxmon_mGeneralErrorPlots[NFlx]->GetYaxis()->SetTitleSize(0.045);
  }

  std::vector<MvtxMonDraw::Quality> status;
  status = analyseForError(mvtxmon_LaneStatusOverview[0][NFlx],mvtxmon_LaneStatusOverview[1][NFlx],mvtxmon_LaneStatusOverview[2][NFlx],mvtxmon_mGeneralErrorFile[NFlx]);

  returnCode += PublishHistogram(Pad[padID],1,mvtxmon_LaneStatusOverview[0][NFlx],"lcolz");
  DrawPave(status, 0);
  returnCode += PublishHistogram(Pad[padID],2,mvtxmon_LaneStatusOverview[1][NFlx],"lcolz");
  DrawPave(status, 3);
  returnCode += PublishHistogram(Pad[padID],3,mvtxmon_LaneStatusOverview[2][NFlx],"lcolz");
  DrawPave(status, 6);
  returnCode += PublishHistogram(Pad[padID],4,mvtxmon_mGeneralOccupancy[NFlx]);
  returnCode += PublishHistogram(Pad[padID],4,mvtxmon_mGeneralOccupancy[NFlx],"colz same");
  returnCode += PublishHistogram(Pad[padID],5,mvtxmon_mGeneralErrorPlots[NFlx]);
  returnCode += PublishHistogram(Pad[padID],6,mvtxmon_mGeneralErrorFile[NFlx],"lcol");

  TPaveText *ptt5 = new TPaveText(.1,.85,.9,.95,"blNDC");
  ptt5->SetTextSize(0.04);
  ptt5->SetFillColor(0);
  ptt5->SetLineColor(0);
  ptt5->SetBorderSize(1);
  ptt5->AddText(" FLX 0   FLX 1   FLX 2   FLX 3   FLX 4   FLX 5");
  ptt5->Draw();
 
  Pad[padID]->cd(7);

  TPaveText *bulbRed = new TPaveText(0.1, 0.82, 0.9, 0.92, "NDC");
  bulbRed->SetName("BulbRed");
  formatPaveText(bulbRed, 0.05, kRed, 12, "#color[2]{It's Burning! Call The Expert!!!}");

  TPaveText *bulbYellow = new TPaveText(0.1, 0.82, 0.9, 0.92, "NDC");
  bulbYellow->SetName("BulbYellow");
  formatPaveText(bulbYellow, 0.05, kRed, 12, "#color[808]{Did you break it? Call The Expert!}");

  TPaveText *bulbGreen = new TPaveText(0.1, 0.82, 0.9, 0.92, "NDC");
  bulbGreen->SetName("BulbGreen");
  formatPaveText(bulbGreen, 0.05, kGreen, 12, "#color[418]{Chill, Have Some Coffee}");

  TPaveText *bulb = new TPaveText(0.05, 0.05, 0.95, 0.95, "NDC");
  
  
  TPaveText *ptt4 = new TPaveText(.1,.07,.9,.58,"blNDC");
  ptt4->SetTextSize(0.04);
  ptt4->SetFillColor(0);
  ptt4->SetTextAlign(22);
  ptt4->SetLineColor(0);
  ptt4->SetBorderSize(1); 
  ptt4->AddText("Alarms:");


  if(status.at(9) == Quality::Good && status.at(10) == Quality::Good && status.at(11) == Quality::Good){
    bulb->SetFillColor(kGreen);
    bulb->Draw();
    bulbGreen->Draw("same");
  }

  if(status.at(9) == Quality::Medium || status.at(10) == Quality::Medium || status.at(11) == Quality::Medium){
    if(status.at(9) == Quality::Medium) ptt4->AddText("#color[808]{QA Layer 0 Medium}");
    if(status.at(10) == Quality::Medium) ptt4->AddText("#color[808]{QA Layer 1 Medium}");
    if(status.at(11) == Quality::Medium) ptt4->AddText("#color[808]{QA Layer 2 Medium}");
    bulb->SetFillColor(kYellow);
    bulb->Draw();
    bulbYellow->Draw("same");
  }

  if(status.at(9) == Quality::Bad || status.at(10) == Quality::Bad || status.at(11) == Quality::Bad ||status.at(12) == Quality::Bad||status.at(13) == Quality::Bad||status.at(14) == Quality::Bad||status.at(15) == Quality::Bad||status.at(16) == Quality::Bad||status.at(17) == Quality::Bad ||(bitsetAND & 0x3F) != 0x3F){
    if(status.at(9) == Quality::Bad) ptt4->AddText("#color[2]{QA Layer 0 Bad}");
    if(status.at(10) == Quality::Bad) ptt4->AddText("#color[2]{QA Layer 1 Bad}");
    if(status.at(11) == Quality::Bad) ptt4->AddText("#color[2]{QA Layer 2 Bad}");
    if(status.at(12) == Quality::Bad) ptt4->AddText("#color[2]{Felix 0 Decoder Errors}");
    if(status.at(13) == Quality::Bad) ptt4->AddText("#color[2]{Felix 1 Decoder Errors}");
    if(status.at(14) == Quality::Bad) ptt4->AddText("#color[2]{Felix 2 Decoder Errors}");
    if(status.at(15) == Quality::Bad) ptt4->AddText("#color[2]{Felix 3 Decoder Errors}");
    if(status.at(16) == Quality::Bad) ptt4->AddText("#color[2]{Felix 4 Decoder Errors}");
    if(status.at(17) == Quality::Bad) ptt4->AddText("#color[2]{Felix 5 Decoder Errors}");
    if((bitsetAND & 0x3F) != 0x3F) ptt4->AddText("#color[2]{Some Servers are Offline or in Error}");
    bulb->SetFillColor(kRed);
    bulb->Draw();
    bulbRed->Draw("same");
  }

  

  ptt4->Draw();

  TPaveText *tlegend = new TPaveText(.2,.6,.8,.8,"blNDC");
  tlegend->SetTextSize(0.04);
  tlegend->SetFillColor(0);
  tlegend->SetTextAlign(22);
  tlegend->SetLineColor(0);
  tlegend->SetBorderSize(1);
  tlegend->AddText("Legend:");
  tlegend->AddText("#color[418]{Green = QA OK}");
  tlegend->AddText("#color[808]{Yellow = QA Medium}");
  tlegend->AddText("#color[2]{Red = QA Bad}");
  tlegend->Draw();
 
  Pad[padID]->cd(8);
  TPaveText *pt = new TPaveText(.05,.1,.95,.8);
  pt->AddText("Online Monitoring Server Status");
  for (int iFelix = 0; iFelix < NFlx; iFelix++){
    std::string serverStatus = "Felix " + std::to_string(iFelix);
    if((bitsetOR & (1<<iFelix)) == 1<<iFelix && (bitsetAND & (1<<iFelix)) == 1<<iFelix) serverStatus += " #color[418]{ONLINE} ";
    if((bitsetOR & (1<<iFelix)) == 0 && (bitsetAND & (1<<iFelix)) == 0) serverStatus += " #color[2]{OFFLINE}";
    if((bitsetOR & (1<<iFelix)) == 1<<iFelix && (bitsetAND & (1<<iFelix)) == 0) serverStatus += " #color[2]{ERROR}";
    pt->AddText(serverStatus.c_str());
  }
  pt->Draw();
  PublishStatistics(TC[canvasID],cl);
  TC[canvasID]->SetEditable(false);
  return returnCode < 0 ? -1 : 0;
}

int MvtxMonDraw::DrawFEE(const std::string & /* what */)
{
  OnlMonClient *cl = OnlMonClient::instance();
  const int canvasID = 2;
  const int padID = 2;

  TH2I* mTriggerVsFeeId[NFlx+1] = {nullptr}; 
  TH1I* mTrigger[NFlx+1] ={nullptr}; 
  //TH2I* mLaneInfo[NFlx+1] = {nullptr}; 
  TH2I* mLaneStatus[3][NFlx+1] = {nullptr}; 
  TH2I* mLaneStatusCumulative[3][NFlx+1] = {nullptr}; 
  TH1I* mLaneStatusSummary[3][NFlx+1] = {nullptr}; 
  TH1I* mLaneStatusSummaryIB[NFlx+1] = {nullptr}; 

  for (int iFelix = 0; iFelix <NFlx; iFelix++){
    mTriggerVsFeeId[iFelix] = dynamic_cast<TH2I*>(cl->getHisto(Form("MVTXMON_%d",iFelix),"MVTXMON_FEE_TriggerVsFeeid"));
    mTrigger[iFelix] = dynamic_cast<TH1I*>(cl->getHisto(Form("MVTXMON_%d",iFelix),"MVTXMON_FEE_TriggerFlag"));
    //mLaneInfo[iFelix] = dynamic_cast<TH2I*>(cl->getHisto(Form("MVTXMON_%d",iFelix),"MVTXMON/FEE/LaneInfo"));
    mLaneStatusSummaryIB[iFelix]= dynamic_cast<TH1I*>(cl->getHisto(Form("MVTXMON_%d",iFelix),"MVTXMON_LaneStatusSummary_LaneStatusSummary"));
    for (int i = 0; i < 3; i++) {
      mLaneStatus[i][iFelix] = dynamic_cast<TH2I*>(cl->getHisto(Form("MVTXMON_%d",iFelix),Form("MVTXMON_LaneStatus_laneStatusFlag%s", mLaneStatusFlag[i].c_str())));
      mLaneStatusCumulative[i][iFelix] = dynamic_cast<TH2I*>(cl->getHisto(Form("MVTXMON_%d",iFelix),Form("MVTXMON_LaneStatus_laneStatusFlagCumulative%s", mLaneStatusFlag[i].c_str())));
      mLaneStatusSummary[i][iFelix] = dynamic_cast<TH1I*>(cl->getHisto(Form("MVTXMON_%d",iFelix),Form("MVTXMON_LaneStatusSummary_LaneStatusSummaryL%i", i)));
    }
  }

  for (int i = 0; i < 3; i++) {
    MergeServers<TH2I*>(mLaneStatus[i]);
    MergeServers<TH2I*>(mLaneStatusCumulative[i]);
    MergeServers<TH1I*>(mLaneStatusSummary[i]);
  }
  MergeServers<TH2I*>(mTriggerVsFeeId);
  MergeServers<TH1I*>(mTrigger);
  //MergeServers<TH2I*>(mLaneInfo);
  MergeServers<TH1I*>(mLaneStatusSummaryIB);

  if (!gROOT->FindObject("MvtxMon_FEE"))
  {
    MakeCanvas("MvtxMon_FEE");
  }

  TC[canvasID]->SetEditable(true);
  TC[canvasID]->Clear("D");
  Pad[padID]->Divide(4,3);

  for (int i = 0; i < 3; i++) {
    if(mLaneStatus[i][NFlx]){
      mLaneStatus[i][NFlx]->GetXaxis()->SetTitleSize(0.055);
      mLaneStatus[i][NFlx]->GetXaxis()->SetTitleOffset(0.85);
      mLaneStatus[i][NFlx]->GetYaxis()->SetTitleSize(0.06);
      mLaneStatus[i][NFlx]->GetXaxis()->SetLabelSize(0.05);
      mLaneStatus[i][NFlx]->GetYaxis()->SetLabelSize(0.05);
    }
    if(mLaneStatusCumulative[i][NFlx]){
      mLaneStatusCumulative[i][NFlx]->GetXaxis()->SetTitleSize(0.055);
      mLaneStatusCumulative[i][NFlx]->GetXaxis()->SetTitleOffset(0.85);
      mLaneStatusCumulative[i][NFlx]->GetYaxis()->SetTitleSize(0.06);
      mLaneStatusCumulative[i][NFlx]->GetXaxis()->SetLabelSize(0.05);
      mLaneStatusCumulative[i][NFlx]->GetYaxis()->SetLabelSize(0.05);
    }
    if(mLaneStatusSummary[i][NFlx]){
      mLaneStatusSummary[i][NFlx]->GetXaxis()->SetLabelSize(0.07);
      mLaneStatusSummary[i][NFlx]->GetYaxis()->SetTitleSize(0.06);
      mLaneStatusSummary[i][NFlx]->GetYaxis()->SetLabelSize(0.05);
      mLaneStatusSummary[i][NFlx]->GetYaxis()->SetTitleOffset(0.9);
    }
  }
  if(mLaneStatusSummaryIB[NFlx]){
    mLaneStatusSummaryIB[NFlx]->GetXaxis()->SetLabelSize(0.07);
    mLaneStatusSummaryIB[NFlx]->GetYaxis()->SetTitleSize(0.06);
    mLaneStatusSummaryIB[NFlx]->GetYaxis()->SetLabelSize(0.05);
    mLaneStatusSummaryIB[NFlx]->GetYaxis()->SetTitleOffset(0.9);
  }
  if(mTriggerVsFeeId[NFlx]){
    mTriggerVsFeeId[NFlx]->GetXaxis()->SetTitleSize(0.055);
    mTriggerVsFeeId[NFlx]->GetXaxis()->SetTitleOffset(0.75);
    mTriggerVsFeeId[NFlx]->GetYaxis()->SetTitleOffset(1.2);
    mTriggerVsFeeId[NFlx]->GetYaxis()->SetTitleSize(0.06);
    mTriggerVsFeeId[NFlx]->GetYaxis()->SetLabelSize(0.06);
    mTriggerVsFeeId[NFlx]->GetXaxis()->SetLabelSize(0.05);
  }
  if(mTrigger[NFlx]){
    mTrigger[NFlx]->GetXaxis()->SetTitleSize(0.055);
    mTrigger[NFlx]->GetXaxis()->SetTitleOffset(0.85);
    mTrigger[NFlx]->GetYaxis()->SetTitleSize(0.06);
    mTrigger[NFlx]->GetYaxis()->SetLabelSize(0.06);
    mTrigger[NFlx]->GetXaxis()->SetLabelSize(0.05);
    mTrigger[NFlx]->GetYaxis()->SetTitleOffset(0.85);
  }

  TPaveText *tlayer[3] = {nullptr};
  
  for (int i = 0; i < 3 ; i++){
    double shift[3] = {0,0.25,0.55};
    tlayer[i] = new TPaveText(.14+shift[i],.87,.24+shift[i],.93,"blNDC");
    tlayer[i]->SetTextSize(0.05);
    tlayer[i]->SetFillColor(0);
    tlayer[i]->SetTextAlign(22);
    tlayer[i]->SetLineColor(0);
    tlayer[i]->SetBorderSize(1);
    tlayer[i]->AddText(Form("Layer %d", i));
  }

  Pad[padID]->cd(3)->SetTopMargin(0.16);
  Pad[padID]->cd(4)->SetTopMargin(0.16);
  Pad[padID]->cd(4)->SetTopMargin(0.145);
  Pad[padID]->cd(7)->SetTopMargin(0.16);
  Pad[padID]->cd(8)->SetTopMargin(0.16);
  Pad[padID]->cd(11)->SetTopMargin(0.16);
  Pad[padID]->cd(12)->SetTopMargin(0.16);

  int returnCode = 0;
   Pad[padID]->cd(1)->SetLeftMargin(0.16);
  returnCode += PublishHistogram(Pad[padID],1,mTriggerVsFeeId[NFlx],"lcol");
  returnCode += PublishHistogram(Pad[padID],5,mTrigger[NFlx]);
  //returnCode += PublishHistogram(Pad[9],3,mLaneInfo[NFlx]);
  returnCode += PublishHistogram(Pad[padID],3,mLaneStatus[0][NFlx]);
  for(int i = 0;i<3;i++)tlayer[i]->Draw();
  returnCode += PublishHistogram(Pad[padID],7,mLaneStatus[1][NFlx]);
  for(int i = 0;i<3;i++)tlayer[i]->Draw();
  returnCode += PublishHistogram(Pad[padID],11,mLaneStatus[2][NFlx]);
  for(int i = 0;i<3;i++)tlayer[i]->Draw();
  returnCode += PublishHistogram(Pad[padID],4,mLaneStatusCumulative[0][NFlx]);
  for(int i = 0;i<3;i++)tlayer[i]->Draw();
  returnCode += PublishHistogram(Pad[padID],8,mLaneStatusCumulative[1][NFlx]);
  for(int i = 0;i<3;i++)tlayer[i]->Draw();
  returnCode += PublishHistogram(Pad[padID],12,mLaneStatusCumulative[2][NFlx]);
  for(int i = 0;i<3;i++)tlayer[i]->Draw();
  returnCode += PublishHistogram(Pad[padID],2,mLaneStatusSummary[0][NFlx]);
  returnCode += PublishHistogram(Pad[padID],6,mLaneStatusSummary[1][NFlx]);
  returnCode += PublishHistogram(Pad[padID],10,mLaneStatusSummary[2][NFlx]);
  returnCode += PublishHistogram(Pad[padID],9,mLaneStatusSummaryIB[NFlx]);
  
  PublishStatistics(TC[canvasID],cl);
  TC[canvasID]->SetEditable(false);
  return returnCode < 0 ? -1 : 0;
}

int MvtxMonDraw::DrawOCC(const std::string & /* what */)
{
  OnlMonClient *cl = OnlMonClient::instance();
  const int canvasID = 3;
  const int padID = 3;

  TH1D* hOccupancyPlot[3][NFlx+1] = {nullptr};
  TH2D* hChipStaveOccupancy[3][NFlx+1] = {nullptr};
  TH1D *mvtxmon_ChipStave1D[NFlx+1] = {nullptr};
  TH1D *mvtxmon_ChipFiredHis[NFlx+1] = {nullptr};
  TH1I *hChipStrobes[NFlx+1] = {nullptr};
  TH1I *hChipL1[NFlx+1] = {nullptr};

  for (int aLayer = 0; aLayer < 3; aLayer++) {
    for (int iFelix = 0; iFelix <NFlx; iFelix++){
      hOccupancyPlot[aLayer][iFelix] = dynamic_cast<TH1D*>(cl->getHisto(Form("MVTXMON_%d",iFelix),Form("MVTXMON_Occupancy_Layer%dOccupancy", aLayer)));
      hChipStaveOccupancy[aLayer][iFelix] = dynamic_cast<TH2D*>(cl->getHisto(Form("MVTXMON_%d",iFelix),Form("MVTXMON_Occupancy_Layer%d_Layer%dChipStave", aLayer, aLayer)));
    }
  }
 
  for (int iFelix = 0; iFelix <NFlx; iFelix++){
    mvtxmon_ChipStave1D[iFelix] = dynamic_cast<TH1D*>(cl->getHisto(Form("MVTXMON_%d",iFelix),"MVTXMON_RawData_ChipStave1D"));
    mvtxmon_ChipFiredHis[iFelix] = dynamic_cast<TH1D*>(cl->getHisto(Form("MVTXMON_%d",iFelix),"MVTXMON_RawData_ChipFiredHis"));
    hChipStrobes[iFelix] = dynamic_cast<TH1I*>(cl->getHisto(Form("MVTXMON_%d",iFelix),"hChipStrobes"));
    hChipL1[iFelix] = dynamic_cast<TH1I*>(cl->getHisto(Form("MVTXMON_%d",iFelix),"hChipL1"));
  }

  for (int i = 0; i < 3; i++) {
    MergeServers<TH1D*>(hOccupancyPlot[i]);
    MergeServers<TH2D*>(hChipStaveOccupancy[i]); 
  }

  TH2D *hChipStaveOccupancy_low[3] =  {nullptr};

  for (int i = 0; i < 3; i++) {
    if(hOccupancyPlot[i][NFlx]){
      hOccupancyPlot[i][NFlx]->GetXaxis()->SetLabelSize(0.05);  
      hOccupancyPlot[i][NFlx]->GetXaxis()->SetTitleSize(0.05);
      hOccupancyPlot[i][NFlx]->GetXaxis()->SetTitleOffset(0.85);
      hOccupancyPlot[i][NFlx]->GetYaxis()->SetLabelSize(0.045);
      hOccupancyPlot[i][NFlx]->GetYaxis()->SetTitleSize(0.05);
      hOccupancyPlot[i][NFlx]->GetYaxis()->SetTitleOffset(1);
      hChipStaveOccupancy[i][NFlx]->SetMinimum(0);
    }  
    if(hChipStaveOccupancy[i][NFlx]){
      hChipStaveOccupancy[i][NFlx]->GetXaxis()->SetLabelSize(0.05);  
      hChipStaveOccupancy[i][NFlx]->GetXaxis()->SetTitleSize(0.05);
      hChipStaveOccupancy[i][NFlx]->GetXaxis()->SetTitleOffset(0.85);
      hChipStaveOccupancy[i][NFlx]->GetYaxis()->SetLabelSize(0.06);
      hChipStaveOccupancy[i][NFlx]->GetYaxis()->SetTitleSize(0.05);
      hChipStaveOccupancy[i][NFlx]->GetYaxis()->SetTitleOffset(0.7);

      hChipStaveOccupancy_low[i] = (TH2D*)(hChipStaveOccupancy[i][NFlx])->Clone(Form("hChipStaveOccupancy_low_%d",i));
      TString tmp = "Chip Occupancy < 1\%";
      hChipStaveOccupancy_low[i]->SetTitle(tmp + Form(" of average occupancy in layer %d",i));
    }
  }

   

  MergeServers<TH1D*>(mvtxmon_ChipStave1D);
  MergeServers<TH1D*>(mvtxmon_ChipFiredHis);
  MergeServers<TH1I*>(hChipStrobes); 
  MergeServers<TH1I*>(hChipL1); 

  if(mvtxmon_ChipStave1D[NFlx]){
    mvtxmon_ChipStave1D[NFlx]->GetXaxis()->SetLabelSize(0.05);  
    mvtxmon_ChipStave1D[NFlx]->GetXaxis()->SetTitleSize(0.055);
    mvtxmon_ChipStave1D[NFlx]->GetXaxis()->SetTitleOffset(0.75);
    mvtxmon_ChipStave1D[NFlx]->GetYaxis()->SetLabelSize(0.045);
    mvtxmon_ChipStave1D[NFlx]->GetYaxis()->SetTitleSize(0.06);
    mvtxmon_ChipStave1D[NFlx]->GetYaxis()->SetTitleOffset(1.2);
  }

  if(mvtxmon_ChipFiredHis[NFlx]){
    mvtxmon_ChipFiredHis[NFlx]->GetXaxis()->SetLabelSize(0.05);  
    mvtxmon_ChipFiredHis[NFlx]->GetXaxis()->SetTitleSize(0.055);
    mvtxmon_ChipFiredHis[NFlx]->GetXaxis()->SetTitleOffset(0.75);
    mvtxmon_ChipFiredHis[NFlx]->GetYaxis()->SetLabelSize(0.06);
    mvtxmon_ChipFiredHis[NFlx]->GetYaxis()->SetTitleSize(0.06);
    mvtxmon_ChipFiredHis[NFlx]->GetYaxis()->SetTitleOffset(1.2);
  }

 double avr_occ[3] = {0};

  for (int iLayer = 0; iLayer < 3; iLayer++) {
    for (int iStave = 0; iStave < NStaves[iLayer]; iStave++) {
      for (int iChip = 0; iChip < 9; iChip++) {      
        double occ = hChipStaveOccupancy[iLayer][NFlx]->GetBinContent(iChip+1,iStave+1)/(hChipStrobes[NFlx]->GetBinContent((StaveBoundary[iLayer]+iStave)*9+iChip+1)*1024*512);
        if(occ > 10e-50)hChipStaveOccupancy[iLayer][NFlx]->SetBinContent(iChip+1,iStave+1,occ);
        if(occ > 10e-50)avr_occ[iLayer]+=occ;
      }
    }
  }

for (int iLayer = 0; iLayer < 3; iLayer++) {
  avr_occ[iLayer]/=(NStaves[iLayer]*9);
    for (int iStave = 0; iStave < NStaves[iLayer]; iStave++) {
      for (int iChip = 0; iChip < 9; iChip++) {      
        double occ = hChipStaveOccupancy[iLayer][NFlx]->GetBinContent(iChip+1,iStave+1);
        if(occ <0.01* avr_occ[iLayer]) hChipStaveOccupancy_low[iLayer]->SetBinContent(iChip+1,iStave+1,1);
	else hChipStaveOccupancy_low[iLayer]->SetBinContent(iChip+1,iStave+1,0);
      }
    }
  }

  

  if (!gROOT->FindObject("MvtxMon_OCC"))
  {
    MakeCanvas("MvtxMon_OCC");


  }

  TC[canvasID]->SetEditable(true);
  TC[canvasID]->Clear("D");
  Pad[padID]->Divide(4,3);
  Pad[padID]->cd(1)->SetLogy();
  Pad[padID]->cd(2)->SetLogy();
  Pad[padID]->cd(3)->SetLogy();


  hChipStaveOccupancy[0][NFlx]->GetZaxis()->SetTitle("");
  hChipStaveOccupancy[1][NFlx]->GetZaxis()->SetTitle("");
  hChipStaveOccupancy[2][NFlx]->GetZaxis()->SetTitle("");

  hChipStrobes[NFlx]->SetTitle("Chip Strobes (L1 triggers) vs Chip*Stave");
  hChipStrobes[NFlx]->GetYaxis()->SetTitle("Number of strobes");


  int returnCode = 0;
  returnCode += PublishHistogram(Pad[padID],1,hOccupancyPlot[0][NFlx]);
  returnCode += PublishHistogram(Pad[padID],2,hOccupancyPlot[1][NFlx]);
  returnCode += PublishHistogram(Pad[padID],3,hOccupancyPlot[2][NFlx]);
  returnCode += PublishHistogram(Pad[padID],5,hChipStaveOccupancy[0][NFlx],"colz");
  returnCode += PublishHistogram(Pad[padID],6,hChipStaveOccupancy[1][NFlx],"colz");
  returnCode += PublishHistogram(Pad[padID],7,hChipStaveOccupancy[2][NFlx],"colz");
  returnCode += PublishHistogram(Pad[padID],9,hChipStaveOccupancy_low[0],"col");
  returnCode += PublishHistogram(Pad[padID],10,hChipStaveOccupancy_low[1],"col");
  returnCode += PublishHistogram(Pad[padID],11,hChipStaveOccupancy_low[2],"col");
  returnCode += PublishHistogram(Pad[padID],4,mvtxmon_ChipStave1D[NFlx]);
  returnCode += PublishHistogram(Pad[padID],8,mvtxmon_ChipFiredHis[NFlx]);
  returnCode += PublishHistogram(Pad[padID],12,hChipStrobes[NFlx]);

//for( int i = 1; i<= hChipL1[NFlx]->GetNbinsX();i++) std::cout<<hChipL1[NFlx]->GetBinContent(i)<<" ";

   Float_t rightmax = 2*hChipL1[NFlx]->GetMaximum();
   Float_t scale = hChipStrobes[NFlx]->GetMaximum()/rightmax;
   hChipL1[NFlx]->SetLineColor(kRed);
   hChipL1[NFlx]->Scale(scale);
   //std::cout<<rightmax<<" "<<hChipStrobes[NFlx]->GetMaximum()<<" "<<scale<<std::endl;
   returnCode += PublishHistogram(Pad[padID],12,hChipL1[NFlx],"same hist");

   TGaxis *axis = new TGaxis(48*9,hChipStrobes[NFlx]->GetMinimum(), 48*9, hChipStrobes[NFlx]->GetMaximum(),0,rightmax,510,"+L");
   axis->SetLineColor(kRed);
   axis->SetLabelColor(kRed);
   axis->SetLabelFont(42);
   axis->SetTitleFont(42);
   axis->SetTitleColor(kRed);
   axis->SetTitle("Number of L1 triggers");
   axis->Draw();

  TPaveText *pt = new TPaveText(6,12,8.5,13);
  pt->SetBorderSize(0);
  pt->SetFillStyle(0);
  TText *t1 = pt->AddText("Not saved in .root");
  t1->SetTextColor(kRed);
  Pad[padID]->cd(9);
  pt->Draw();
  Pad[padID]->cd(10);
  pt->Draw();
  Pad[padID]->cd(11);
  pt->Draw();

  PublishStatistics(TC[canvasID],cl);
  TC[canvasID]->SetEditable(false);
  return returnCode < 0 ? -1 : 0;
}

int MvtxMonDraw::DrawFHR(const std::string & /* what */)
{
  OnlMonClient *cl = OnlMonClient::instance();
  const int canvasID = 4;
  const int padID = 4;

  TH2D* mDeadChipPos[3][NFlx+1] = {nullptr};
  TH2D* mAliveChipPos[3][NFlx+1] = {nullptr};
  //TH2D* mChipStaveOccupancy[3][NFlx];
  TH1D* mOccupancyPlot[3][NFlx+1] = {nullptr};
  
  TH2I* mErrorVsFeeid[NFlx+1] = {nullptr}; 
  TH2Poly* mGeneralOccupancy[NFlx+1] = {nullptr}; 
  TH2Poly* mGeneralNoisyPixel[NFlx+1] = {nullptr}; 
  TH2D* mTotalDeadChipPos[NFlx+1] = {nullptr}; 
  TH2D* mTotalAliveChipPos[NFlx+1] = {nullptr}; 
  TH1D *mvtxmon_EvtHitChip[NFlx+1] = {nullptr}; 
  TH1D *mvtxmon_EvtHitDis[NFlx+1] = {nullptr}; 

  TH1I *mRCDAQevt[NFlx+1] = {nullptr};

  TH2D* mChipStaveNoisy[3][NFlx+1]; 

  int nFLX[3]= {0};

  for (int iFelix = 0; iFelix <NFlx; iFelix++){
    mErrorVsFeeid[iFelix] =  dynamic_cast<TH2I*>(cl->getHisto(Form("MVTXMON_%d",iFelix),"MVTXMON_General_ErrorVsFeeid"));
    mGeneralOccupancy[iFelix] = dynamic_cast<TH2Poly*>(cl->getHisto(Form("MVTXMON_%d",iFelix),"MVTXMON_General_Occupancy"));
    mGeneralNoisyPixel[iFelix] = dynamic_cast<TH2Poly*>(cl->getHisto(Form("MVTXMON_%d",iFelix),"MVTXMON_General_Noisy_Pixel"));
    mTotalDeadChipPos[iFelix] = dynamic_cast<TH2D*>(cl->getHisto(Form("MVTXMON_%d",iFelix),"MVTXMON_Occupancy_TotalDeadChipPos"));
    mTotalAliveChipPos[iFelix] = dynamic_cast<TH2D*>(cl->getHisto(Form("MVTXMON_%d",iFelix),"MVTXMON_Occupancy_TotalAliveChipPos"));
    mvtxmon_EvtHitChip[iFelix] = dynamic_cast<TH1D*>(cl->getHisto(Form("MVTXMON_%d",iFelix),"MVTXMON_RawData_EvtHitChip"));
    mvtxmon_EvtHitDis[iFelix] = dynamic_cast<TH1D*>(cl->getHisto(Form("MVTXMON_%d",iFelix),"MVTXMON_RawData_EvtHitDis"));
    mRCDAQevt[iFelix] = dynamic_cast<TH1I*>(cl->getHisto(Form("MVTXMON_%d",iFelix),"RCDAQ_evt"));
    for (int mLayer = 0; mLayer < 3; mLayer++) {
      mDeadChipPos[mLayer][iFelix] =  dynamic_cast<TH2D*>(cl->getHisto(Form("MVTXMON_%d",iFelix),Form("MVTXMON_Occupancy_Layer%d_Layer%dDeadChipPos", mLayer, mLayer)));
      if(mDeadChipPos[mLayer][iFelix]) nFLX[mLayer]++;
      mAliveChipPos[mLayer][iFelix] =  dynamic_cast<TH2D*>(cl->getHisto(Form("MVTXMON_%d",iFelix),Form("MVTXMON_Occupancy_Layer%d_Layer%dAliveChipPos", mLayer, mLayer)));
     // mChipStaveOccupancy[mLayer][iFelix] =  dynamic_cast<TH2D*>(cl->getHisto(Form("MVTXMON/Occupancy/Layer%d/Layer%dChipStaveC", mLayer, mLayer)));
      mOccupancyPlot[mLayer][iFelix] =  dynamic_cast<TH1D*>(cl->getHisto(Form("MVTXMON_%d",iFelix),Form("MVTXMON_Occupancy_Layer%dOccupancy_LOG", mLayer)));

      mChipStaveNoisy[mLayer][iFelix] =  dynamic_cast<TH2D*>(cl->getHisto(Form("MVTXMON_%d",iFelix),Form("MVTXMON_Noisy_Layer%d_ChipStave", mLayer)));
    }
  }

  for (int iFelix = 0; iFelix < NFlx; iFelix++) {
    for (int mLayer = 0; mLayer < 3; mLayer++) {
	if(mAliveChipPos[mLayer][iFelix] && mRCDAQevt[iFelix]) mAliveChipPos[mLayer][iFelix]->Scale(1./mRCDAQevt[iFelix]->Integral());       

     if(mChipStaveNoisy[mLayer][iFelix] && mRCDAQevt[iFelix]) mChipStaveNoisy[mLayer][iFelix]->Scale(1./mRCDAQevt[iFelix]->Integral());       
    }
    if(mGeneralNoisyPixel[iFelix] && mRCDAQevt[iFelix]) mGeneralNoisyPixel[iFelix]->Scale(1./mRCDAQevt[iFelix]->Integral());
  }

  for (int mLayer = 0; mLayer < 3; mLayer++) {
    MergeServers<TH2D*>(mDeadChipPos[mLayer]);
    MergeServers<TH2D*>(mAliveChipPos[mLayer]);
    MergeServers<TH1D*>(mOccupancyPlot[mLayer]);
    MergeServers<TH2D*>(mChipStaveNoisy[mLayer]);
    if(mDeadChipPos[mLayer][NFlx]){
      mDeadChipPos[mLayer][NFlx]->SetMinimum(0);
      mDeadChipPos[mLayer][NFlx]->SetMaximum(1);
    }
    /*if(mAliveChipPos[mLayer][NFlx]){
      mAliveChipPos[mLayer][NFlx]->SetMinimum(0);
      mAliveChipPos[mLayer][NFlx]->SetMaximum(1);
    }*/
  }



  MergeServers<TH2I*>(mErrorVsFeeid);
  MergeServers<TH2Poly*>(mGeneralOccupancy);
  MergeServers<TH2Poly*>(mGeneralNoisyPixel);
  MergeServers<TH2D*>(mTotalDeadChipPos);
  MergeServers<TH2D*>(mTotalAliveChipPos);
  MergeServers<TH1D*>(mvtxmon_EvtHitChip);
  MergeServers<TH1D*>(mvtxmon_EvtHitDis);
  MergeServers<TH1I*>(mRCDAQevt);

  if(mTotalDeadChipPos[NFlx]){
    mTotalDeadChipPos[NFlx]->SetMinimum(0);
    mTotalDeadChipPos[NFlx]->SetMaximum(1);
  }
 /* if(mTotalAliveChipPos[NFlx]){
    mTotalAliveChipPos[NFlx]->SetMinimum(0);
    mTotalAliveChipPos[NFlx]->SetMaximum(1);
  }*/

   for (int mLayer = 0; mLayer < 3; mLayer++) {
    if(mDeadChipPos[mLayer][NFlx]){
        for (int binx = 0; binx < mDeadChipPos[mLayer][NFlx]->GetNbinsX(); binx++) {
       for (int biny = 0; biny < mDeadChipPos[mLayer][NFlx]->GetNbinsY(); biny++) {

      mDeadChipPos[mLayer][NFlx]->SetBinContent(binx+1,biny+1,mDeadChipPos[mLayer][NFlx]->GetBinContent(binx+1,biny+1)-nFLX[mLayer]+1);
}
}
    }     
   }



  if (!gROOT->FindObject("MvtxMon_FHR"))
  {
    MakeCanvas("MvtxMon_FHR");
  }

  TC[canvasID]->SetEditable(true);
  TC[canvasID]->Clear("D");
  Pad[padID]->Divide(5,3);

  int returnCode = 0;
 // returnCode += PublishHistogram(Pad[padID],1,mErrorVsFeeid[0]);
 // returnCode += PublishHistogram(Pad[padID],6,mGeneralOccupancy[NFlx],"COLZ");
  returnCode += PublishHistogram(Pad[padID],5,mGeneralNoisyPixel[NFlx],"COLZ");
  returnCode += PublishHistogram(Pad[padID],1,mDeadChipPos[0][NFlx],"COL");
  returnCode += PublishHistogram(Pad[padID],6,mDeadChipPos[1][NFlx],"COL");
  returnCode += PublishHistogram(Pad[padID],11,mDeadChipPos[2][NFlx],"COL");
  returnCode += PublishHistogram(Pad[padID],2,mAliveChipPos[0][NFlx],"COLZ");
  returnCode += PublishHistogram(Pad[padID],7,mAliveChipPos[1][NFlx],"COLZ");
  returnCode += PublishHistogram(Pad[padID],12,mAliveChipPos[2][NFlx],"COLZ");
 // returnCode += PublishHistogram<TH2D*>(TC[canvasID],10,mChipStaveOccupancy[0][0]);
 // returnCode += PublishHistogram<TH2D*>(TC[canvasID],11,mChipStaveOccupancy[1][0]);
//  returnCode += PublishHistogram<TH2D*>(TC[canvasID],12,mChipStaveOccupancy[2][0]);
  returnCode += PublishHistogram(Pad[padID],3,mChipStaveNoisy[0][NFlx],"COLZ");
  returnCode += PublishHistogram(Pad[padID],8,mChipStaveNoisy[1][NFlx],"COLZ");
  returnCode += PublishHistogram(Pad[padID],13,mChipStaveNoisy[2][NFlx],"COLZ");
  returnCode += PublishHistogram(Pad[padID],4,mOccupancyPlot[0][NFlx]);
  returnCode += PublishHistogram(Pad[padID],9,mOccupancyPlot[1][NFlx]);
  returnCode += PublishHistogram(Pad[padID],14,mOccupancyPlot[2][NFlx]);  

  PublishHistogram(Pad[padID],10,mvtxmon_EvtHitChip[NFlx]);
  PublishHistogram(Pad[padID],15,mvtxmon_EvtHitDis[NFlx]);

  //returnCode += PublishHistogram(Pad[11],16,mTotalDeadChipPos[NFlx],"COL");
  //returnCode += PublishHistogram(Pad[11],17,mTotalAliveChipPos[NFlx],"COL");

  
  PublishStatistics(TC[canvasID],cl);
  TC[canvasID]->SetEditable(false);
  return returnCode < 0 ? -1 : 0;
}

int MvtxMonDraw::SavePlot(const std::string &what, const std::string &type)
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

int MvtxMonDraw::MakeHtml(const std::string &what)
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
  //      << "</TITLE></HEAD>" << std::endl;
  // out2 << "<P>Some status output would go here." << std::endl;
  // out2.close();
  // cl->SaveLogFile(*this);
  return 0;
}

int MvtxMonDraw::DrawHistory(const std::string & /* what */)
{

  const int canvasID = 5;
  const int padID = 5;

  int iret[NFlx] = {0,0,0,0,0,0};
  // you need to provide the following vectors
  // which are filled from the db
  std::vector<float> var[NFlx];
  std::vector<float> varerr[NFlx];
  std::vector<time_t> timestamp[NFlx];
  std::vector<int> runnumber[NFlx];
  std::string varname = "n_events";
  // this sets the time range from whihc values should be returned
  time_t begin = 0;            // begin of time (1.1.1970)
  time_t end = time(nullptr);  // current time (right NOW)

  for (int iFelix = 0; iFelix <NFlx; iFelix++){
    iret[iFelix] = dbvars[iFelix]->GetVar(begin, end, varname, timestamp[iFelix], runnumber[iFelix], var[iFelix], varerr[iFelix]);
    if (iret[iFelix]){
      std::cout << __PRETTY_FUNCTION__ << " Error in db access: FELIX "<<iFelix << std::endl;
    //return iret;
    }
  }


  if (!gROOT->FindObject("MvtxMon3"))
  {
    MakeCanvas("MvtxMon3");
  }
  TC[canvasID]->SetEditable(true);
  TC[canvasID]->Clear("D");
  Pad[padID]->Divide(3,2);

  for (int iFelix = 0; iFelix <NFlx; iFelix++){
    // timestamps come sorted in ascending order
    float *x = new float[var[iFelix].size()];
    float *y = new float[var[iFelix].size()];
    float *ex = new float[var[iFelix].size()];
    float *ey = new float[var[iFelix].size()];
    int n = var[iFelix].size();
    for (unsigned int i = 0; i < var[iFelix].size(); i++)
    {
      //       std::cout << "timestamp: " << ctime(&timestamp[i])
      // 	   << ", run: " << runnumber[i]
      // 	   << ", var: " << var[i]
      // 	   << ", varerr: " << varerr[i]
      // 	   << std::endl;
      x[i] = timestamp[iFelix][i] - TimeOffsetTicks;
      y[i] = var[iFelix][i];
      ex[i] = 0;
      ey[i] = varerr[iFelix][i];
    }
    Pad[padID]->cd(iFelix+1);
    if (gr[iFelix])
    {
      delete gr[iFelix];
    }
    gr[iFelix] = new TGraphErrors(n, x, y, ex, ey);
    gr[iFelix]->SetMarkerColor(4);
    gr[iFelix]->SetMarkerStyle(21);
    gr[iFelix]->Draw("ALP");
    gr[iFelix]->GetXaxis()->SetTimeDisplay(1);
    gr[iFelix]->GetXaxis()->SetLabelSize(0.03);
    // the x axis labeling looks like crap
    // please help me with this, the SetNdivisions
    // don't do the trick
    gr[iFelix]->GetXaxis()->SetNdivisions(-1006);
    gr[iFelix]->GetXaxis()->SetTimeOffset(TimeOffsetTicks);
    gr[iFelix]->GetXaxis()->SetTimeFormat("%Y/%m/%d %H:%M");
    delete[] x;
    delete[] y;
    delete[] ex;
    delete[] ey;
  }

 /* varname = "mvtxmoncount";
  iret = dbvars->GetVar(begin, end, varname, timestamp, runnumber, var, varerr);
  if (iret)
  {
    std::cout << __PRETTY_FUNCTION__ << " Error in db access" << std::endl;
    return iret;
  }
  x = new float[var.size()];
  y = new float[var.size()];
  ex = new float[var.size()];
  ey = new float[var.size()];
  n = var.size();
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
  Pad[6]->cd();
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
  delete[] x;
  delete[] y;
  delete[] ex;
  delete[] ey;*/

  TC[canvasID]->Update();
  TC[canvasID]->SetEditable(false);
  return 0;
}

//template <typename T>
int MvtxMonDraw::PublishHistogram(TCanvas *c, int pad, TH1 *h, const char* opt ){
  if(c && pad !=0){
    c->cd(pad);
    //std::cout<<"ups"<<std::endl;
  }
  if(h){  
    h->DrawCopy(opt);
    return 0;
  }
  else{
    //DrawDeadServer(transparent[0]);
    return -1;
  }
}

//template <typename T>
int MvtxMonDraw::PublishHistogram(TPad *p, int pad, TH1 *h, const char* opt ){
  if(p && pad !=0){
    p->cd(pad);
    TCanvas *c = nullptr;
    return PublishHistogram(c,0,h,opt);
  }
  else{
    return -1;
  }
}

//template <typename T>
int MvtxMonDraw::PublishHistogram(TPad *p, TH1 *h, const char* opt ){
  if(p){
    p->cd();
    TCanvas *c = nullptr;
    return PublishHistogram(c,0,h,opt);
  }

  else{
    return -1;
  }
}

template <typename T>
int MvtxMonDraw::MergeServers(T *h){
  bool cloned = false;
  int bitset = 0;
  for (int iFelix = 0; iFelix < NFlx; iFelix++){
    if (cloned == false){
      if(h[iFelix]){
        h[NFlx] = dynamic_cast<T>(h[iFelix]->Clone());
        bitset |= (1U << (iFelix));
        cloned = true;
      }
      else continue;
    }
    else {
      if(h[iFelix]){
        h[NFlx]->Add(h[iFelix],1.);
        bitset |= (1U << (iFelix));
      }
    }
  }
  return bitset;
}

void MvtxMonDraw::PublishStatistics(TCanvas *c,OnlMonClient *cl){
  c->cd();
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
  c->Update();
  c->Show();
}

void MvtxMonDraw::formatPaveText(TPaveText* aPT, float aTextSize, Color_t aTextColor, short aTextAlign, const char* aText)
{
  aPT->SetTextSize(aTextSize);
  aPT->SetTextAlign(aTextAlign);
  aPT->SetFillColor(0);
  aPT->SetTextAlign(22);
  aPT->SetTextColor(aTextColor);
  aPT->AddText(aText);
}

std::vector<MvtxMonDraw::Quality> MvtxMonDraw::analyseForError(TH2Poly* over1, TH2Poly* over2,TH2Poly* over3, TH1* decErr){
      std::vector<Quality> result((NFlags+1)*NLAYERS+6);
      bool badStave[NFlags][NSTAVE];
      for (int iflag = 0; iflag < NFlags; iflag++) {
        for (int ilayer = 0; ilayer < NLAYERS; ilayer++) {
          result.at((3*iflag)+ilayer) = Quality::Good;
}
}
      //count dead staves
      for (int iflag = 0; iflag < NFlags; iflag++) {
        for (int ilayer = 0; ilayer < NLAYERS; ilayer++) {
          int countStave = 0;
          for (int ibin = StaveBoundary[ilayer] + 1; ibin <= StaveBoundary[ilayer + 1]; ++ibin) {
            badStave[iflag][ibin-1] = false;
            double bincontent = 0;
            if(iflag==0 && over1) bincontent = over1->GetBinContent(ibin);
            if(iflag==1 && over2) bincontent = over2->GetBinContent(ibin);
            if(iflag==2 && over3) bincontent = over3->GetBinContent(ibin);
            if (bincontent/*hp[iflag][NFlx]->GetBinContent(ibin)*/ > maxbadchips / 9.) {
              badStave[iflag][ibin-1] = true;
              //std::cout<<"bad stave"<<std::endl;
              countStave++;
              result.at((3*iflag)+ilayer) = Quality::Medium;
            }
          }
          if (countStave > 0.25 * NStaves[ilayer]) {
            result.at((3*iflag)+ilayer) = Quality::Bad;
          }
        }
      }
      
      // check if more than 25% staves with error, warning or fault
      for (int ilayer = 0; ilayer < NLAYERS; ilayer++) {
        int countStave = 0;
        result.at(ilayer+9) = Quality::Good;
        for (int ibin = StaveBoundary[ilayer] + 1; ibin <= StaveBoundary[ilayer + 1]; ++ibin) {
            if(badStave[0][ibin-1] || badStave[1][ibin-1] || badStave[2][ibin-1]){
              result.at(ilayer+9) = Quality::Medium;
              countStave++;
            }
          }
          if (countStave > 0.25 * NStaves[ilayer]) {
            result.at(ilayer+9) = Quality::Bad;
          }
        }

      for (int iflx = 0; iflx < 6; iflx++) {
        if(decErr){
          if(decErr->Integral(6*iflx+1,6*iflx+7)>0) result.at(12+iflx) = Quality::Bad;
          else result.at(12+iflx) = Quality::Good;
        }
        else result.at(12+iflx) = Quality::Good;
      }
      
      return result;



     /* if (mo->getName() == Form("LaneStatus/laneStatusFlag%s", mLaneStatusFlag[iflag].c_str())) {
        result = Quality::Good;
        auto* h = dynamic_cast<TH2I*>(mo->getObject());
        if (h->GetMaximum() > 0) {
          result.set(Quality::Bad);
        }
      }*/
      //if (mo->getName() == Form("LaneStatus/laneStatusOverviewFlag%s", mLaneStatusFlag[iflag].c_str())) {
     /*   result.at(iflag) = Quality.Good;
        auto* hp = dynamic_cast<TH2Poly*>(mo->getObject());
        badStave = false;
        // Initialization of metaData for IB, ML, OL
        for (int ilayer = 0; ilayer < NLAYERS; ilayer++) {
          int countStave = 0;
          badStaveCount = false;
          for (int ibin = StaveBoundary[ilayer] + 1; ibin <= StaveBoundary[ilayer + 1]; ++ibin) {
              // Check if there are staves in the IB with lane in Bad (bins are filled with %)
              if (hp->GetBinContent(ibin) > maxbadchips / 9.) {
                badStave = true;
                result.at(iflag) = Quality.Medium;
                countStave++;
              }
          } // end loop bins (staves)
          // Initialize metadata for the 7 layers
          result.addMetadata(Form("Layer%d", ilayer), "good");
          // Check if there are more than 25% staves in Bad per layer
          if (countStave > 0.25 * NStaves[ilayer]) {
            badStaveCount = true;
            result.updateMetadata(Form("Layer%d", ilayer), "bad");
          }
        } // end loop over layers
        if (badStave) {
          result.at(iflag) = Quality.Medium;
        }
        if (badStaveCount) {
          result.at(iflag) = Quality.Bad;
        }
     // } // end lanestatusOverview
    }*/

}

void MvtxMonDraw::DrawPave(std::vector<MvtxMonDraw::Quality> status, int position, const char*/* what*/){
  TPaveText *pt = new TPaveText(.4,.4,.6,.6,"blNDC");
  pt->SetTextSize(0.04);
  pt->SetFillColor(0);
  pt->SetLineColor(0);
  pt->SetBorderSize(1);
  if(status.at(position) == Quality::Good && status.at(position+1) == Quality::Good && status.at(position+2) == Quality::Good) pt->AddText("#color[418]{QA OK}");
  if(status.at(position) == Quality::Medium) pt->AddText("#color[808]{QA Layer 0 Medium}");
  if(status.at(position+1) == Quality::Medium) pt->AddText("#color[808]{QA Layer 1 Medium}");
  if(status.at(position+2) == Quality::Medium) pt->AddText("#color[808]{QA Layer 2 Medium}");
  if(status.at(position) == Quality::Bad) pt->AddText("#color[2]{QA Layer 0 Bad}");
  if(status.at(position+1) == Quality::Bad) pt->AddText("#color[2]{QA Layer 1 Bad}");
  if(status.at(position+2) == Quality::Bad) pt->AddText("#color[2]{QA Layer 2 Bad}");
  pt->Draw();

}

time_t MvtxMonDraw::getTime()
{
  OnlMonClient *cl = OnlMonClient::instance();
  time_t currtime = cl->EventTime("CURRENT");
  return currtime;
}
