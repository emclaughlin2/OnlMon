#include "CemcMonDraw.h"

#include <onlmon/OnlMonClient.h>

#include <TAxis.h>  // for TAxis
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TH2D.h>
#include <TPad.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TText.h>
#include <TStyle.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TLegend.h>
#include <TProfile.h>

#include <cstring>  // for memset
#include <ctime>
#include <fstream>
#include <iostream>  // for operator<<, basic_ostream, basic_os...
#include <sstream>
#include <vector>  // for vector

CemcMonDraw::CemcMonDraw(const std::string &name)
: OnlMonDraw(name)
{
  
  return;
}

int CemcMonDraw::Init()
{
  std::cout << "initializing" << std::endl;

  cemcStyle = new TStyle("cemcStyle","cemcStyle");

  Int_t font=42; // Helvetica
  cemcStyle->SetLabelFont(font,"x");
  cemcStyle->SetTitleFont(font,"x");
  cemcStyle->SetLabelFont(font,"y");
  cemcStyle->SetTitleFont(font,"y");
  cemcStyle->SetLabelFont(font,"z");
  cemcStyle->SetTitleFont(font,"z");
  cemcStyle->SetOptStat(0);
  cemcStyle->SetPadTickX(1);
  cemcStyle->SetPadTickY(1);

  gROOT->SetStyle("cemcStyle");
  gROOT->ForceStyle();

  return 0;
}


int CemcMonDraw::MakeCanvas(const std::string &name)
{
  OnlMonClient *cl = OnlMonClient::instance();
  int xsize = cl->GetDisplaySizeX();
  int ysize = cl->GetDisplaySizeY();
  if (name == "CemcMon1")
    {
      // xpos (-1) negative: do not draw menu bar
      TC[0] = new TCanvas(name.c_str(), "CemcMon Tower Hits",  xsize/3 , 0, xsize/3, ysize*0.9);
      // root is pathetic, whenever a new TCanvas is created root piles up
      // 6kb worth of X11 events which need to be cleared with
      // gSystem->ProcessEvents(), otherwise your process will grow and
      // grow and grow but will not show a definitely lost memory leak
      gSystem->ProcessEvents();
      Pad[0] = new TPad("cemcpad0","hit map",0.,0.2,1.,1.);
      Pad[0]->Draw();
      // this one is used to plot the run number on the canvas
      transparent[0] = new TPad("transparent0", "this does not show", 0, 0, 1, 1);
      transparent[0]->SetFillStyle(4000);
      transparent[0]->Draw();
    
      // warning 
      warning[0] = new TPad("warning0", "hot tower warnings", 0, 0, 1, 0.2);
      warning[0]->SetFillStyle(4000);
      warning[0]->Draw();
      TC[0]->SetEditable(0);
    }
  else if (name == "CemcMon2")
    {
      // xpos negative: do not draw menu bar
      TC[1] = new TCanvas(name.c_str(), "CemcMon2 Packet Information",  2*xsize /3, 0,2*xsize / 3, ysize*0.9);
      gSystem->ProcessEvents();
      Pad[1] = new TPad("cemcpad1", "packet event check", 0.0, 0.6, 1.0/2, 0.95, 0);
      Pad[2] = new TPad("cemcpad2", "packet size", 0.0, 0.3, 1.0/2, 0.6, 0);
      Pad[3] = new TPad("cemcpad3", "packet channels", 0.0, 0.0, 1.0/2, 0.3, 0);
      Pad[1]->Draw();
      Pad[2]->Draw();
      Pad[3]->Draw();
      // this one is used to plot the run number on the canvas
      transparent[1] = new TPad("transparent1", "this does not show", 0, 0, 1., 1);
      transparent[1]->SetFillStyle(4000);
      transparent[1]->Draw();
      
      //packet warnings
      warning[1] = new TPad("warning1","packet warnings",0.5,0,1,1);
      warning[1] -> SetFillStyle(4000);
      warning[1] -> Draw();
      TC[1]->SetEditable(0);
    }
  else if (name == "CemcMon3")
    {
      TC[2] = new TCanvas(name.c_str(), "CemcMon Waveform Processing", xsize/3 , 0, xsize/3, ysize*0.9);
      gSystem->ProcessEvents();
      Pad[4] = new TPad("cemcpad4", "who needs this?", 0.0, 0.6, 1.0, 0.95, 0);
      Pad[5] = new TPad("cemcpad5", "who needs this?", 0.0, 0.3, 1.0, 0.6, 0);
      Pad[6] = new TPad("cemcpad6", "who needs this?", 0.0, 0.0, 1.0, 0.3, 0);
      Pad[4]->Draw();
      Pad[5]->Draw();
      Pad[6]->Draw();
      // this one is used to plot the run number on the canvas
      transparent[2] = new TPad("transparent3", "this does not show", 0, 0, 1, 1);
      transparent[2]->SetFillStyle(4000);
      transparent[2]->Draw();
      TC[2]->SetEditable(0);
    }
  else if (name == "CemcMon4")
    {
      TC[3] = new TCanvas(name.c_str(), "CemcMon Waveform Processing Expert", xsize / 3, 0, xsize / 3, ysize*0.9);
      gSystem->ProcessEvents();
      Pad[7] = new TPad("cemcpad7", "who needs this?", 0.0, 0.6, 1.0, 0.95, 0);
      Pad[8] = new TPad("cemcpad8", "who needs this?", 0.0, 0.3, 1.0, 0.6, 0);
      Pad[9] = new TPad("cemcpad9", "who needs this?", 0.0, 0.0, 1.0, 0.3, 0);
      Pad[7]->Draw();
      Pad[8]->Draw();
      Pad[9]->Draw();
      // this one is used to plot the run number on the canvas
      transparent[3] = new TPad("transparent3", "this does not show", 0, 0, 1, 1);
      transparent[3]->SetFillStyle(4000);
      transparent[3]->Draw();
      TC[3]->SetEditable(0);
    }
  return 0;
}



int CemcMonDraw::Draw(const std::string &what)
{
  int iret = 0;
  int idraw = 0;
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
  if (what == "ALL" || what == "THIRD")
    {
      iret += DrawThird(what);
      idraw++;
    }
  if (what == "ALL" || what == "FOURTH")
    {
      iret += DrawFourth(what);
      idraw++;
    }
  if (!idraw)
    {
      std::cout << __PRETTY_FUNCTION__ << " Unimplemented Drawing option: " << what << std::endl;
      iret = -1;
    }
  return iret;
}



int CemcMonDraw::DrawFirst(const std::string & /* what */)
{
  OnlMonClient *cl = OnlMonClient::instance();
  //watch the absolute insanity as we merge all these 
  //histograms from across seven different machines
  TH2* hist1[nSEBs];
  const int nHists = 4;
  int start[nHists];
  start[0] = -1;
  for(int i = 0; i < nSEBs; i++)
    {
      hist1[i] = (TH2*)cl->getHisto(Form("CEMCMON_%d",i),"h2_cemc_rm");
      if(hist1[i] && start[0] == -1) start[0] = i;
      if(start[0] > -1 && hist1[i])
	{
	  hist1[i] -> SetName(Form("h2_cemc_rm_%d",i));
	  if(i != start[0])hist1[start[0]] -> Add(hist1[i],1);
	}
    }
  

  TH2* h2_cemc_mean[nSEBs];
  start[1] = -1;
  for(int i = 0; i < nSEBs; i++)
    {
      h2_cemc_mean[i] = (TH2*)cl->getHisto(Form("CEMCMON_%d",i),"h2_cemc_mean");
      if( h2_cemc_mean[i] && start[1] == -1) start[1] = i;
      if(start[1] > -1 && h2_cemc_mean[i])
	{
	  h2_cemc_mean[i] -> SetName(Form("h2_cemc_mean_%d",i));
	  if(i != start[1])h2_cemc_mean[start[1]] -> Add(h2_cemc_mean[i],1);
	}
    }					 

  TH1* h_event[nSEBs];
  TH1* h_eventSource[nSEBs];
  start[2] = -1;
  float maxEvent = 1;
  for(int i = 0; i < nSEBs; i++)
    {
      h_eventSource[i] = (TH1*)cl->getHisto(Form("CEMCMON_%d",i),"h1_event");
      if(h_eventSource[i] && start[2] == -1) start[2] = i;
      if(start[2] > -1 && h_eventSource[i])
	{
	  h_event[i] = (TH1*)h_eventSource[i] -> Clone();
	 
	  h_event[i] -> SetName(Form("h1_event_%d",i));
	  if(/*i != start[2]*/h_event[i] -> GetEntries() > maxEvent)
	    {
	      //h_event[start[2]] -> Add(h_event[i],1);
	      maxEvent = h_event[i] -> GetEntries();
	    }
	}
    }
  
  if (start[0] < 0)
    {
      DrawDeadServer(transparent[0]);
      TC[0]->SetEditable(0);
      return -1;
    }

  //h_event[start[2]] -> Scale(1./divisor);
  
  TH1* adcCount[nSEBs];
  start[3] = -1;
  for(int i = 0; i < nSEBs; i++)
    {
      adcCount[i] = (TH1*)cl->getHisto(Form("CEMCMON_%i",i),"h1_cemc_adc");
      
      if(adcCount[i] && start[3] == -1) start[3] = i;
      if(start[3] > -1 && adcCount[i])
	{
	  adcCount[i] -> SetName(Form("h1_cemc_adc_%d",i));
	  if(i != start[3])adcCount[start[3]] -> Add(adcCount[i],1);
	}
    }
  
  if (!gROOT->FindObject("CemcMon1"))
    {
      MakeCanvas("CemcMon1");
    }

 
  if(maxEvent > 0)
    {
      
      if(maxEvent > templateDepth)h2_cemc_mean[start[1]]->Scale(1./templateDepth);
      else h2_cemc_mean[start[1]]->Scale(1./maxEvent);
      
	
      //hist1[start[0]] -> Divide(h2_cemc_mean[start[1]]);
      //hist1[start[0]] -> Scale(1./hist1[start[0]] -> GetMean());
      //hist1[start[0]] -> Scale(1/10.);
      if(adcCount[start[3]]->GetMean())h2_cemc_mean[start[1]]->Scale(1./adcCount[start[3]]->GetMean());
  
      if(adcCount[start[3]]->GetMean())hist1[start[0]] -> Scale(1./adcCount[start[3]]->GetMean());
  
      for(int i = 0; i < nTowersEta; i++)
      	{
      	  for(int j = 0; j < nTowersPhi; j++)
      	    {
      	      if(i < 8) continue;
      	      if(h2_cemc_mean[start[1]]-> GetBinContent(i+1, j+1) < 0.75 && hist1[start[0]] -> GetBinContent(i+1, j+1) < 0.75) hist1[start[0]]->SetBinContent(i+1, j+1, h2_cemc_mean[start[1]] -> GetBinContent(i+1, j+1));
      	      else hist1[start[0]] -> SetBinContent(i+1, j+1, hist1[start[0]]->GetBinContent(i+1,j+1)/h2_cemc_mean[start[1]]->GetBinContent(i+1,j+1));
      	    }
      	} 
    }
  else
    {
      for(int i = 0; i < nTowersEta; i++)
	{
	  for(int j = 0; j < nTowersPhi; j++)
	    {
	      if(i < 8) continue;
	     
	      else hist1[start[0]]->SetBinContent(i+1, j+1, 0);
	    }
	}
    }
  //hist1 = (TH2D*)hist1 -> Divide(h2_cemc_mean);
  //TH1 *temp = (TH1*)h2_cemc_mean -> Clone();
  //temp -> GetXaxis() -> SetRangeUser(8,temp -> GetNbinsX());//accounts for the fact that we're missing
  //adc boards on the south side. 

  
  

  TC[0]->SetEditable(1);
  TC[0]->Clear("D");
  Pad[0]->cd();
 
  hist1[start[0]]->GetXaxis()->SetTitle("eta index");
  hist1[start[0]]->GetYaxis()->SetTitle("phi index");
  hist1[start[0]]->GetZaxis()->SetTitle("Tower Running Mean/ Histogram Running Mean");
  hist1[start[0]]->GetXaxis()->CenterTitle();
  hist1[start[0]]->GetYaxis()->CenterTitle();
  hist1[start[0]]->GetZaxis()->CenterTitle();
  hist1[start[0]]->GetXaxis()->SetNdivisions(12,kFALSE);
  hist1[start[0]]->GetYaxis()->SetNdivisions(32,kFALSE);

  float tsize = 0.03;
  hist1[start[0]]->GetXaxis()->SetLabelSize(tsize);
  hist1[start[0]]->GetYaxis()->SetLabelSize(tsize);
  hist1[start[0]]->GetYaxis()->SetTitleOffset(1.4);
  hist1[start[0]]->GetZaxis()->SetLabelSize(tsize);
  hist1[start[0]]->GetXaxis()->SetTitleSize(tsize);
  hist1[start[0]]->GetYaxis()->SetTitleSize(tsize);
  hist1[start[0]]->GetXaxis()->SetTickLength(0.02);

  //hist1->GetZaxis()->SetRangeUser(0,2);
  
  TLine *line_sector[32];
  for(int i_line=0;i_line<32;i_line++)
    {
      line_sector[i_line] = new TLine(0,(i_line+1)*8,96,(i_line+1)*8);
      line_sector[i_line]->SetLineColor(1);
      line_sector[i_line]->SetLineWidth(1);
      line_sector[i_line]->SetLineStyle(1);
    }


  const int numVertDiv = 12;
  int dEI = 96/numVertDiv;
  TLine *l_board[numVertDiv-1];
  for(int il=1; il<numVertDiv; il++){
    l_board[il-1] = new TLine(dEI*il,0,dEI*il,256);
    l_board[il-1]->SetLineColor(1);
    l_board[il-1]->SetLineWidth(1);
    l_board[il-1]->SetLineStyle(1);
    if(il==6) l_board[il-1]->SetLineWidth(2);
  }
  
  gPad->SetTopMargin(0.08);
  gPad->SetBottomMargin(0.07);
  gPad->SetLeftMargin(0.08);
  gPad->SetRightMargin(0.12);
  
  

  // modify palette to black, green, and red
  Int_t palette[3] = {kGray+2,8,2};
  cemcStyle->SetPalette(3,palette);
  gROOT->SetStyle("cemcStyle");
  gROOT->ForceStyle();
  gStyle->SetPalette(3,palette);
  double_t levels[4] = {0,0.75,1.25,2};
  hist1[start[0]] -> GetZaxis() -> SetRangeUser(0,2);
  hist1[start[0]] ->SetContour(4,levels);
  gStyle -> SetOptStat(0);
  hist1[start[0]]->DrawCopy("colz");
  //h2_cemc_mean[start[1]]->DrawCopy("colz");
  for(int i_line=0;i_line<32;i_line++) line_sector[i_line]->Draw();
  for(int il=0; il<numVertDiv-1; il++) l_board[il]->Draw();

  FindHotTower(warning[0],hist1[start[0]]);
  TText PrintRun;
  PrintRun.SetTextFont(62);
  PrintRun.SetTextSize(0.04);
  PrintRun.SetNDC();          // set to normalized coordinates
  PrintRun.SetTextAlign(23);  // center/top alignment
  std::ostringstream runnostream;
  std::ostringstream runnostream2;
  std::string runstring;
  time_t evttime = cl->EventTime("CURRENT");
  // fill run number and event time into string
  runnostream  << "Running mean of Tower Hits Normalized by All Towers " ;
  runnostream2 << "Run " << cl->RunNumber() << ", Time: " << ctime(&evttime);
  transparent[0]->cd();
  runstring = runnostream.str();
  PrintRun.DrawText(0.5, 0.99, runstring.c_str());
  runstring = runnostream2.str();
  PrintRun.DrawText(0.5, 0.966, runstring.c_str());

  TC[0]->Update();
  TC[0]->Show();
  TC[0]->SetEditable(0);
  if(save)TC[0] -> SaveAs("plots/towerHits.pdf");
  return 0;
}



int CemcMonDraw::DrawSecond(const std::string & /* what */)
{
  
  OnlMonClient *cl = OnlMonClient::instance();
  gStyle -> SetOptStat(0);
  TH1* h_event[nSEBs];
  TH1* h_eventSource[nSEBs];
  int start[4];
  start[0] = -1;
  float maxEvent = -1;
  for(int i = 0; i < nSEBs; i++)
    {
      h_eventSource[i] = (TH1*) cl->getHisto(Form("CEMCMON_%d",i),"h1_event");
      if(h_eventSource[i] && start[0] == -1) start[0] = i;
      if(start[0] > -1 && h_eventSource[i])
	{
	  h_event[i] = (TH1*)h_eventSource[i] -> Clone();
	  h_event[i] -> SetName(Form("h1_event_%d",i));
	  if(/*i != start[0] &&*/ h_event[i] -> GetEntries() > maxEvent)
	    {
	      //h_event[start[0]] -> Add(h_event[i],1);
	      maxEvent = h_event[i] -> GetEntries();
	    }
	}
    }
  
  TH1* h1_packet_number[nSEBs];
  start[1] = -1;
  for(int i = 0; i < nSEBs; i++)
    {
      h1_packet_number[i] = (TH1*) cl->getHisto(Form("CEMCMON_%d",i),"h1_packet_number");
      if(h1_packet_number[i] && start[1] == -1) start[1] = i;
      if(start[1] > -1 && h1_packet_number[i])
	{
	  h1_packet_number[i] -> SetName(Form("h1_papcket_number_%d",i));
	  if(i != start[1])h1_packet_number[start[1]] -> Add(h1_packet_number[i],1);
	}
    }

  TH1* h1_packet_length[nSEBs];
  start[2] = -1;
  for(int i = 0; i < nSEBs; i++)
    {
      h1_packet_length[i] = (TH1*) cl->getHisto(Form("CEMCMON_%d",i),"h1_packet_length");
      if(h1_packet_length[i] && start[2] == -1) start[2] = i;
      if(start[2] > -1 && h1_packet_length[i])
	{
	  h1_packet_length[i] -> SetName(Form("h1_papcket_length_%d",i));
	  if(i != start[2])h1_packet_length[start[2]] -> Add(h1_packet_length[i],1);
	}
    }

  
  TH1* h1_packet_chans[nSEBs];
  start[3] = -1;
  for(int i = 0; i < nSEBs; i++)
    {
      h1_packet_chans[i] = (TH1*) cl->getHisto(Form("CEMCMON_%d",i),"h1_packet_chans");
      if(h1_packet_chans[i] && start[3] == -1) start[3] = i;
      if(start[3] > -1 && h1_packet_chans[i])
	{
	  h1_packet_chans[i] -> SetName(Form("h1_papcket_chans_%d",i));
	  if(i != start[3])h1_packet_chans[start[3]] -> Add(h1_packet_chans[i],1);
	}
    }
  
    if (start[0] < 0  || start[1] < 0 || start[2] < 0)
    {
      DrawDeadServer(transparent[1]);
      TC[1]->SetEditable(0);
      return -1;
    }
  
    //h_event[start[0]] -> Scale(1./divisor);
  
  if(maxEvent > 0)
    {
      // h1_packet_number[start[1]] -> Scale(1./h_event[start[0]] -> GetBinContent(1));
      // h1_packet_length[start[2]] -> Scale(1./h_event[start[0]] -> GetBinContent(1));
      // h1_packet_chans[start[3]] -> Scale(1./h_event[start[0]] -> GetBinContent(1));
      h1_packet_number[start[1]] -> Scale(1./maxEvent);
      h1_packet_length[start[2]] -> Scale(1./maxEvent);
      h1_packet_chans[start[3]] -> Scale(1./maxEvent);
    }
  
  if (!gROOT->FindObject("CemcMon2"))
    {
      MakeCanvas("CemcMon2");
    }
  
  TC[1]->SetEditable(1);
  TC[1]->Clear("D");
  
  TLine *one = new TLine(6000.5,1,6128.5,1);
  one -> SetLineStyle(7);
  
  TLine *goodSize = new TLine(6000.5,5981,6128.5,5981);
  goodSize -> SetLineStyle(7);
  
  TLine *goodChans = new TLine(6000.5,192,6128.5,192);
  goodChans -> SetLineStyle(7);

  float param = 0.95;
  //float param = 0.99;
  
  TLegend *leg = new TLegend(0.3,0.70,0.95,0.90);
  leg -> SetFillStyle(0);
  leg -> SetBorderSize(0);
  
  TLine *warnLineOne = new TLine(6000.5,param*1,6128.5,param*1);
  warnLineOne -> SetLineStyle(7);
  warnLineOne -> SetLineColor(2);
  
  leg -> AddEntry(warnLineOne,Form("%g%% Threshold",param*100),"l");
  
  TLine *warnLineOneS = new TLine(6000.5,param*1,6128.5,param*1);
  warnLineOneS -> SetLineStyle(10);
  warnLineOneS -> SetLineColor(2);

  leg -> AddEntry(warnLineOneS,Form("%g%% Threshold, High Eta, South",100*param),"l");

  TLine *warnLineSize = new TLine(6000.5,param*5981.,6128.5,param*5981.);
  warnLineSize -> SetLineStyle(7);
  warnLineSize -> SetLineColor(2);

  TLine *warnLineSizeS = new TLine(6000.5,param*3991.,6128.5,param*3991.);
  warnLineSizeS -> SetLineStyle(10);
  warnLineSizeS -> SetLineColor(2);

  TLine *warnLineChans = new TLine(6000.5,param*192.,6128.5,param*192.);
  warnLineChans -> SetLineStyle(7);
  warnLineChans -> SetLineColor(2);

  TLine *warnLineChansS = new TLine(6000.5,param*128.,6128.5,param*128.);
  warnLineChansS -> SetLineStyle(10);
  warnLineChansS -> SetLineColor(2);

  Pad[1]->cd();
  float tsize = 0.08;
  h1_packet_number[start[1]] -> GetYaxis() -> SetRangeUser(0.0,1.3);
  std::vector<std::vector<int>> badPackets;
  badPackets.push_back(getBadPackets(h1_packet_number[start[1]],0,param));
 
  h1_packet_number[start[1]]->GetXaxis()->SetNdivisions(510,kTRUE);
  h1_packet_number[start[1]]->GetXaxis()->SetTitle("packet #");
  h1_packet_number[start[1]]->GetYaxis()->SetTitle("% Of Events Present");
  h1_packet_number[start[1]]->GetXaxis()->SetLabelSize(tsize-0.01);
  h1_packet_number[start[1]]->GetYaxis()->SetLabelSize(tsize-0.01);
  h1_packet_number[start[1]]->GetXaxis()->SetTitleSize(tsize-0.01);
  h1_packet_number[start[1]]->GetYaxis()->SetTitleSize(tsize-0.01);
  h1_packet_number[start[1]]->GetXaxis()->SetTitleOffset(1);
  gPad->SetBottomMargin(0.16);
  gPad->SetLeftMargin(0.16);
  gPad->SetRightMargin(0.05);
  gPad->SetLeftMargin(0.15);
  gPad->SetTicky();
  gPad->SetTickx();
  h1_packet_number[start[1]] -> DrawCopy("hist");
  one -> Draw("same");
  warnLineOne -> Draw("same");

  Pad[2]->cd();
  h1_packet_length[start[2]] -> GetYaxis() -> SetRangeUser(0,6500);
  badPackets.push_back(getBadPackets(h1_packet_length[start[2]],1,param));
 
  h1_packet_length[start[2]]->GetXaxis()->SetNdivisions(510,kTRUE);
  h1_packet_length[start[2]]->GetXaxis()->SetTitle("packet #");
  h1_packet_length[start[2]]->GetYaxis()->SetTitle("Average Packet Size");
  h1_packet_length[start[2]]->GetXaxis()->SetLabelSize(tsize-.01);
  h1_packet_length[start[2]]->GetYaxis()->SetLabelSize(tsize);
  h1_packet_length[start[2]]->GetXaxis()->SetTitleSize(tsize-.01);
  h1_packet_length[start[2]]->GetYaxis()->SetTitleSize(tsize);
  h1_packet_length[start[2]]->GetXaxis()->SetTitleOffset(1);
  h1_packet_length[start[2]]->GetYaxis()->SetTitleOffset(0.9);
  h1_packet_length[start[2]]->GetYaxis()->SetRangeUser(0,10000);
  gPad->SetBottomMargin(0.16);
  gPad->SetLeftMargin(0.16);
  gPad->SetRightMargin(0.05);
  gPad->SetLeftMargin(0.15);
  gPad->SetTicky();
  gPad->SetTickx();
  h1_packet_length[start[2]] -> DrawCopy("hist");
  goodSize -> Draw("same");
  warnLineSize -> Draw("same");
  warnLineSizeS -> Draw("same");
  leg -> Draw("same");

 
  Pad[3]->cd();
  h1_packet_chans[start[3]] -> GetYaxis() -> SetRangeUser(0,212);
  badPackets.push_back(getBadPackets(h1_packet_chans[start[3]],2,param));
  h1_packet_chans[start[3]]->GetXaxis()->SetNdivisions(510,kTRUE);
  h1_packet_chans[start[3]]->GetXaxis()->SetTitle("packet #");
  h1_packet_chans[start[3]]->GetYaxis()->SetTitle("Average # of Channels");
  h1_packet_chans[start[3]]->GetXaxis()->SetLabelSize(tsize-.01);
  h1_packet_chans[start[3]]->GetYaxis()->SetLabelSize(tsize);
  h1_packet_chans[start[3]]->GetXaxis()->SetTitleSize(tsize-.01);
  h1_packet_chans[start[3]]->GetYaxis()->SetTitleSize(tsize);
  h1_packet_chans[start[3]]->GetXaxis()->SetTitleOffset(0.8);
  h1_packet_chans[start[3]]->GetYaxis()->SetTitleOffset(0.8);
  gPad->SetBottomMargin(0.16);
  gPad->SetLeftMargin(0.16);
  gPad->SetRightMargin(0.05);
  gPad->SetLeftMargin(0.15);
  gPad->SetTicky();
  gPad->SetTickx();
  h1_packet_chans[start[3]] -> DrawCopy("hist");
  goodChans -> Draw("same");
  warnLineChans -> Draw("same");
  warnLineChansS -> Draw("same");

  warning[1] -> cd();
  TLegend *badPacks = new TLegend(0,1/4.,1,0.9);
  TPaveText *title = new TPaveText(0,0.9,1,0.95);
  title -> AddText("Bad Packets");
  badPacks -> SetNColumns(8);
  badPacks -> SetTextSize(0.03);
  float badboys = 0;
  if(maxEvent > 0)
    {
      for(int i = 0; i < 3; i++)
	{
	  for(int j = 0; j < (int)badPackets[i].size(); j++)
	    {
	      //there's most certainly a better way to do this but it's 5:00 on day 5 of owl shift
	      //just want to prevent a packet showing up multiple times and crowding the screen
	      if(badPackets[i][j] == 0) continue;//need this to prevent seg faulting
	  
	      if(i == 0)
		{
		  badPacks -> AddEntry("",Form("%d",badPackets[i][j]),"");
		  badboys++;
		}
	      else if(i == 1)  
		{
		  if(!(std::count(badPackets[i-1].begin(),badPackets[i-1].end(),badPackets[i][j])))
		    {
		      badPacks -> AddEntry("",Form("%d",badPackets[i].at(j)),"");
		      badboys++;
		    }
		}
	      else if(i == 2)
		{
		  if(!(std::count(badPackets[i-1].begin(),badPackets[i-1].end(),badPackets[i][j]))) 
		    {
		      badPacks -> AddEntry("",Form("%d",badPackets[i].at(j)),"");
		      badboys++;
		    }
		  if(!(std::count(badPackets[i-2].begin(),badPackets[i-2].end(),badPackets[i][j])) && !(std::count(badPackets[i-1].begin(),badPackets[i-1].end(),badPackets[i][j])))
		    {
		      badPacks -> AddEntry("",Form("%d",badPackets[i].at(j)),"");
		      badboys++;
		    }
		}
	    }
	}
    }
  badPacks -> Draw();
  title -> Draw();
  TPaveText *desc = new TPaveText(0,0,1,0.25);
  desc -> AddText(Form("Currently %.2g%% of packets are reporting a problem",badboys/128.*100));
  desc -> AddText("Packets will be reported bad as above for the following reasons:" );
  desc -> AddText(Form("A packet appears in less than %g %% of events",param*100));
  desc -> AddText(Form("A packet is less than %g %% of size 5981",param*100));
  desc -> AddText(Form("A packet sees fewer than %g %% of 192 channels",param*100));
  desc -> Draw();
  TText PrintRun;
  PrintRun.SetTextFont(62);
  PrintRun.SetTextSize(0.02);
  PrintRun.SetNDC();          // set to normalized coordinates
  PrintRun.SetTextAlign(23);  // center/top alignment
  std::ostringstream runnostream;
  std::string runstring;
  std::ostringstream runnostream2;
  time_t evttime = cl->EventTime("CURRENT");
  // fill run number and event time into string
  
  runnostream << "Packet Information";
  runnostream2 << " Run " << cl->RunNumber() << ", Time: " << ctime(&evttime);
  transparent[1]->cd();

  runstring = runnostream.str();
  PrintRun.DrawText(0.5,.99, runstring.c_str());
  
  runstring = runnostream2.str();
  PrintRun.DrawText(0.5, .966, runstring.c_str());
  TC[1]->Update();
  TC[1]->Show();
  TC[1]->SetEditable(0);
  if(save)TC[1] -> SaveAs("plots/packets.pdf");
  return 0;
}



int CemcMonDraw::DrawThird(const std::string & /* what */)
{
  OnlMonClient *cl = OnlMonClient::instance();
 
  TH2* h2_waveform_twrAvg[nSEBs]; 
  int start[3];
  start[0] = -1;
 
  for(int i = 0; i < nSEBs; i++)
    {
      h2_waveform_twrAvg[i] = (TH2*) cl->getHisto(Form("CEMCMON_%d",i),"h2_waveform_twrAvg");
      if(h2_waveform_twrAvg[i] && start[0] == -1) start[0] = i;
      if(start[0] > -1 && h2_waveform_twrAvg[i])
	{
	  h2_waveform_twrAvg[i] -> SetName(Form("h2_waveform_twrAvg_%d",i));
	  h2_waveform_twrAvg[start[0]] -> Add(h2_waveform_twrAvg[i],1);
	}
    }
  
  TH1* h1_waveform_time[nSEBs]; 
  start[1] = -1;
  for(int i = 0; i < nSEBs; i++)
    {
      h1_waveform_time[i] = (TH1*) cl->getHisto(Form("CEMCMON_%d",i),"h1_waveform_time");
      if(h1_waveform_time[i] && start[1] == -1) start[1] = i;
      if(start[1] > -1 && h1_waveform_time[i])
	{
	  h1_waveform_time[i] -> SetName(Form("h1_waveform_time_%d",i));
	  h1_waveform_time[start[1]] -> Add(h1_waveform_time[i],1);
	}
    }
  
  TH1* h1_waveform_pedestal[nSEBs]; 
  start[2] = -1;
  for(int i = 0; i < nSEBs; i++)
    {
      h1_waveform_pedestal[i] = (TH1*) cl->getHisto(Form("CEMCMON_%d",i),"h1_waveform_pedestal");

      if(h1_waveform_pedestal[i] && start[2] == -1) start[2] = i;
      if(start[2] > -1 && h1_waveform_pedestal[i])
	{
	  h1_waveform_pedestal[i] -> SetName(Form("h1_waveform_pedestal_%d",i));
	  h1_waveform_pedestal[start[2]] -> Add(h1_waveform_pedestal[i],1);
	}
    }

  if (!gROOT->FindObject("CemcMon3"))
    {
      MakeCanvas("CemcMon3");
    }
  
  TC[2]->SetEditable(1);
  TC[2]->Clear("D");
  Pad[4]->cd();
  if (start[0] < 0 || start[1] < 0|| start[2] < 0)
    {
      DrawDeadServer(transparent[2]);
      TC[2]->SetEditable(0);
      return -1;
    }

  gStyle->SetTitleFontSize(0.03);
  float ymaxp = h2_waveform_twrAvg[start[0]]->ProfileX()->GetMaximum();
  h2_waveform_twrAvg[start[0]] -> GetYaxis() -> SetRangeUser(0,ymaxp*10);

  float tsize = 0.06;
  float tsize2 = 0.08;

  h2_waveform_twrAvg[start[0]]->GetXaxis()->SetNdivisions(16);
  h2_waveform_twrAvg[start[0]]->GetXaxis()->SetTitle("sample #");
  h2_waveform_twrAvg[start[0]]->GetYaxis()->SetTitle("Waveform ADC in latest event [ADC]");
  h2_waveform_twrAvg[start[0]]->GetXaxis()->SetLabelSize(tsize-.01);
  h2_waveform_twrAvg[start[0]]->GetYaxis()->SetLabelSize(tsize);
  h2_waveform_twrAvg[start[0]]->GetXaxis()->SetTitleSize(tsize-.01);
  h2_waveform_twrAvg[start[0]]->GetYaxis()->SetTitleSize(tsize);
  h2_waveform_twrAvg[start[0]]->GetXaxis()->SetTitleOffset(1.);
  h2_waveform_twrAvg[start[0]]->GetYaxis()->SetTitleOffset(1.2);
  TLine *windowLow1 = new TLine(4,0,4,ymaxp*10);
  TLine *windowHigh1 = new TLine(10,0,10,ymaxp*10);
  gStyle -> SetOptStat(0);
  gPad->SetBottomMargin(0.16);
  gPad->SetLeftMargin(0.16);
  gPad->SetRightMargin(0.05);
  gPad->SetLeftMargin(0.15);
  gPad->SetTicky();
  gPad->SetTickx();
  h2_waveform_twrAvg[start[0]]->DrawCopy("colz");
  windowLow1 -> Draw("same");
  windowHigh1 -> Draw("same");
  gStyle -> SetPalette(57);
  gPad -> SetLogz();

  TText PrintRun;
  PrintRun.SetTextFont(62);
  PrintRun.SetTextSize(0.03);
  PrintRun.SetNDC();          // set to normalized coordinates
  PrintRun.SetTextAlign(23);  // center/top alignment
  std::ostringstream runnostream;
  std::ostringstream runnostream2;
  std::string runstring;
  time_t evttime = cl->EventTime("CURRENT");
  // fill run number and event time into string
  runnostream << "Waveform fitting";
  runnostream2 <<"Run " << cl->RunNumber() << ", Time: " << ctime(&evttime);
  
  transparent[2]->cd();

  runstring = runnostream.str();
  PrintRun.DrawText(0.5, 0.99, runstring.c_str());

  runstring = runnostream2.str();
  PrintRun.DrawText(0.5, 0.966, runstring.c_str());
  
  Pad[5]->cd();

  gStyle->SetTitleFontSize(0.06);

  h1_waveform_time[start[1]]->GetXaxis()->SetNdivisions(16);
  h1_waveform_time[start[1]]->GetXaxis()->SetTitle("waveform peak position [sample #]");
  h1_waveform_time[start[1]]->GetYaxis()->SetTitle("Fraction of Towers");
  h1_waveform_time[start[1]]->GetXaxis()->SetLabelSize(tsize2);
  h1_waveform_time[start[1]]->GetYaxis()->SetLabelSize(tsize2);
  h1_waveform_time[start[1]]->GetXaxis()->SetTitleSize(tsize2);
  h1_waveform_time[start[1]]->GetYaxis()->SetTitleSize(tsize2);
  h1_waveform_time[start[1]]->GetXaxis()->SetTitleOffset(1.0);
  h1_waveform_time[start[1]]->GetYaxis()->SetTitleOffset(.85);
  
  if(h1_waveform_time[start[1]]->GetEntries())h1_waveform_time[start[1]] -> Scale(1./h1_waveform_time[start[1]]->GetEntries());
  TLine *windowLow2 = new TLine(4,0,4,h1_waveform_time[start[1]]->GetMaximum());
  TLine *windowHigh2 = new TLine(10,0,10,h1_waveform_time[start[1]]->GetMaximum());
  gPad->SetTopMargin(0.06);
  gPad->SetBottomMargin(0.18);
  gPad->SetRightMargin(0.05);
  gStyle -> SetOptStat(0);
  gPad->SetLeftMargin(0.15);
  gPad->SetTicky();
  gPad->SetTickx();
  h1_waveform_time[start[1]]->DrawCopy("hist");
  windowLow2 -> Draw("same");
  windowHigh2 -> Draw("same");

  Pad[6]->cd();

  gStyle->SetTitleFontSize(0.06);

  h1_waveform_pedestal[start[2]]->GetXaxis()->SetNdivisions(8);
  h1_waveform_pedestal[start[2]]->GetXaxis()->SetTitle("ADC Pedestal");
  h1_waveform_pedestal[start[2]]->GetYaxis()->SetTitle("Fraction of Towers");
  h1_waveform_pedestal[start[2]]->GetXaxis()->SetLabelSize(tsize2);
  h1_waveform_pedestal[start[2]]->GetYaxis()->SetLabelSize(tsize2);
  h1_waveform_pedestal[start[2]]->GetXaxis()->SetTitleSize(tsize2);
  h1_waveform_pedestal[start[2]]->GetYaxis()->SetTitleSize(tsize2);
  h1_waveform_pedestal[start[2]]->GetXaxis()->SetTitleOffset(1);
  h1_waveform_pedestal[start[2]]->GetYaxis()->SetTitleOffset(0.85);
  if(h1_waveform_pedestal[start[2]]->GetEntries())h1_waveform_pedestal[start[2]] -> Scale(1./h1_waveform_pedestal[start[2]]->GetEntries());
  //h1_waveform_pedestal -> GetXaxis() -> SetRangeUser(1000,2000);
  gPad->SetTopMargin(0.06);
  gPad->SetBottomMargin(0.18);
  gPad->SetRightMargin(0.05);
  gStyle -> SetOptStat(0);
  gPad->SetLeftMargin(0.15);
  gPad->SetTicky();
  gPad->SetTickx();
  h1_waveform_pedestal[start[2]]->DrawCopy("hist");

  TC[2]->Update();
  TC[2]->Show();
  TC[2]->SetEditable(0);
  if(save)TC[2] -> SaveAs("plots/waveform.pdf");
  return 0;
}

int CemcMonDraw::DrawFourth(const std::string & /* what */)
{

  OnlMonClient *cl = OnlMonClient::instance();
  
  TH1 *h_waveform_sigDiff[nSEBs]; 
  int start[3];
  start[0] = -1;
  for(int i = 0; i < nSEBs; i++)
    {
      h_waveform_sigDiff[i] = (TH1*) cl->getHisto(Form("CEMCMON_%d",i),"h1_fitting_sigDiff");
      if(h_waveform_sigDiff[i] && start[0] == -1) start[0] = i;
      if(start[0] > -1 && h_waveform_sigDiff[i])
	{
	  h_waveform_sigDiff[i] -> SetName(Form("h_fitting_sigDiff_%d",i));
	  h_waveform_sigDiff[start[0]] -> Add(h_waveform_sigDiff[i],1);
	}
    }

  TH1 *h_waveform_pedDiff[nSEBs]; 
  start[1] = -1;
  for(int i = 0; i < nSEBs; i++)
    {
      h_waveform_pedDiff[i] = (TH1*) cl->getHisto(Form("CEMCMON_%d",i),"h1_fitting_pedDiff");
      if(h_waveform_pedDiff[i] && start[1] == -1) start[1] = i;
      if(start[1] > -1 && h_waveform_pedDiff[i])
	{
	  h_waveform_pedDiff[i] -> SetName(Form("h_fitting_pedDiff_%d",i));
	  h_waveform_pedDiff[start[1]] -> Add(h_waveform_pedDiff[i],1);
	}
    }
  
  TH1 *h_waveform_timeDiff[nSEBs];
  start[2] = -1;
  for(int i = 0; i < nSEBs; i++)
    {
      h_waveform_timeDiff[i] = (TH1*) cl->getHisto(Form("CEMCMON_%d",i),"h1_fitting_timeDiff");
      if(h_waveform_timeDiff[i] && start[2] == -1) start[2] = i;
      if(start[2] > -1 && h_waveform_timeDiff[i])
	{
	  h_waveform_timeDiff[i] -> SetName(Form("h_fitting_timeDiff_%d",i));
	  h_waveform_timeDiff[start[2]] -> Add(h_waveform_timeDiff[i],1);
	}
    }
  
  if(!gROOT -> FindObject("CemcMon4"))
    {
      MakeCanvas("CemcMon4");
    }
  
  TC[3]->SetEditable(1);
  TC[3]->Clear("D");
  if(start[0] < 0|| start[1] < 0 || start[2] < 0)
    {
      DrawDeadServer(transparent[3]);
      TC[3]->SetEditable(0);
      return -1;
    }
  
  gStyle->SetTitleFontSize(0.03);
    
  
  TText PrintRun;
  PrintRun.SetTextFont(62);
  PrintRun.SetTextSize(0.03);
  PrintRun.SetNDC();          // set to normalized coordinates
  PrintRun.SetTextAlign(23);  // center/top alignment
  std::ostringstream runnostream;
  std::ostringstream runnostream2;
  std::string runstring;
  time_t evttime = cl->EventTime("CURRENT");
  // fill run number and event time into string
  runnostream << "Waveform Template vs. Fast Fitting";
  runnostream2 << "Run " << cl->RunNumber() << ", Time: " << ctime(&evttime);
  
  transparent[3]->cd();

  runstring = runnostream.str();
  PrintRun.DrawText(0.5, 0.99, runstring.c_str());
  
  
  runstring = runnostream2.str();
  PrintRun.DrawText(0.5, 0.966, runstring.c_str());
  
  
  Pad[7]->cd();

  gStyle->SetTitleFontSize(0.06);
  
  float tsize2 = 0.08;
  h_waveform_sigDiff[start[0]]->GetXaxis()->SetNdivisions(16);
  h_waveform_sigDiff[start[0]]->GetXaxis()->SetTitle("Fast ADC/Template ADC");
  h_waveform_sigDiff[start[0]]->GetYaxis()->SetTitle("Fraction of Towers");
  h_waveform_sigDiff[start[0]]->GetXaxis()->SetLabelSize(tsize2-.01);
  h_waveform_sigDiff[start[0]]->GetYaxis()->SetLabelSize(tsize2-.01);
  h_waveform_sigDiff[start[0]]->GetXaxis()->SetTitleSize(tsize2-.01);
  h_waveform_sigDiff[start[0]]->GetYaxis()->SetTitleSize(tsize2-.01);
  h_waveform_sigDiff[start[0]]->GetXaxis()->SetTitleOffset(.9);
  h_waveform_sigDiff[start[0]]->GetYaxis()->SetTitleOffset(1.1);
  if(h_waveform_sigDiff[start[0]] -> GetEntries())h_waveform_sigDiff[start[0]] -> Scale(1./h_waveform_sigDiff[start[0]] -> GetEntries());
  gPad->SetTopMargin(0.06);
  gPad->SetBottomMargin(0.18);
  gPad->SetRightMargin(0.05);
  gPad->SetLeftMargin(0.15);
  gStyle->SetOptStat(0);
  gPad->SetTicky();
  gPad->SetTickx();
  h_waveform_sigDiff[start[0]]->DrawCopy("hist");

  Pad[8]->cd();

  gStyle->SetTitleFontSize(0.06);
  
  h_waveform_pedDiff[start[1]]->GetXaxis()->SetNdivisions(16);
  h_waveform_pedDiff[start[1]]->GetXaxis()->SetTitle("Fast Pedestal/Template Pedestal");
  h_waveform_pedDiff[start[1]]->GetYaxis()->SetTitle("Fraction of Towers");
  h_waveform_pedDiff[start[1]]->GetXaxis()->SetLabelSize(tsize2);
  h_waveform_pedDiff[start[1]]->GetYaxis()->SetLabelSize(tsize2);
  h_waveform_pedDiff[start[1]]->GetXaxis()->SetTitleSize(tsize2);
  h_waveform_pedDiff[start[1]]->GetYaxis()->SetTitleSize(tsize2);
  h_waveform_pedDiff[start[1]]->GetXaxis()->SetTitleOffset(0.9);
  h_waveform_pedDiff[start[1]]->GetYaxis()->SetTitleOffset(0.9);
  if(h_waveform_pedDiff[start[1]] -> GetEntries())h_waveform_pedDiff[start[1]] -> Scale(1./ h_waveform_pedDiff[start[1]] -> GetEntries());
  gPad->SetTopMargin(0.06);
  gPad->SetBottomMargin(0.18);
  gPad->SetRightMargin(0.05);
  gPad->SetLeftMargin(0.15);
  gStyle->SetOptStat(0);
  gPad->SetTicky();
  gPad->SetTickx();
  h_waveform_pedDiff[start[1]]->DrawCopy("hist");

  Pad[9]->cd();

  h_waveform_timeDiff[start[2]]->GetXaxis()->SetNdivisions(16);
  h_waveform_timeDiff[start[2]]->GetXaxis()->SetTitle("Fast Peak Time - Template Peak Time");
  h_waveform_timeDiff[start[2]]->GetYaxis()->SetTitle("Fraction of Towers");
  h_waveform_timeDiff[start[2]]->GetXaxis()->SetLabelSize(tsize2);
  h_waveform_timeDiff[start[2]]->GetYaxis()->SetLabelSize(tsize2);
  h_waveform_timeDiff[start[2]]->GetXaxis()->SetTitleSize(tsize2);
  h_waveform_timeDiff[start[2]]->GetYaxis()->SetTitleSize(tsize2);
  h_waveform_timeDiff[start[2]]->GetXaxis()->SetTitleOffset(0.9);
  h_waveform_timeDiff[start[2]]->GetYaxis()->SetTitleOffset(0.9);
  if(h_waveform_timeDiff[start[2]]->GetEntries())h_waveform_timeDiff[start[2]] -> Scale(1./ h_waveform_timeDiff[start[2]]->GetEntries());
  gPad->SetTopMargin(0.06);
  gPad->SetBottomMargin(0.18);
  gPad->SetRightMargin(0.05);
  gPad->SetLeftMargin(0.15);
  gPad->SetTicky();
  gPad->SetTickx();
  h_waveform_timeDiff[start[2]]->DrawCopy("hist");


  TC[3]->Update();
  TC[3]->Show();
  TC[3]->SetEditable(0);
  gStyle->SetOptStat(0);
  if(save)TC[3] -> SaveAs("plots/waveformExpert.pdf");
  return 0;
}


int CemcMonDraw::FindHotTower(TPad *warningpad,TH2* hhit){
  float nhott = 0;
  float ndeadt = 0;
  float hot_threshold  = 1.25;
  float dead_threshold = 0.75;
  float nTowerTotal = 24576.-2048.;//-2048 to account for the non-functioning towers at the edge of the south
  for(int ieta=0; ieta<nTowersEta; ieta++){
    for(int iphi=0; iphi<nTowersPhi; iphi++){
      
      if(ieta < 8) continue;
      if(hhit -> GetBinContent(ieta+1, iphi+1) == 0)continue;
      double nhit = hhit->GetBinContent(ieta+1, iphi+1);
	
      if(nhit > hot_threshold)	nhott++;
      
	
      if(nhit < dead_threshold)	ndeadt++;
      
    }
  }

  //draw warning here
  warningpad->cd();
  TPaveText *dead = new TPaveText(0.01,0.7,0.33,1);
  dead -> SetFillColor(kGray+2);
  dead -> SetTextColor(kWhite);
  dead -> AddText(Form("Cold towers: %.3g%%",100*ndeadt/nTowerTotal));
  // if(100*ndeadt/nTowerTotal > 2.5)
  //   {
  //     dead -> AddText("");
  //   }
  TPaveText *good = new TPaveText(0.33,0.7,0.66,1);
  good -> SetFillColor(kGreen+1);
  good -> AddText(Form("Good towers: %.3g%%",100*(nTowerTotal-ndeadt-nhott)/nTowerTotal));
  TPaveText *hot = new TPaveText(0.66,0.7,1,1);
  hot -> SetFillColor(kRed-9);
  hot -> AddText(Form("Hot towers: %.3g%%",100*nhott/nTowerTotal));
  // if(100*nhott/nTowerTotal > 3.5)
  //   {
  //     hot -> AddText("");
  //   }
  TPaveText *warn = new TPaveText(0.01,0.1,1,0.7);
  warn -> SetTextSize(0.1);
  warn -> AddText("Helpful Numbers: 1 Box = Interface Board (IB)");
  warn -> AddText("3 Boxes (horiz) = 1 Packet; 6 Boxes (horiz) = 1 Sector");
  //warn -> AddText("Inform Expert if Number of Dead Towers Exceeds: 2.5%");
  //warn -> AddText("Inform Expert if Number of Hot Towers Exceeds: 3.5%");
  warn -> AddText("For now, watch for entire IB's or sectors going dead");
  dead -> Draw();
  good -> Draw();

  hot -> Draw();
  warn -> Draw();
  warningpad->Update();
  return 0;
}



int CemcMonDraw::SavePlot(const std::string &what, const std::string &type)
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

int CemcMonDraw::MakeHtml(const std::string &what)
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


std::vector<int> CemcMonDraw::getBadPackets(TH1 *hist, int what, float cutoff)
{
  float params[3] = {1., 5981., 192.};
  float params2[3] = {1., 3991., 128.};
  //float cutoff = 0.75;
  std::vector<int> badpacks = {0};
  
  for(int i = 1; i < hist -> GetNbinsX(); i++)
    {
      if((hist -> GetBinContent(i) < params[what]*cutoff ||  hist -> GetBinContent(i) > params[what]) &&  !(((6000 + i - 2.)/4.) == floor(((6000 + i - 2)/4))) ) badpacks.push_back(i+6000);
      
      else if((hist -> GetBinContent(i) < params2[what]*cutoff || hist -> GetBinContent(i) > params2[what]) &&  (((6000 + i - 2.)/4.) == floor(((6000 + i - 2)/4))) ) badpacks.push_back(i+6000);
    }
  
  return badpacks;
}

time_t CemcMonDraw::getTime()
{
  OnlMonClient *cl = OnlMonClient::instance();
  time_t currtime = cl->EventTime("CURRENT");
  return currtime;
}
