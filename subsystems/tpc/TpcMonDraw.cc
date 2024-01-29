#include "TpcMonDraw.h"

#include <onlmon/OnlMonClient.h>

#include <TAxis.h>  // for TAxis
#include <TCanvas.h>
#include <TDatime.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TPad.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TText.h>
#include <TMath.h>
#include <TPaveLabel.h>
#include <TStyle.h>
#include <TString.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TEllipse.h>

#include <cstring>  // for memset
#include <ctime>
#include <fstream>
#include <iostream>  // for operator<<, basic_ostream, basic_os...
#include <sstream>
#include <vector>  // for vector

#include <cmath>
#include <cstdio>  // for printf
#include <string>  // for allocator, string, char_traits

TpcMonDraw::TpcMonDraw(const std::string &name)
  : OnlMonDraw(name)
{
  return;
}

int TpcMonDraw::Init()
{
  return 0;
}

int TpcMonDraw::MakeCanvas(const std::string &name)
{
  OnlMonClient *cl = OnlMonClient::instance();
  int xsize = cl->GetDisplaySizeX();
  int ysize = cl->GetDisplaySizeY();
  if (name == "TpcMon1")
  {
    // xpos (-1) negative: do not draw menu bar
    TC[0] = new TCanvas(name.c_str(), "TpcMon Example Monitor", -1, 0, xsize / 2, ysize);
    // root is pathetic, whenever a new TCanvas is created root piles up
    // 6kb worth of X11 events which need to be cleared with
    // gSystem->ProcessEvents(), otherwise your process will grow and
    // grow and grow but will not show a definitely lost memory leak
    gSystem->ProcessEvents();
    Pad[0] = new TPad("tpcpad1", "who needs this?", 0.1, 0.5, 0.9, 0.9, 0);
    Pad[1] = new TPad("tpcpad2", "who needs this?", 0.1, 0.05, 0.9, 0.45, 0);
    Pad[0]->Draw();
    Pad[1]->Draw();
    // this one is used to plot the run number on the canvas
    transparent[0] = new TPad("transparent0", "this does not show", 0, 0, 1, 1);
    transparent[0]->SetFillStyle(4000);
    transparent[0]->Draw();
    TC[0]->SetEditable(false);
  }
  else if (name == "TpcMon2")
  {
    // xpos negative: do not draw menu bar
    TC[1] = new TCanvas(name.c_str(), "TpcMon2 Example Monitor", -xsize / 2, 0, xsize / 2, ysize);
    gSystem->ProcessEvents();
    Pad[2] = new TPad("tpcpad3", "who needs this?", 0.1, 0.5, 0.9, 0.9, 0);
    Pad[3] = new TPad("tpcpad4", "who needs this?", 0.1, 0.05, 0.9, 0.45, 0);
    Pad[2]->Draw();
    Pad[3]->Draw();
    // this one is used to plot the run number on the canvas
    transparent[1] = new TPad("transparent1", "this does not show", 0, 0, 1, 1);
    transparent[1]->SetFillStyle(4000);
    transparent[1]->Draw();
    TC[1]->SetEditable(false);
  }
  else if (name == "TpcMon3")
  {
    TC[2] = new TCanvas(name.c_str(), "TpcMon3 Example Monitor", xsize / 2, 0, xsize / 2, ysize);
    gSystem->ProcessEvents();
    Pad[4] = new TPad("tpcpad5", "who needs this?", 0.1, 0.5, 0.9, 0.9, 0);
    Pad[5] = new TPad("tpcpad6", "who needs this?", 0.1, 0.05, 0.9, 0.45, 0);
    Pad[4]->Draw();
    Pad[5]->Draw();
    // this one is used to plot the run number on the canvas
    //        transparent[2] = new TPad("transparent2", "this does not show", 0, 0, 1, 1);
    //        transparent[2]->SetFillStyle(4000);
    //        transparent[2]->Draw();
    //      TC[2]->SetEditable(0);
  }
  else if (name == "TPCModules")
  {
    TC[3] = new TCanvas(name.c_str(), "ADC Count by GEM Example", 1350, 700);
    gSystem->ProcessEvents();
    TC[3]->Divide(2,1);
    //gStyle->SetPalette(57); //kBird CVD friendly
    // this one is used to plot the run number on the canvas
    transparent[2] = new TPad("transparent2", "this does not show", 0, 0, 1, 1);
    transparent[2]->SetFillStyle(4000);
    transparent[2]->Draw();
    TC[3]->SetEditable(false);
  }
  else if (name == "TPCSampleSize")
  {
    TC[4] = new TCanvas(name.c_str(), "TPC Sample Size Distribution in Events", -1, 0, xsize , ysize );
    gSystem->ProcessEvents();
    TC[4]->Divide(4,7);
    transparent[3] = new TPad("transparent3", "this does not show", 0, 0, 1, 1);
    transparent[3]->SetFillStyle(4000);
    transparent[3]->Draw();
    TC[4]->SetEditable(false);
  }
  else if (name == "TPCCheckSumError")
  {
    TC[5] = new TCanvas(name.c_str(), "TPC CheckSumError Probability in Events",-1, 0, xsize , ysize );
    gSystem->ProcessEvents();
    TC[5]->Divide(4,7);
    transparent[4] = new TPad("transparent4", "this does not show", 0, 0, 1, 1);
    transparent[4]->SetFillStyle(4000);
    transparent[4]->Draw();
    TC[5]->SetEditable(false);
  }
  else if (name == "TPCADCSample")
  {
    TC[6] = new TCanvas(name.c_str(), "TPC ADC vs Sample in Whole Sector",-1, 0, xsize , ysize);
    gSystem->ProcessEvents();
    //gStyle->SetPalette(57); //kBird CVD friendly
    TC[6]->Divide(4,7);
    transparent[5] = new TPad("transparent5", "this does not show", 0, 0, 1, 1);
    transparent[5]->SetFillStyle(4000);
    transparent[5]->Draw();
    TC[6]->SetEditable(false);
  }
  else if (name == "TPCMaxADCModule")
  {
    TC[7] = new TCanvas(name.c_str(), "(MAX ADC - pedestal) in SLIDING WINDOW for each Module in Sector", -1, 0, xsize , ysize);
    gSystem->ProcessEvents();
    //gStyle->SetPalette(57); //kBird CVD friendly
    TC[7]->Divide(4,7);
    transparent[6] = new TPad("transparent6", "this does not show", 0, 0, 1, 1);
    transparent[6]->SetFillStyle(4000);
    transparent[6]->Draw();
    TC[7]->SetEditable(false);
  }
  else if (name == "TPCRawADC1D")
  {
    TC[8] = new TCanvas(name.c_str(), "TPC RAW ADC 1D distribution", -1, 0, xsize , ysize);
    gSystem->ProcessEvents();
    //gStyle->SetPalette(57); //kBird CVD friendly
    TC[8]->Divide(4,7);
    transparent[7] = new TPad("transparent7", "this does not show", 0, 0, 1, 1);
    transparent[7]->SetFillStyle(4000);
    transparent[7]->Draw();
    TC[8]->SetEditable(false);
  }
  else if (name == "TPCMaxADC1D")
  {
    TC[9] = new TCanvas(name.c_str(), "(MAX ADC - pedestal) in SLIDING WINDOW 1D distribution",-1, 0, xsize , ysize);
    gSystem->ProcessEvents();
    //gStyle->SetPalette(57); //kBird CVD friendly
    TC[9]->Divide(4,7);
    transparent[8] = new TPad("transparent8", "this does not show", 0, 0, 1, 1);
    transparent[8]->SetFillStyle(4000);
    transparent[8]->Draw();
    TC[9]->SetEditable(false);
  }
  else if (name == "TPCClusterXY")
  {
    TC[10] = new TCanvas(name.c_str(), "(MAX ADC - pedestal)> 20 ADC for NS and SS, WEIGHTED", 1350, 700);
    gSystem->ProcessEvents();
    //gStyle->SetPalette(57); //kBird CVD friendly
    TC[10]->Divide(2,1);
    // this one is used to plot the run number on the canvas
    transparent[9] = new TPad("transparent9", "this does not show", 0, 0, 1, 1);
    transparent[9]->SetFillStyle(4000);
    transparent[9]->Draw();
    TC[10]->SetEditable(false);
  }
  else if (name == "TPCClusterXY_unw")
  {
    TC[11] = new TCanvas(name.c_str(), "(MAX ADC - pedestal)> 20 ADC for NS and SS, UNWEIGHTED", 1350, 700);
    gSystem->ProcessEvents();
    //gStyle->SetPalette(57); //kBird CVD friendly
    TC[11]->Divide(2,1);
    // this one is used to plot the run number on the canvas
    transparent[10] = new TPad("transparent10", "this does not show", 0, 0, 1, 1);
    transparent[10]->SetFillStyle(4000);
    transparent[10]->Draw();
    TC[11]->SetEditable(false);
  }
  else if (name == "TPCADCSamplelarge")
  {
    TC[12] = new TCanvas(name.c_str(), "TPC ADC vs Large Sample in Whole Sector",-1, 0, xsize , ysize);
    gSystem->ProcessEvents();
    //gStyle->SetPalette(57); //kBird CVD friendly
    TC[12]->Divide(4,7);
    transparent[11] = new TPad("transparent11", "this does not show", 0, 0, 1, 1);
    transparent[11]->SetFillStyle(4000);
    transparent[11]->Draw();
    TC[12]->SetEditable(false);
  }

  else if (name == "TPCClusterZY")
  {
    TC[13] = new TCanvas(name.c_str(), "(MAX ADC - pedestal)> 20 ADC for NS and SS, WEIGHTED", 1350, 700);
    gSystem->ProcessEvents();
    //gStyle->SetPalette(57); //kBird CVD friendly
    TC[13]->Divide(1,1);
    // this one is used to plot the run number on the canvas
    transparent[12] = new TPad("transparent12", "this does not show", 0, 0, 1, 1);
    transparent[12]->SetFillStyle(4000);
    transparent[12]->Draw();
    TC[13]->SetEditable(false);
  }

  else if (name == "TPCClusterZY_unw")
  {
    TC[14] = new TCanvas(name.c_str(), "(MAX ADC - pedestal)> 20 ADC for NS and SS, UNWEIGHTED", 1350, 700);
    gSystem->ProcessEvents();
    //gStyle->SetPalette(57); //kBird CVD friendly
    TC[14]->Divide(1,1);
    // this one is used to plot the run number on the canvas
    transparent[13] = new TPad("transparent13", "this does not show", 0, 0, 1, 1);
    transparent[13]->SetFillStyle(4000);
    transparent[13]->Draw();
    TC[14]->SetEditable(false);
  }
 
  
  return 0;
}

int TpcMonDraw::Draw(const std::string &what)
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
  if (what == "ALL" || what == "TPCMODULE")
  {
    iret += DrawTPCModules(what);
    idraw++;
  }
  if (what == "ALL" || what == "TPCSAMPLESIZE")
  {
    iret += DrawTPCSampleSize(what);
    idraw++;
  }
  if (what == "ALL" || what == "TPCCHECKSUMERROR")
  {
    iret += DrawTPCCheckSum(what);
    idraw++;
  }
  if (what == "ALL" || what == "TPCADCVSSAMPLE")
  {
    iret += DrawTPCADCSample(what);
    idraw++;
  }
  if (what == "ALL" || what == "TPCMAXADCMODULE")
  {
    iret += DrawTPCMaxADCModule(what);
    idraw++;
  }
  if (what == "ALL" || what == "TPCRAWADC1D")
  {
    iret += DrawTPCRawADC1D(what);
    idraw++;
  }
  if (what == "ALL" || what == "TPCMAXADC1D")
  {
    iret += DrawTPCMaxADC1D(what);
    idraw++;
  }
  if (what == "ALL" || what == "TPCCLUSTERSXYWEIGTHED")
  {
    iret += DrawTPCXYclusters(what);
    idraw++;
  }
  if (what == "ALL" || what == "TPCCLUSTERSXYUNWEIGTHED")
  {
    iret += DrawTPCXYclusters_unweighted(what);
    idraw++;
  }
  if (what == "ALL" || what == "TPCADCVSSAMPLELARGE")
  {
    iret += DrawTPCADCSampleLarge(what);
    idraw++;
  }
  if (what == "ALL" || what == "TPCCLUSTERSZYWEIGTHED")
  {
    iret += DrawTPCZYclusters(what);
    idraw++;
  }
  if (what == "ALL" || what == "TPCCLUSTERSZYUNWEIGTHED")
  {
    iret += DrawTPCZYclusters_unweighted(what);
    idraw++;
  }
  if (!idraw)
  {
    std::cout << __PRETTY_FUNCTION__ << " Unimplemented Drawing option: " << what << std::endl;
    iret = -1;
  }
  return iret;
}

int TpcMonDraw::DrawFirst(const std::string & /* what */)
{
  OnlMonClient *cl = OnlMonClient::instance();
  TH1 *tpcmon_hist1 =  cl->getHisto("TPCMON_0","tpcmon_hist1");
  TH1 *tpcmon_hist2 =  cl->getHisto("TPCMON_0","tpcmon_hist1");

  if (!gROOT->FindObject("TpcMon1"))
  {
    MakeCanvas("TpcMon1");
  }
  TC[0]->SetEditable(true);
  TC[0]->Clear("D");
  Pad[0]->cd();
  if (tpcmon_hist1)
  {
    tpcmon_hist1->DrawCopy();
  }
  else
  {
    DrawDeadServer(transparent[0]);
    TC[0]->SetEditable(false);
    return -1;
  }
  Pad[1]->cd();
  if (tpcmon_hist2)
  {
    tpcmon_hist2->DrawCopy();
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
  TC[0]->SetEditable(false);
  return 0;
}

int TpcMonDraw::DrawSecond(const std::string & /* what */)
{
  OnlMonClient *cl = OnlMonClient::instance();
  TH1 *tpcmon_hist1 =  cl->getHisto("TPCMON_0","tpcmon_hist2");
  TH1 *tpcmon_hist2 =  cl->getHisto("TPCMON_0","tpcmon_hist2");
  if (!gROOT->FindObject("TpcMon2"))
  {
    MakeCanvas("TpcMon2");
  }
  TC[1]->SetEditable(true);
  TC[1]->Clear("D");
  Pad[2]->cd();
  if (tpcmon_hist1)
  {
    tpcmon_hist1->DrawCopy();
  }
  else
  {
    DrawDeadServer(transparent[1]);
    TC[1]->SetEditable(false);
    return -1;
  }
  Pad[3]->cd();
  if (tpcmon_hist2)
  {
    tpcmon_hist2->DrawCopy();
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
  TC[1]->SetEditable(false);
  return 0;
}

int TpcMonDraw::DrawTPCModules(const std::string & /* what */)
{
  OnlMonClient *cl = OnlMonClient::instance();

  TH2 *tpcmon_NSIDEADC[24] = {nullptr};
  TH2 *tpcmon_SSIDEADC[24] = {nullptr};

  char TPCMON_STR[100];
  // TPC ADC pie chart
  for( int i=0; i<24; i++ ) 
  {
    //const TString TPCMON_STR( Form( "TPCMON_%i", i ) );
    sprintf(TPCMON_STR,"TPCMON_%i",i);
    tpcmon_NSIDEADC[i] = (TH2*) cl->getHisto(TPCMON_STR,"NorthSideADC");
    tpcmon_SSIDEADC[i] = (TH2*) cl->getHisto(TPCMON_STR,"SouthSideADC");
  }


  //TH2 *tpcmon_NSIDEADC1 = (TH2*) cl->getHisto("TPCMON_0","NorthSideADC");
  //TH2 *tpcmon_SSIDEADC1 = (TH2*) cl->getHisto("TPCMON_0","SouthSideADC");

  if (!gROOT->FindObject("TPCModules"))
  {
    MakeCanvas("TPCModules");
  }

  dummy_his1 = new TH2F("dummy_his1", "ADC Counts North Side", 100, -1.5, 1.5, 100, -1.5, 1.5); //dummy histos for titles
  dummy_his2 = new TH2F("dummy_his2", "ADC Counts South Side", 100, -1.5, 1.5, 100, -1.5, 1.5);

  //labels
  NS18 = new TPaveLabel( 1.046586,-0.1938999,1.407997,0.2144871, "18" );
  NS17 = new TPaveLabel( 0.962076,0.4382608,1.323487,0.8466479 , "17" );
  NS16 = new TPaveLabel( 0.4801947,0.8802139,0.8416056,1.288601 , "16" );
  NS15 = new TPaveLabel( -0.1823921,1.011681,0.1790189,1.425662, "15" );
  NS14 = new TPaveLabel( -0.8449788,0.8690253,-0.4835679,1.288601 , "14" );
  NS13 = new TPaveLabel( -1.30879,0.441058,-0.9473786,0.8550394 , "13" );
  NS12 = new TPaveLabel( -1.411009,-0.2050886,-1.049598,0.2144871, "12" );
  NS23 = new TPaveLabel( -1.302585,-0.7757116,-0.9471979,-0.3561359 , "23" );
  NS22 = new TPaveLabel( -0.8449788,-1.309971,-0.4835679,-0.8848013 , "22" );
  NS21 = new TPaveLabel( -0.1823921,-1.426557,0.1790189,-1.006982 , "21" );
  NS20 = new TPaveLabel( 0.4801947,-1.309076,0.8416056,-0.8839062 , "20" );
  NS19 = new TPaveLabel( 0.9622567,-0.7785088,1.323668,-0.3533387 , "19" );

  SS00 = new TPaveLabel( 1.046586,-0.1938999,1.407997,0.2144871, "00" );
  SS01 = new TPaveLabel( 0.962076,0.4382608,1.323487,0.8466479 , "01" );
  SS02 = new TPaveLabel( 0.4801947,0.8802139,0.8416056,1.288601 , "02" );
  SS03 = new TPaveLabel( -0.1823921,1.011681,0.1790189,1.425662, "03" );
  SS04 = new TPaveLabel( -0.8449788,0.8690253,-0.4835679,1.288601 , "04" );
  SS05 = new TPaveLabel( -1.30879,0.441058,-0.9473786,0.8550394 , "05" );
  SS06 = new TPaveLabel( -1.411009,-0.2050886,-1.049598,0.2144871, "06" );
  SS07 = new TPaveLabel( -1.302585,-0.7757116,-0.9471979,-0.3561359 , "07" );
  SS08 = new TPaveLabel( -0.8449788,-1.309971,-0.4835679,-0.8848013 , "08" );
  SS09 = new TPaveLabel( -0.1823921,-1.426557,0.1790189,-1.006982 , "09" );
  SS10 = new TPaveLabel( 0.4801947,-1.309076,0.8416056,-0.8839062 , "10" );
  SS11 = new TPaveLabel( 0.9622567,-0.7785088,1.323668,-0.3533387 , "11" );

  NS18->SetFillColor(0);
  NS17->SetFillColor(0);
  NS16->SetFillColor(0);
  NS15->SetFillColor(0);
  NS14->SetFillColor(0);
  NS13->SetFillColor(0);
  NS12->SetFillColor(0);
  NS23->SetFillColor(0);
  NS22->SetFillColor(0);
  NS21->SetFillColor(0);
  NS20->SetFillColor(0);
  NS19->SetFillColor(0);

  SS00->SetFillColor(0);
  SS01->SetFillColor(0);
  SS02->SetFillColor(0);
  SS03->SetFillColor(0);
  SS04->SetFillColor(0);
  SS05->SetFillColor(0);
  SS06->SetFillColor(0);
  SS07->SetFillColor(0);
  SS08->SetFillColor(0);
  SS09->SetFillColor(0);
  SS10->SetFillColor(0);
  SS11->SetFillColor(0);

  TC[3]->SetEditable(true);
  TC[3]->Clear("D");
  TC[3]->cd(1);
  gPad->SetTopMargin(0.15);
  gStyle->SetOptStat(0);
  dummy_his1->Draw("colpolzsame");
 
  float NS_max = 0;
  for( int i=0; i<12; i++ )
  {
    if( tpcmon_NSIDEADC[i] ){
    TC[3]->cd(1);
    tpcmon_NSIDEADC[i] -> DrawCopy("colpolzsame");
    if( tpcmon_NSIDEADC[i]->GetBinContent(tpcmon_NSIDEADC[i]->GetMaximumBin()) > NS_max)
    {
      NS_max = tpcmon_NSIDEADC[i]->GetBinContent(tpcmon_NSIDEADC[i]->GetMaximumBin());
      dummy_his1->SetMaximum( NS_max );
    }
    gStyle->SetPalette(57); //kBird CVD friendly
    }

  }
  TC[3]->Update();
  TC[3]->cd(1);
  SS00->Draw("same");
  SS01->Draw("same");
  SS02->Draw("same");
  SS03->Draw("same");
  SS04->Draw("same");
  SS05->Draw("same");
  SS06->Draw("same");
  SS07->Draw("same");
  SS08->Draw("same");
  SS09->Draw("same");
  SS10->Draw("same");
  SS11->Draw("same");

  TC[3]->cd(2);
  gPad->SetTopMargin(0.15);
  gStyle->SetOptStat(0);
  dummy_his2->Draw("colpolzsame");

  float SS_max = 0;
  for( int i=0; i<12; i++ )
  {
    if( tpcmon_SSIDEADC[i+12] ){
    //std::cout<<"tpcmon_SSIDEADC i: "<< i+12 <<std::endl;
    TC[3]->cd(2);
    tpcmon_SSIDEADC[i+12] -> DrawCopy("colpolzsame");

    if ( tpcmon_SSIDEADC[i+12]->GetBinContent(tpcmon_SSIDEADC[i+12]->GetMaximumBin()) > SS_max)
    {
      SS_max = tpcmon_SSIDEADC[i+12]->GetBinContent(tpcmon_SSIDEADC[i+12]->GetMaximumBin());
      dummy_his2->SetMaximum( SS_max );
    }

    gStyle->SetPalette(57); //kBird CVD friendly
    }
  }
  TC[3]->Update();

  TC[3]->cd(2);
  NS18->Draw("same");
  NS17->Draw("same");
  NS16->Draw("same");
  NS15->Draw("same");
  NS14->Draw("same");
  NS13->Draw("same");
  NS12->Draw("same");
  NS23->Draw("same");
  NS22->Draw("same");
  NS21->Draw("same");
  NS20->Draw("same");
  NS19->Draw("same");


  TText PrintRun;
  PrintRun.SetTextFont(62);
  PrintRun.SetTextSize(0.04);
  PrintRun.SetNDC();          // set to normalized coordinates
  PrintRun.SetTextAlign(23);  // center/top alignment
  std::ostringstream runnostream;
  std::string runstring;
  time_t evttime = cl->EventTime("CURRENT");
  // fill run number and event time into string
  runnostream << ThisName << "_SUMADC Run " << cl->RunNumber()
              << ", Time: " << ctime(&evttime);
  runstring = runnostream.str();
  transparent[2]->cd();
  PrintRun.DrawText(0.5, 0.91, runstring.c_str());

  //turn off stats box
  //dummy_his1->SetStats(0);
  //dummy_his2->SetStats(0);
  
  TC[3]->Update();



  TC[3]->Show();
  TC[3]->SetEditable(false);
  
  return 0;
}

int TpcMonDraw::DrawTPCSampleSize(const std::string & /* what */)
{
  OnlMonClient *cl = OnlMonClient::instance();

  TH1 *tpcmon_samplesizedist[24] = {nullptr};

  char TPCMON_STR[100];
  for( int i=0; i<24; i++ ) 
  {
    //const TString TPCMON_STR( Form( "TPCMON_%i", i ) );
    sprintf(TPCMON_STR,"TPCMON_%i",i);
    tpcmon_samplesizedist[i] = (TH1*) cl->getHisto(TPCMON_STR,"sample_size_hist");
  }

  if (!gROOT->FindObject("TPCSampleSize"))
  {
    MakeCanvas("TPCSampleSize");
  }

  TC[4]->SetEditable(true);
  TC[4]->Clear("D");
  for( int i=0; i<24; i++ ) 
  {
    if( tpcmon_samplesizedist[i] )
    {
      TC[4]->cd(i+5);
      tpcmon_samplesizedist[i]->DrawCopy("");
      gPad->SetLogx(kTRUE);
      gPad->SetLogy(kTRUE);
    }
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
  runnostream << ThisName << "_SampleSize Run " << cl->RunNumber()
              << ", Time: " << ctime(&evttime);
  runstring = runnostream.str();
  transparent[3]->cd();
  PrintRun.DrawText(0.5, 0.91, runstring.c_str());

  TC[4]->Update();
  TC[4]->SetLogx();
  TC[4]->Show();
  TC[4]->SetEditable(false);

  return 0;
}

int TpcMonDraw::DrawTPCCheckSum(const std::string & /* what */)
{
  OnlMonClient *cl = OnlMonClient::instance();
  
  TH1 *tpcmon_checksumerror[24] = {nullptr};
  TH1 *tpcmon_checksums[24] = {nullptr};

  char TPCMON_STR[100];
  for( int i=0; i<24; i++ ) 
  {
    //const TString TPCMON_STR( Form( "TPCMON_%i", i ) );
    sprintf(TPCMON_STR,"TPCMON_%i",i);
    tpcmon_checksumerror[i] = (TH1*) cl->getHisto(TPCMON_STR,"Check_Sum_Error");
    tpcmon_checksums[i] = (TH1*) cl->getHisto(TPCMON_STR,"Check_Sums");
  }

  if (!gROOT->FindObject("TPCCheckSumError"))
  {
    MakeCanvas("TPCCheckSumError");
  }

  TC[5]->SetEditable(true);
  TC[5]->Clear("D");
  TC[5]->cd(1);

  for( int i=0; i<24; i++ )
  {
    if( tpcmon_checksumerror[i] && tpcmon_checksums[i] )
    {
      TC[5]->cd(i+5);

      tpcmon_checksumerror[i]->Divide(tpcmon_checksums[i]);
      tpcmon_checksumerror[i]->GetYaxis()->SetRangeUser(0.0001,1);
      tpcmon_checksumerror[i]->DrawCopy("HIST");
    }
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
  runnostream << ThisName << "_CheckSumError Run " << cl->RunNumber()
              << ", Time: " << ctime(&evttime);
  runstring = runnostream.str();
  transparent[4]->cd();
  PrintRun.DrawText(0.5, 0.91, runstring.c_str());

  TC[5]->Update();
  //TC[5]->SetLogy();
  TC[5]->Show();
  TC[5]->SetEditable(false);


  return 0;
}

int TpcMonDraw::DrawTPCADCSample(const std::string & /* what */)
{
  OnlMonClient *cl = OnlMonClient::instance();

  TH2 *tpcmon_ADCSAMPLE[24] = {nullptr};

  char TPCMON_STR[100];
  for( int i=0; i<24; i++ ) 
  {
    //const TString TPCMON_STR( Form( "TPCMON_%i", i ) );
    sprintf(TPCMON_STR,"TPCMON_%i",i);
    tpcmon_ADCSAMPLE[i] = (TH2*) cl->getHisto(TPCMON_STR,"ADC_vs_SAMPLE");
  }


  if (!gROOT->FindObject("TPCADCSample"))
  {
    MakeCanvas("TPCADCSample");
  }  

  TC[6]->SetEditable(true);
  TC[6]->Clear("D");

  for( int i=0; i<24; i++ )
  {
    if( tpcmon_ADCSAMPLE[i] )
    {
      TC[6]->cd(i+5);
      gStyle->SetPalette(57); //kBird CVD friendly
      gPad->SetLogz(kTRUE);
      tpcmon_ADCSAMPLE[i] -> DrawCopy("colz");
    }
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
  runnostream << ThisName << "_ADC_vs_SAMPLE Run " << cl->RunNumber()
              << ", Time: " << ctime(&evttime);
  runstring = runnostream.str();
  transparent[5]->cd();
  PrintRun.DrawText(0.5, 0.91, runstring.c_str());


  TC[6]->Update();
  TC[6]->Show();
  TC[6]->SetEditable(false);

  return 0;
}

int TpcMonDraw::DrawTPCMaxADCModule(const std::string & /* what */)
{
  OnlMonClient *cl = OnlMonClient::instance();

  TH2 *tpcmon_MAXADC_MODULE[24] = {nullptr};

  char TPCMON_STR[100];
  for( int i=0; i<24; i++ ) 
  {
    //const TString TPCMON_STR( Form( "TPCMON_%i", i ) );
    sprintf(TPCMON_STR,"TPCMON_%i",i);
    tpcmon_MAXADC_MODULE[i] = (TH2*) cl->getHisto(TPCMON_STR,"MAXADC");
  }


  if (!gROOT->FindObject("TPCMaxADCModule"))
  {
    MakeCanvas("TPCMaxADCModule");
  }  

  TC[7]->SetEditable(true);
  TC[7]->Clear("D");

  for( int i=0; i<24; i++ )
  {
    if( tpcmon_MAXADC_MODULE[i] )
    {
      TC[7]->cd(i+5);
      gStyle->SetPalette(57); //kBird CVD friendly
      gPad->SetLogz(kTRUE);
      tpcmon_MAXADC_MODULE[i] -> DrawCopy("colz");
    }
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
  runnostream << ThisName << "_MAXADC2D Run " << cl->RunNumber()
              << ", Time: " << ctime(&evttime);
  runstring = runnostream.str();
  transparent[6]->cd();
  PrintRun.DrawText(0.5, 0.91, runstring.c_str());



  TC[7]->Update();
  TC[7]->Show();
  TC[7]->SetEditable(false);

  return 0;
}

int TpcMonDraw::DrawTPCRawADC1D(const std::string & /* what */)
{
  //std::cout<<"Made it inside DrawTPCRawADC1D"<<std::endl;
  OnlMonClient *cl = OnlMonClient::instance();

  TH1 *tpcmon_RAWADC1D[24][3] = {nullptr};

  char TPCMON_STR[100];
  for( int i=0; i<24; i++ ) 
  {
    //const TString TPCMON_STR( Form( "TPCMON_%i", i ) );
    sprintf(TPCMON_STR,"TPCMON_%i",i);
    tpcmon_RAWADC1D[i][0] = (TH1*) cl->getHisto(TPCMON_STR,"RAWADC_1D_R1");
    tpcmon_RAWADC1D[i][1] = (TH1*) cl->getHisto(TPCMON_STR,"RAWADC_1D_R2");
    tpcmon_RAWADC1D[i][2] = (TH1*) cl->getHisto(TPCMON_STR,"RAWADC_1D_R3");
  }


  if (!gROOT->FindObject("TPCRawADC1D"))
  {
    MakeCanvas("TPCRawADC1D");
  }  

  TC[8]->SetEditable(true);
  TC[8]->Clear("D");

  auto legend = new TLegend(0.7,0.65,0.98,0.95);
  bool draw_leg = 0;

  for( int i=0; i<24; i++ )
  {
    TC[8]->cd(i+5);
    for( int j = 2; j>-1; j-- )
    {
      if( tpcmon_RAWADC1D[i][j] )
      {
        if(j == 2){tpcmon_RAWADC1D[i][j] -> DrawCopy("HIST");}
        else      {tpcmon_RAWADC1D[i][j] -> DrawCopy("HISTsame");} //assumes that R3 will always exist and is most entries
      }
    }
    gPad->Update();
    gPad->SetLogy(kTRUE);  

    if(draw_leg == 0 && tpcmon_RAWADC1D[i][0] && tpcmon_RAWADC1D[i][1] && tpcmon_RAWADC1D[i][2]) //if you have not drawn the legend yet, draw it BUT ONLY ONCE
    {
      legend->AddEntry(tpcmon_RAWADC1D[i][0], "R1");
      legend->AddEntry(tpcmon_RAWADC1D[i][1], "R2");
      legend->AddEntry(tpcmon_RAWADC1D[i][2], "R3");
      TC[8]->cd(i+5);
      legend->Draw();
      draw_leg = 1; //don't draw it again
    }
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
  runnostream << ThisName << "_RAWADC Run " << cl->RunNumber()
              << ", Time: " << ctime(&evttime);
  runstring = runnostream.str();
  transparent[7]->cd();
  PrintRun.DrawText(0.5, 0.91, runstring.c_str());

  TC[8]->Update();
  TC[8]->Show();
  TC[8]->SetEditable(false);

  return 0;
}

int TpcMonDraw::DrawTPCMaxADC1D(const std::string & /* what */)
{
  //std::cout<<"Made it inside DrawTPCMaxADC1D"<<std::endl;
  OnlMonClient *cl = OnlMonClient::instance();

  TH1 *tpcmon_MAXADC1D[24][3] = {nullptr};

  char TPCMON_STR[100];
  for( int i=0; i<24; i++ ) 
  {
    //const TString TPCMON_STR( Form( "TPCMON_%i", i ) );
    sprintf(TPCMON_STR,"TPCMON_%i",i);
    tpcmon_MAXADC1D[i][0] = (TH1*) cl->getHisto(TPCMON_STR,"MAXADC_1D_R1");
    tpcmon_MAXADC1D[i][1] = (TH1*) cl->getHisto(TPCMON_STR,"MAXADC_1D_R2");
    tpcmon_MAXADC1D[i][2] = (TH1*) cl->getHisto(TPCMON_STR,"MAXADC_1D_R3");
  }


  if (!gROOT->FindObject("TPCMaxADC1D"))
  {
    MakeCanvas("TPCMaxADC1D");
  }  

  TC[9]->SetEditable(true);
  TC[9]->Clear("D");

  auto legend = new TLegend(0.7,0.65,0.98,0.95);
  bool draw_leg = 0;

  for( int i=0; i<24; i++ )
  {
    TC[9]->cd(i+5);
    for( int j = 2; j>-1; j-- )
    {
      if( tpcmon_MAXADC1D[i][j] )
      {
        if(j == 2){tpcmon_MAXADC1D[i][j] -> DrawCopy("HIST");}
        else      {tpcmon_MAXADC1D[i][j] -> DrawCopy("HISTsame");} //assumes that R3 will always exist and is max
      }
    }
    gPad->Update();
    gPad->SetLogy(kTRUE);
     
    if(draw_leg == 0 && tpcmon_MAXADC1D[i][0] && tpcmon_MAXADC1D[i][1] && tpcmon_MAXADC1D[i][2]) //if you have not drawn the legend yet, draw it BUT ONLY ONCE
    {
      TC[9]->cd(i+5);
      legend->AddEntry(tpcmon_MAXADC1D[i][0], "R1");
      legend->AddEntry(tpcmon_MAXADC1D[i][1], "R2");
      legend->AddEntry(tpcmon_MAXADC1D[i][2], "R3");
      legend->Draw();
      draw_leg = 1; //don't draw it again
    }
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
  runnostream << ThisName << "_LOCALMAXADC-PEDESTAL Run " << cl->RunNumber()
              << ", Time: " << ctime(&evttime);
  runstring = runnostream.str();
  transparent[8]->cd();
  PrintRun.DrawText(0.5, 0.91, runstring.c_str());


  TC[9]->Update();
  TC[9]->Show();
  TC[9]->SetEditable(false);
  
  return 0;
}
int TpcMonDraw::DrawTPCXYclusters(const std::string & /* what */)
{
  OnlMonClient *cl = OnlMonClient::instance();

  TH2 *tpcmon_NSTPC_clusXY[24][3] = {nullptr};
  TH2 *tpcmon_SSTPC_clusXY[24][3] = {nullptr};

  dummy_his1_XY = new TH2F("dummy_his1_XY", "(ADC-Pedestal) > 20 ADC North Side, WEIGHTED", 400, -800, 800, 400, -800, 800); //dummy histos for titles
  dummy_his2_XY = new TH2F("dummy_his2_XY", "(ADC-Pedestal) > 20 ADC South Side, WEIGHTED", 400, -800, 800, 400, -800, 800);

  dummy_his1_XY->SetXTitle("X [mm]");
  dummy_his1_XY->SetYTitle("Y [mm]");
  dummy_his1_XY->GetYaxis()->SetTitleSize(0.02);

  dummy_his2_XY->SetXTitle("-X [mm]"); //SS x is flipped from global coordinates
  dummy_his2_XY->SetYTitle("Y [mm]");
  dummy_his2_XY->GetYaxis()->SetTitleSize(0.02);

  //the lines are for the sector boundaries
  Double_t sec_gap_inner = (2*M_PI - 0.5024*12.0)/12.0;

  Double_t sec_gap_outer = (2*M_PI - 0.5097*12.0)/12.0;

  Double_t sec_gap = (sec_gap_inner + sec_gap_outer)/2.0;

  Double_t sec_phi = (0.5024 + 0.5097)/2.0;

  TLine *lines[12];

  for(int ln=0;ln<12;ln++)
  {
    lines[ln] = new TLine(311.05*cos((sec_phi+sec_gap)/2.0+ln*(sec_phi+sec_gap)),311.05*sin((sec_phi+sec_gap)/2.0+ln*(sec_phi+sec_gap)),759.11*cos((sec_phi+sec_gap)/2.0+ln*(sec_phi+sec_gap)),759.11*sin((sec_phi+sec_gap)/2.0+ln*(sec_phi+sec_gap)));
  }

  TEllipse *e1 = new TEllipse(0.0,0.0,311.05,311.05);
  TEllipse *e2 = new TEllipse(0.0,0.0,(402.49+411.53)/2.0,(402.49+411.53)/2.0);
  TEllipse *e3 = new TEllipse(0.0,0.0,(583.67+574.75)/2.0,(583.67+574.75)/2.0);
  TEllipse *e4 = new TEllipse(0.0,0.0,759.11,759.11);
  //__________________

  char TPCMON_STR[100];
  for( int i=0; i<24; i++ ) 
  {
    //const TString TPCMON_STR( Form( "TPCMON_%i", i ) );
    sprintf(TPCMON_STR,"TPCMON_%i",i);
    tpcmon_NSTPC_clusXY[i][0] = (TH2*) cl->getHisto(TPCMON_STR,"NorthSideADC_clusterXY_R1");
    tpcmon_NSTPC_clusXY[i][1] = (TH2*) cl->getHisto(TPCMON_STR,"NorthSideADC_clusterXY_R2");
    tpcmon_NSTPC_clusXY[i][2] = (TH2*) cl->getHisto(TPCMON_STR,"NorthSideADC_clusterXY_R3");

    tpcmon_SSTPC_clusXY[i][0] = (TH2*) cl->getHisto(TPCMON_STR,"SouthSideADC_clusterXY_R1");
    tpcmon_SSTPC_clusXY[i][1] = (TH2*) cl->getHisto(TPCMON_STR,"SouthSideADC_clusterXY_R2");
    tpcmon_SSTPC_clusXY[i][2] = (TH2*) cl->getHisto(TPCMON_STR,"SouthSideADC_clusterXY_R3");
  }

  if (!gROOT->FindObject("TPCClusterXY"))
  {
    MakeCanvas("TPCClusterXY");
  }  

  TC[10]->SetEditable(true);
  TC[10]->Clear("D");

  TText PrintRun;
  PrintRun.SetTextFont(62);
  PrintRun.SetTextSize(0.04);
  PrintRun.SetNDC();          // set to normalized coordinates
  PrintRun.SetTextAlign(23);  // center/top alignment
  std::ostringstream runnostream;
  std::string runstring;
  time_t evttime = cl->EventTime("CURRENT");
  // fill run number and event time into string
  runnostream << ThisName << "_ADC-Pedestal>20 ADC WEIGHTED, Run" << cl->RunNumber()
              << ", Time: " << ctime(&evttime);
  runstring = runnostream.str();
  transparent[9]->cd();
  PrintRun.DrawText(0.5, 0.91, runstring.c_str());

  TC[10]->cd(1);
  gStyle->SetOptStat(kFALSE);
  gPad->SetTopMargin(0.15);
  gPad->SetLogz(kTRUE);
  dummy_his1_XY->Draw("colzsame");

  float NS_max = 0;
  for( int i=0; i<12; i++ )
  {
    for( int j=0; j<3; j++ )
    {
      if( tpcmon_NSTPC_clusXY[i][j] )
      {
        TC[10]->cd(1);
        tpcmon_NSTPC_clusXY[i][j] -> Draw("colzsame");
        //gStyle->SetLogz(kTRUE);
        if ( tpcmon_NSTPC_clusXY[i][j]->GetBinContent(tpcmon_NSTPC_clusXY[i][j]->GetMaximumBin()) > NS_max)
        {
          NS_max = tpcmon_NSTPC_clusXY[i][j]->GetBinContent(tpcmon_NSTPC_clusXY[i][j]->GetMaximumBin());
          dummy_his1_XY->SetMaximum( NS_max );
        }
        gStyle->SetPalette(57); //kBird CVD friendly
      }

    }
  }
  TC[10]->cd(1);
  e1->SetFillStyle(0);
  e2->SetFillStyle(0);
  e3->SetFillStyle(0);
  e4->SetFillStyle(0);

  e1->Draw("same");
  e2->Draw("same");
  e3->Draw("same");
  e4->Draw("same");
  for(int ln2=0;ln2<12;ln2++)
  {
    lines[ln2]->Draw("same"); 
  }
  TC[10]->Update();

  TC[10]->cd(2);
  gStyle->SetOptStat(kFALSE);
  gPad->SetTopMargin(0.15);
  gPad->SetLogz(kTRUE);
  dummy_his2_XY->Draw("colzsame");

  float SS_max = 0;
  for( int i=0; i<12; i++ )
  {
    for( int j=0; j<3; j++ )
    {
      if( tpcmon_SSTPC_clusXY[i+12][j] )
      {
        //std::cout<<"South Side Custer XY i: "<< i+12 <<", j: "<<j<<std::endl;
        TC[10]->cd(2);
        tpcmon_SSTPC_clusXY[i+12][j] -> Draw("colzsame");
        //gStyle->SetLogz(kTRUE);
        if ( tpcmon_SSTPC_clusXY[i+12][j]->GetBinContent(tpcmon_SSTPC_clusXY[i+12][j]->GetMaximumBin()) > SS_max)
        {
          SS_max = tpcmon_SSTPC_clusXY[i+12][j]->GetBinContent(tpcmon_SSTPC_clusXY[i+12][j]->GetMaximumBin());
          dummy_his2_XY->SetMaximum( SS_max );
        }
        gStyle->SetPalette(57); //kBird CVD friendly
      }
    }

  }
  TC[10]->cd(2);
  e1->SetFillStyle(0);
  e2->SetFillStyle(0);
  e3->SetFillStyle(0);
  e4->SetFillStyle(0);

  e1->Draw("same");
  e2->Draw("same");
  e3->Draw("same");
  e4->Draw("same");
  for(int ln2=0;ln2<12;ln2++)
  {
    lines[ln2]->Draw("same"); 
  }

  TC[10]->Update();
  TC[10]->Show();
  TC[10]->SetEditable(false);

  return 0;
}

int TpcMonDraw::DrawTPCXYclusters_unweighted(const std::string & /* what */)
{
  OnlMonClient *cl = OnlMonClient::instance();

  TH2 *tpcmon_NSTPC_clusXY[24][3] = {nullptr};
  TH2 *tpcmon_SSTPC_clusXY[24][3] = {nullptr};

  dummy_his1_XY_unw = new TH2F("dummy_his1_XY_unw", "(ADC-Pedestal) > 20 ADC North Side, UNWEIGHTED", 400, -800, 800, 400, -800, 800); //dummy histos for titles
  dummy_his2_XY_unw = new TH2F("dummy_his2_XY_unw", "(ADC-Pedestal) > 20 ADC South Side, UNWEIGHTED", 400, -800, 800, 400, -800, 800);

  dummy_his1_XY_unw->SetXTitle("X [mm]");
  dummy_his1_XY_unw->SetYTitle("Y [mm]");
  dummy_his1_XY_unw->GetYaxis()->SetTitleSize(0.02);

  dummy_his2_XY_unw->SetXTitle("-X [mm]"); //SS x is flipped from global coordinates
  dummy_his2_XY_unw->SetYTitle("Y [mm]");
  dummy_his2_XY_unw->GetYaxis()->SetTitleSize(0.02);

  //the lines are for the sector boundaries
  Double_t sec_gap_inner = (2*M_PI - 0.5024*12.0)/12.0;

  Double_t sec_gap_outer = (2*M_PI - 0.5097*12.0)/12.0;

  Double_t sec_gap = (sec_gap_inner + sec_gap_outer)/2.0;

  Double_t sec_phi = (0.5024 + 0.5097)/2.0;

  TLine *lines[12];

  for(int ln=0;ln<12;ln++)
  {
    lines[ln] = new TLine(311.05*cos((sec_phi+sec_gap)/2.0+ln*(sec_phi+sec_gap)),311.05*sin((sec_phi+sec_gap)/2.0+ln*(sec_phi+sec_gap)),759.11*cos((sec_phi+sec_gap)/2.0+ln*(sec_phi+sec_gap)),759.11*sin((sec_phi+sec_gap)/2.0+ln*(sec_phi+sec_gap)));
  }

  TEllipse *e1 = new TEllipse(0.0,0.0,311.05,311.05);
  TEllipse *e2 = new TEllipse(0.0,0.0,(402.49+411.53)/2.0,(402.49+411.53)/2.0);
  TEllipse *e3 = new TEllipse(0.0,0.0,(583.67+574.75)/2.0,(583.67+574.75)/2.0);
  TEllipse *e4 = new TEllipse(0.0,0.0,759.11,759.11);
  //__________________

  char TPCMON_STR[100];
  for( int i=0; i<24; i++ ) 
  {
    //const TString TPCMON_STR( Form( "TPCMON_%i", i ) );
    sprintf(TPCMON_STR,"TPCMON_%i",i);
    tpcmon_NSTPC_clusXY[i][0] = (TH2*) cl->getHisto(TPCMON_STR,"NorthSideADC_clusterXY_R1_unw");
    tpcmon_NSTPC_clusXY[i][1] = (TH2*) cl->getHisto(TPCMON_STR,"NorthSideADC_clusterXY_R2_unw");
    tpcmon_NSTPC_clusXY[i][2] = (TH2*) cl->getHisto(TPCMON_STR,"NorthSideADC_clusterXY_R3_unw");

    tpcmon_SSTPC_clusXY[i][0] = (TH2*) cl->getHisto(TPCMON_STR,"SouthSideADC_clusterXY_R1_unw");
    tpcmon_SSTPC_clusXY[i][1] = (TH2*) cl->getHisto(TPCMON_STR,"SouthSideADC_clusterXY_R2_unw");
    tpcmon_SSTPC_clusXY[i][2] = (TH2*) cl->getHisto(TPCMON_STR,"SouthSideADC_clusterXY_R3_unw");
  }

  if (!gROOT->FindObject("TPCClusterXY_unw"))
  {
    MakeCanvas("TPCClusterXY_unw");
  }  

  TC[11]->SetEditable(true);
  TC[11]->Clear("D");

  TText PrintRun;
  PrintRun.SetTextFont(62);
  PrintRun.SetTextSize(0.04);
  PrintRun.SetNDC();          // set to normalized coordinates
  PrintRun.SetTextAlign(23);  // center/top alignment
  std::ostringstream runnostream;
  std::string runstring;
  time_t evttime = cl->EventTime("CURRENT");
  // fill run number and event time into string
  runnostream << ThisName << "_ADC-Pedestal>20, UNWEIGHTED, Run " << cl->RunNumber()
              << ", Time: " << ctime(&evttime);
  runstring = runnostream.str();
  transparent[10]->cd();
  PrintRun.DrawText(0.5, 0.91, runstring.c_str());

  TC[11]->cd(1);
  gStyle->SetOptStat(kFALSE);
  gPad->SetTopMargin(0.15);
  //gPad->SetLogz(kTRUE);
  dummy_his1_XY_unw->Draw("colzsame");

  float NS_max = 0;
  for( int i=0; i<12; i++ )
  {
    for( int j=0; j<3; j++ )
    {
      if( tpcmon_NSTPC_clusXY[i][j] )
      {
        TC[11]->cd(1);
        tpcmon_NSTPC_clusXY[i][j] -> Draw("colzsame");
        //gStyle->SetLogz(kTRUE);
        if ( tpcmon_NSTPC_clusXY[i][j]->GetBinContent(tpcmon_NSTPC_clusXY[i][j]->GetMaximumBin()) > NS_max)
        {
          NS_max = tpcmon_NSTPC_clusXY[i][j]->GetBinContent(tpcmon_NSTPC_clusXY[i][j]->GetMaximumBin());
          dummy_his1_XY_unw->SetMaximum( NS_max );
        }
        gStyle->SetPalette(57); //kBird CVD friendly
      }
    }

  }
  TC[11]->cd(1);
  e1->SetFillStyle(0);
  e2->SetFillStyle(0);
  e3->SetFillStyle(0);
  e4->SetFillStyle(0);

  e1->Draw("same");
  e2->Draw("same");
  e3->Draw("same");
  e4->Draw("same");
  for(int ln2=0;ln2<12;ln2++)
  {
    lines[ln2]->Draw("same"); 
  }
  TC[11]->Update();

  TC[11]->cd(2);
  gStyle->SetOptStat(kFALSE);
  gPad->SetTopMargin(0.15);
  dummy_his2_XY_unw->Draw("colzsame");
  //gPad->SetLogz(kTRUE);

  float SS_max = 0;
  for( int i=0; i<12; i++ )
  {
    for( int j=0; j<3; j++ )
    {
      if( tpcmon_SSTPC_clusXY[i+12][j] )
      {
        //std::cout<<"South Side Custer XY i: "<< i+12 <<", j: "<<j<<std::endl;
        TC[11]->cd(2);
        tpcmon_SSTPC_clusXY[i+12][j] -> Draw("colzsame");
        //gStyle->SetLogz(kTRUE);
        if ( tpcmon_SSTPC_clusXY[i+12][j]->GetBinContent(tpcmon_SSTPC_clusXY[i+12][j]->GetMaximumBin()) > SS_max)
        {
          SS_max = tpcmon_SSTPC_clusXY[i+12][j]->GetBinContent(tpcmon_SSTPC_clusXY[i+12][j]->GetMaximumBin());
          dummy_his2_XY_unw->SetMaximum( SS_max );
        }
        gStyle->SetPalette(57); //kBird CVD friendly
      }
    }

  }
  TC[11]->cd(2);
  e1->SetFillStyle(0);
  e2->SetFillStyle(0);
  e3->SetFillStyle(0);
  e4->SetFillStyle(0);

  e1->Draw("same");
  e2->Draw("same");
  e3->Draw("same");
  e4->Draw("same");
  for(int ln2=0;ln2<12;ln2++)
  {
    lines[ln2]->Draw("same"); 
  }
  TC[11]->Update();

  TC[11]->Show();
  TC[11]->SetEditable(false);

  return 0;
}

int TpcMonDraw::DrawTPCADCSampleLarge(const std::string & /* what */)
{
  OnlMonClient *cl = OnlMonClient::instance();

  TH2 *tpcmon_ADCSAMPLE_large[24] = {nullptr};

  char TPCMON_STR[100];
  for( int i=0; i<24; i++ ) 
  {
    //const TString TPCMON_STR( Form( "TPCMON_%i", i ) );
    sprintf(TPCMON_STR,"TPCMON_%i",i);
    tpcmon_ADCSAMPLE_large[i] = (TH2*) cl->getHisto(TPCMON_STR,"ADC_vs_SAMPLE_large");
  }


  if (!gROOT->FindObject("TPCADCSamplelarge"))
  {
    MakeCanvas("TPCADCSamplelarge");
  }  

  TC[12]->SetEditable(true);
  TC[12]->Clear("D");

  for( int i=0; i<24; i++ )
  {
    if( tpcmon_ADCSAMPLE_large[i] )
    {
      TC[12]->cd(i+5);
      gStyle->SetPalette(57); //kBird CVD friendly
      gPad->SetLogz(kTRUE);
      tpcmon_ADCSAMPLE_large[i] -> DrawCopy("colz");
    }
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
  runnostream << ThisName << "_ADC_vs_SAMPLE_large Run " << cl->RunNumber()
              << ", Time: " << ctime(&evttime);
  runstring = runnostream.str();
  transparent[11]->cd();
  PrintRun.DrawText(0.5, 0.91, runstring.c_str());


  TC[12]->Update();
  TC[12]->Show();
  TC[12]->SetEditable(false);

  return 0;
}



int TpcMonDraw::DrawTPCZYclusters(const std::string & /* what */)
{

  OnlMonClient *cl = OnlMonClient::instance();

  TH2 *tpcmon_NSTPC_clusZY[24] = {nullptr};
  TH2 *tpcmon_SSTPC_clusZY[24] = {nullptr};

  dummy_his1_ZY = new TH2F("dummy_his1_ZY", "(ADC-Pedestal) > 20 ADC, WEIGHTED", 515, -1030, 1030, 400, -800, 800); //dummy histos for titles
  dummy_his1_ZY->SetXTitle("Z [mm]");
  dummy_his1_ZY->SetYTitle("Y [mm]");

  char TPCMON_STR[100];
  for( int i=0; i<24; i++ ) 
  {
    //const TString TPCMON_STR( Form( "TPCMON_%i", i ) );
    sprintf(TPCMON_STR,"TPCMON_%i",i);
    tpcmon_NSTPC_clusZY[i] = (TH2*) cl->getHisto(TPCMON_STR,"NorthSideADC_clusterZY");

    tpcmon_SSTPC_clusZY[i] = (TH2*) cl->getHisto(TPCMON_STR,"SouthSideADC_clusterZY");
  }

  if (!gROOT->FindObject("TPCClusterZY"))
  {
    MakeCanvas("TPCClusterZY");
  }  

  TC[13]->SetEditable(true);
  TC[13]->Clear("D");

  TText PrintRun;
  PrintRun.SetTextFont(62);
  PrintRun.SetTextSize(0.04);
  PrintRun.SetNDC();          // set to normalized coordinates
  PrintRun.SetTextAlign(23);  // center/top alignment
  std::ostringstream runnostream;
  std::string runstring;
  time_t evttime = cl->EventTime("CURRENT");
  // fill run number and event time into string
  runnostream << ThisName << "_ADC-Pedestal>20 ADC, WEIGHTED, Run " << cl->RunNumber()
              << ", Time: " << ctime(&evttime);
  runstring = runnostream.str();
  transparent[12]->cd();
  PrintRun.DrawText(0.5, 0.91, runstring.c_str());

  TC[13]->cd(1);
  gStyle->SetOptStat(kFALSE);
  gPad->SetTopMargin(0.15);
  gPad->SetLogz(kTRUE);
  dummy_his1_ZY->Draw("colzsame");


  float max = 0;
  for( int i=0; i<12; i++ )
  {
    if( tpcmon_NSTPC_clusZY[i] )
    {
      TC[13]->cd(1);
      tpcmon_NSTPC_clusZY[i] -> Draw("colzsame");
      if ( tpcmon_NSTPC_clusZY[i]->GetBinContent(tpcmon_NSTPC_clusZY[i]->GetMaximumBin()) > max ) 
      {
        max = tpcmon_NSTPC_clusZY[i]->GetBinContent(tpcmon_NSTPC_clusZY[i]->GetMaximumBin());
        dummy_his1_ZY->SetMaximum( max );
      }
      gStyle->SetPalette(57); //kBird CVD friendly
    }

  }
  TC[13]->Update();

  for( int i=0; i<12; i++ )
  {
    if( tpcmon_SSTPC_clusZY[i+12] )
    {
      TC[13]->cd(1);
      tpcmon_SSTPC_clusZY[i+12] -> Draw("colzsame");
      if ( tpcmon_SSTPC_clusZY[i+12]->GetBinContent(tpcmon_SSTPC_clusZY[i+12]->GetMaximumBin()) > max ) 
      {
        max = tpcmon_SSTPC_clusZY[i+12]->GetBinContent(tpcmon_SSTPC_clusZY[i+12]->GetMaximumBin());
        dummy_his1_ZY->SetMaximum( max );
      }
      //gStyle->SetLogz(kTRUE);
      gStyle->SetPalette(57); //kBird CVD friendly
    }

  }

  TC[13]->Update();
  TC[13]->Show();
  TC[13]->SetEditable(false);



  return 0;
}


int TpcMonDraw::DrawTPCZYclusters_unweighted(const std::string & /* what */)
{


  OnlMonClient *cl = OnlMonClient::instance();

  TH2 *tpcmon_NSTPC_clusZY_unw[24] = {nullptr};
  TH2 *tpcmon_SSTPC_clusZY_unw[24] = {nullptr};

  dummy_his1_ZY_unw = new TH2F("dummy_his1_ZY_unw", "(ADC-Pedestal) > 20 ADC, UNWEIGHTED", 515, -1030, 1030, 400, -800, 800); //dummy histos for titles
  dummy_his1_ZY_unw->SetXTitle("Z [mm]");
  dummy_his1_ZY_unw->SetYTitle("Y [mm]");

  char TPCMON_STR[100];
  for( int i=0; i<24; i++ ) 
  {
    //const TString TPCMON_STR( Form( "TPCMON_%i", i ) );
    sprintf(TPCMON_STR,"TPCMON_%i",i);
    tpcmon_NSTPC_clusZY_unw[i] = (TH2*) cl->getHisto(TPCMON_STR,"NorthSideADC_clusterZY_unw");

    tpcmon_SSTPC_clusZY_unw[i] = (TH2*) cl->getHisto(TPCMON_STR,"SouthSideADC_clusterZY_unw");
  }

  if (!gROOT->FindObject("TPCClusterZY_unw"))
  {
    MakeCanvas("TPCClusterZY_unw");
  }  

  TC[14]->SetEditable(true);
  TC[14]->Clear("D");

  TText PrintRun;
  PrintRun.SetTextFont(62);
  PrintRun.SetTextSize(0.04);
  PrintRun.SetNDC();          // set to normalized coordinates
  PrintRun.SetTextAlign(23);  // center/top alignment
  std::ostringstream runnostream;
  std::string runstring;
  time_t evttime = cl->EventTime("CURRENT");
  // fill run number and event time into string
  runnostream << ThisName << "_ADC-Pedestal>20 ADC, UNWEIGHTED, Run " << cl->RunNumber()
              << ", Time: " << ctime(&evttime);
  runstring = runnostream.str();
  transparent[13]->cd();
  PrintRun.DrawText(0.5, 0.91, runstring.c_str());

  TC[14]->cd(1);
  gStyle->SetOptStat(kFALSE);
  gPad->SetTopMargin(0.15);
  gPad->SetLogz(kTRUE);
  dummy_his1_ZY_unw->Draw("colzsame");

  float max = 0;
  for( int i=0; i<12; i++ )
  {
    if( tpcmon_NSTPC_clusZY_unw[i] )
    {
      TC[14]->cd(1);
      tpcmon_NSTPC_clusZY_unw[i] -> Draw("colzsame");
      if ( tpcmon_NSTPC_clusZY_unw[i]->GetBinContent(tpcmon_NSTPC_clusZY_unw[i]->GetMaximumBin()) > max ) 
      {
        max = tpcmon_NSTPC_clusZY_unw[i]->GetBinContent(tpcmon_NSTPC_clusZY_unw[i]->GetMaximumBin());
        dummy_his1_ZY->SetMaximum( max );
      }
      //gStyle->SetLogz(kTRUE);
      gStyle->SetPalette(57); //kBird CVD friendly
    }

  }
  TC[14]->Update();

  for( int i=0; i<12; i++ )
  {
    if( tpcmon_SSTPC_clusZY_unw[i+12] )
    {
      TC[14]->cd(1);
      tpcmon_SSTPC_clusZY_unw[i+12] -> Draw("colzsame");
      if ( tpcmon_SSTPC_clusZY_unw[i+12]->GetBinContent(tpcmon_SSTPC_clusZY_unw[i+12]->GetMaximumBin()) > max ) 
      {
        max = tpcmon_SSTPC_clusZY_unw[i+12]->GetBinContent(tpcmon_SSTPC_clusZY_unw[i+12]->GetMaximumBin());
        dummy_his1_ZY->SetMaximum( max );
      }
      //gStyle->SetLogz(kTRUE);
      gStyle->SetPalette(57); //kBird CVD friendly
    }

  }
  TC[14]->Update();

  TC[14]->Show();
  TC[14]->SetEditable(false);



  return 0;
}


int TpcMonDraw::SavePlot(const std::string &what, const std::string &type)
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

int TpcMonDraw::MakeHtml(const std::string &what)
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

time_t TpcMonDraw::getTime()
{
  OnlMonClient *cl = OnlMonClient::instance();
  time_t currtime = 0;
  int i = 0;
  while (currtime == 0 && i <= 23)
  {
    std::string servername = "TPCMON_" + std::to_string(i);
    currtime = cl->EventTime(servername,"CURRENT");
    i++;
  }
  return currtime;
}
