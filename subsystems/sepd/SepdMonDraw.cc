#include "SepdMonDraw.h"

#include <onlmon/OnlMonClient.h>

#include <TAxis.h>  // for TAxis
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TProfile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TText.h>

#include <cstring>  // for memset
#include <ctime>
#include <fstream>
#include <iostream>  // for operator<<, basic_ostream, basic_os...
#include <sstream>
#include <vector>    // for vector

SepdMonDraw::SepdMonDraw(const std::string &name)
  : OnlMonDraw(name)
{
  return;
}

int SepdMonDraw::Init()
{
  return 0;
}

int SepdMonDraw::MakeCanvas(const std::string &name)
{
  OnlMonClient *cl = OnlMonClient::instance();
  int xsize = cl->GetDisplaySizeX();
  int ysize = cl->GetDisplaySizeY();
  if (name == "SepdMon1")
  {
    // xpos (-1) negative: do not draw menu bar
    TC[0] = new TCanvas(name.c_str(), "SepdMon Example Monitor", 1200, 600);
    // root is pathetic, whenever a new TCanvas is created root piles up
    // 6kb worth of X11 events which need to be cleared with
    // gSystem->ProcessEvents(), otherwise your process will grow and
    // grow and grow but will not show a definitely lost memory leak
    gSystem->ProcessEvents();
    Pad[0] = new TPad("sepdpad0", "Left", 0., 0., 0.5, 1);
    Pad[1] = new TPad("sepdpad1", "Right", 0.5, 0., 1, 1);
    Pad[0]->Draw();
    Pad[1]->Draw();
    // this one is used to plot the run number on the canvas
    transparent[0] = new TPad("transparent0", "this does not show", 0, 0, 1, 1);
    transparent[0]->SetFillStyle(4000);
    transparent[0]->Draw();
    TC[0]->SetEditable(false);
  }
  else if (name == "SepdMon2")
  {
    // xpos negative: do not draw menu bar
    TC[1] = new TCanvas(name.c_str(), "SepdMon2 Example Monitor", 1200, 600);
    gSystem->ProcessEvents();
    Pad[2] = new TPad("sepdpad2", "Left", 0., 0., 0.5, 1);
    Pad[3] = new TPad("sepdpad3", "Right", 0.5, 0., 1, 1);
    Pad[2]->Draw();
    Pad[3]->Draw();
    // this one is used to plot the run number on the canvas
    transparent[1] = new TPad("transparent1", "this does not show", 0, 0, 1, 1);
    transparent[1]->SetFillStyle(4000);
    transparent[1]->Draw();
    TC[1]->SetEditable(false);
  }
  /*
  else if (name == "SepdMon3")
  {
  }
  */
  else if (name == "SepdMon4")
  {
    TC[3] = new TCanvas(name.c_str(), "SepdMon3 Waveform Info", xsize / 3, 0, xsize / 3, ysize * 0.9);
    gSystem->ProcessEvents();
    Pad[6] = new TPad("hcalpad6", "who needs this?", 0.0, 0.6, 1.0, 0.95, 0);
    Pad[7] = new TPad("hcalpad7", "who needs this?", 0.0, 0.3, 1.0, 0.6, 0);
    Pad[8] = new TPad("hcalpad8", "who needs this?", 0.0, 0.0, 1.0, 0.3, 0);
    Pad[6]->Draw();
    Pad[7]->Draw();
    Pad[8]->Draw();
    // this one is used to plot the run number on the canvas
    transparent[3] = new TPad("transparent3", "this does not show", 0, 0, 1, 1);
    transparent[3]->SetFillStyle(4000);
    transparent[3]->Draw();
    TC[3]->SetEditable(0);
  }
  else if (name == "SepdMon5")
  {
    // xpos negative: do not draw menu bar
    TC[4] = new TCanvas(name.c_str(), "HCAL Packet Information", 2 * xsize / 3, 0, 2 * xsize / 3, ysize * 0.9);
    gSystem->ProcessEvents();
    Pad[10] = new TPad("pad10", "packet event check", 0.0, 0.6, 1.0 / 2, 0.95, 0);
    Pad[11] = new TPad("pad11", "packet size", 0.0, 0.3, 1.0 / 2, 0.6, 0);
    Pad[12] = new TPad("pad12", "packet channels", 0.0, 0.0, 1.0 / 2, 0.3, 0);
    Pad[13] = new TPad("pad13", "event number offset", 0.5, 0.6, 1.0, 0.95, 0);
    // pad 14 and 15 side by side from left to right for correlation
    // Pad[14] = new TPad("pad14", "correlation0", 0.5, 0.3, 0.75, 0.6, 0);
    // Pad[15] = new TPad("pad15", "correlation1", 0.75, 0.3, 1.0, 0.6, 0);

    Pad[10]->Draw();
    Pad[11]->Draw();
    Pad[12]->Draw();
    Pad[13]->Draw();
    //Pad[14]->Draw();
    //Pad[15]->Draw();
    //  this one is used to plot the run number on the canvas
    transparent[4] = new TPad("transparent1", "this does not show", 0, 0, 1., 1);
    transparent[4]->SetFillStyle(4000);
    transparent[4]->Draw();

    // packet warnings
    warning[1] = new TPad("warning1", "packet warnings", 0.5, 0, 1, 0.2);
    warning[1]->SetFillStyle(4000);
    warning[1]->Draw();
    TC[4]->SetEditable(0);
  }
  return 0;
}

int SepdMonDraw::Draw(const std::string &what)
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
  /*
  if (what == "ALL" || what == "THIRD")
  {
    iret += DrawThird(what);
    idraw++;
  }
  */
  if (what == "ALL" || what == "FOURTH")
  {
    iret += DrawFourth(what);
    idraw++;
  }
  if (what == "ALL" || what == "FIFTH")
  {
    iret += DrawFifth(what);
    idraw++;
  }

  /*
    if (what == "ALL" || what == "HISTORY")
    {
      iret += DrawHistory(what);
      idraw++;
    }
    */
  if (!idraw)
  {
    std::cout << __PRETTY_FUNCTION__ << " Unimplemented Drawing option: " << what << std::endl;
    iret = -1;
  }
  return iret;
}

int SepdMonDraw::DrawFirst(const std::string & /* what */)
{
  OnlMonClient *cl = OnlMonClient::instance();
  TH2 *h_ADC0_s = (TH2 *) cl->getHisto("SEPDMON_0", "h_ADC0_s");
  TH2 *h_ADC0_n = (TH2 *) cl->getHisto("SEPDMON_0", "h_ADC0_n");
  TH2 *h_ADC_s = (TH2 *) cl->getHisto("SEPDMON_0", "h_ADC_s");
  TH2 *h_ADC_n = (TH2 *) cl->getHisto("SEPDMON_0", "h_ADC_n");
  TH2 *h_hits0_s = (TH2 *) cl->getHisto("SEPDMON_0", "h_hits0_s");
  TH2 *h_hits0_n = (TH2 *) cl->getHisto("SEPDMON_0", "h_hits0_n");
  TH2 *h_hits_s = (TH2 *) cl->getHisto("SEPDMON_0", "h_hits_s");
  TH2 *h_hits_n = (TH2 *) cl->getHisto("SEPDMON_0", "h_hits_n");
  time_t evttime = cl->EventTime("CURRENT");

  if (!gROOT->FindObject("SepdMon1"))
  {
    MakeCanvas("SepdMon1");
  }
  TC[0]->SetEditable(true);
  TC[0]->Clear("D");

  if (!h_ADC0_s)
  {
    DrawDeadServer(transparent[0]);
    TC[0]->SetEditable(false);
    return -1;
  }

  int nbinsx0 = h_ADC0_s->GetNbinsX();
  int nbinsy0 = h_ADC0_s->GetNbinsY();
  int nbinsx = h_ADC_s->GetNbinsX();
  int nbinsy = h_ADC_s->GetNbinsY();

  for (int ibx = 0; ibx < nbinsx0; ibx++)
  {
    for (int iby = 0; iby < nbinsy0; iby++)
    {
      double con = h_ADC0_s->GetBinContent(ibx + 1, iby + 1);
      double div = h_hits0_s->GetBinContent(ibx + 1, iby + 1);
      h_ADC0_s->SetBinContent(ibx + 1, iby + 1, con / div);
      con = h_ADC0_n->GetBinContent(ibx + 1, iby + 1);
      div = h_hits0_n->GetBinContent(ibx + 1, iby + 1);
      h_ADC0_n->SetBinContent(ibx + 1, iby + 1, con / div);
    }
  }

  for (int ibx = 0; ibx < nbinsx; ibx++)
  {
    for (int iby = 0; iby < nbinsy; iby++)
    {
      double con = h_ADC_s->GetBinContent(ibx + 1, iby + 1);
      double div = h_hits_s->GetBinContent(ibx + 1, iby + 1);
      h_ADC_s->SetBinContent(ibx + 1, iby + 1, con / div);
      con = h_ADC_n->GetBinContent(ibx + 1, iby + 1);
      div = h_hits_n->GetBinContent(ibx + 1, iby + 1);
      h_ADC_n->SetBinContent(ibx + 1, iby + 1, con / div);
    }
  }
  gStyle->SetOptStat(0);
  Pad[0]->cd();
  h_ADC0_s->Draw("COLZPOL");
  h_ADC_s->Draw("COLZPOL same");
  Pad[1]->cd();
  h_ADC0_n->Draw("COLZPOL");
  h_ADC_n->Draw("COLZPOL same");

  TText PrintRun;
  PrintRun.SetTextFont(62);
  PrintRun.SetTextSize(0.04);
  PrintRun.SetNDC();          // set to normalized coordinates
  PrintRun.SetTextAlign(23);  // center/top alignment
  std::ostringstream runnostream;
  std::string runstring;
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

int SepdMonDraw::DrawSecond(const std::string & /* what */)
{
  OnlMonClient *cl = OnlMonClient::instance();
  TH2 *h_hits0_s = (TH2 *) cl->getHisto("SEPDMON_0", "h_hits0_s");
  TH2 *h_hits0_n = (TH2 *) cl->getHisto("SEPDMON_0", "h_hits0_n");
  TH2 *h_hits_s = (TH2 *) cl->getHisto("SEPDMON_0", "h_hits_s");
  TH2 *h_hits_n = (TH2 *) cl->getHisto("SEPDMON_0", "h_hits_n");

  TH1 *h_event = cl->getHisto("SEPDMON_0", "h_event");
  time_t evttime = cl->EventTime("CURRENT");

  if (!gROOT->FindObject("SepdMon2"))
  {
    MakeCanvas("SepdMon2");
  }
  if (!h_hits0_s)
  {
    DrawDeadServer(transparent[0]);
    TC[0]->SetEditable(false);
    return -1;
  }
  int nevt = h_event->GetEntries();
  h_hits0_s->Scale(1. / nevt);
  h_hits_s->Scale(1. / nevt);
  h_hits0_n->Scale(1. / nevt);
  h_hits_n->Scale(1. / nevt);
  gStyle->SetOptStat(0);
  TC[1]->SetEditable(true);
  TC[1]->Clear("D");
  Pad[2]->cd();
  h_hits0_s->Draw("COLZPOL");
  h_hits_s->Draw("COLZPOL same");
  Pad[3]->cd();
  h_hits0_n->Draw("COLZPOL");
  h_hits_n->Draw("COLZPOL same");
  TText PrintRun;
  PrintRun.SetTextFont(62);
  PrintRun.SetTextSize(0.04);
  PrintRun.SetNDC();          // set to normalized coordinates
  PrintRun.SetTextAlign(23);  // center/top alignment
  std::ostringstream runnostream;
  std::string runstring;
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

int SepdMonDraw::DrawThird(const std::string & /* what */)
{
  OnlMonClient *cl = OnlMonClient::instance();
  TH2 *h_ADC_corr = (TH2 *) cl->getHisto("SEPDMON_0", "h_ADC_corr");
  TH2 *h_hits_corr = (TH2 *) cl->getHisto("SEPDMON_0", "h_hits_corr");
  time_t evttime = cl->EventTime("CURRENT");
  if (!gROOT->FindObject("SepdMon3"))
  {
    MakeCanvas("SepdMon3");
  }
  if (!h_ADC_corr)
  {
    DrawDeadServer(transparent[0]);
    TC[0]->SetEditable(false);
    return -1;
  }
  TC[2]->SetEditable(true);
  TC[2]->Clear("D");
  Pad[4]->cd();
  h_ADC_corr->Draw("COLZ");
  Pad[5]->cd();
  h_hits_corr->Draw("COLZ");
  TText PrintRun;
  PrintRun.SetTextFont(62);
  PrintRun.SetTextSize(0.04);
  PrintRun.SetNDC();          // set to normalized coordinates
  PrintRun.SetTextAlign(23);  // center/top alignment
  std::ostringstream runnostream;
  std::string runstring;
  // fill run number and event time into string
  runnostream << ThisName << "_2 Run " << cl->RunNumber()
              << ", Time: " << ctime(&evttime);
  runstring = runnostream.str();
  transparent[1]->cd();
  PrintRun.DrawText(0.5, 1., runstring.c_str());
  TC[2]->Update();
  TC[2]->Show();
  TC[2]->SetEditable(false);
  return 0;
}

int SepdMonDraw::DrawFourth(const std::string & /* what */)
{
  OnlMonClient *cl = OnlMonClient::instance();

  TH1F *h_waveform_time = (TH1F *) cl->getHisto("SEPDMON_0", "h1_waveform_time");
  TH1F *h_waveform_pedestal = (TH1F *) cl->getHisto("SEPDMON_0", "h1_waveform_pedestal");
  TH2F *h2_sepd_waveform = (TH2F *) cl->getHisto("SEPDMON_0", "h2_sepd_waveform");

  if (!gROOT->FindObject("SepdMon4"))
  {
    MakeCanvas("SepdMon4");
  }

  TC[3]->SetEditable(1);
  TC[3]->Clear("D");
  Pad[6]->cd();
  if (!h2_sepd_waveform || !h_waveform_time || !h_waveform_pedestal)
  {
    //cout which one is not found
    if (!h2_sepd_waveform) std::cout << "h2_sepd_waveform not found" << std::endl;
    if (!h_waveform_time) std::cout << "h_waveform_time not found" << std::endl;
    if (!h_waveform_pedestal) std::cout << "h_waveform_pedestal not found" << std::endl;
    DrawDeadServer(transparent[3]);
    TC[3]->SetEditable(0);
    return -1;
  }

  Pad[6]->cd();
  gStyle->SetTitleFontSize(0.03);
  float ymaxp = h2_sepd_waveform->ProfileX()->GetMaximum();
  h2_sepd_waveform->GetYaxis()->SetRangeUser(0, ymaxp * 20);

  h2_sepd_waveform->Draw("colz");

  float tsize = 0.06;
  h2_sepd_waveform->GetXaxis()->SetNdivisions(510, kTRUE);
  h2_sepd_waveform->GetXaxis()->SetTitle("Sample #");
  h2_sepd_waveform->GetYaxis()->SetTitle("Waveform [ADC]");
  h2_sepd_waveform->GetXaxis()->SetLabelSize(tsize);
  h2_sepd_waveform->GetYaxis()->SetLabelSize(tsize);
  h2_sepd_waveform->GetXaxis()->SetTitleSize(tsize);
  h2_sepd_waveform->GetYaxis()->SetTitleSize(tsize);
  h2_sepd_waveform->GetXaxis()->SetTitleOffset(1.2);
  h2_sepd_waveform->GetYaxis()->SetTitleOffset(0.75);
  gPad->SetLogz();
  gPad->SetBottomMargin(0.16);
  gPad->SetLeftMargin(0.2);
  gPad->SetRightMargin(0.05);
  gPad->SetLeftMargin(0.15);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(57);
  gPad->SetTicky();
  gPad->SetTickx();

  TText PrintRun;
  PrintRun.SetTextFont(62);
  PrintRun.SetTextSize(0.03);
  PrintRun.SetNDC();          // set to normalized coordinates
  PrintRun.SetTextAlign(23);  // center/top alignment
  std::ostringstream runnostream;
  std::string runstring;
  time_t evttime = getTime();
  // fill run number and event time into string
  runnostream << ThisName << ": Pulse fitting, Run" << cl->RunNumber()
              << ", Time: " << ctime(&evttime);
  runstring = runnostream.str();
  transparent[3]->cd();
  PrintRun.DrawText(0.5, 0.99, runstring.c_str());

  Pad[7]->cd();

  gStyle->SetTitleFontSize(0.06);

  float tsize2 = 0.08;
  h_waveform_time->Draw("hist");
  h_waveform_time->GetXaxis()->SetNdivisions(510, kTRUE);
  h_waveform_time->GetXaxis()->SetTitle("Sample #");
  h_waveform_time->GetYaxis()->SetTitle("Towers");
  h_waveform_time->GetXaxis()->SetLabelSize(tsize2);
  h_waveform_time->GetYaxis()->SetLabelSize(tsize2);
  h_waveform_time->GetXaxis()->SetTitleSize(tsize2);
  h_waveform_time->GetYaxis()->SetTitleSize(tsize2);
  h_waveform_time->GetXaxis()->SetTitleOffset(1.0);
  h_waveform_time->GetYaxis()->SetTitleOffset(0.85);
  gPad->SetTopMargin(0.06);
  gPad->SetBottomMargin(0.18);
  gPad->SetRightMargin(0.05);
  gPad->SetLeftMargin(0.2);
  gStyle->SetOptStat(0);
  gPad->SetTicky();
  gPad->SetTickx();

  Pad[8]->cd();

  gStyle->SetTitleFontSize(0.06);

  h_waveform_pedestal->Draw("hist");
  h_waveform_pedestal->GetXaxis()->SetNdivisions(510, kTRUE);
  h_waveform_pedestal->GetXaxis()->SetTitle("ADC Pedistal");
  h_waveform_pedestal->GetYaxis()->SetTitle("Towers");
  h_waveform_pedestal->GetXaxis()->SetLabelSize(tsize2);
  h_waveform_pedestal->GetYaxis()->SetLabelSize(tsize2);
  h_waveform_pedestal->GetXaxis()->SetTitleSize(tsize2);
  h_waveform_pedestal->GetYaxis()->SetTitleSize(tsize2);
  h_waveform_pedestal->GetXaxis()->SetTitleOffset(0.9);
  h_waveform_pedestal->GetYaxis()->SetTitleOffset(0.85);
  gPad->SetTopMargin(0.06);
  gPad->SetBottomMargin(0.18);
  gPad->SetRightMargin(0.05);
  gPad->SetLeftMargin(0.2);
  gStyle->SetOptStat(0);
  gPad->SetTicky();
  gPad->SetTickx();

  TC[3]->Update();
  TC[3]->Show();
  TC[3]->SetEditable(0);

  return 0;
}
int SepdMonDraw::DrawFifth(const std::string & /* what */)
{
  OnlMonClient *cl = OnlMonClient::instance();

  TH1F *h1_packet_number = (TH1F *) cl->getHisto("SEPDMON_0", "h1_packet_number");
  TH1F *h1_packet_length = (TH1F *) cl->getHisto("SEPDMON_0", "h1_packet_length");
  TH1F *h1_packet_chans = (TH1F *) cl->getHisto("SEPDMON_0", "h1_packet_chans");
  TH1F *h1_packet_event = (TH1F *) cl->getHisto("SEPDMON_0", "h1_packet_event");

  if (!gROOT->FindObject("SepdMon5"))
  {
    MakeCanvas("SepdMon5");
  }
  TC[4]->SetEditable(1);
  TC[4]->Clear("D");

  if (!h1_packet_number || !h1_packet_length || !h1_packet_chans || !h1_packet_event)
  {
    // print out which is not found
    if (!h1_packet_number) std::cout << "h1_packet_number not found" << std::endl;
    if (!h1_packet_length) std::cout << "h1_packet_length not found" << std::endl;
    if (!h1_packet_chans) std::cout << "h1_packet_chans not found" << std::endl;
    if (!h1_packet_event) std::cout << "h1_packet_event not found" << std::endl;

    DrawDeadServer(transparent[5]);
    TC[5]->SetEditable(0);
    return -1;
  }

  // int maxbin = h1_packet_event->GetMaximumBin();
  int maxy = h1_packet_event->GetMaximum();
  // substract all none zero bin by maxy
  for (int i = 1; i <= h1_packet_event->GetNbinsX(); i++)
  {
    if (h1_packet_event->GetBinContent(i) != 0)
    {
      h1_packet_event->SetBinContent(i, h1_packet_event->GetBinContent(i) - maxy);
    }
  }

  // find the x range for h1_packet_number
  double xmin = h1_packet_number->GetXaxis()->GetXmin();
  double xmax = h1_packet_number->GetXaxis()->GetXmax();

  TLine *one = new TLine(xmin, 1, xmax, 1);
  one->SetLineStyle(7);

  TLine *goodSize = new TLine(xmin, 3991, xmax, 3991);
  goodSize->SetLineStyle(7);

  TLine *goodChans = new TLine(xmin, 128, xmax, 128);
  goodChans->SetLineStyle(7);

  float param = 0.95;

  TLegend *leg = new TLegend(0.3, 0.16, 0.95, 0.4);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);

  TLine *warnLineOne = new TLine(xmin, param * 1, xmax, param * 1);
  warnLineOne->SetLineStyle(7);
  warnLineOne->SetLineColor(2);

  leg->AddEntry(warnLineOne, "95% Threshold", "l");

  TLine *warnLineSize = new TLine(xmin, param * 3991., xmax, param * 3991.);
  warnLineSize->SetLineStyle(7);
  warnLineSize->SetLineColor(2);

  TLine *warnLineChans = new TLine(xmin, param * 128., xmax, param * 128.);
  warnLineChans->SetLineStyle(7);
  warnLineChans->SetLineColor(2);

  Pad[10]->cd();
  float tsize = 0.08;
  h1_packet_number->GetYaxis()->SetRangeUser(0.0, 1.3);
  h1_packet_number->Draw("hist");
  one->Draw("same");
  warnLineOne->Draw("same");
  h1_packet_number->GetXaxis()->SetNdivisions(510, kTRUE);
  h1_packet_number->GetXaxis()->SetTitle("Packet #");
  h1_packet_number->GetYaxis()->SetTitle("% Of Events Present");
  h1_packet_number->GetXaxis()->SetLabelSize(tsize - 0.01);
  h1_packet_number->GetYaxis()->SetLabelSize(tsize - 0.01);
  h1_packet_number->GetXaxis()->SetTitleSize(tsize - 0.01);
  h1_packet_number->GetYaxis()->SetTitleSize(tsize - 0.01);
  h1_packet_number->GetXaxis()->SetTitleOffset(1);
  gPad->SetBottomMargin(0.16);
  gPad->SetLeftMargin(0.16);
  gPad->SetRightMargin(0.05);
  gPad->SetLeftMargin(0.15);
  gStyle->SetOptStat(0);
  gPad->SetTicky();
  gPad->SetTickx();

  Pad[11]->cd();
  h1_packet_length->Draw("hist");
  h1_packet_length->GetYaxis()->SetRangeUser(0, 6500);
  goodSize->Draw("same");
  warnLineSize->Draw("same");
  leg->Draw("same");
  h1_packet_length->GetXaxis()->SetNdivisions(510, kTRUE);
  h1_packet_length->GetXaxis()->SetTitle("Packet #");
  h1_packet_length->GetYaxis()->SetTitle("Average Packet Size");
  h1_packet_length->GetXaxis()->SetLabelSize(tsize - .01);
  h1_packet_length->GetYaxis()->SetLabelSize(tsize);
  h1_packet_length->GetXaxis()->SetTitleSize(tsize - .01);
  h1_packet_length->GetYaxis()->SetTitleSize(tsize);
  h1_packet_length->GetXaxis()->SetTitleOffset(1);
  h1_packet_length->GetYaxis()->SetTitleOffset(0.8);
  gPad->SetBottomMargin(0.16);
  gPad->SetLeftMargin(0.16);
  gPad->SetRightMargin(0.05);
  gPad->SetLeftMargin(0.15);
  gStyle->SetOptStat(0);
  gPad->SetTicky();
  gPad->SetTickx();

  Pad[12]->cd();
  h1_packet_chans->Draw("hist");
  h1_packet_chans->GetYaxis()->SetRangeUser(0, 212);
  goodChans->Draw("same");
  warnLineChans->Draw("same");
  h1_packet_chans->GetXaxis()->SetNdivisions(510, kTRUE);
  h1_packet_chans->GetXaxis()->SetTitle("Packet #");
  h1_packet_chans->GetYaxis()->SetTitle("Average # of Channels");
  h1_packet_chans->GetXaxis()->SetLabelSize(tsize - .01);
  h1_packet_chans->GetYaxis()->SetLabelSize(tsize);
  h1_packet_chans->GetXaxis()->SetTitleSize(tsize - .01);
  h1_packet_chans->GetYaxis()->SetTitleSize(tsize);
  h1_packet_chans->GetXaxis()->SetTitleOffset(0.8);
  h1_packet_chans->GetYaxis()->SetTitleOffset(0.8);
  gPad->SetBottomMargin(0.16);
  gPad->SetLeftMargin(0.16);
  gPad->SetRightMargin(0.05);
  gPad->SetLeftMargin(0.15);
  gStyle->SetOptStat(0);
  gPad->SetTicky();
  gPad->SetTickx();

  Pad[13]->cd();
  h1_packet_event->Draw("hist");
  double ymax = h1_packet_event->GetMaximum();
  double ymin = h1_packet_event->GetMinimum();

  h1_packet_event->GetYaxis()->SetRangeUser(ymin - 0.3 * (ymax - ymin + 30), ymax + 0.3 * (ymax - ymin + 30));
  h1_packet_event->GetXaxis()->SetTitle("Packet #");
  h1_packet_event->GetYaxis()->SetTitle("clock offset");
  h1_packet_event->GetXaxis()->SetLabelSize(tsize - .01);
  h1_packet_event->GetYaxis()->SetLabelSize(tsize);
  h1_packet_event->GetXaxis()->SetTitleSize(tsize - .01);
  h1_packet_event->GetYaxis()->SetTitleSize(tsize);
  h1_packet_event->GetXaxis()->SetTitleOffset(0.8);
  h1_packet_event->GetYaxis()->SetTitleOffset(1.2);
  gPad->SetBottomMargin(0.16);
  gPad->SetLeftMargin(0.2);
  gPad->SetRightMargin(0.05);
  gPad->SetLeftMargin(0.2);
  gStyle->SetOptStat(0);
  gPad->SetTicky();
  gPad->SetTickx();


  // draw the warning here:
  
  warning[1]->cd();
  
  std::vector<int> badPackets;
  std::vector<std::string> whatswrong;
  for (int i = 1; i <= 6; i++)
  {
    bool missing = false;
    bool badnumber = false;
    bool badlength = false;
    bool badchans = false;
    if (h1_packet_number->GetBinContent(i) == 0)
    {
      missing = true;
    }
    if (h1_packet_number->GetBinContent(i) < param)
    {
      badnumber = true;
    }
    if (h1_packet_length->GetBinContent(i) < param * 3991.)
    {
      badlength = true;
    }
    if (h1_packet_chans->GetBinContent(i) < param * 128.)
    {
      badchans = true;
    }
    if (badnumber || badlength || badchans || missing)
    {
      badPackets.push_back((int) h1_packet_number->GetBinCenter(i));
      std::string reason = "";
      if (missing)
      {
        reason += "packet lost! ";
      }
      else
      {
        if (badnumber)
        {
          reason += "some events are missing, ";
        }
        if (badlength)
        {
          reason += "too short, ";
        }
        if (badchans)
        {
          reason += "too few channels, ";
        }
        // remove the last two characters
        reason = reason.substr(0, reason.size() - 2);
        reason += ".";
      }
      whatswrong.push_back(reason);
    }
  }
  bool mismatch = false;
  for (int i = 1; i <= h1_packet_event->GetNbinsX(); i++)
  {
    if (h1_packet_event->GetBinContent(i) != 0)
    {
      mismatch = true;
    }
  }
  // draw the mismatch warning on the pad
  TText mismatchWarn;
  mismatchWarn.SetTextFont(62);
  mismatchWarn.SetTextSize(0.06);
  mismatchWarn.SetTextColor(2);
  mismatchWarn.SetNDC();
  mismatchWarn.SetTextAlign(23);
  if (mismatch)
  {
    mismatchWarn.DrawText(0.5, 0.95, "Packet misaligned!");
  }

  mismatchWarn.SetTextColor(1);
  mismatchWarn.SetTextSize(0.05);
 
  // draw the bad packet warning here
  TText PacketWarn;
  PacketWarn.SetTextFont(62);
  PacketWarn.SetTextSize(0.04);
  PacketWarn.SetTextColor(1);
  PacketWarn.SetNDC();
  PacketWarn.SetTextAlign(23);
  PacketWarn.DrawText(0.5, 0.75, "Bad Packets:");
  for (int i = 0; i < (int) badPackets.size(); i++)
  {
    PacketWarn.DrawText(0.5, 0.7 - 0.05 * i, Form("%i: %s", badPackets[i], whatswrong[i].c_str()));
  }

  TText PrintRun;
  PrintRun.SetTextFont(62);
  PrintRun.SetTextSize(0.02);
  PrintRun.SetNDC();          // set to normalized coordinates
  PrintRun.SetTextAlign(23);  // center/top alignment
  std::ostringstream runnostream;
  std::string runstring;
  std::ostringstream runnostream2;
  time_t evttime = getTime();
  // fill run number and event time into string

  runnostream << "Packet Information";
  runnostream2 << " Run " << cl->RunNumber() << ", Time: " << ctime(&evttime);
  transparent[4]->cd();

  runstring = runnostream.str();
  PrintRun.DrawText(0.5, .99, runstring.c_str());

  runstring = runnostream2.str();
  PrintRun.DrawText(0.5, .966, runstring.c_str());
  TC[4]->Update();
  TC[4]->Show();
  TC[4]->SetEditable(0);
  // if (save) TC[5]->SaveAs("plots/packets.pdf");
  return 0;
}

int SepdMonDraw::SavePlot(const std::string &what, const std::string &type)
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

int SepdMonDraw::MakeHtml(const std::string &what)
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

time_t SepdMonDraw::getTime()
{
  OnlMonClient *cl = OnlMonClient::instance();
  time_t currtime = 0;

  currtime = cl->EventTime("CURRENT");
  return currtime;
}
