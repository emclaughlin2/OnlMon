#include "LL1MonDraw.h"

#include <onlmon/OnlMonClient.h>

#include <TAxis.h>  // for TAxis
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TPad.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TText.h>
#include <TStyle.h>

#include <cstring>  // for memset
#include <ctime>
#include <fstream>
#include <iostream>  // for operator<<, basic_ostream, basic_os...
#include <sstream>
#include <vector>  // for vector

LL1MonDraw::LL1MonDraw(const std::string &name)
  : OnlMonDraw(name)
{
  return;
}

int LL1MonDraw::Init()
{
  return 0;
}

int LL1MonDraw::MakeCanvas(const std::string &name)
{
  OnlMonClient *cl = OnlMonClient::instance();
  int xsize = cl->GetDisplaySizeX();
  int ysize = cl->GetDisplaySizeY();
  gStyle->SetOptStat(0);
  if (name == "LL1Mon1")
  {
    // xpos (-1) negative: do not draw menu bar
    TC[0] = new TCanvas(name.c_str(), "ll1Mon Monitor", -1, 0, xsize/2, ysize/2);
    TC[0]->SetTicks(1,1);
    // root is pathetic, whenever a new TCanvas is created root piles up
    // 6kb worth of X11 events which need to be cleared with
    // gSystem->ProcessEvents(), otherwise your process will grow and
    // grow and grow but will not show a definitely lost memory leak
    gSystem->ProcessEvents();
    for(int ipad=0; ipad<nPad1; ipad++){
      Pad[ipad] = new TPad(Form("ll1pad%d",ipad),Form("ll1pad%d",ipad),padx1[ipad],pady1[ipad],padx2[ipad],pady2[ipad],0);
      Pad[ipad]->Draw();
      Pad[ipad]->SetTicks(1,1);
    }
    transparent[0] = new TPad("transparent0", "this does not show", 0, 0, 1, 1);
    transparent[0]->SetFillStyle(4000);
    transparent[0]->Draw();
    TC[0]->SetEditable(false);
  }
  else if (name == "LL1Mon2")
  {
    // xpos negative: do not draw menu bar
    gStyle->SetOptStat(0);
    TC[1] = new TCanvas(name.c_str(), "ll1Mon2 Monitor", -xsize / 2, 0, xsize / 2, ysize/2);
    TC[1]->SetTicks(1,1);
    gSystem->ProcessEvents();
    // this one is used to plot the run number on the canvas
    transparent[1] = new TPad("transparent1", "this does not show", 0, 0, 1, 1);
    transparent[1]->SetFillStyle(4000);
    transparent[1]->Draw();
    TC[1]->SetEditable(false);
  }
  else if (name == "LL1Mon3")
  {
    TC[2] = new TCanvas(name.c_str(), "ll1Mon3 Monitor", xsize / 2, 0, xsize / 2, ysize/2);
    TC[2]->SetTicks(1,1);
    gSystem->ProcessEvents();
    // this one is used to plot the run number on the canvas
            transparent[2] = new TPad("transparent2", "this does not show", 0, 0, 1, 1);
            transparent[2]->SetFillStyle(4000);
            transparent[2]->Draw();
          TC[2]->SetEditable(0);
  }
  return 0;
}

int LL1MonDraw::Draw(const std::string &what)
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
  if (!idraw)
  {
    std::cout << __PRETTY_FUNCTION__ << " Unimplemented Drawing option: " << what << std::endl;
    iret = -1;
  }
  return iret;
}

int LL1MonDraw::DrawFirst(const std::string & /* what */)
{
  OnlMonClient *cl = OnlMonClient::instance();
  TH1 *h_nhit_n1 = cl->getHisto("LL1MON_0","h_nhit_n1");
  TH1 *h_nhit_n2 = cl->getHisto("LL1MON_0","h_nhit_n2");
  TH1 *h_nhit_s1 = cl->getHisto("LL1MON_0","h_nhit_s1");
  TH1 *h_nhit_s2 = cl->getHisto("LL1MON_0","h_nhit_s2");
  time_t evttime = cl->EventTime("CURRENT");
  if (!gROOT->FindObject("LL1Mon1"))
  {
    MakeCanvas("LL1Mon1");
  }
  TC[0]->SetEditable(true);
  TC[0]->Clear("D");
  Pad[0]->cd();
  if (h_nhit_n1)
  {
    h_nhit_n1->DrawCopy();
  }
  else
  {
    DrawDeadServer(transparent[0]);
    TC[0]->SetEditable(false);
    return -1;
  }
  Pad[1]->cd();
  h_nhit_n2->DrawCopy();
  Pad[2]->cd();
  h_nhit_s1->DrawCopy();
  Pad[3]->cd();
  h_nhit_s2->DrawCopy();

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
  gStyle->SetOptStat(0);
  TC[0]->Update();
  TC[0]->SetTicks(1,1);
  TC[0]->Show();
  TC[0]->SetEditable(false);
  return 0;
}

int LL1MonDraw::DrawSecond(const std::string & /* what */)
{
  OnlMonClient *cl = OnlMonClient::instance();
  TH2 *h_nhit_corr= (TH2D*) cl->getHisto("LL1MON_0","h_nhit_corr");
  time_t evttime = cl->EventTime("CURRENT");
  if (!gROOT->FindObject("LL1Mon2"))
  {
    MakeCanvas("LL1Mon2");
  }
  TC[1]->SetEditable(true);
  TC[1]->Clear("D");
  if (h_nhit_corr)
  {
    h_nhit_corr->Draw("colz");
  }
  else
  {
    DrawDeadServer(transparent[1]);
    TC[1]->SetEditable(false);
    return -1;
  }
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
  TC[1]->SetTicks(1,1);
  TC[1]->Update();
  TC[1]->Show();
  TC[1]->SetEditable(false);
  return 0;
}

int LL1MonDraw::DrawThird(const std::string & /* what */)
{
  OnlMonClient *cl = OnlMonClient::instance();
  TH2 *h_line_up= (TH2*) cl->getHisto("LL1MON_0","h_line_up");
  time_t evttime = cl->EventTime("CURRENT");
  if (!gROOT->FindObject("LL1Mon3"))
  {
    MakeCanvas("LL1Mon3");
  }
  TC[2]->SetEditable(true);
  TC[2]->Clear("D");
  gStyle->SetOptStat(0);
  if (h_line_up)
  {
    h_line_up->Draw("colz");
  }
  else
  {
    DrawDeadServer(transparent[2]);
    TC[2]->SetEditable(false);
    return -1;
  }
  TText PrintRun;
  PrintRun.SetTextFont(62);
  PrintRun.SetTextSize(0.04);
  PrintRun.SetNDC();          // set to normalized coordinates
  PrintRun.SetTextAlign(23);  // center/top alignment
  std::ostringstream runnostream;
  std::string runstring;
  // fill run number and event time into string
  runnostream << ThisName << "_3 Run " << cl->RunNumber()
              << ", Time: " << ctime(&evttime);
  runstring = runnostream.str();
  transparent[2]->cd();
  PrintRun.DrawText(0.5, 1., runstring.c_str());
  TC[2]->Update();
  TC[2]->SetTicks(1,1);
  TC[2]->Show();
  TC[2]->SetEditable(false);
  return 0;
}

int LL1MonDraw::SavePlot(const std::string &what, const std::string &type)
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

int LL1MonDraw::MakeHtml(const std::string &what)
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
