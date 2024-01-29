#include "ZdcMonDraw.h"

#include <onlmon/OnlMonClient.h>
#include <onlmon/OnlMonDB.h>

#include <TAxis.h>  // for TAxis
#include <TCanvas.h>
#include <TDatime.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TPad.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TText.h>

#include <cstring>  // for memset
#include <ctime>
#include <fstream>
#include <iostream>  // for operator<<, basic_ostream, basic_os...
#include <sstream>
#include <vector>  // for vector

ZdcMonDraw::ZdcMonDraw(const std::string &name)
  : OnlMonDraw(name)
{
  return;
}

int ZdcMonDraw::Init()
{
  return 0;
}

int ZdcMonDraw::MakeCanvas(const std::string &name)
{
  OnlMonClient *cl = OnlMonClient::instance();
  int xsize = cl->GetDisplaySizeX();
  int ysize = cl->GetDisplaySizeY();
  if (name == "ZdcMon1")
  {
    // xpos (-1) negative: do not draw menu bar
    TC[0] = new TCanvas(name.c_str(), "ZdcMon Example Monitor", -1, 0, xsize / 2, ysize);
    // root is pathetic, whenever a new TCanvas is created root piles up
    // 6kb worth of X11 events which need to be cleared with
    // gSystem->ProcessEvents(), otherwise your process will grow and
    // grow and grow but will not show a definitely lost memory leak
    gSystem->ProcessEvents();
    Pad[0] = new TPad("zdcpad1", "who needs this?", 0.1, 0.5, 0.9, 0.9, 0);
    Pad[1] = new TPad("zdcpad2", "who needs this?", 0.1, 0.05, 0.9, 0.45, 0);
    Pad[0]->Draw();
    Pad[1]->Draw();
    // this one is used to plot the run number on the canvas
    transparent[0] = new TPad("transparent0", "this does not show", 0, 0, 1, 1);
    transparent[0]->SetFillStyle(4000);
    transparent[0]->Draw();
    TC[0]->SetEditable(false);
  }
  else if (name == "ZdcMon2")
  {
    // xpos negative: do not draw menu bar
    TC[1] = new TCanvas(name.c_str(), "ZdcMon2 Example Monitor", -xsize / 2, 0, xsize / 2, ysize);
    gSystem->ProcessEvents();
    Pad[2] = new TPad("zdcpad3", "who needs this?", 0.1, 0.5, 0.9, 0.9, 0);
    Pad[3] = new TPad("zdcpad4", "who needs this?", 0.1, 0.05, 0.9, 0.45, 0);
    Pad[2]->Draw();
    Pad[3]->Draw();
    // this one is used to plot the run number on the canvas
    transparent[1] = new TPad("transparent1", "this does not show", 0, 0, 1, 1);
    transparent[1]->SetFillStyle(4000);
    transparent[1]->Draw();
    TC[1]->SetEditable(false);
  }
  return 0;
}

int ZdcMonDraw::Draw(const std::string &what)
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
  if (!idraw)
  {
    std::cout << __PRETTY_FUNCTION__ << " Unimplemented Drawing option: " << what << std::endl;
    iret = -1;
  }
  return iret;
}

int ZdcMonDraw::DrawFirst(const std::string & /* what */)
{
  OnlMonClient *cl = OnlMonClient::instance();
  TH1 *zdc_adc_south = cl->getHisto("ZDCMON_0","zdc_adc_south");
  TH1 *zdc_adc_north = cl->getHisto("ZDCMON_0","zdc_adc_north");
  if (!gROOT->FindObject("ZdcMon1"))
  {
    MakeCanvas("ZdcMon1");
  }
  TC[0]->SetEditable(true);
  TC[0]->Clear("D");
  Pad[0]->cd();
  gPad->SetLogy();
  gPad->SetLogx();
  if (zdc_adc_south)
  {
      zdc_adc_south->Scale(1/zdc_adc_south->Integral(), "width");
      zdc_adc_south->DrawCopy();
  }
  else
  {
    DrawDeadServer(transparent[0]);
    TC[0]->SetEditable(false);
    return -1;
  }
  Pad[1]->cd();
  gPad->SetLogy();
  gPad->SetLogx();
  if (zdc_adc_north)
  {
      zdc_adc_north->Scale(1/zdc_adc_north->Integral(), "width");
      zdc_adc_north->DrawCopy();
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

int ZdcMonDraw::DrawSecond(const std::string & /* what */)
{
  OnlMonClient *cl = OnlMonClient::instance();
  TH1 *zdc_adc_south = cl->getHisto("ZDCMON_0","zdc_adc_south");
  TH1 *zdc_adc_north = cl->getHisto("ZDCMON_0","zdc_adc_north");
  if (!gROOT->FindObject("ZdcMon2"))
  {
    MakeCanvas("ZdcMon2");
  }
  TC[1]->SetEditable(true);
  TC[1]->Clear("D");
  Pad[2]->cd();
  if (zdc_adc_south)
  {
    zdc_adc_south->DrawCopy();
  }
  else
  {
    DrawDeadServer(transparent[1]);
    TC[1]->SetEditable(false);
    return -1;
  }
  Pad[3]->cd();
  if (zdc_adc_north)
  {
     zdc_adc_north->DrawCopy();
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

int ZdcMonDraw::SavePlot(const std::string &what, const std::string &type)
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

int ZdcMonDraw::MakeHtml(const std::string &what)
{
  int iret = Draw(what);
  if (iret)  // on error no html output please
  {
    return iret;
  }

  OnlMonClient *cl = OnlMonClient::instance();

  // Register the 1st canvas png file to the menu and produces the png file.
  std::string pngfile = cl->htmlRegisterPage(*this, "First Canvas", "1", "png");
  cl->CanvasToPng(TC[0], pngfile);

  // idem for 2nd canvas.
  pngfile = cl->htmlRegisterPage(*this, "Second Canvas", "2", "png");
  cl->CanvasToPng(TC[1], pngfile);
  // Now register also EXPERTS html pages, under the EXPERTS subfolder.

  std::string logfile = cl->htmlRegisterPage(*this, "EXPERTS/Log", "log", "html");
  std::ofstream out(logfile.c_str());
  out << "<HTML><HEAD><TITLE>Log file for run " << cl->RunNumber()
      << "</TITLE></HEAD>" << std::endl;
  out << "<P>Some log file output would go here." << std::endl;
  out.close();

  std::string status = cl->htmlRegisterPage(*this, "EXPERTS/Status", "status", "html");
  std::ofstream out2(status.c_str());
  out2 << "<HTML><HEAD><TITLE>Status file for run " << cl->RunNumber()
       << "</TITLE></HEAD>" << std::endl;
  out2 << "<P>Some status output would go here." << std::endl;
  out2.close();
  cl->SaveLogFile(*this);
  return 0;
}
