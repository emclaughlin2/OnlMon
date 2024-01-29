#include "PktSizeMonDraw.h"
#include "PktSizeCommon.h"
#include "PktSizeDBodbc.h"

#include <onlmon/OnlMonClient.h>
#include <onlmon/RunDBodbc.h>

#include <phool/phool.h>

#include <TCanvas.h>
#include <TDatime.h>
#include <TGraph.h>
#include <TH2.h>
#include <TLine.h>
#include <TMarker.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TText.h>

#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <sstream>

// make history for currently noisy packets, not packets whose size has
// changed significantly with respect to previous runs
#define CURRENTNOISY

static int bytelimit = 220;
static float MAXSIZEDISP = 1000.;
static float MAXPKTDISP = 26000.;
static float GRANTXTOFFSET = 20.;

PktSizeMonDraw::PktSizeMonDraw(const std::string &name)
  : OnlMonDraw(name)
{
  memset(transparent, 0, sizeof(transparent));
  memset(TC, 0, sizeof(TC));
  memset(Pad, 0, sizeof(Pad));
  tm = new TMarker();
  tm->SetMarkerStyle(8);
  tm->SetMarkerColor(1); 
  lastrun = 0;
  return;
}

PktSizeMonDraw::~PktSizeMonDraw()
{
  packetmap.clear();
  knownbig.clear();
  activepackets.clear();
  delete rd;
  delete db;
  for ( int q = 0; q < 12; q++)
  {
  delete Frames[q];
  }
  return;
}

int PktSizeMonDraw::Init()
{
  PktSizeCommon::fillgranules(granulepacketlimits);
  db = new PktSizeDBodbc(ThisName);
  rd = new RunDBodbc();
  return 0;
}

int PktSizeMonDraw::MakeCanvas(const char *name)
{
  OnlMonClient *cl = OnlMonClient::instance();
  int xsize = cl->GetDisplaySizeX();
  int ysize = cl->GetDisplaySizeY();
  if (!strcmp(name, "PktSizeMon0"))
  {
    // xpos (-1) negative: do not draw menu bar
    TC[0] = new TCanvas(name, "Packet Size Monitor", -1, 0, xsize / 2, ysize);
    // gSystem->ProcessEvents();
    // Pad[0] = new TPad("pktpad1", "who needs this?", 0.1, 0.05, 0.9, 0.9, 0);
    // Pad[0]->Draw();
    // // this one is used to plot the run number on the canvas
    // transparent[0] = new TPad("transparent0", "this does not show", 0, 0, 1, 1);
    // transparent[0]->SetFillStyle(4000);
    // transparent[0]->Draw();
  }
  else if (!strcmp(name, "PktSizeMon1"))
  {
    // xpos negative: do not draw menu bar
    TC[1] = new TCanvas(name, "Packet Size History", -xsize / 2, 0, xsize / 2, ysize);
    gSystem->ProcessEvents();
    Pad[1] = new TPad("pktpad1", "who needs this?", 0.1, 0.05, 0.9, 0.9, 0);
    Pad[1]->Draw();
    // this one is used to plot the run number on the canvas
    transparent[1] = new TPad("transparent1", "this does not show", 0, 0, 1, 1);
    transparent[1]->SetFillStyle(4000);
    transparent[1]->Draw();
  }
  return 0;
}

int PktSizeMonDraw::Draw(const std::string &what)
{
  int iret = 0;
  int idraw = 0;
  if (what == "ALL" || what == "FIRST")
  {
    iret += DrawFirst(what);
    idraw++;
  }
  if (what == "HISTORY")
  {
    iret += DrawHistory(what);
    idraw++;
  }
  if (!idraw)
  {
    std::cout << PHWHERE << " Unimplemented Drawing option: " << what << std::endl;
    iret = -1;
  }
  return iret;
}

int PktSizeMonDraw::DrawFirst(const std::string & /*what*/)
{
  for ( int w = 0; w < 12; w++)
  {
  delete Frames[w];
  }
  OnlMonClient *cl = OnlMonClient::instance();
  if (!gROOT->FindObject("PktSizeMon0"))
  {
    MakeCanvas("PktSizeMon0");
  }
  //  TC[0]->Clear("D");
  TH1 *pktsize_hist = cl->getHisto("PKTSIZEMON_0", "pktsize_hist");
  if (!pktsize_hist)
  {
    DrawDeadServer(transparent[0]);
    TC[0]->Update();
    return -1;
  }
  TPad* Pads[12];
  Pads[0] = new TPad("", "", 0, 0.66, 0.25, 1, 0);
  Pads[1] = new TPad("", "", 0.25, 0.66, 0.5, 1, 0);   
  Pads[2] = new TPad("", "", 0.5, 0.66, 0.75, 1, 0);  
  Pads[3] = new TPad("", "", 0.75, 0.66, 1, 1, 0);   
  Pads[4] = new TPad("", "", 0, 0.33, 0.25, 0.66, 0);  
  Pads[5] = new TPad("", "", 0.25, 0.33, 0.5, 0.66, 0);   
  Pads[6] = new TPad("", "", 0.5, 0.33, 0.75, 0.66, 0);   
  Pads[7] = new TPad("", "", 0.75, 0.33, 1, 0.66, 0);   
  Pads[8] = new TPad("", "", 0, 0, .25, 0.33, 0);   
  Pads[9] = new TPad("", "", 0.25, 0, 0.5, 0.33, 0);   
  Pads[10] = new TPad("", "", 0.5, 0, 0.75, 0.33, 0);   
  Pads[11] = new TPad("", "", 0.75, 0, 1, 0.33, 0);
  
  Frames[0] = new TH2F("h2", "MBD", 2, 1000, 1003, 2, 0, 10000);
  Frames[1] = new TH2F("h2", "MVTX", 2, 2000, 2200, 2, 0, 10000);
  Frames[2] = new TH2F("h2", "INTT", 2, 3000, 3200, 2, 0, 10000);
  Frames[3] = new TH2F("h2", "TPC", 2, 4230, 4232, 2, 0, 10000);
  Frames[4] = new TH2F("h2", "TPOT", 2, 5000, 5200, 2, 0, 10000);
  Frames[5] = new TH2F("h2", "CEMC", 2, 6000, 6126, 2, 0, 10000);
  Frames[6] = new TH2F("h2", "IHCAL", 2, 7000, 7008, 2, 0, 10000);
  Frames[7] = new TH2F("h2", "OHCAL", 2, 8000, 8008, 2, 0, 10000);
  Frames[8] = new TH2F("h2", "SEPD", 2, 9000, 9002, 2, 0, 10000);
  Frames[9] = new TH2F("h2", "ZDC", 2, 12000, 12002, 2, 0, 10000);
  Frames[10] = new TH2F("h2", "Local lvl.1", 2, 13999, 14001, 2, 0, 10000);
  Frames[11] = new TH2F("h2", "Global lvl.1", 2, 14000, 14002, 2, 0, 10000);

  int n = 0;
  while (n < 12)
  {
  Pads[n]->SetLeftMargin(.125);
  Pads[n]->Draw();
  n = n+1;
  }
  int x = 0;
  while (x < 12)
  {
  Pads[x]->cd();
  Frames[x]->GetXaxis()->SetTitle("Packet ID");
  Frames[x]->GetXaxis()->CenterTitle(true);
  Frames[x]->GetYaxis()->SetTitle("Packet size (bytes)");
  Frames[x]->GetYaxis()->CenterTitle(true); 
  Frames[x]->SetStats(0);
  Frames[x]->DrawClone();
  x = x+1;
  }
 
  std::map <int, double> packet;
  for (int i = 1; i <= pktsize_hist->GetNbinsX(); i++)
  {
  int a = pktsize_hist->GetBinError(i);
  double b = pktsize_hist->GetBinContent(i);
  packet[a] = b;  
  if (a > 999 && a < 1003)
  {
  Pads[0]-> cd();
  tm->DrawMarker(a, packet[a]); 
  }
  else if (a < 1000 || (a > 1002 && a < 2000))
    {
    std::cout<< "unknown packet id: " << (a) << '\n';
    }
  if (a > 1999 && a < 2201)
  {
  Pads[1]-> cd();
  tm->DrawMarker(a, packet[a]);
  } 
  else if (a > 2200 && a < 3000)
    {
    std::cout<< "unknown packet id: " << (a) << '\n';
    }
  if (a > 2999 && a < 3201)
  {
  Pads[2]-> cd();  
  tm->DrawMarker(a, packet[a]);
  }
  else if (a > 3200 && a < 4231)
    {
    std::cout<< "unknown packet id: " << (a) << '\n';
    }
  if (a == 4231)
  {
  Pads[3]-> cd();
  tm->DrawMarker(a, packet[a]);
  } 
  else if (a > 4231 && a < 5000)
    {
    std::cout<< "unknown packet id: " << (a) << '\n';
    }
  if (a > 4999 && a < 5201)
  {
  Pads[4]-> cd();
  tm->DrawMarker(a, packet[a]);
  } 
  else if (a > 5200 && a < 6000)
    {
    std::cout<< "unknown packet id: " << (a) << '\n';
    } 
  if (a > 5999 && a < 6127)
  {
  Pads[5]-> cd();
  tm->DrawMarker(a, packet[a]);
  } 
  else if (a > 6126 && a < 7000)
    {
    std::cout<< "unknown packet id: " << (a) << '\n';
    }  
  if (a > 6999 && a < 7009)
  {
  Pads[6]-> cd(); 
  tm->DrawMarker(a, packet[a]); 
  }
  else if ( a > 7008 && a < 8000)
    {
    std::cout<< "unknown packet id: " << (a) << '\n';
    }
  if (a > 7999 && a < 8009)
  {
  Pads[7]-> cd();
  tm->DrawMarker(a, packet[a]);
  } 
  else if (a > 8008 && a < 9001)
    {
    std::cout<< "unknown packet id: " << (a) << '\n';
    }
  if (a == 9001)
  {
  Pads[8]-> cd();  
  tm->DrawMarker(a, packet[a]);
  }
  else if (a > 9001 && a < 12001)
    {
    std::cout<< "unknown packet id: " << (a) << '\n';
    }
  if (a == 12001)
  {
  Pads[9]-> cd();
  tm->DrawMarker(a, packet[a]);
  }   
  else if (a > 12001 && a < 14000)
    {
    std::cout<< "unknown packet id: " << (a) << '\n';
    }
  if (a == 14000)
  {
  Pads[10]-> cd();
  tm->DrawMarker(a, packet[a]);
  } 
  if (a == 14001)
  {
  Pads[11]-> cd();
  tm->DrawMarker(a, packet[a]);
  } 
  else if (a > 14001)
    {
    std::cout<< "unknown packet id: " << (a) << '\n';
    }
  }
  return 0;
}

int PktSizeMonDraw::DrawOldFirst(const std::string & /*what*/)
{
  OnlMonClient *cl = OnlMonClient::instance();
  if (!gROOT->FindObject("PktSizeMon0"))
  {
    MakeCanvas("PktSizeMon0");
  }
  TC[0]->Clear("D");
  TH1 *pktsize_hist = cl->getHisto("PKTSIZEMON_0", "pktsize_hist");
  if (!pktsize_hist)
  {
    DrawDeadServer(transparent[0]);
    TC[0]->Update();
    return -1;
  }
  FillPacketMap(pktsize_hist);
  TText PrintRun;
  PrintRun.SetTextFont(62);
  PrintRun.SetTextSize(0.03);
  PrintRun.SetNDC();          // set to normalized coordinates
  PrintRun.SetTextAlign(23);  // center/top alignment
  std::ostringstream runnostream;
  time_t evttime = cl->EventTime("CURRENT");
  // fill run number and event time into string
  int runnumber = cl->RunNumber();
  runnostream << "Packet Size Display Run " << runnumber
              << ", Time: " << ctime(&evttime);
  //  transparent[0]->cd();
  PrintRun.DrawText(0.5, 0.98, runnostream.str().c_str());
  runnostream.str("");
  runnostream << "Based on " << pktsize_hist->GetBinContent(0) << " Events";
  PrintRun.DrawText(0.5, 0.94, runnostream.str().c_str());
  Pad[0]->cd();
  TH2 *htmp = new TH2F("pktsize", "", 2, 0, MAXPKTDISP, 2, 0, MAXSIZEDISP);
  htmp->SetStats(kFALSE);
  htmp->GetXaxis()->SetNoExponent();
  htmp->GetXaxis()->SetTitle("Packet Id");
  htmp->DrawCopy();
  delete htmp;
  TLine tl;
  tl.DrawLine(0, bytelimit, MAXPKTDISP, bytelimit);
  TMarker marker;
  TText outtxt;
  outtxt.SetTextFont(62);
  outtxt.SetTextSize(0.02);
  outtxt.SetTextColor(1);
  outtxt.SetTextAlign(13);
  TLine ta;
  TMarker ArrowHead;
  ArrowHead.SetMarkerSize(1.5);
  ArrowHead.SetMarkerColor(2);
  ArrowHead.SetMarkerStyle(26);

  ta.SetLineColor(2);
  ta.SetLineStyle(1);
  std::map<int, std::map<unsigned int, float> >::const_iterator iter;
  std::map<unsigned int, float>::const_iterator piter;
  iter = packetmap.find(runnumber);
  float ypos = MAXSIZEDISP - 100.;
  int foundnoisypacket = 0;
  if (iter != packetmap.end())
  {
    marker.SetMarkerStyle(22);
    marker.SetMarkerSize(2);
    marker.SetMarkerColor(2);
    for (piter = iter->second.begin(); piter != iter->second.end(); ++piter)
    {
      if (piter->second > bytelimit)
      {
        if (!IsKnownBig(piter->first, piter->second))
        {
          float linetopend = 9;
          marker.DrawMarker((double) (piter->first), piter->second);
          std::ostringstream pktid;
          pktid << std::setprecision(3) << piter->first << "(" << piter->second << " w)";
          int xpos = ((piter->first) / 1000) * 1000;
          if (xpos == 21000)
          {
            xpos = ((piter->first) / 100) * 100;
            if (xpos == 21300)
            {
              linetopend = 50;
            }
          }
          else if (xpos == 7000)
          {
            xpos = ((piter->first) / 100) * 100;
            if (xpos == 7100)
            {
              linetopend = 50;
            }
          }
          else if (xpos == 24000)
          {
            xpos = ((piter->first) / 100) * 100;
            if (xpos == 24100)
            {
              linetopend = 50;
            }
          }
          float arrowoffset = linetopend + 6;
          ta.DrawLine(xpos, 0, xpos, MAXSIZEDISP + linetopend);
          ArrowHead.DrawMarker(xpos, MAXSIZEDISP + arrowoffset);
          outtxt.DrawText(xpos, ypos, pktid.str().c_str());
          ypos -= 20;
          foundnoisypacket = 1;
        }
      }
    }
  }
  if (!foundnoisypacket)
  {
    outtxt.SetNDC();
    outtxt.SetTextSize(0.05);
    outtxt.SetTextColor(3);
    outtxt.SetTextAlign(22);
    outtxt.DrawText(0.5, 0.7, "Congratulations");
    outtxt.DrawText(0.5, 0.5, "No");
    outtxt.DrawText(0.5, 0.3, "Noisy Packets");
  }
  TText granlabel;
  granlabel.SetTextAlign(12);
  granlabel.SetTextAngle(90);
  std::map<std::string, std::pair<unsigned int, unsigned int> >::const_iterator graniter;
  for (graniter = granulepacketlimits.begin(); graniter != granulepacketlimits.end(); ++graniter)
  {
    granlabel.SetTextSize(0.03);
    int xpos = ((graniter->second.first) / 1000) * 1000;
    if (xpos != 14000)
    {
      if (xpos == 24000)
      {
        granlabel.SetTextSize(0.015);
        xpos = ((graniter->second.first) / 100) * 100;
        if (xpos == 24000)
        {
          granlabel.DrawText(xpos, MAXSIZEDISP + 20, "VTXP");
        }
        else
        {
          granlabel.DrawText(xpos, MAXSIZEDISP + 70, "VTXS");
        }
      }
      else if (xpos == 7000)
      {
        xpos = ((graniter->second.first) / 100) * 100;
        granlabel.SetTextSize(0.015);
        if (xpos == 7000)
        {
          granlabel.DrawText(xpos, MAXSIZEDISP + 20, "TOFE");
        }
        else
        {
          granlabel.DrawText(xpos, MAXSIZEDISP + 70, "TOFW");
        }
      }
      else if (xpos == 21000)
      {
        xpos = ((graniter->second.first) / 100) * 100;
        granlabel.SetTextSize(0.015);
        if (xpos == 21100)
        {
          granlabel.DrawText(xpos, MAXSIZEDISP + 20, "MPC");
        }
        else
        {
          granlabel.DrawText(xpos, MAXSIZEDISP + 70, "MPCEX");
        }
      }
      else
      {
        granlabel.DrawText(xpos, MAXSIZEDISP + GRANTXTOFFSET, (graniter->first).c_str());
      }
    }
  }
  TC[0]->Update();
  return 0;
}

int PktSizeMonDraw::FillPacketMap(const TH1 *pktsize_hist)
{
  OnlMonClient *cl = OnlMonClient::instance();
  int runnumber = cl->RunNumber();
  std::map<int, std::map<unsigned int, float> >::iterator iter;
  iter = packetmap.find(runnumber);
  if (iter != packetmap.end())
  {
    iter->second.clear();
  }
  else
  {
    std::map<unsigned int, float> newmap;
    packetmap[runnumber] = newmap;
    iter = packetmap.find(runnumber);
  }
  for (int i = 1; i <= pktsize_hist->GetNbinsX(); i++)
  {
    unsigned int packetid = (unsigned int) pktsize_hist->GetBinError(i);
    iter->second[packetid] = pktsize_hist->GetBinContent(i);
  }
  ExtractActivePackets(iter->second);  // make set with all currently active packets
  return 0;
}

int PktSizeMonDraw::SavePlot(const std::string &what, const std::string &type)
{
  OnlMonClient *cl = OnlMonClient::instance();
  int iret = Draw(what);
  if (iret)  // on error no ps files please
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

int PktSizeMonDraw::MakeHtml(const std::string &what)
{
  // bbc
  AddKnownBig(1002, 225);
  AddKnownBig(1003, 225);
  // mpc
  for (int i = 21101; i <= 21106; i++)
  {
    AddKnownBig(i, 300);
  }
  // mpc-ex
  // north
  for (int i = 21301; i <= 21308; i++)
  {
    AddKnownBig(i, 600);
  }
  // south
  for (int i = 21351; i <= 21358; i++)
  {
    AddKnownBig(i, 600);
  }

  // zdc
  AddKnownBig(13001, 175);
  // vtxp
  for (int i = 24001; i < 24061; i++)
  {
    AddKnownBig(i, 400);
  }
  // vtxs
  for (int i = 24101; i < 24141; i++)
  {
    AddKnownBig(i, 400);
  }
  int iret = Draw(what);
  if (iret)  // on error no html output please
  {
    return iret;
  }

  OnlMonClient *cl = OnlMonClient::instance();

  // Register the 1st canvas png file to the menu and produces the png file.
  std::string pngfile = cl->htmlRegisterPage(*this, "PktSize", "1", "png");
  cl->CanvasToPng(TC[0], pngfile);
  pngfile = cl->htmlRegisterPage(*this, "History", "2", "png");
  cl->CanvasToPng(TC[1], pngfile);

  cl->SaveLogFile(*this);

  return 0;
}

int PktSizeMonDraw::DrawHistory(const std::string & /* what */)
{
  if (!gROOT->FindObject("PktSizeMon1"))
  {
    MakeCanvas("PktSizeMon1");
  }
  OnlMonClient *cl = OnlMonClient::instance();
  TH1 *pktsize_hist = cl->getHisto("PKTSIZEMON_0", "pktsize_hist");
  if (!pktsize_hist)
  {
    DrawDeadServer(transparent[1]);
    TC[1]->Update();
    return -1;
  }
  // prevent runnumber=0/-1 from startup to screw things up
  if (cl->RunNumber() <= 0)
  {
    return 0;
  }
  if (!lastrun)
  {
    lastrun = cl->RunNumber();
    loadpreviousruns();
  }
  if (lastrun != cl->RunNumber())
  {
    // plot only physics run history
    // if the last run was physics, read it from DB and
    // enter into list
    // if not, remove the packet entries from this run from internal list
    if (rd->RunType(lastrun) == "PHYSICS")
    {
      ReplaceRunFromDB(lastrun);
      CleanOldRuns(20);
    }
    else
    {
      RemoveRun(lastrun);
    }
    lastrun = cl->RunNumber();
  }
  MakeNoisyCandidates();
  int firstrun = *(runlist.begin());
  // somehow I see run=0, no idea where this comes from
  // this led to infinite (or very long) loops in the erase
  //   if ( firstrun < 300000)
  //     {
  //       Print("RUNS");
  //       runlist.erase(runlist.begin());
  //       firstrun = *(runlist.begin());
  //     }
  TText PrintRun;
  PrintRun.SetTextFont(62);
  PrintRun.SetTextSize(0.03);
  PrintRun.SetNDC();          // set to normalized coordinates
  PrintRun.SetTextAlign(23);  // center/top alignment
  int runnumber = cl->RunNumber();
  std::ostringstream runnostream;
  time_t evttime = cl->EventTime("CURRENT");
  runnostream << "Packet Size History Run " << runnumber
              << ", Time: " << ctime(&evttime);
  transparent[1]->cd();
  PrintRun.DrawText(0.5, 0.98, runnostream.str().c_str());
  Pad[1]->cd();
  TH2 *htmp = new TH2F("noisepkts", "", 2, firstrun - 2, lastrun + 2, 2, 0, MAXSIZEDISP);
  htmp->SetStats(kFALSE);
  htmp->GetXaxis()->SetNoExponent();
  htmp->GetXaxis()->SetTitle("Run");
  htmp->GetXaxis()->SetLabelSize(0.025);
  htmp->DrawCopy();
  delete htmp;
  TGraph gr;
  TMarker tr;
  TText tx;
  std::set<unsigned int>::const_iterator iter;
  int icnt = 0;
  for (iter = noisypackets.begin(); iter != noisypackets.end(); ++iter)
  {
    PlotNoisy(gr, tr, tx, *iter, icnt);
    icnt++;
  }
  TLine tl;
  tl.DrawLine(firstrun - 2, bytelimit, lastrun + 2, bytelimit);
  TC[1]->Update();

  return 0;
}

int PktSizeMonDraw::loadpreviousruns(const int nruns)
{
  std::set<int>::const_iterator iter;
  rd->GetRunNumbers(runlist, "PHYSICS", nruns, lastrun);

  for (iter = runlist.begin(); iter != runlist.end(); ++iter)
  {
    std::map<unsigned int, float> pkts;
    FillRunPacketList(pkts, *iter);
    packetmap[*iter] = pkts;
  }
  return 0;
}

int PktSizeMonDraw::AddKnownBig(const unsigned int packetid, const float maxsize)
{
  knownbig[packetid] = maxsize;
  return 0;
}

int PktSizeMonDraw::IsKnownBig(const unsigned int packetid, const float size)
{
  std::map<unsigned int, float>::const_iterator iter = knownbig.find(packetid);
  if (iter != knownbig.end())
  {
    if (iter->second > size)
    {
      return 1;
    }
  }
  return 0;
}

int PktSizeMonDraw::ReplaceRunFromDB(const int runno)
{
  std::map<unsigned int, float> pkts;
  FillRunPacketList(pkts, runno);
  packetmap[runno] = pkts;
  return 0;
}

int PktSizeMonDraw::FillRunPacketList(std::map<unsigned int, float> &pkts, const int runnumber)
{
  std::map<std::string, std::pair<unsigned int, unsigned int> >::const_iterator graniter;
  for (graniter = granulepacketlimits.begin(); graniter != granulepacketlimits.end(); ++graniter)
  {
    db->GetPacketContent(pkts, runnumber, graniter->first);
  }
  return 0;
}

int PktSizeMonDraw::CleanOldRuns(const unsigned int maxrun)
{
  if (verbosity > 0)
  {
    std::cout << "Size before cleanup: " << packetmap.size() << std::endl;
  }
  while (packetmap.size() > maxrun + 1)  // current run is in this map -> maxrun+1
  {
    int runno = (packetmap.begin())->first;
    runlist.erase(runno);
    packetmap.erase(packetmap.begin());
  }
  if (verbosity > 0)
  {
    std::cout << "Size after cleanup: " << packetmap.size() << std::endl;
  }
  return 0;
}

int PktSizeMonDraw::RemoveRun(const int runno)
{
  std::map<int, std::map<unsigned int, float> >::iterator iter;
  iter = packetmap.find(runno);
  if (iter != packetmap.end())
  {
    packetmap.erase(iter);
    runlist.erase(runno);
  }
  else
  {
    std::cout << "Could not find run " << runno << " in list of runs" << std::endl;
  }
  return 0;
}

void PktSizeMonDraw::Print(const std::string &what) const
{
  if (what == "PACKETS" || what == "ALL")
  {
    std::map<int, std::map<unsigned int, float> >::const_iterator iter;
    std::map<unsigned int, float>::const_iterator piter;
    for (iter = packetmap.begin(); iter != packetmap.end(); ++iter)
    {
      for (piter = iter->second.begin(); piter != iter->second.end(); ++piter)
      {
        std::cout << "Run " << iter->first
                  << ", packetid: " << piter->first
                  << ", size: " << piter->second
                  << std::endl;
      }
    }
  }
  if (what == "RUNS" || what == "ALL")
  {
    std::set<int>::const_iterator riter;
    std::cout << "List of runs in Run set: " << std::endl;
    for (riter = runlist.begin(); riter != runlist.end(); ++riter)
    {
      std::cout << "Run " << *riter << std::endl;
    }
  }
  return;
}

int PktSizeMonDraw::ExtractActivePackets(const std::map<unsigned int, float> &packetsize)
{
  std::map<unsigned int, float>::const_iterator piter;
  for (piter = packetsize.begin(); piter != packetsize.end(); ++piter)
  {
    activepackets.insert(piter->first);
  }
  return 0;
}

int PktSizeMonDraw::MakeNoisyCandidates()
{
  noisypackets.clear();
#ifdef CURRENTNOISY
  OnlMonClient *cl = OnlMonClient::instance();
  int runnumber = cl->RunNumber();
  std::map<int, std::map<unsigned int, float> >::const_iterator iter;
  std::map<unsigned int, float>::const_iterator piter;
  iter = packetmap.find(runnumber);
  if (iter != packetmap.end())
  {
    for (piter = iter->second.begin(); piter != iter->second.end(); ++piter)
    {
      if (piter->second > bytelimit)
      {
        if (!IsKnownBig(piter->first, piter->second))
        {
          noisypackets.insert(piter->first);
        }
      }
    }
  }

#else
  std::set<unsigned int>::const_iterator piter;
  std::map<int, std::map<unsigned int, float> >::const_iterator siter;
  std::map<unsigned int, float>::const_iterator psizeiter;
  for (piter = activepackets.begin(); piter != activepackets.end(); ++piter)
  {
    float size = -1;
    for (siter = packetmap.begin(); siter != packetmap.end(); ++siter)
    {
      psizeiter = siter->second.find(*piter);
      if (psizeiter != siter->second.end())
      {
        if (size < 0)
        {
          size = psizeiter->second;
        }
        else
        {
          if (psizeiter->second > 200 && fabs(psizeiter->second - size) > size / 5.)
          {
            if (verbosity > 0)
            {
              std::cout << "Adding noisy packet " << *piter
                        << " old size " << size
                        << " current size " << psizeiter->second
                        << " current Run " << siter->first
                        << std::endl;
            }
            noisypackets.insert(*piter);
          }
          size = psizeiter->second;
        }
      }
    }
  }
#endif
  return noisypackets.size();
}

int PktSizeMonDraw::PlotNoisy(TGraph &gr, TMarker &tr, TText &tx, const unsigned int ipkt, const unsigned int icnt)
{
  int imarker = 20;
  int icol = 1;
  std::map<int, std::map<unsigned int, float> >::const_iterator siter;
  std::map<unsigned int, float>::const_iterator psizeiter;
  std::vector<double> runno, meansize;

  for (siter = packetmap.begin(); siter != packetmap.end(); ++siter)
  {
    psizeiter = siter->second.find(ipkt);
    if (psizeiter != siter->second.end())
    {
      runno.push_back((double) siter->first);
      meansize.push_back(psizeiter->second);
    }
  }
  double *x = new double[runno.size()];
  double *y = new double[runno.size()];
  for (unsigned int i = 0; i < runno.size(); i++)
  {
    x[i] = runno[i];
    y[i] = meansize[i];
  }
  icol += icnt - (icnt / 9) * 9;
  imarker += icnt / 9;
  gr.SetMarkerStyle(imarker);
  gr.SetMarkerColor(icol);
  gr.SetLineColor(icol);
  gr.SetMarkerSize(1);
  gr.DrawGraph(runno.size(), x, y, "LP");
  tr.SetNDC();
  tr.SetMarkerStyle(imarker);
  tr.SetMarkerColor(icol);
  tr.DrawMarker(x[0], (MAXSIZEDISP - 100) - 20 * icnt);
  std::ostringstream pktstring;
  pktstring << ipkt;
  tx.SetTextColor(icol);
  tx.SetTextSize(0.03);
  double txtXcoord = x[0] + (lastrun - x[0]) / 70;
  tx.SetTextAlign(12);
  tx.DrawText(txtXcoord, (MAXSIZEDISP - 100) - 20 * icnt, pktstring.str().c_str());
  delete[] x;
  delete[] y;
  return 0;
}
