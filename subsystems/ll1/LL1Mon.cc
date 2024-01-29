// use #include "" only for your local include and put
// those in the first line(s) before any #include <>
// otherwise you are asking for weird behavior
// (more info - check the difference in include path search when using "" versus <>)

#include "LL1Mon.h"
#include "ll1setup.cc"

#include <onlmon/OnlMon.h>  // for OnlMon
#include <onlmon/OnlMonDB.h>
#include <onlmon/OnlMonServer.h>

#include <Event/msg_profile.h>

#include <TH1.h>
#include <TH2.h>
#include <TRandom.h>

#include <cmath>
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

LL1Mon::LL1Mon(const std::string &name)
  : OnlMon(name)
{
  // leave ctor fairly empty, its hard to debug if code crashes already
  // during a new Bbcll1Mon()
  return;
}

LL1Mon::~LL1Mon()
{
  // you can delete NULL pointers it results in a NOOP (No Operation)
  return;
}

int LL1Mon::Init()
{
  gRandom->SetSeed(rand());
  // read our calibrations from LL1MonData.dat
  const char *LL1calib = getenv("LL1CALIB");
  if (!LL1calib)
  {
    std::cout << "LL1CALIB environment variable not set" << std::endl;
    exit(1);
  }
  std::string fullfile = std::string(LL1calib) + "/" + "LL1MonData.dat";
  std::ifstream calib(fullfile);
  calib.close();
  // use printf for stuff which should go the screen but not into the message
  // system (all couts are redirected)
  printf("doing the Init\n");
  h_line_up = new TH2D("h_line_up",";Sample;Channel", 20, -0.5, 19.5, 60, -0.5, 59.5);
  h_nhit_corr = new TH2D("h_nhit_corr",";N_{hit}^{north};N_{hit}^{south}",nhitbins, binstart, binend, nhitbins, binstart, binend);
  //h_nhit_n_corr = new TH2D("h_nhit_n_corr",";N_{hit}^{N1};N_{hit}^{N2}",nhitbins/2 + 1,binstart,binend2,nhitbins/2+1,binstart,binend2);
  //h_nhit_s_corr = new TH2D("h_nhit_s_corr",";N_{hit}^{S1};N_{hit}^{S2}",nhitbins/2 + 1,binstart,binend2, nhitbins/2+1,binstart,binend2);
  h_nhit_n1 = new TH1D("h_nhit_n1",";N_{hit}^{N1};counts", nhitbins/2 + 1, binstart,binend2);
  h_nhit_n2 = new TH1D("h_nhit_n2",";N_{hit}^{N2};counts", nhitbins/2 + 1, binstart,binend2);
  h_nhit_s1 = new TH1D("h_nhit_s1",";N_{hit}^{S1};counts", nhitbins/2 + 1, binstart,binend2);
  h_nhit_s2 = new TH1D("h_nhit_s2",";N_{hit}^{S2};counts", nhitbins/2 + 1, binstart,binend2);
  //h_hit_check = new TH2D("h_hit_check",";N_{hit}^{north} sample index;N_{hit}^{south} sample index",3,8,11,3,8,11);
  //h_time_diff = new TH1D("h_time_diff",";time_diff;counts", 50, -0.5, 300);

  OnlMonServer *se = OnlMonServer::instance();
  // register histograms with server otherwise client won't get them
  se->registerHisto(this, h_line_up); 
  se->registerHisto(this, h_nhit_corr);  
  se->registerHisto(this, h_nhit_n1);  
  se->registerHisto(this, h_nhit_n2);  
  se->registerHisto(this, h_nhit_s1);  
  se->registerHisto(this, h_nhit_s2);  
  Reset();
  return 0;
}

int LL1Mon::BeginRun(const int /* runno */)
{
  // if you need to read calibrations on a run by run basis
  // this is the place to do it
  return 0;
}

int LL1Mon::process_event(Event * evt )
{
  evtcnt++;
  //   int ibd = 0;
  LL1HEADER *ll1h = new LL1HEADER();
  ll1setup(evt,ll1h);

  Packet *p902 = evt->getPacket(920);
  if (p902)
  {
    thresh = p902->iValue(0);
    delete p902;
  }

  int id=9;
  if(ll1h->nhit_n[9] >= thresh && ll1h->nhit_s[9]>=thresh){
    id=9;
  }
  else if (ll1h->nhit_n[ll1h->idxhitn] >= thresh && ll1h->nhit_s[ll1h->idxhitn] >= thresh){
    id = ll1h->idxhitn;
  }
  else if (ll1h->nhit_n[ll1h->idxhits] >= thresh && ll1h->nhit_s[ll1h->idxhits] >= thresh){
    id = ll1h->idxhits;
  }

  h_nhit_corr->Fill(ll1h->nhit_n[id], ll1h->nhit_s[id]);
  h_nhit_n1->Fill(ll1h->nhit_n1[id]);
  h_nhit_n2->Fill(ll1h->nhit_n2[id]);
  h_nhit_s1->Fill(ll1h->nhit_s1[id]);
  h_nhit_s2->Fill(ll1h->nhit_s2[id]);
  for (int is = 0; is < nSamples; is++)
  { 
    for (int ic = 0; ic < nChannels; ic++)
    {
      h_line_up->Fill(is, 60 - ic, ll1h->channel[ic][is]);
    }
    for (int ic = 0; ic < 8; ic++)
    {
      h_line_up->Fill(is, 8 - ic, ll1h->triggerwords[ic][is]);
    }
  }
  delete ll1h;

  return 0;
}

int LL1Mon::Reset()
{
  // reset our internal counters
  evtcnt = 0;
  idummy = 0;
  return 0;
}

