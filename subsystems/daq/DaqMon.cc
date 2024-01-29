// use #include "" only for your local include and put
// those in the first line(s) before any #include <>
// otherwise you are asking for weird behavior
// (more info - check the difference in include path search when using "" versus <>)

#include "DaqMon.h"

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

DaqMon::DaqMon(const std::string &name)
  : OnlMon(name)
{
  // leave ctor fairly empty, its hard to debug if code crashes already
  // during a new DaqMon()
  return;
}

DaqMon::~DaqMon()
{
  // you can delete NULL pointers it results in a NOOP (No Operation)
  return;
}

int DaqMon::Init()
{
  gRandom->SetSeed(rand());
  // read our calibrations from DaqMonData.dat
  const char *daqcalib = getenv("DAQCALIB");
  if (!daqcalib)
  {
    std::cout << "DAQCALIB environment variable not set" << std::endl;
    exit(1);
  }
  std::string fullfile = std::string(daqcalib) + "/" + "DaqMonData.dat";
  std::ifstream calib(fullfile);
  calib.close();
  // use printf for stuff which should go the screen but not into the message
  // system (all couts are redirected)
  printf("doing the Init\n");
  daqhist1 = new TH1F("daqmon_hist1", "test 1d histo", 101, 0., 100.);
  daqhist2 = new TH2F("daqmon_hist2", "test 2d histo", 101, 0., 100., 101, 0., 100.);
  OnlMonServer *se = OnlMonServer::instance();
  // register histograms with server otherwise client won't get them
  se->registerHisto(this, daqhist1);  // uses the TH1->GetName() as key
  se->registerHisto(this, daqhist2);
  Reset();
  return 0;
}

int DaqMon::BeginRun(const int /* runno */)
{
  // if you need to read calibrations on a run by run basis
  // this is the place to do it
  return 0;
}

int DaqMon::process_event(Event * /* evt */)
{
  evtcnt++;
  // get temporary pointers to histograms
  // one can do in principle directly se->getHisto("daqhist1")->Fill()
  // but the search in the histogram Map is somewhat expensive and slows
  // things down if you make more than one operation on a histogram
  daqhist1->Fill(gRandom->Gaus(50,10));
  daqhist2->Fill(gRandom->Gaus(50,10), gRandom->Gaus(50,10), 1.);

  return 0;
}

int DaqMon::Reset()
{
  // reset our internal counters
  evtcnt = 0;
  idummy = 0;
  return 0;
}

