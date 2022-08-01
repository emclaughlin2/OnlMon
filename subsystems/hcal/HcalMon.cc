// use #include "" only for your local include and put
// those in the first line(s) before any #include <>
// otherwise you are asking for weird behavior
// (more info - check the difference in include path search when using "" versus <>)

#include "HcalMon.h"

#include <onlmon/OnlMon.h>  // for OnlMon
#include <onlmon/OnlMonDB.h>
#include <onlmon/OnlMonServer.h>

#include <Event/Event.h>
#include <Event/msg_profile.h>

#include <TH1.h>
#include <TH2.h>

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

HcalMon::HcalMon(const std::string &name)
  : OnlMon(name)
{
  // leave ctor fairly empty, its hard to debug if code crashes already
  // during a new HcalMon()
  return;
}

HcalMon::~HcalMon()
{
  // you can delete NULL pointers it results in a NOOP (No Operation)
  delete dbvars;
  return;
}

const int hcal_etabin[] = 
  { 1,1,0,0,3,3,2,2,
    5,5,4,4,7,7,6,6,
    9,9,8,8,11,11,10,10,
    13,13,12,12,15,15,14,14,
    17,17,16,16,19,19,18,18,
    21,21,20,20,23,23,22,22 };

const int hcal_phybin[] =
  { 1,0,1,0,1,0,1,0,
    1,0,1,0,1,0,1,0,
    1,0,1,0,1,0,1,0,
    1,0,1,0,1,0,1,0,
    1,0,1,0,1,0,1,0,
    1,0,1,0,1,0,1,0 };

const int depth = 5;
const int x_range = 10;

int HcalMon::Init()
{
  // read our calibrations from HcalMonData.dat
  std::string fullfile = std::string(getenv("HCALCALIB")) + "/" + "HcalMonData.dat";
  std::ifstream calib(fullfile);
  calib.close();
  // use printf for stuff which should go the screen but not into the message
  // system (all couts are redirected)
  printf("doing the Init\n");

  h2_hcal_hits = new TH2F("h2_hcal_hits","",24,0,24,64,0,64);
  h2_hcal_mean = new TH2F("h2_hcal_mean","",24,0,24,64,0,64);
  h2_hcal_history = new TH2F("h2_hcal_history","",1536,-0.5,1535.5,10,0.5,10.5);

  OnlMonServer *se = OnlMonServer::instance();
  // register histograms with server otherwise client won't get them
  se->registerHisto(this, h2_hcal_hits);
  se->registerHisto(this, h2_hcal_mean);
  se->registerHisto(this, h2_hcal_history);

  dbvars = new OnlMonDB(ThisName);  // use monitor name for db table name
  DBVarInit();
  Reset();
  return 0;
}

double HcalMon::getSignal(Packet *p, const int channel)
{
  double baseline = 0;
  for (int s = 0;  s < 3; s++) {
      baseline += p->iValue(s,channel);
  }
  baseline /= 3.;

  double signal = 0;
  float x = 0;
  for (int s = 3;  s < p->iValue(0,"SAMPLES"); s++) {
      x++;
      signal += p->iValue(s,channel) - baseline;
  }
  signal /= x;

  return signal;
}


int HcalMon::BeginRun(const int /* runno */)
{
  // if you need to read calibrations on a run by run basis
  // this is the place to do it
  return 0;
}

int HcalMon::process_event(Event *e /* evt */)
{
  evtcnt++;
  OnlMonServer *se = OnlMonServer::instance();
  // using ONLMONBBCLL1 makes this trigger selection configurable from the outside
  // e.g. if the BBCLL1 has problems or if it changes its name
  if (!se->Trigger("ONLMONBBCLL1"))
  {
    std::ostringstream msg;
    msg << "Processing Event " << evtcnt
        << ", Trigger : 0x" << std::hex << se->Trigger()
        << std::dec;
    // severity levels and id's for message sources can be found in
    // $ONLINE_MAIN/include/msg_profile.h
    // The last argument is a message type. Messages of the same type
    // are throttled together, so distinct messages should get distinct
    // message types
    se->send_message(this, MSG_SOURCE_UNSPECIFIED, MSG_SEV_INFORMATIONAL, msg.str(), TRGMESSAGE);
  }

  for (int packet = 8001; packet < 8033; packet++) {
    Packet *p = e->getPacket(packet);
    if (p) {
      int ip = p->getIdentifier() - 8001;

      for (int c = 0; c < 48; c++) {
        float signal =  getSignal(p,c);
        double phi_bin =  (2 * ip + hcal_phybin[c] + 0.5) ;
        double eta_bin = hcal_etabin[c] + 0.5;
        int ih = int(64 * (eta_bin - 0.5) + (phi_bin - 0.5)); // 1D ieta+iphi indexing

        // fill hits histogram
        if (signal > 100) {
          h2_hcal_hits->Fill(eta_bin,phi_bin);
        }

        // fill running mean histogram
        if (evtcnt == 1) h2_hcal_mean->SetBinContent(h2_hcal_mean->FindBin(eta_bin,phi_bin),0.0);
        float running_mean = 0.0;
        float old_mean = h2_hcal_mean->GetBinContent(h2_hcal_mean->FindBin(eta_bin,phi_bin));
        if (evtcnt <= depth) {
          running_mean = old_mean * (evtcnt - 1) + signal;
        } else {
          running_mean = old_mean * (depth - 1) + signal;
        }
        h2_hcal_mean->SetBinContent(h2_hcal_mean->FindBin(eta_bin,phi_bin),running_mean);
          
        // fill time history histogram
        if (evtcnt <= x_range) {
          h2_hcal_history->SetBinContent(h2_hcal_history->FindBin(ih,evtcnt),signal);
          h2_hcal_history->SetTitle(TString::Format("%d",evtcnt));
        } else {
          for (int ix = 0; ix < 10; ix++) {
            // move bin content over from index i to index i-1
            h2_hcal_history->SetBinContent(h2_hcal_history->FindBin(ih,ix),h2_hcal_history->GetBinContent(h2_hcal_history->FindBin(ih,ix+1)));
            h2_hcal_history->SetTitle(TString::Format("%d",evtcnt));
          }
          h2_hcal_history->SetBinContent(h2_hcal_history->FindBin(ih,10),signal);
        }

      }// channel loop

      delete p;

    }// if packet good

  }// packet loop
  
  if (idummy++ > 10)
  {
    if (dbvars)
    {
      dbvars->SetVar("hcalmoncount", (float) evtcnt, 0.1 * evtcnt, (float) evtcnt);
      dbvars->SetVar("hcalmondummy", sin((double) evtcnt), cos((double) se->Trigger()), (float) evtcnt);
      dbvars->SetVar("hcalmonnew", (float) se->Trigger(), 10000. / se->CurrentTicks(), (float) evtcnt);
      dbvars->DBcommit();
    }
    std::ostringstream msg;
    msg << "Filling Histos";
    se->send_message(this, MSG_SOURCE_UNSPECIFIED, MSG_SEV_INFORMATIONAL, msg.str(), FILLMESSAGE);
    idummy = 0;
  }
  
  return 0;
}

int HcalMon::Reset()
{
  // reset our internal counters
  evtcnt = 0;
  idummy = 0;
  return 0;
}

int HcalMon::DBVarInit()
{
  // variable names are not case sensitive
  
  std::string varname;
  varname = "hcalmoncount";
  dbvars->registerVar(varname);
  varname = "hcalmondummy";
  dbvars->registerVar(varname);
  varname = "hcalmonnew";
  dbvars->registerVar(varname);
  if (verbosity > 0)
  {
    dbvars->Print();
  }
  dbvars->DBInit();
  return 0;
}
