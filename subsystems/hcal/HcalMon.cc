// use #include "" only for your local include and put
// those in the first line(s) before any #include <>
// otherwise you are asking for weird behavior
// (more info - check the difference in include path search when using "" versus <>)

#include "HcalMon.h"
#include "pseudoRunningMean.h"

#include <onlmon/OnlMon.h>  // for OnlMon
#include <onlmon/OnlMonDB.h>
#include <onlmon/OnlMonServer.h>

#include <Event/Event.h>
#include <Event/EventTypes.h>
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

const int depth = 50;
const int hist_range = 10;
const int n_packet = 32;
const int n_channel = 48;

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

  // make the per-packet running mean objects
  // 32 packets and 48 channels for hcal detectors
  for ( int i = 0; i < n_packet; i++) {
    rm_vector.push_back( new pseudoRunningMean(n_channel,depth));
  }

  return 0;
}

int HcalMon::BeginRun(const int /* runno */)
{
  // if you need to read calibrations on a run by run basis
  // this is the place to do it
  std::vector<runningMean*>::iterator rm_it;
  for ( rm_it = rm_vector.begin(); rm_it != rm_vector.end(); ++rm_it) {
    (*rm_it)->Reset();
  }

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
  float norm = 0;
  for (int s = 3;  s < p->iValue(0,"SAMPLES"); s++) {
      norm++;
      signal += p->iValue(s,channel) - baseline;
  }
  signal /= norm;

  return signal;
}

int HcalMon::process_event(Event *e /* evt */)
{
  //evtcnt++;
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
  
  if ( e->getEvtType() == BEGRUNEVENT) { // see what kind of run this is, LED or Physics
    Packet *p961 = e->getPacket(961);
    if ( p961) { // this is only printing a message
      p961->dump();
      delete p961;
    }
    Packet *p962 = e->getPacket(962);
    if ( p962) {   // we extract the flag 0 = Physics, 1= LED, more can be defined
      switch (p962->iValue(0) ) {
        case 0:
          runtypestr = "Physics";
          break;
        case 1:
          runtypestr = "LED";
          break;
        default:
          runtypestr = "Unknown";
          break;
        }
      delete p962;
    }
    return 0; // ends the process_event here if it is a begin run event
  }

  // set event count not by number of loops through process_event but by the event sequence number 
  // this accounts for processing begin run events 
  if ( e->getEvtType() == 1) {
    h2_hcal_hits->Reset();
    h2_hcal_mean->Reset();
  }
  evtcnt = e->getEvtSequence();
  h2_hcal_hits->SetTitle(TString::Format("Hcal Hits Run %d Event %d RunType %s",e->getRunNumber(), e->getEvtSequence(), runtypestr.c_str()));
  h2_hcal_mean->SetTitle(TString::Format("Hcal Running Mean Run %d Event %d RunType %s",e->getRunNumber(), e->getEvtSequence(), runtypestr.c_str()));

  for (int packet = 8001; packet < 8033; packet++) {
    Packet *p = e->getPacket(packet);
    if (p) {
      int ip = p->getIdentifier() - 8001; // packet number indexed [0,32)
      double signalvector[n_channel] = {0.};
      for (int c = 0; c < n_channel; c++) {
        signalvector[c] = getSignal(p,c);
      }
      rm_vector[ip]->Add(signalvector);

      for (int c = 0; c < n_channel; c++) {
        double signal =  getSignal(p,c);
        double phi_bin =  (2 * ip + hcal_phybin[c] + 0.5) ;
        double eta_bin = hcal_etabin[c] + 0.5;
        int ih = int(64 * (eta_bin - 0.5) + (phi_bin - 0.5)); // 1D ieta+iphi indexing

        // fill hits histogram
        if (signal > 100) h2_hcal_hits->Fill(eta_bin,phi_bin);

        // fill running mean histogram
        h2_hcal_mean->SetBinContent(h2_hcal_mean->FindBin(eta_bin,phi_bin),rm_vector[ip]->getMean(c));

        // fill time history histogram
        if (evtcnt <= hist_range) {
          h2_hcal_history->SetBinContent(h2_hcal_history->FindBin(ih,evtcnt),rm_vector[ip]->getMean(c));
          h2_hcal_history->SetTitle(TString::Format("%d",evtcnt));
        } else {
          for (int ihr = 1; ihr < hist_range; ihr++) {
            // move bin content over from index i to index i-1 for events > time history range
            h2_hcal_history->SetBinContent(h2_hcal_history->FindBin(ih,ihr),h2_hcal_history->GetBinContent(h2_hcal_history->FindBin(ih,ihr+1)));
          }
          h2_hcal_history->SetTitle(TString::Format("%d",evtcnt));
          h2_hcal_history->SetBinContent(h2_hcal_history->FindBin(ih,hist_range),rm_vector[ip]->getMean(c));
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
      // dbvars->SetVar("hcalmonnew", (float) se->Trigger(), 10000. / se->CurrentTicks(), (float) evtcnt);
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
  // varname = "hcalmonval_0_63";
  // dbvars->registerVar(varname);
  if (verbosity > 0)
  {
    dbvars->Print();
  }
  dbvars->DBInit();
  return 0;
}
