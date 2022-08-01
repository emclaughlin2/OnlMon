// use #include "" only for your local include and put
// those in the first line(s) before any #include <>
// otherwise you are asking for weird behavior
// (more info - check the difference in include path search when using "" versus <>)

#include "MbdMon.h"
#include "MbdMonDefs.h"

#include <onlmon/OnlMon.h>  // for OnlMon
#include <onlmon/OnlMonDB.h>
#include <onlmon/OnlMonServer.h>

#include <Event/msg_profile.h>

#include <TH1.h>
#include <TH2.h>
#include <TRandom3.h>
#include <TString.h>

#include <cmath>
#include <cstdio>  // for printf
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>  // for allocator, string, char_traits

using namespace mbd_onlmon;

enum
{
  TRGMESSAGE = 1,
  FILLMESSAGE = 2
};

MbdMon::MbdMon(const std::string &name)
  : OnlMon(name)
{
  // leave ctor fairly empty, its hard to debug if code crashes already
  // during a new MbdMon()
  return;
}

MbdMon::~MbdMon()
{
  // you can delete NULL pointers it results in a NOOP (No Operation)
  delete dbvars;
  return;
}

int MbdMon::Init()
{
  // use printf for stuff which should go the screen but not into the message
  // system (all couts are redirected)
  printf("MbdMon::Init()\n");

  trand3 = new TRandom3(0);

  // read our calibrations from MbdMonData.dat
  std::string fullfile = std::string(getenv("MBDCALIB")) + "/" + "MbdMonData.dat";
  std::ifstream calib(fullfile);
  calib.close();

  // Book Histograms
  //mbdhist1 = new TH1F("mbdmon_hist1", "test 1d histo", 101, 0., 100.);
  //mbdhist2 = new TH2F("mbdmon_hist2", "test 2d histo", 101, 0., 100., 101, 0., 100.);
  std::ostringstream name, title;
 
  // TDC Distribution ----------------------------------------------------

  name << "bbc_tdc" ;
  title << "BBC Raw TDC Distribution" ;
  bbc_tdc = new TH2F(name.str().c_str(), title.str().c_str(),
      nPMT_MBD, -.5, nPMT_MBD - .5,
      nBIN_TDC, 0, tdc_max_overflow * TDC_CONVERSION_FACTOR );
  name.str("");
  title.str("");

  // TDC Overflow Deviation ----------------------------------------------
  name << "bbc_tdc_overflow" ;
  title << "MBD/BBC TDC Overflow Deviation" ;
  bbc_tdc_overflow = new TH2F(name.str().c_str(), title.str().c_str(),
      nPMT_MBD, -.5, nPMT_MBD - .5,
      int(VIEW_OVERFLOW_MAX - VIEW_OVERFLOW_MIN + 1),
      VIEW_OVERFLOW_MIN - .5, VIEW_OVERFLOW_MAX + .5 );
  name.str("");
  title.str("");


  // TDC Overflow Distribution for each PMT ------------------------------
  for ( int ipmt = 0 ; ipmt < nPMT_MBD ; ipmt++ )
  {
    name << "bbc_tdc_overflow_" << std::setw(3) << std::setfill('0') << ipmt ;
    title << "MBD/BBC TDC Overflow Deviation of #" << std::setw(3) << std::setfill('0') << ipmt ;
    bbc_tdc_overflow_each[ipmt] = new TH1F(name.str().c_str(), title.str().c_str(),
        int(VIEW_OVERFLOW_MAX - VIEW_OVERFLOW_MIN + 1),
        VIEW_OVERFLOW_MIN, VIEW_OVERFLOW_MAX );
    name.str("");
    title.str("");
  }

  // ADC Distribution --------------------------------------------------------

  bbc_adc = new TH2F("bbc_adc", "MBD/BBC ADC(Charge) Distribution", nPMT_MBD, -.5, nPMT_MBD - .5, nBIN_ADC, 0, MAX_ADC_MIP );

  /*
     for ( int trig = 0 ; trig < nTRIGGER ; trig++ )
     {
  // nHit ----------------------------------------------------------------

  name << "bbc_nhit_" << TRIGGER_str[trig] ;
  title << "BBC nHIT by " << TRIGGER_str[trig] ;
  bbc_nhit[trig] = new TH1D(name.str().c_str(), title.str().c_str(),
  nPMT_MBD, -.5, nPMT_MBD - .5 );
  name.str("");
  title.str("");
  }
  */

  bbc_tdc_armhittime = new TH2F("bbc_tdc_armhittime", "Arm-Hit-Time Correlation of North and South BBC",
      64, min_armhittime, max_armhittime,
      64, min_armhittime, max_armhittime );
  bbc_tdc_armhittime->GetXaxis()->SetTitle("South[ns]");
  bbc_tdc_armhittime->GetYaxis()->SetTitle("North[ns]");

  bbc_zvertex = new TH1F("bbc_zvertex", "BBC ZVertex", 128, min_zvertex, max_zvertex );
  bbc_zvertex->GetXaxis()->SetTitle("BBC Raw ZVertex [cm]");
  bbc_zvertex->GetYaxis()->SetTitle("Number of Event");
  bbc_zvertex->GetXaxis()->SetTitleSize( 0.05);
  bbc_zvertex->GetYaxis()->SetTitleSize( 0.05);
  bbc_zvertex->GetXaxis()->SetTitleOffset(0.70);
  bbc_zvertex->GetYaxis()->SetTitleOffset(1.75);

  bbc_zvertex_bbll1 = new TH1F("bbc_zvertex_bbll1", "BBC ZVertex triggered by BBLL1",
      zvtnbin, min_zvertex, max_zvertex );
  bbc_zvertex_bbll1->Sumw2();
  bbc_zvertex_bbll1->GetXaxis()->SetTitle("ZVertex [cm]");
  bbc_zvertex_bbll1->GetYaxis()->SetTitle("Number of Event");
  bbc_zvertex_bbll1->GetXaxis()->SetTitleSize( 0.05);
  bbc_zvertex_bbll1->GetYaxis()->SetTitleSize( 0.05);
  bbc_zvertex_bbll1->GetXaxis()->SetTitleOffset(0.70);
  bbc_zvertex_bbll1->GetYaxis()->SetTitleOffset(1.75);

  /*
     bbc_zvertex_zdc = new TH1F("bbc_zvertex_zdc",
     "BBC ZVertex triggered by ZDC",
     zvtnbin/2, min_zvertex, max_zvertex );
     bbc_zvertex_zdc->Sumw2();
     bbc_zvertex_zdc->GetXaxis()->SetTitle("ZVertex [cm]");
     bbc_zvertex_zdc->GetYaxis()->SetTitle("Number of Event");
     bbc_zvertex_zdc->GetXaxis()->SetTitleSize( 0.05);
     bbc_zvertex_zdc->GetYaxis()->SetTitleSize( 0.05);
     bbc_zvertex_zdc->GetXaxis()->SetTitleOffset(0.70);
     bbc_zvertex_zdc->GetYaxis()->SetTitleOffset(1.75);

     bbc_zvertex_zdc_scale3 = new TH1F("bbc_zvertex_zdc_scale3",
     "BBC ZVertex triggered by ZDC (scale 3 times up)",
     zvtnbin, min_zvertex, max_zvertex );
     bbc_zvertex_zdc_scale3->Sumw2();
     bbc_zvertex_zdc_scale3->GetXaxis()->SetTitle("ZVertex [cm]");
     bbc_zvertex_zdc_scale3->GetYaxis()->SetTitle("Number of Event");
     bbc_zvertex_zdc_scale3->GetXaxis()->SetTitleSize( 0.05);
     bbc_zvertex_zdc_scale3->GetYaxis()->SetTitleSize( 0.05);
     bbc_zvertex_zdc_scale3->GetXaxis()->SetTitleOffset(0.70);
     bbc_zvertex_zdc_scale3->GetYaxis()->SetTitleOffset(1.75);
  */

  bbc_zvertex_bbll1_novtx = new TH1F("bbc_zvertex_bbll1_novtx", "BBC ZVertex triggered by BBLL1(noVtxCut)",
      zvtnbin/2, min_zvertex, max_zvertex );
  bbc_zvertex_bbll1_novtx->Sumw2();
  bbc_zvertex_bbll1_novtx->GetXaxis()->SetTitle("ZVertex [cm]");
  bbc_zvertex_bbll1_novtx->GetYaxis()->SetTitle("Number of Event");
  bbc_zvertex_bbll1_novtx->GetXaxis()->SetTitleSize( 0.05);
  bbc_zvertex_bbll1_novtx->GetYaxis()->SetTitleSize( 0.05);
  bbc_zvertex_bbll1_novtx->GetXaxis()->SetTitleOffset(0.70);
  bbc_zvertex_bbll1_novtx->GetYaxis()->SetTitleOffset(1.75);

  bbc_zvertex_bbll1_narrowvtx = new TH1F("bbc_zvertex_bbll1_narrowvtx",//RUN11 AuAu
      "BBC ZVertex triggered by BBLL1(narrowVtxCut)",
      zvtnbin, min_zvertex, max_zvertex );
  bbc_zvertex_bbll1_narrowvtx->Sumw2();
  bbc_zvertex_bbll1_narrowvtx->GetXaxis()->SetTitle("ZVertex [cm]");
  bbc_zvertex_bbll1_narrowvtx->GetYaxis()->SetTitle("Number of Event");
  bbc_zvertex_bbll1_narrowvtx->GetXaxis()->SetTitleSize( 0.05);
  bbc_zvertex_bbll1_narrowvtx->GetYaxis()->SetTitleSize( 0.05);
  bbc_zvertex_bbll1_narrowvtx->GetXaxis()->SetTitleOffset(0.70);
  bbc_zvertex_bbll1_narrowvtx->GetYaxis()->SetTitleOffset(1.75);


  /*
     bbc_zvertex_bbll1_zdc = new TH1F("bbc_zvertex_bbll1_zdc",
     "BBC ZVertex triggered by BBLL1&ZDCNS",
  //bbc_zvertex_bbll1_zdc = new TH1F("bbc_zvertex_bbll1_zdc",
  //    "BBC ZVertex triggered by BBLL1&ZDCLL1wide",
  zvtnbin, min_zvertex, max_zvertex );
  bbc_zvertex_bbll1_zdc->Sumw2();
  bbc_zvertex_bbll1_zdc->GetXaxis()->SetTitle("ZVertex [cm]");
  bbc_zvertex_bbll1_zdc->GetYaxis()->SetTitle("Number of Event");
  bbc_zvertex_bbll1_zdc->GetXaxis()->SetTitleSize( 0.05);
  bbc_zvertex_bbll1_zdc->GetYaxis()->SetTitleSize( 0.05);
  bbc_zvertex_bbll1_zdc->GetXaxis()->SetTitleOffset(0.70);
  bbc_zvertex_bbll1_zdc->GetYaxis()->SetTitleOffset(1.75);
  */

  bbc_nevent_counter = new TH1F("bbc_nevent_counter",
      "The nEvent Counter bin1:Total Event bin2:Collision Event bin3:Laser Event",
      16, 0, 16 );

  //bbc_tzero_zvtx = new TH2F("bbc_tzero_zvtx",
  //    "TimeZero vs ZVertex", 100, -200, 200, 110, -11, 11 );
  bbc_tzero_zvtx = new TH2F("bbc_tzero_zvtx", "TimeZero vs ZVertex", 100, -200, 200, 110, -6, 16 );
  bbc_tzero_zvtx->SetXTitle("ZVertex[cm]");
  bbc_tzero_zvtx->SetYTitle("TimeZero[ns]");

  bbc_avr_hittime = new TH1F("bbc_avr_hittime", "BBC Average Hittime", 128, 0, 24 );
  bbc_south_hittime = new TH1F("bbc_south_hittime", "BBC South Hittime", 128, 0, 24 );
  bbc_north_hittime = new TH1F("bbc_north_hittime", "BBC North Hittime", 128, 0, 24 );
  bbc_south_chargesum = new TH1F("bbc_south_chargesum", "BBC South ChargeSum [MIP]", 128, 0, MAX_CHARGE_SUM );
  bbc_north_chargesum = new TH1F("bbc_north_chargesum", "BBC North ChargeSum [MIP]", 128, 0, MAX_CHARGE_SUM );
  bbc_avr_hittime->Sumw2();
  bbc_avr_hittime->GetXaxis()->SetTitle("Avr HitTime [ns]");
  bbc_avr_hittime->GetYaxis()->SetTitle("Number of Event");
  bbc_avr_hittime->GetXaxis()->SetTitleSize( 0.05);
  bbc_avr_hittime->GetYaxis()->SetTitleSize( 0.05);
  bbc_avr_hittime->GetXaxis()->SetTitleOffset(0.70);
  bbc_avr_hittime->GetYaxis()->SetTitleOffset(1.75);


  bbc_south_hittime->GetXaxis()->SetTitle("South HitTime [ns]");
  bbc_south_hittime->GetYaxis()->SetTitle("Number of Event");
  bbc_south_hittime->GetXaxis()->SetTitleSize( 0.05);
  bbc_south_hittime->GetYaxis()->SetTitleSize( 0.05);
  bbc_south_hittime->GetXaxis()->SetTitleOffset(0.70);
  bbc_south_hittime->GetYaxis()->SetTitleOffset(1.75);

  bbc_north_hittime->GetXaxis()->SetTitle("North HitTime [ns]");
  bbc_north_hittime->GetYaxis()->SetTitle("Number of Event");
  bbc_north_hittime->GetXaxis()->SetTitleSize( 0.05);
  bbc_north_hittime->GetYaxis()->SetTitleSize( 0.05);
  bbc_north_hittime->GetXaxis()->SetTitleOffset(0.70);
  bbc_north_hittime->GetYaxis()->SetTitleOffset(1.75);

  bbc_north_chargesum->SetTitle("BBC ChargeSum [MIP]");
  bbc_north_chargesum->GetXaxis()->SetTitle("ChargeSum [MIP]");
  //bbc_north_chargesum->GetXaxis()->SetTitle("North ChargeSum [MIP]");
  bbc_north_chargesum->GetYaxis()->SetTitle("Number of Event");
  bbc_north_chargesum->GetXaxis()->SetTitleSize( 0.05);
  bbc_north_chargesum->GetYaxis()->SetTitleSize( 0.05);
  bbc_north_chargesum->GetXaxis()->SetTitleOffset(0.70);
  bbc_north_chargesum->GetYaxis()->SetTitleOffset(1.75);

  //bbc_south_chargesum->GetXaxis()->SetTitle("South ChargeSum [MIP]");
  bbc_south_chargesum->GetYaxis()->SetTitle("Number of Event");
  bbc_south_chargesum->GetXaxis()->SetTitleSize( 0.05);
  bbc_south_chargesum->GetYaxis()->SetTitleSize( 0.05);
  bbc_south_chargesum->GetXaxis()->SetTitleOffset(0.70);
  bbc_south_chargesum->GetYaxis()->SetTitleOffset(1.75);

  // scale down factor for each trigger
  bbc_prescale_hist = new TH1F("bbc_prescale_hist", "", 100, 0, 100);

  // register histograms with server otherwise client won't get them
  OnlMonServer *se = OnlMonServer::instance();
  //se->registerHisto(this, mbdhist1);  // uses the TH1->GetName() as key
  //se->registerHisto(this, mbdhist2);

  se->registerHisto( this, bbc_adc );
  se->registerHisto( this, bbc_tdc );
  se->registerHisto( this, bbc_tdc_overflow );
  for ( int ipmt = 0 ; ipmt < nPMT_MBD ; ipmt++ )
  {
    se->registerHisto( this, bbc_tdc_overflow_each[ipmt] );
  }

  se->registerHisto( this, bbc_tdc_armhittime );
  se->registerHisto( this, bbc_zvertex );
  se->registerHisto( this, bbc_zvertex_bbll1 );
  //registHist( se->registerHisto( this, bbc_zvertex_zdc );
  //registHist( se->registerHisto( this, bbc_zvertex_zdc_scale3 );
  se->registerHisto( this, bbc_zvertex_bbll1_novtx );
  se->registerHisto( this, bbc_zvertex_bbll1_narrowvtx );
  //registHist( se->registerHisto( this, bbc_zvertex_bbll1_zdc );
  se->registerHisto( this, bbc_nevent_counter );
  se->registerHisto( this, bbc_tzero_zvtx );
  se->registerHisto( this, bbc_prescale_hist );
  se->registerHisto( this, bbc_avr_hittime );
  se->registerHisto( this, bbc_north_hittime );
  se->registerHisto( this, bbc_south_hittime );
  se->registerHisto( this, bbc_north_chargesum );
  se->registerHisto( this, bbc_south_chargesum );

  dbvars = new OnlMonDB(ThisName);  // use monitor name for db table name
  DBVarInit();
  Reset();

  return 0;
}

int MbdMon::BeginRun(const int /* runno */)
{
  // if you need to read calibrations on a run by run basis
  // this is the place to do it
  return 0;
}

int MbdMon::process_event(Event * /* evt */)
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


  // get temporary pointers to histograms
  // one can do in principle directly se->getHisto("mbdhist1")->Fill()
  // but the search in the histogram Map is somewhat expensive and slows
  // things down if you make more than one operation on a histogram
  double zvtx = trand3->Gaus(0, 10.0);
  double t0 = trand3->Gaus(0, 10/C);

  bbc_adc->Fill( 10, 1000 );
  bbc_zvertex->Fill( zvtx );
  bbc_tzero_zvtx->Fill( zvtx, t0 );

  if (idummy++ > 10)
  {
    if (dbvars)
    {
      dbvars->SetVar("mbdmoncount", (float) evtcnt, 0.1 * evtcnt, (float) evtcnt);
      dbvars->SetVar("mbdmondummy", sin((double) evtcnt), cos((double) se->Trigger()), (float) evtcnt);
      dbvars->SetVar("mbdmonnew", (float) se->Trigger(), 10000. / se->CurrentTicks(), (float) evtcnt);
      dbvars->DBcommit();
    }
    std::ostringstream msg;
    msg << "Filling Histos";
    se->send_message(this, MSG_SOURCE_UNSPECIFIED, MSG_SEV_INFORMATIONAL, msg.str(), FILLMESSAGE);
    idummy = 0;
  }
  return 0;
}

int MbdMon::Reset()
{
  // reset our internal counters
  evtcnt = 0;
  idummy = 0;
  return 0;
}

int MbdMon::DBVarInit()
{
  // variable names are not case sensitive
  std::string varname;
  varname = "mbdmoncount";
  dbvars->registerVar(varname);
  varname = "mbdmondummy";
  dbvars->registerVar(varname);
  varname = "mbdmonnew";
  dbvars->registerVar(varname);
  if (verbosity > 0)
  {
    dbvars->Print();
  }
  dbvars->DBInit();
  return 0;
}
