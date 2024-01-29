// use #include "" only for your local include and put
// those in the first line(s) before any #include <>
// otherwise you are asking for weird behavior
// (more info - check the difference in include path search when using "" versus <>)

#include "BbcMon.h"
#include "BbcMonDefs.h"

#include <onlmon/OnlMon.h>
#include <onlmon/OnlMonDB.h>
#include <onlmon/OnlMonServer.h>

#include <Event/msg_profile.h>

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/packet.h>

#include "OnlBbcEvent.h"
#include "OnlBbcSig.h"

#include <bbc/BbcGeomV1.h>

#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TH2Poly.h>
#include <TString.h>
#include <TGraphErrors.h>

#include <cmath>
#include <cstdio>   // for printf
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>   // for allocator, string, char_traits

enum
{
  TRGMESSAGE = 1,
  FILLMESSAGE = 2
};

BbcMon::BbcMon(const std::string &name)
  : OnlMon(name)
{
  // leave ctor fairly empty, its hard to debug if code crashes already
  // during a new BbcMon()
  return;
}

BbcMon::~BbcMon()
{
  delete bevt;
  delete bbcgeom;
  return;
}

int BbcMon::Init()
{
  // use printf for stuff which should go the screen but not into the message
  // system (all couts are redirected)
  printf("BbcMon::Init()\n");

  bevt = new OnlBbcEvent();

  /*
  // read our calibrations from BbcMonData.dat
  const char *bbccalib = getenv("BBCCALIB");
  if (!bbccalib)
  {
    std::cout << "BBCCALIB environment variable not set" << std::endl;
    exit(1);
  }
  std::string gainfile = std::string(bbccalib) + "/" + "bbc_mip.calib";
  Read_Charge_Calib( gainfile );
  */

  // Book Histograms
  std::ostringstream name, title;

  // TDC Distribution ----------------------------------------------------

  name << "bbc_tdc";
  title << "BBC Raw TDC Distribution";
  bbc_tdc = new TH2F(name.str().c_str(), title.str().c_str(),
                     nPMT_BBC, -.5, nPMT_BBC - .5,
                     bbc_onlmon::nBIN_TDC, 0, bbc_onlmon::tdc_max_overflow * bbc_onlmon::TDC_CONVERSION_FACTOR);
  name.str("");
  title.str("");

  // TDC Overflow Deviation ----------------------------------------------
  name << "bbc_tdc_overflow";
  title << "BBC/MBD TDC Overflow Deviation";
  bbc_tdc_overflow = new TH2F(name.str().c_str(), title.str().c_str(),
                              nPMT_BBC, -.5, nPMT_BBC - .5,
                              int(bbc_onlmon::VIEW_OVERFLOW_MAX - bbc_onlmon::VIEW_OVERFLOW_MIN + 1),
                              bbc_onlmon::VIEW_OVERFLOW_MIN - .5, bbc_onlmon::VIEW_OVERFLOW_MAX + .5);
  name.str("");
  title.str("");

  // TDC Overflow Distribution for each PMT ------------------------------
  for (int ipmt = 0; ipmt < nPMT_BBC; ipmt++)
  {
    name << "bbc_tdc_overflow_" << std::setw(3) << std::setfill('0') << ipmt;
    title << "BBC/MBD TDC Overflow Deviation of #" << std::setw(3) << std::setfill('0') << ipmt;
    bbc_tdc_overflow_each[ipmt] = new TH1F(name.str().c_str(), title.str().c_str(),
                                           int(bbc_onlmon::VIEW_OVERFLOW_MAX - bbc_onlmon::VIEW_OVERFLOW_MIN + 1),
                                           bbc_onlmon::VIEW_OVERFLOW_MIN, bbc_onlmon::VIEW_OVERFLOW_MAX);
    name.str("");
    title.str("");
  }

  // ADC Distribution --------------------------------------------------------

  bbc_adc = new TH2F("bbc_adc", "BBC/MBD ADC(Charge) Distribution", nPMT_BBC, -.5, nPMT_BBC - .5, bbc_onlmon::nBIN_ADC, 0, bbc_onlmon::MAX_ADC_MIP);

  bbc_tdc_armhittime = new TH2F("bbc_tdc_armhittime", "Arm-Hit-Time Correlation of North and South BBC/MBD",
                                64, bbc_onlmon::min_armhittime, bbc_onlmon::max_armhittime,
                                64, bbc_onlmon::min_armhittime, bbc_onlmon::max_armhittime);
  bbc_tdc_armhittime->GetXaxis()->SetTitle("South [ns]");
  bbc_tdc_armhittime->GetYaxis()->SetTitle("North [ns]");

  // Vertex Distributions --------------------------------------------------------

  bbc_zvertex = new TH1F("bbc_zvertex", "BBC/MBD ZVertex (all trigs, wide)", 128, bbc_onlmon::min_zvertex, bbc_onlmon::max_zvertex);
  bbc_zvertex->GetXaxis()->SetTitle("ZVertex [cm]");
  bbc_zvertex->GetYaxis()->SetTitle("Number of Event");
  bbc_zvertex->GetXaxis()->SetTitleSize(0.05);
  bbc_zvertex->GetYaxis()->SetTitleSize(0.05);
  bbc_zvertex->GetXaxis()->SetTitleOffset(0.70);
  bbc_zvertex->GetYaxis()->SetTitleOffset(1.75);

  bbc_zvertex_bbll1 = new TH1F("bbc_zvertex_bbll1", "BBC/MBD ZVertex (All triggers)",
                               bbc_onlmon::zvtnbin, bbc_onlmon::min_zvertex, bbc_onlmon::max_zvertex);
  bbc_zvertex_bbll1->Sumw2();
  bbc_zvertex_bbll1->GetXaxis()->SetTitle("ZVertex [cm]");
  bbc_zvertex_bbll1->GetYaxis()->SetTitle("Number of Event");
  bbc_zvertex_bbll1->GetXaxis()->SetTitleSize(0.05);
  bbc_zvertex_bbll1->GetYaxis()->SetTitleSize(0.05);
  bbc_zvertex_bbll1->GetXaxis()->SetTitleOffset(0.70);
  bbc_zvertex_bbll1->GetYaxis()->SetTitleOffset(1.75);

  bbc_zvertex_short = new TH1F("bbc_zvertex_short", "BBC/MBD ZVertex (All triggers), short time scale",
                               bbc_onlmon::zvtnbin, bbc_onlmon::min_zvertex, bbc_onlmon::max_zvertex);
  bbc_zvertex_short->Sumw2();
  bbc_zvertex_short->GetXaxis()->SetTitle("ZVertex [cm]");
  bbc_zvertex_short->GetYaxis()->SetTitle("Number of Event");
  bbc_zvertex_short->GetXaxis()->SetTitleSize(0.05);
  bbc_zvertex_short->GetYaxis()->SetTitleSize(0.05);
  bbc_zvertex_short->GetXaxis()->SetTitleOffset(0.70);
  bbc_zvertex_short->GetYaxis()->SetTitleOffset(1.75);

  f_zvtx = new TF1("f_zvtx", "gaus", -30., 30.);
  bbc_nevent_counter = new TH1F("bbc_nevent_counter",
                                "The nEvent Counter bin1:Total Event bin2:Collision Event bin3:Laser Event",
                                16, 0, 16);

  // bbc_tzero_zvtx = new TH2F("bbc_tzero_zvtx",
  //     "TimeZero vs ZVertex", 100, -200, 200, 110, -11, 11 );
  bbc_tzero_zvtx = new TH2F("bbc_tzero_zvtx", "TimeZero vs ZVertex", 100, -200, 200, 110, -16, 16);
  bbc_tzero_zvtx->SetXTitle("ZVertex [cm]");
  bbc_tzero_zvtx->SetYTitle("TimeZero[ns]");

  bbc_avr_hittime = new TH1F("bbc_avr_hittime", "BBC/MBD Average Hittime", 128, 0, 24);
  bbc_avr_hittime->Sumw2();
  bbc_avr_hittime->GetXaxis()->SetTitle("Avr HitTime [ns]");
  bbc_avr_hittime->GetYaxis()->SetTitle("Number of Event");
  bbc_avr_hittime->GetXaxis()->SetTitleSize(0.05);
  bbc_avr_hittime->GetYaxis()->SetTitleSize(0.05);
  bbc_avr_hittime->GetXaxis()->SetTitleOffset(0.70);
  bbc_avr_hittime->GetYaxis()->SetTitleOffset(1.75);

  // should make plots of the raw tdc for time channels
  //bbc_south_hittime = new TH1F("bbc_south_hittime", "BBC South Hittime", 164, -100, 16300);
  //bbc_north_hittime = new TH1F("bbc_north_hittime", "BBC North Hittime", 164, -100, 16300);
  bbc_south_hittime = new TH1F("bbc_south_hittime", "BBC/MBD South Hittime", 150, -15, 15);
  bbc_north_hittime = new TH1F("bbc_north_hittime", "BBC/MBD North Hittime", 150, -15, 15);

  bbc_south_hittime->GetXaxis()->SetTitle("South HitTime [ns]");
  bbc_south_hittime->GetYaxis()->SetTitle("Number of Event");
  bbc_south_hittime->GetXaxis()->SetTitleSize(0.05);
  bbc_south_hittime->GetYaxis()->SetTitleSize(0.05);
  bbc_south_hittime->GetXaxis()->SetTitleOffset(0.70);
  bbc_south_hittime->GetYaxis()->SetTitleOffset(1.75);

  bbc_north_hittime->GetXaxis()->SetTitle("North HitTime [ns]");
  bbc_north_hittime->GetYaxis()->SetTitle("Number of Event");
  bbc_north_hittime->GetXaxis()->SetTitleSize(0.05);
  bbc_north_hittime->GetYaxis()->SetTitleSize(0.05);
  bbc_north_hittime->GetXaxis()->SetTitleOffset(0.70);
  bbc_north_hittime->GetYaxis()->SetTitleOffset(1.75);

  //bbc_south_chargesum = new TH1F("bbc_south_chargesum", "BBC South ChargeSum [MIP]", 128, 0, bbc_onlmon::MAX_CHARGE_SUM);
  //bbc_north_chargesum = new TH1F("bbc_north_chargesum", "BBC North ChargeSum [MIP]", 128, 0, bbc_onlmon::MAX_CHARGE_SUM);
  bbc_south_chargesum = new TH1F("bbc_south_chargesum", "BBC South ChargeSum [AU]", 128, 0, bbc_onlmon::MAX_CHARGE_SUM);
  bbc_north_chargesum = new TH1F("bbc_north_chargesum", "BBC North ChargeSum [AU]", 128, 0, bbc_onlmon::MAX_CHARGE_SUM);

  bbc_north_chargesum->SetTitle("BBC/MBD ChargeSum [MIP]");
  bbc_north_chargesum->GetXaxis()->SetTitle("ChargeSum [MIP]");
  // bbc_north_chargesum->GetXaxis()->SetTitle("North ChargeSum [MIP]");
  bbc_north_chargesum->GetYaxis()->SetTitle("Number of Event");
  bbc_north_chargesum->GetXaxis()->SetTitleSize(0.05);
  bbc_north_chargesum->GetYaxis()->SetTitleSize(0.05);
  bbc_north_chargesum->GetXaxis()->SetTitleOffset(0.70);
  bbc_north_chargesum->GetYaxis()->SetTitleOffset(1.75);

  // bbc_south_chargesum->GetXaxis()->SetTitle("South ChargeSum [MIP]");
  bbc_south_chargesum->GetYaxis()->SetTitle("Number of Event");
  bbc_south_chargesum->GetXaxis()->SetTitleSize(0.05);
  bbc_south_chargesum->GetYaxis()->SetTitleSize(0.05);
  bbc_south_chargesum->GetXaxis()->SetTitleOffset(0.70);
  bbc_south_chargesum->GetYaxis()->SetTitleOffset(1.75);

  // scale down factor for each trigger
  bbc_prescale_hist = new TH1F("bbc_prescale_hist", "", 100, 0, 100);

  // waveforms
  bbc_time_wave = new TH2F("bbc_time_wave","time waveforms by ch",bbc_onlmon::BBC_NSAMPLES,-0.5,bbc_onlmon::BBC_NSAMPLES-0.5,128,0,128);
  bbc_charge_wave = new TH2F("bbc_charge_wave","charge waveforms by ch",bbc_onlmon::BBC_NSAMPLES,-0.5,bbc_onlmon::BBC_NSAMPLES-0.5,128,0,128);

  bbc_time_wave->GetXaxis()->SetTitle("Sample");
  bbc_time_wave->GetYaxis()->SetTitle("Ch");
  bbc_time_wave->GetXaxis()->SetTitleSize(0.05);
  bbc_time_wave->GetYaxis()->SetTitleSize(0.05);
  bbc_time_wave->GetXaxis()->SetTitleOffset(0.70);
  bbc_time_wave->GetYaxis()->SetTitleOffset(1.75);
  bbc_charge_wave->GetXaxis()->SetTitle("Sample");
  bbc_charge_wave->GetYaxis()->SetTitle("Ch");
  bbc_charge_wave->GetXaxis()->SetTitleSize(0.05);
  bbc_charge_wave->GetYaxis()->SetTitleSize(0.05);
  bbc_charge_wave->GetXaxis()->SetTitleOffset(0.70);
  bbc_charge_wave->GetYaxis()->SetTitleOffset(1.75);

  // hitmaps
  bbc_south_hitmap = new TH2Poly();
  bbc_south_hitmap->SetName( "bbc_south_hitmap" );
  bbc_south_hitmap->SetTitle( "BBC/MBD South Hitmap" );
  bbc_south_hitmap->GetXaxis()->SetTitle("x (cm)");
  bbc_south_hitmap->GetYaxis()->SetTitle("y (cm)");
  bbc_south_hitmap->GetXaxis()->SetTitleSize(0.05);
  bbc_south_hitmap->GetYaxis()->SetTitleSize(0.05);
  bbc_south_hitmap->GetXaxis()->SetTitleOffset(0.70);
  bbc_south_hitmap->GetYaxis()->SetTitleOffset(0.70);
  bbc_south_hitmap->SetMinimum(0.1);

  bbc_north_hitmap = new TH2Poly();
  bbc_north_hitmap->SetName( "bbc_north_hitmap" );
  bbc_north_hitmap->SetTitle( "BBC/MBD North Hitmap" );
  bbc_north_hitmap->GetXaxis()->SetTitle("x (cm)");
  bbc_north_hitmap->GetYaxis()->SetTitle("y (cm)");
  bbc_north_hitmap->GetXaxis()->SetTitleSize(0.05);
  bbc_north_hitmap->GetYaxis()->SetTitleSize(0.05);
  bbc_north_hitmap->GetXaxis()->SetTitleOffset(0.70);
  bbc_north_hitmap->GetYaxis()->SetTitleOffset(0.70);
  bbc_north_hitmap->SetMinimum(0.1);

  // Get the detector geometry
  bbcgeom = new BbcGeomV1{};
  Double_t x[6];    // x,y location of the 6 points of the BBC hexagonal PMT's, in cm
  Double_t y[6];
  for (int ipmt=0; ipmt<128; ipmt++)
  {
    float xcent = bbcgeom->get_x(ipmt);
    float ycent = bbcgeom->get_y(ipmt);
    int arm = bbcgeom->get_arm(ipmt);
    std::cout << ipmt << "\t" << xcent << "\t" << ycent << std::endl;

    // create hexagon
    x[0] = xcent - 0.8; // in cm
    y[0] = ycent + 1.4;
    x[1] = xcent + 0.8;
    y[1] = ycent + 1.4;
    x[2] = xcent + 1.6;
    y[2] = ycent;
    x[3] = xcent + 0.8;
    y[3] = ycent - 1.4;
    x[4] = xcent - 0.8;
    y[4] = ycent - 1.4;
    x[5] = xcent - 1.6;
    y[5] = ycent;

    if ( arm==0 )
    {
      bbc_south_hitmap->AddBin(6,x,y);
    }
    else if ( arm==1 )
    {
      bbc_north_hitmap->AddBin(6,x,y);
    }

  }


  // register histograms with server otherwise client won't get them
  OnlMonServer *se = OnlMonServer::instance();

  se->registerHisto(this, bbc_adc);
  se->registerHisto(this, bbc_tdc);
  se->registerHisto(this, bbc_tdc_overflow);
  for (auto &ipmt : bbc_tdc_overflow_each)
  {
    se->registerHisto(this, ipmt);
  }

  se->registerHisto(this, bbc_tdc_armhittime);
  se->registerHisto(this, bbc_zvertex);
  se->registerHisto(this, bbc_zvertex_bbll1);
  se->registerHisto(this, bbc_nevent_counter);
  se->registerHisto(this, bbc_tzero_zvtx);
  se->registerHisto(this, bbc_prescale_hist);
  se->registerHisto(this, bbc_avr_hittime);
  se->registerHisto(this, bbc_north_hittime);
  se->registerHisto(this, bbc_south_hittime);
  se->registerHisto(this, bbc_north_chargesum);
  se->registerHisto(this, bbc_south_chargesum);
  se->registerHisto(this, bbc_north_hitmap);
  se->registerHisto(this, bbc_south_hitmap);
  /*
  se->registerHisto(this, bbc_tmax[0]);
  se->registerHisto(this, bbc_tmax[1]);
  */
  se->registerHisto(this, bbc_time_wave);
  se->registerHisto(this, bbc_charge_wave);

  /*
  dbvars = new OnlMonDB(ThisName);  // use monitor name for db table name
  DBVarInit();
  */
  Reset();

  return 0;
}

int BbcMon::BeginRun(const int runno)
{
  // if you need to read calibrations on a run by run basis
  // this is the place to do it
  std::cout << "BbcMon::BeginRun(), run " << runno << std::endl;
  Reset();
  bevt->InitRun();

  return 0;
}

int BbcMon::EndRun(const int /* runno */)
{
  // This does nothing for now, but can put summary info here for debugging

  return 0;
}


int BbcMon::process_event(Event *evt)
{
  if (evt->getEvtType() != DATAEVENT)
  {
    return 0;
  }

  evtcnt++;
  bbc_nevent_counter->Fill(0);

  [[maybe_unused]] OnlMonServer *se = OnlMonServer::instance();

  Packet *p[2];
  p[0] = evt->getPacket(1001);
  p[1] = evt->getPacket(1002);

  // Check that we have both MBD/BBC packets
  if ( !p[0] || !p[1] )
  {
    se = OnlMonServer::instance();
    std::ostringstream msg;
    msg << "BBC/MBD packet not found" ;
    se->send_message(this,MSG_SOURCE_BBC,MSG_SEV_WARNING, msg.str(),1);
    msg.str("");
    bbc_nevent_counter->Fill(3);  // bad event, missing packets

    delete p[0];
    delete p[1];

    return 0;
  }

  // Check that both MBD/BBC packets have good checksums
  if ( (p[0]->iValue(0,"EVENCHECKSUMOK") == 0) || (p[0]->iValue(0,"ODDCHECKSUMOK") == 0) ||
       (p[1]->iValue(0,"EVENCHECKSUMOK") == 0) || (p[1]->iValue(0,"ODDCHECKSUMOK") == 0) )
  {
    se = OnlMonServer::instance();
    std::ostringstream msg;
    msg << "BBC/MBD packets have bad checksum(s)" ;
    se->send_message(this,MSG_SOURCE_BBC,MSG_SEV_WARNING, msg.str(),1);
    msg.str("");
    bbc_nevent_counter->Fill(4);  // bad event, missing packets

    delete p[0];
    delete p[1];

    return 0;
  }

  int f_evt = evt->getEvtSequence();

  // calculate BBC
  bevt->Clear();
  bevt->setRawData(evt);
  bevt->calculate();

  if ( bevt->calib_is_done() == 0 ) 
  {
    delete p[0];
    delete p[1];
    return 0;
  }

  bbc_nevent_counter->Fill(1);

  double zvtx = bevt->get_bbcz();
  double t0 = bevt->get_t0();
  double qsum[2] = {0,0};
  qsum[0] = bevt->getChargeSum(0);
  qsum[1] = bevt->getChargeSum(1);

  static int counter = 0;
  evtcnt++;
  if ( counter < 10 )
  {
    std::cout << "zt\t" << f_evt << "\t" << zvtx << "\t" << t0 << std::endl;
    counter++;
  }

  // vertex and t0
  bbc_zvertex->Fill(zvtx);
  bbc_zvertex_bbll1->Fill(zvtx);
  bbc_tzero_zvtx->Fill(zvtx, t0);
  //should this be the mean hit time, not per ipmt?
  //bbc_south_hittime->Fill( tq );
  //bbc_north_hittime->Fill( tq );

  bbc_zvertex_short->Fill(zvtx);

  int n_goodevt = bbc_nevent_counter->GetBinContent(2);
  if ( n_goodevt%1000 == 0 )
  {
    f_zvtx->SetRange( -75., 75. );
    f_zvtx->SetParameters( 250, 0., 10 );
    bbc_zvertex_short->Fit( f_zvtx, "RNQ");

    // Report z-vertex mean and width
    Double_t mean = f_zvtx->GetParameter(1);
    Double_t rms = f_zvtx->GetParameter(2);
    // we should do a check of a good fit here (skip for now)

    std::ostringstream msg;
    msg << "BBC/MBD zvertex mean/width: " << mean << " " << rms;
    se->send_message(this,MSG_SOURCE_BBC,MSG_SEV_INFORMATIONAL, msg.str(),1);
    std::cout << "BBC/MBD zvtx mean/width: " << mean << " " << rms << std::endl;

    bbc_zvertex_short->Reset();
  }

  for (int ipmt=0; ipmt<128; ipmt++)
  {
    float q = bevt->getQ(ipmt);
    bbc_adc->Fill(ipmt, q);

    //std::cout << f_evt << "\tipmt " << ipmt << "\t" << q << "\t";
    if ( q>0.5 )
    {
      float tq = bevt->getTQ(ipmt);
      if (ipmt<64 )
        bbc_south_hittime->Fill( tq );
      else
        bbc_north_hittime->Fill( tq );

      //std::cout << tq;
    }
    //std::cout << std::endl;
  }

  // charge
  bbc_south_chargesum->Fill(qsum[0]);
  bbc_north_chargesum->Fill(qsum[1]);

  // raw waveforms
  for (int ipmt=0; ipmt<128; ipmt++)
  {
    int tch = (ipmt/8)*16 + ipmt%8;
    OnlBbcSig *bbcsig = bevt->GetSig( tch );
    TGraphErrors *gwave = bbcsig->GetGraph();
    Int_t n = gwave->GetN();
    Double_t *x = gwave->GetX();
    Double_t *y = gwave->GetY();
    for (int isamp=0; isamp<n; isamp++)
    {
      bbc_time_wave->Fill( x[isamp], ipmt, y[isamp] );
    }

    int qch = tch + 8;
    bbcsig = bevt->GetSig( qch );
    gwave = bbcsig->GetGraph();
    n = gwave->GetN();
    x = gwave->GetX();
    y = gwave->GetY();
    for (int isamp=0; isamp<n; isamp++)
    {
      bbc_charge_wave->Fill( x[isamp], ipmt, y[isamp] );
    }

    // hit map
    float xcent = bbcgeom->get_x(ipmt);
    float ycent = bbcgeom->get_y(ipmt);
    int arm = bbcgeom->get_arm(ipmt);
    float q = bevt->getQ(ipmt);

    if ( arm==0 )
    {
      bbc_south_hitmap->Fill( xcent, ycent, q );
    }
    else if ( arm==1 )
    {
      bbc_north_hitmap->Fill( xcent, ycent, q );
    }

  }

  delete p[0];
  delete p[1];

  return 0;
}

int BbcMon::Reset()
{
  // reset our internal counters
  evtcnt = 0;
  //idummy = 0;

  bbc_adc->Reset();
  bbc_tdc->Reset();
  bbc_tdc_overflow->Reset();
  bbc_tdc_armhittime->Reset();
  bbc_nevent_counter->Reset();
  bbc_zvertex->Reset();
  bbc_zvertex_bbll1->Reset();
  bbc_zvertex_short->Reset();
  bbc_tzero_zvtx->Reset();
  bbc_avr_hittime->Reset();
  bbc_south_hittime->Reset();
  bbc_north_hittime->Reset();
  bbc_south_chargesum->Reset();
  bbc_north_chargesum->Reset();
  bbc_prescale_hist->Reset();
  bbc_time_wave->Reset();
  bbc_charge_wave->Reset();

  return 0;
}

int BbcMon::DBVarInit()
{
  // variable names are not case sensitive
  /*
  std::string varname;
  varname = "bbcmoncount";
  dbvars->registerVar(varname);
  varname = "bbcmondummy";
  dbvars->registerVar(varname);
  varname = "bbcmonnew";
  dbvars->registerVar(varname);
  if (verbosity > 0)
  {
    dbvars->Print();
  }
  dbvars->DBInit();
  */
  return 0;
}

