#ifndef __ONLBBCEVENT_H__
#define __ONLBBCEVENT_H__

#include "OnlBbcSig.h"

class PHCompositeNode;
class Event;
class Packet;
class TF1;
class TCanvas;

class OnlBbcEvent
{
public:
  OnlBbcEvent ();
  ~OnlBbcEvent ();

  bool setRawData(Event *);
  int calculate();
  int InitRun();
  void Clear();

  int calib_is_done() { return calib_done; }

  float get_bbcz() { return f_bbcz; }
  float get_t0() { return f_bbct0; }

  int getEventNumber (void) const
  {
    return (EventNumber);
  }

  int getAdc (int PmtIndx) const
  {
    //return (Adc[PmtIndx]);
    return (f_pmtq[PmtIndx]);
  }

  int getTdc0 (int PmtIndx) const
  {
    return (f_pmtt0[PmtIndx]);
  }

  int getTdc1 (int PmtIndx) const
  {
    return (f_pmtt1[PmtIndx]);
  }

  float getQ (int PmtIndx) const
  {
    //return (Adc[PmtIndx]);
    return (f_pmtq[PmtIndx]);
  }

  // Time from Timing Channel
  float getTT (int PmtIndx) const
  {
    return (f_pmtt0[PmtIndx]);
  }

  // Time from Charge Channel
  float getTQ (int PmtIndx) const
  {
    return (f_pmtt1[PmtIndx]);
  }

  /*
  int getCutVal ()
  {
    return (calib.getCutVal ());
  }

  BbcTime_t getEventTime (void) const
  {
    return (0);
  }

  float getTrueAdc (int PmtIndx) const
  {
    return (TrueAdc[PmtIndx]);
  }

  float getCharge (int PmtIndx) const
  {
    return (Charge[PmtIndx]);
  }
  float getHitTime0 (int PmtIndx) const
  {
    return (HitTime0[PmtIndx]);
  }
  float getHitTime1 (int PmtIndx) const
  {
    return (HitTime1[PmtIndx]);
  }
  int isHit (int PmtIndx) const
  {
    return (iHit[PmtIndx]);
  }
 
  int isArmHitPmt (int PmtIndx) const
  {
    return (armHitPmt[PmtIndx]);
  }

  int getnHitPmt (Bbc::ArmType arm) const
  {
    return (nHitPmt[arm]);
  }
  */

  float getChargeSum (int arm) const
  {
    return (ChargeSum[arm]);
  }
  /*
  float getArmHitTime (Bbc::ArmType arm) const
  {
    return (ArmHitTime[arm]);
  }
  float getArmHitTimeError (Bbc::ArmType arm) const
  {
    return (ArmHitTimeError[arm]);
  }
  */

  float getZVertex (void) const
  {
    return (ZVertex);
  }
  float getZVertexError (void) const
  {
    return (ZVertexError);
  }
  float getTimeZero (void) const
  {
    return (TimeZero);
  }

  /*
  float getTimeZeroError (void) const
  {
    return (TimeZeroError);
  }

  void setEventNumber (int givenEventNumber)
  {
    EventNumber = givenEventNumber;
  }
  void setAdc (int givenAdc, int PmtIndx)
  {
    Adc[PmtIndx] = givenAdc;
  }
  void setTdc0 (int givenTdc0, int PmtIndx)
  {
    Tdc0[PmtIndx] = givenTdc0;
  }
  void setTdc1 (int givenTdc1, int PmtIndx)
  {
    Tdc1[PmtIndx] = givenTdc1;
  }

  int printAdcTdc (void) const;
  int printPmtChTime (void) const;

  float TimeLagOfTransitTime (int PmtIndx, float ZVertex2) const;
  */

  OnlBbcSig* GetSig(const int ipmt) { return &bbcsig[ipmt]; }


private:
  int Read_Charge_Calib( const char *calfname );
  int Read_TQ_T0_Offsets( const char *calfname );
  float gaincorr[128] {};         // gain corrections
  float tq_t0_offsets[128] {};    // t0 offsets in charge channels
  float bz_offset{0.};

  static const int NCH = 256;
  static const int NSAMPLES = 31;
  static const int NCHPERPKT = 128;

  int verbose;
  int EventNumber;
  int calib_done;
  Packet *p[2] {nullptr,nullptr} ;

  //raw data
  Float_t f_adc[NCH][NSAMPLES];     // raw waveform
  Float_t f_samp[NCH][NSAMPLES];    // raw waveform
  Float_t f_ampl[NCH];              // amplitude
  //Float_t f_t0[NCH];                // time

  std::vector<OnlBbcSig> bbcsig;
  
  static const int BBC_N_PMT = 128;

  // converted (corrected) data
  int iHit[BBC_N_PMT];
  int armHitPmt[BBC_N_PMT];

  // Adc value w/ pedestal subtracted
  //float TrueAdc[BBC_N_PMT];		

  float Charge[BBC_N_PMT];
  float HitTime0[BBC_N_PMT];
  float HitTime1[BBC_N_PMT];

  // End product data
  int nHitPmt[2];
  float ChargeSum[2];
  float ArmHitTime[2];
  float ArmHitTimeError[2];
  float ZVertex;
  float ZVertexError;
  float TimeZero;
  float TimeZeroError;

  Float_t f_pmtq[BBC_N_PMT]{};   // npe in each arm
  Float_t f_pmtt0[BBC_N_PMT]{};  // time in each arm
  Float_t f_pmtt1[BBC_N_PMT]{};  // time in each arm
  Short_t f_bbcn[2]{};     // num hits for each arm (north and south)
  Float_t f_bbcq[2]{};     // total charge (currently npe) in each arm
  Float_t f_bbct[2]{};     // time in arm
  Float_t f_bbcte[2]{};    // earliest hit time in arm
  Float_t f_bbcz = NAN;    // z-vertex
  Float_t f_bbct0 = NAN;   // start time
  TH1 *hevt_bbct[2]{};  // time in each bbc, per event
  TF1 *gaussian = nullptr;
  Float_t _tres = NAN;  // time resolution of one channel

  TH2 *h2_tmax[2] = {};  // [0 == time ch, 1 == chg ch], max sample in evt vs ch

  float TRIG_SAMP[16];  // [board]
  float bbc_tq_t0_offsets[BBC_N_PMT] = {};

  TCanvas *ac;
  /*
  int calcEndProduct (void);
  int calcArmProduct (Bbc::ArmType arm);
  int calcPmtProduct (int PmtIndx);

  // must exec after calcArmProduct
  int calcFlag ();	
  */

};

#endif /* __ONLBBCEVENT_H__ */
