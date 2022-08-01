#ifndef MBD_MBDMON_H
#define MBD_MBDMON_H

#include <onlmon/OnlMon.h>

class Event;
class OnlMonDB;
class TH1;
class TH2;
class TRandom3;

class MbdMon : public OnlMon
{
 public:
  MbdMon(const std::string &name = "MBDMON");
  virtual ~MbdMon();

  int process_event(Event *evt);
  int Init();                     // Done once at start of server
  int BeginRun(const int runno);  // Done at every start of run
  int Reset();

 protected:
  int DBVarInit();
  int evtcnt = 0;
  int idummy = 0;
  OnlMonDB *dbvars = nullptr;
  
  TH2 *bbc_adc;
  TH2 *bbc_tdc;
  TH2 *bbc_tdc_overflow;
  TH1 *bbc_tdc_overflow_each[128];    // should be [nPMT_BBC], need to fix
  //TH1 *bbc_nhit[nTRIGGER];

  TH2 *bbc_tdc_armhittime;
  TH1 *bbc_nevent_counter;

  TH1 *bbc_zvertex;
  TH1 *bbc_zvertex_bbll1;
  //TH1 *bbc_zvertex_zdc;
  //TH1 *bbc_zvertex_zdc_scale3;
  TH1 *bbc_zvertex_bbll1_novtx;
  TH1 *bbc_zvertex_bbll1_narrowvtx;//Run11 pp
  //TH1 *bbc_zvertex_bbll1_zdc;

  TH2 *bbc_tzero_zvtx;
  TH1 *bbc_avr_hittime;
  TH1 *bbc_south_hittime;
  TH1 *bbc_north_hittime;
  TH1 *bbc_south_chargesum;
  TH1 *bbc_north_chargesum;
  TH1 *bbc_prescale_hist; 

  const double C = 29.9792458; // speed of light [cm/ns]

  TRandom3 *trand3;
};

#endif /* MBD_MBDMON_H */
