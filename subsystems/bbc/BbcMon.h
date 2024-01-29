#ifndef BBC_BBCMON_H
#define BBC_BBCMON_H

#include <onlmon/OnlMon.h>

class Event;
//class OnlMonDB;
class TH1;
class TH2;
class TF1;
class TH2Poly;
class OnlBbcEvent;
class BbcGeom;

class BbcMon : public OnlMon
{
public:
  explicit BbcMon(const std::string &name = "BBCMON");
  virtual ~BbcMon();

  int process_event(Event *evt) override;
  int Init() override;                     // Done once at start of server
  int BeginRun(const int runno) override;  // Done at every start of run
  int EndRun(const int runno) override;    // Done at the end of every run
  int Reset() override;

protected:
  int DBVarInit();

  OnlBbcEvent *bevt { nullptr };

  int evtcnt = 0;
  //OnlMonDB *dbvars = nullptr;

  BbcGeom *bbcgeom = nullptr;   // contains positions of BBC PMTs

  TH2 *bbc_adc = nullptr;
  TH2 *bbc_tdc = nullptr;
  TH2 *bbc_tdc_overflow = nullptr;
  TH1 *bbc_tdc_overflow_each[128] = {};  // should be [nPMT_BBC], need to fix
  // TH1 *bbc_nhit[nTRIGGER];

  TH2 *bbc_tdc_armhittime = nullptr;
  TH1 *bbc_nevent_counter = nullptr;

  TH1 *bbc_zvertex = nullptr;
  TH1 *bbc_zvertex_bbll1 = nullptr;
  TH1 *bbc_zvertex_short {nullptr};   // Used for short time scales
  TF1 *f_zvtx {nullptr};              // Used for fitting central vertex peak
  //TH1 *bbc_zvertex_bbll1_novtx = nullptr;
  //TH1 *bbc_zvertex_bbll1_narrowvtx = nullptr;  // Run11 pp

  TH2 *bbc_tzero_zvtx = nullptr;
  TH1 *bbc_avr_hittime = nullptr;
  TH1 *bbc_south_hittime = nullptr;
  TH1 *bbc_north_hittime = nullptr;
  TH1 *bbc_south_chargesum = nullptr;
  TH1 *bbc_north_chargesum = nullptr;
  TH1 *bbc_prescale_hist = nullptr;

  //TH2 *bbc_tmax[2] = {};  // [0 == time ch, 1 == chg ch], max sample in evt vs ch
  TH2 *bbc_time_wave = nullptr;      // raw waveforms, time ch
  TH2 *bbc_charge_wave = nullptr;    // raw waveforms, charge ch
  //TH2 *bbc_trange = nullptr;       // subtracted tdc at maxsamp vs ch

  TH2Poly *bbc_south_hitmap = nullptr;    // hitmap
  TH2Poly *bbc_north_hitmap = nullptr;    // hitmap

};

#endif /* BBC_BBCMON_H */
