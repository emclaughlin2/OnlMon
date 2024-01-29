#ifndef ZDC_ZDCMON_H
#define ZDC_ZDCMON_H

#include <onlmon/OnlMon.h>

#include <onlmon/OnlMon.h>
#include <cmath>
#include <vector>

class CaloWaveformFitting;
class TowerInfoContainer;
class Event;
class TH1;
class TH2;
class Packet;

class ZdcMon : public OnlMon
{
 public:
  ZdcMon(const std::string &name = "ZDCMON");
  virtual ~ZdcMon();

  int process_event(Event *evt);
  int Init();
  int BeginRun(const int runno);
  int Reset();

 protected:
  std::vector<float> anaWaveformFast(Packet *p, const int channel);
  CaloWaveformFitting *WaveformProcessingFast = nullptr;

  int evtcnt = 0;
  int idummy = 0;
  TH1 *zdc_adc_north = nullptr;
  TH1 *zdc_adc_south = nullptr;

};

#endif /* ZDC_ZDCMON_H */
