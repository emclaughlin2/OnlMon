#ifndef HCAL_HCALMON_H
#define HCAL_HCALMON_H

#include <onlmon/OnlMon.h>

class Event;
class OnlMonDB;
class TH1;
class TH2;
class Packet;

class HcalMon : public OnlMon
{
 public:
  HcalMon(const std::string &name = "HCALMON");
  virtual ~HcalMon();

  int process_event(Event *evt);
  int Init();
  int BeginRun(const int runno);
  int Reset();
  double getSignal(Packet *p, const int channel);

 protected:
  int DBVarInit();
  int evtcnt = 0;
  int idummy = 0;
  OnlMonDB *dbvars = nullptr;
  TH2 *h2_hcal_hits = nullptr;
  TH2 *h2_hcal_mean = nullptr;
  TH2 *h2_hcal_history = nullptr;
};

#endif /* HCAL_HCALMON_H */
