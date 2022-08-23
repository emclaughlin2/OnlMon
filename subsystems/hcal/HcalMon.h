#ifndef HCAL_HCALMON_H
#define HCAL_HCALMON_H

#include <onlmon/OnlMon.h>
#include <vector>

class Event;
class OnlMonDB;
class TH1;
class TH2;
class Packet;
class runningMean;

class HcalMon : public OnlMon
{
 public:
  HcalMon(const std::string &name = "HCALMON");
  virtual ~HcalMon();

  int process_event(Event *evt);
  int Init();
  int BeginRun(const int runno);
  int Reset();

 protected:
  double getSignal(Packet *p, const int channel);
  int DBVarInit();
  int evtcnt = 0;
  int idummy = 0;
  bool debug = false;
  OnlMonDB *dbvars = nullptr;
  TH2 *h2_hcal_hits = nullptr;
  TH2 *h2_hcal_mean = nullptr;
  TH2 *h2_hcal_history = nullptr;

  std::string runtypestr = "Unknown";
  std::vector<runningMean*> rm_vector;
};

#endif /* HCAL_HCALMON_H */
