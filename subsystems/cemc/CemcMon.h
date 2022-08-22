#ifndef CEMC_CEMCMON_H
#define CEMC_CEMCMON_H


#include <onlmon/OnlMon.h>

#include <vector>

class Event;
class OnlMonDB;
class TH1;
class TH2;
class Packet;
class runningMean;

class CemcMon : public OnlMon
{
 public:
  CemcMon(const std::string &name = "CEMCMON", const std::string &id ="");
  virtual ~CemcMon();

  int process_event(Event *evt);
  int Init();
  int BeginRun(const int runno);
  int Reset();

 protected:
  double getSignal(Packet *p, const int channel);

  int DBVarInit();
  int evtcnt = 0;
  int idummy = 0;
  OnlMonDB *dbvars = nullptr;
  TH2 *cemc_occupancy = nullptr;
  TH2 *cemc_runningmean = nullptr;
  TH1 *cemc_signal = nullptr;

  std::string runtypestr = "Unknown";
  std::string id_string;

  std::vector<runningMean*> rm_vector; 

};

#endif /* CEMC_CEMCMON_H */
