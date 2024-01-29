#ifndef DAQ_DAQMON_H
#define DAQ_DAQMON_H

#include <onlmon/OnlMon.h>

class Event;
class TH1;
class TH2;

class DaqMon : public OnlMon
{
 public:
  DaqMon(const std::string &name = "DAQMON");
  virtual ~DaqMon();

  int process_event(Event *evt);
  int Init();
  int BeginRun(const int runno);
  int Reset();

 protected:
  int evtcnt = 0;
  int idummy = 0;
  TH1 *daqhist1 = nullptr;
  TH2 *daqhist2 = nullptr;
};

#endif /* DAQ_DAQMON_H */
