#ifndef EXAMPLE_MYMON_H
#define EXAMPLE_MYMON_H

#include <onlmon/OnlMon.h>

class Event;
class OnlMonDB;
class TH1;
class TH2;

class MyMon : public OnlMon
{
 public:
  MyMon(const std::string &name = "MYMON");
  virtual ~MyMon();

  int process_event(Event *evt);
  int Init();
  int BeginRun(const int runno);
  int Reset();

 protected:
  int DBVarInit();
  int evtcnt = 0;
  int idummy = 0;
  OnlMonDB *dbvars = nullptr;
  TH1 *myhist1 = nullptr;
  TH2 *myhist2 = nullptr;
};

#endif /* EXAMPLE_MYMON_H */
