#ifndef ONLMONSERVER_ONLMONSERVER_H
#define ONLMONSERVER_ONLMONSERVER_H

//#define USE_MUTEX

#include "OnlMonBase.h"
#include "OnlMonDefs.h"

#include <pthread.h>
#include <ctime>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

class Event;
class MessageSystem;
class OnlMon;
class OnlMonStatusDB;
class TH1;

class OnlMonServer : public OnlMonBase
{
 public:
  static OnlMonServer *instance();
  ~OnlMonServer() override;

  // delete copy ctor and assignment operator (cppcheck)
  explicit OnlMonServer(const OnlMonServer &) = delete;
  OnlMonServer &operator=(const OnlMonServer &) = delete;

  void registerHisto(const std::string &monitorname, const std::string &hname, TH1 *h1d, const int replace = 0);
  void registerHisto(const OnlMon *monitor, TH1 *h1d);

  void registerCommonHisto(TH1 *h1d);
  TH1 *getHisto(const std::string &subsys, const std::string &hname) const;
  TH1 *getCommonHisto(const std::string &hname) const;
  TH1 *getHisto(const unsigned int ihisto) const;
  const std::string getHistoName(const unsigned int ihisto) const;
  unsigned int nHistos() const { return CommonHistoMap.size(); }
  int RunNumber() const { return runnumber; }
  void RunNumber(const int irun);
  int EventNumber() const { return eventnumber; }
  void EventNumber(const int iev) { eventnumber = iev; }
  int PortNumber() const { return portnumber; }
  void PortNumber(const int i) { portnumber = i; }
  void Print(const std::string &what = "ALL", std::ostream& os = std::cout) const;
  void PrintFile(const std::string &fname) const;

  void InitAll();

  void registerMonitor(OnlMon *Monitor);
  OnlMon *getMonitor(const std::string &name);
  void dumpHistos(const std::string &filename);
  int process_event(Event *);
  int Reset();
  int BeginRun(const int runno);
  int EndRun(const int runno);
  int WriteHistoFile();

  uint64_t CurrentTicks() const { return currentticks; }
  void CurrentTicks(const uint64_t ival) { currentticks = ival; }
  uint64_t BorTicks() const { return borticks; }
  void BorTicks(const uint64_t ival) { borticks = ival; }

  int BadEvents() const { return badevents; }
  void AddBadEvent() { badevents++; }
  void BadEvents(const int ibad) { badevents = ibad; }

  int WriteLogFile(const std::string &name, const std::string &msg) const;

  int IsPacketActive(const unsigned int ipkt);
  // set status if something went wrong

  int SetSubsystemStatus(OnlMon *Monitor, const int status);
  int SetSubsystemRunStatus(OnlMon *Monitor, const int status);
  int LookAtMe(OnlMon *Monitor, const int level, const std::string &message);
  std::string GetRunType() const { return RunType; }

  int send_message(const OnlMon *Monitor, const int msgsource, const int severity, const std::string &err_message, const int msgtype) const;

#ifdef USE_MUTEX
  void GetMutex(pthread_mutex_t &lock) { lock = mutex; }
#endif
  void SetThreadId(const pthread_t &id) { serverthreadid = id; }

  //int LoadActivePackets();
  //  int parse_granuleDef(std::set<std::string> &pcffilelist);
  // void parse_pcffile(const std::string &lfn);
  void clearactivepackets() { activepackets.clear(); }
  //int LoadLL1Packets();
  int isStandaloneRun() const { return standalone; }
  int isCosmicRun() const { return cosmicrun; }

  int run_empty(const int nevents);
  std::map<std::string, std::map<std::string, TH1 *>>::const_iterator monibegin() { return MonitorHistoSet.begin(); }
  std::map<std::string, std::map<std::string, TH1 *>>::const_iterator moniend() { return MonitorHistoSet.end(); }
  std::vector<OnlMon *>::const_iterator monitor_vec_begin() { return MonitorList.begin(); }
  std::vector<OnlMon *>::const_iterator monitor_vec_end() { return MonitorList.end(); }

 private:
  OnlMonServer(const std::string &name = "OnlMonServer");
  int send_message(const int severity, const std::string &err_message, const int msgtype) const;
  int CacheRunDB(const int runno);
  void registerHisto(const std::string &hname, TH1 *h1d, const int replace = 0);

  static OnlMonServer *__instance;
  int runnumber = -1;
  int eventnumber = 0;
  int portnumber = OnlMonDefs::MONIPORT;
  int badevents = 0;
  uint64_t currentticks = 0;
  uint64_t borticks = 0;
  int activepacketsinit = 0;
  unsigned int scaledtrigmask = 0xFFFFFFFF;
  int scaledtrigmask_used = 0;
  int standalone = 0;
  int cosmicrun = 0;
  std::string TriggerConfig = "UNKNOWN";
  std::string RunType = "UNKNOWN";

  TH1 *serverrunning = nullptr;
  OnlMonStatusDB *statusDB = nullptr;
  OnlMonStatusDB *RunStatusDB = nullptr;
  std::map<const std::string, TH1 *> CommonHistoMap;
  std::vector<OnlMon *> MonitorList;
  std::set<unsigned int> activepackets;
  std::map<std::string, MessageSystem *> MsgSystem;
  std::map<std::string, std::map<std::string, TH1 *>> MonitorHistoSet;
  pthread_mutex_t mutex;
  pthread_t serverthreadid = 0;
};

#endif /* __ONLMONSERVER_H */
