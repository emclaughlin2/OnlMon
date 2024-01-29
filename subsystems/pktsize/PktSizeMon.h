#ifndef __PKTSIZEMON_H__
#define __PKTSIZEMON_H__

#include <onlmon/OnlMon.h>

#include <map>
#include <string>

class Event;
class Packet;
class TH1;
class PktSizeDBodbc;

class PktSizeMon : public OnlMon
{
 public:
  PktSizeMon(const std::string &name = "PKTSIZEMON");
  virtual ~PktSizeMon();

  int Init();
  int process_event(Event *e);
  int Reset();
  int EndRun(const int runno);
  void Print(const std::string &what = "ALL");
  int UpdateDB(const int runno = 0);
  //  virtual void Verbosity(const int i);

 protected:
  int putmapinhisto();
  int nevnts;
  TH1 *sizehist = nullptr;
  Packet *plist[10000]{};
  PktSizeDBodbc *db;
  std::map<unsigned int, unsigned int> packetsize;
  std::map<std::string, std::pair<unsigned int, unsigned int> > granulepacketlimits;
  std::map<unsigned int, std::string> dcmgroups;
  std::map<std::string, unsigned int> dcmgroupsize;
  std::map<unsigned int, std::string> fibergroups;
  std::map<std::string, unsigned int> fibergroupsize;
};

#endif /* __PKTSIZEMON_H__ */
