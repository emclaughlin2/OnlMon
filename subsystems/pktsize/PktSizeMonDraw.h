#ifndef __PKTSIZEMONDRAW_H__
#define __PKTSIZEMONDRAW_H__

#include <onlmon/OnlMonDraw.h>
#include <map>
#include <set>

class RunDBodbc;
class PktSizeDBodbc;
class TCanvas;
class TGraph;
class TH1;
class TH2;
class TMarker;
class TPad;
class TText;

class PktSizeMonDraw : public OnlMonDraw
{
 public:
  PktSizeMonDraw(const std::string &name = "PKTSIZEMON");
  virtual ~PktSizeMonDraw();

  int Init();
  int Draw(const std::string &what = "ALL");
  int MakeHtml(const std::string &what = "ALL");
  int SavePlot(const std::string &what = "ALL", const std::string &type = "png");
  int loadpreviousruns(const int nruns = 5);
  int AddKnownBig(const unsigned int packetid, const float maxsize);
  void Print(const std::string &what = "ALL") const;
  int CleanOldRuns(const unsigned int maxrun);
  int RemoveRun(const int runno);

 protected:
  int MakeCanvas(const char *name);
  int DrawFirst(const std::string &what = "ALL");
  int DrawOldFirst(const std::string &what = "ALL");
  int DrawHistory(const std::string &what = "ALL");
  int FillPacketMap(const TH1 *histo);
  int FillRunPacketList(std::map<unsigned int, float> &pkts, const int runnumber);
  int ExtractActivePackets(const std::map<unsigned int, float> &packetsize);
  int MakeNoisyCandidates();
  int IsKnownBig(const unsigned int packetid, const float size);
  int PlotNoisy(TGraph &gr, TMarker &tr, TText &tx, const unsigned int ipkt, const unsigned int icnt);
  int ReplaceRunFromDB(const int runno);

  TCanvas *TC[2]{};
  TPad *transparent[2]{};
  TPad *Pad[2]{};
  TH2* Frames[12]{};
  PktSizeDBodbc *db = nullptr;
  RunDBodbc *rd = nullptr;
  TMarker *tm = nullptr;
  int lastrun;
  std::map<int, std::map<unsigned int, float> > packetmap;
  std::map<std::string, std::pair<unsigned int, unsigned int> > granulepacketlimits;
  std::map<unsigned int, float> knownbig;
  std::set<unsigned int> activepackets;
  std::set<unsigned int> noisypackets;
  std::set<int> runlist;
};

#endif /*__PKTSIZEMONDRAW_H__ */
