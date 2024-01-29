#ifndef ONLMONCLIENT_ONLMONCLIENT_H
#define ONLMONCLIENT_ONLMONCLIENT_H

#include <onlmon/OnlMonBase.h>
#include <onlmon/OnlMonDefs.h>

#include <ctime>
#include <list>
#include <map>
#include <set>
#include <string>
#include <vector>

class ClientHistoList;
class OnlMonDraw;
class OnlMonHtml;
class TCanvas;
class TH1;
class TStyle;

class OnlMonClient : public OnlMonBase
{
 public:
  static OnlMonClient *instance();
  ~OnlMonClient() override;
  using OnlMonBase::Verbosity;
  void Verbosity(const int i) override;

  int UpdateServerHistoMap(const std::string &hname, const std::string &subsys, const std::string &hostname);
  void PutHistoInMap(const std::string &hname, const std::string &subsys, const std::string &hostname, const int port);
  void updateHistoMap(const std::string &subsys, const std::string &hname, TH1 *h1d);
  int requestMonitorList(const std::string &hostname, const int moniport);
  TH1 *getHisto(const std::string &monitor, const std::string &hname);
  OnlMonDraw *getDrawer(const std::string &name);
  int requestHisto(const std::string &what = "ALL", const std::string &hostname = "localhost", const int moniport = OnlMonDefs::MONIPORT);
  int requestHistoList(const std::string &subsys, const std::string &hostname, const int moniport, std::list<std::string> &histolist);
  int requestHistoByName(const std::string &subsystem, const std::string &what = "ALL");
  int requestHistoBySubSystem(const std::string &subsystem, int getall = 0);
  void registerHisto(const std::string &hname, const std::string &subsys);
  void Print(const char *what = "ALL");

  void AddServerHost(const std::string &hostname);
  void registerDrawer(OnlMonDraw *Drawer);
  int ReadHistogramsFromFile(const std::string &filename);
  int Draw(const char *who = "ALL", const char *what = "ALL");
  int MakePS(const char *who = "ALL", const char *what = "ALL");
  int MakeHtml(const char *who = "ALL", const char *what = "ALL");
  int SavePlot(const std::string &who = "ALL", const std::string &what = "ALL");

  std::string htmlRegisterPage(const OnlMonDraw &drawer,
                               const std::string &path,
                               const std::string &basefilename,
                               const std::string &ext);

  void htmlAddMenu(const OnlMonDraw &drawer, const std::string &path,
                   const std::string &relfilename);

  void htmlNamer(const OnlMonDraw &drawer, const std::string &basefilename,
                 const std::string &ext, std::string &fullfilename,
                 std::string &filename);

  int LocateHistogram(const std::string &hname, const std::string &subsys);
  int RunNumber();
  time_t EventTime(const std::string &which);
  time_t EventTime(const std::string &servername, const std::string &which);
  int SendCommand(const char *hostname, const int port, const char *cmd);

  void SetDisplaySizeX(int xsize) { display_sizex = xsize; }
  void SetDisplaySizeY(int ysize) { display_sizey = ysize; }
  int GetDisplaySizeX() { return display_sizex; }
  int GetDisplaySizeY() { return display_sizey; }
  int CanvasToPng(TCanvas *canvas, std::string const &filename);
  int HistoToPng(TH1 *histo, std::string const &pngfilename, const char *drawopt = "", const int statopt = 11);

  int SaveLogFile(const OnlMonDraw &drawer);
  int SetStyleToDefault();
  int isCosmicRun();
  int isStandalone();
  std::string RunType();
  void CacheRunDB(const int runno);
  void FindAllMonitors();
  int FindMonitor(const std::string &name);
  int IsMonitorRunning(const std::string &name);
  std::string ExtractSubsystem(const std::string &filename);

 private:
  OnlMonClient(const std::string &name = "ONLMONCLIENT");
  int DoSomething(const std::string &who, const std::string &what, const std::string &opt);
  void InitAll();

  static OnlMonClient *__instance;
  OnlMonHtml *fHtml = nullptr;
  TH1 *clientrunning = nullptr;
  TStyle *defaultStyle = nullptr;

  int display_sizex = 0;
  int display_sizey = 0;
  int cosmicrun = 0;
  int standalone = 0;
  int cachedrun = 0;

  std::string runtype = "UNKNOWN";
  std::set<std::string> m_MonitorFetchedSet;
  std::map<std::string, std::map<const std::string, ClientHistoList *>> SubsysHisto;
  std::map<std::string, std::pair<std::string, unsigned int>> MonitorHostPorts;
  std::map<const std::string, ClientHistoList *> Histo;
  std::map<const std::string, OnlMonDraw *> DrawerList;
  std::vector<std::string> MonitorHosts;
};

#endif /* ONLMONCLIENT_ONLMONCLIENT_H */
