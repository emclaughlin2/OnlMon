#include "OnlMonClient.h"
#include "ClientHistoList.h"
#include "OnlMonDraw.h"
#include "OnlMonHtml.h"

#include <onlmon/HistoBinDefs.h>
#include <onlmon/OnlMonBase.h>  // for OnlMonBase
#include <onlmon/OnlMonDefs.h>

#include <MessageTypes.h>  // for kMESS_STRING, kMESS_OBJECT
#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGClient.h>  // for gClient, TGClient
#include <TGFrame.h>
#include <TH1.h>
#include <TImage.h>
#include <TIterator.h>
#include <TList.h>  // for TList
#include <TMessage.h>
#include <TROOT.h>
#include <TSeqCollection.h>
#include <TSocket.h>
#include <TStyle.h>
#include <TSystem.h>

#include <odbc++/connection.h>
#include <odbc++/drivermanager.h>
#include <odbc++/resultset.h>
#include <odbc++/statement.h>  // for Statement
#include <odbc++/types.h>      // for SQLException

#include <sys/stat.h>
#include <sys/utsname.h>
#include <unistd.h>
#include <uuid/uuid.h>
#include <algorithm>
#include <cstdio>   // for printf, remove
#include <cstdlib>  // for getenv, exit
#include <cstring>  // for strcmp
#include <filesystem>
#include <fstream>
#include <iostream>
#include <list>
#include <sstream>
#include <utility>  // for pair

OnlMonClient *OnlMonClient::__instance = nullptr;

int pinit()
{
  return 0;
}

OnlMonClient *OnlMonClient::instance()
{
  if (__instance)
  {
    return __instance;
  }
  __instance = new OnlMonClient("ONLMONCLIENT");
  return __instance;
}

OnlMonClient::OnlMonClient(const std::string &name)
  : OnlMonBase(name)
{
  defaultStyle = new TStyle();
  SetStyleToDefault();
  InitAll();
}

void OnlMonClient::InitAll()
{
  if (gROOT->FindObject("ServerRunning"))
  {
    std::cout << "Don't run Server and Client in same session, exiting" << std::endl;
    exit(1);
  }
  if (!gClient)  // cannot draw on display, warn and bail out
  {
    const char *env_display = getenv("DISPLAY");
    if (env_display)
    {
      std::string displaystring = env_display;
      std::cout << "Cannot open Display, Display Env Var is set to "
                << displaystring << std::endl;
      if (displaystring.find("unix") != std::string::npos)
      {
        utsname ThisNode;
        uname(&ThisNode);
        std::cout << "presumably the virtual framebuffer is not running on " << ThisNode.nodename
                  << ", check if process /usr/X11R6/bin/Xvfb is alive" << std::endl;
      }
      else if (displaystring.find("localhost") != std::string::npos)
      {
        std::cout << "Check your ssh forwarding" << std::endl;
        std::cout << "your $HOME/.ssh/config has to contain the line" << std::endl;
        std::cout << "ForwardX11 yes" << std::endl;
      }
    }
    else
    {
      std::cout << "Display not set, cannot continue" << std::endl;
    }
    exit(1);
  }
  // we have a working display
  clientrunning = new TH1F("ClientRunning", "ClientRunning", 1, 0, 1);
  if (getenv("ONLMON_HTMLDIR") == nullptr)
  {
    std::cout << "ONLMON_HTMLDIR not set, exiting" << std::endl;
    exit(1);
  }
  fHtml = new OnlMonHtml(getenv("ONLMON_HTMLDIR"));

  TGFrame *rootWin = (TGFrame *) gClient->GetRoot();
  display_sizex = rootWin->GetDefaultWidth();
  display_sizey = rootWin->GetDefaultHeight();
  for (int i = 0; i < kMAXSIGNALS; i++)
  {
    gSystem->IgnoreSignal((ESignals) i);
  }
  return;
}

OnlMonClient::~OnlMonClient()
{
  while (DrawerList.begin() != DrawerList.end())
  {
    if (verbosity > 0)
    {
      std::cout << "deleting " << DrawerList.begin()->first << std::endl;
    }
    delete DrawerList.begin()->second;
    DrawerList.erase(DrawerList.begin());
  }
  while (Histo.begin() != Histo.end())
  {
    delete Histo.begin()->second;
    Histo.erase(Histo.begin());
  }
  delete clientrunning;
  delete fHtml;
  delete defaultStyle;
  // this deletes all open canvases
  TSeqCollection *allCanvases = gROOT->GetListOfCanvases();
  TCanvas *canvas = nullptr;
  while ((canvas = static_cast<TCanvas *>(allCanvases->First())))
  {
    if (verbosity > 0)
    {
      std::cout << "Deleting Canvas " << canvas->GetName() << std::endl;
    }
    delete canvas;
  }
  __instance = nullptr;
  return;
}

void OnlMonClient::registerHisto(const std::string &hname, const std::string &subsys)
{
  auto subsysiter = SubsysHisto.find(subsys);
  if (subsysiter == SubsysHisto.end())
  {
    std::map<const std::string, ClientHistoList *> entry;
    // c++11 map.insert returns pair<iterator,bool>
    subsysiter = SubsysHisto.insert(std::make_pair(subsys, entry)).first;
    if (Verbosity() > 2)
    {
      std::cout << "inserting " << subsys << " into SubsysHisto, readback: " << subsysiter->first << std::endl;
    }
  }
  auto hiter = subsysiter->second.find(hname);
  if (hiter == subsysiter->second.end())
  {
    ClientHistoList *newhisto = new ClientHistoList();
    newhisto->SubSystem(subsys);
    subsysiter->second.insert(std::make_pair(hname, newhisto));
    if (Verbosity() > 0)
    {
      std::cout << "new histogram " << hname << " of subsystem " << subsys << std::endl;
    }
  }
  else
  {
    hiter->second->SubSystem(subsys);
  }
  // register FrameWorkVars histo if we don't have it already
  std::string frameworkhistoname = "FrameWorkVars";
  hiter = subsysiter->second.find(frameworkhistoname);
  if (hiter == subsysiter->second.end())
  {
    ClientHistoList *newhisto = new ClientHistoList();
    newhisto->SubSystem(subsys);
    subsysiter->second.insert(std::make_pair(frameworkhistoname, newhisto));
    if (Verbosity() > 0)
    {
      std::cout << "new histogram " << frameworkhistoname << " of subsystem " << subsys << std::endl;
    }
  }
  else
  {
    hiter->second->SubSystem(subsys);
  }

  return;
}

int OnlMonClient::requestHistoBySubSystem(const std::string &subsys, int getall)
{
  std::string mysubsys = subsys.substr(0,subsys.find('_'));
  for (auto frwrkiter : m_MonitorFetchedSet)
  {
    if (frwrkiter.find(mysubsys) == std::string::npos)
      {
	m_MonitorFetchedSet.clear();
	break;
      }
  }
  int iret = 0;
  std::map<const std::string, ClientHistoList *>::const_iterator histoiter;
  std::map<const std::string, ClientHistoList *>::const_iterator histonewiter;
  if (!IsMonitorRunning(subsys))  // test if saved monitor server is still running
  {
    // monitor is not running - remove its entry if it exists
    auto moniiter = MonitorHostPorts.find(subsys);
    if (moniiter != MonitorHostPorts.end())
    {
      MonitorHostPorts.erase(moniiter);
    }
    // Find monitor, if not found reset ClientHistoList entries for host and port
    if (FindMonitor(subsys) == 0)
    {
      auto subsyshistos = SubsysHisto.find(subsys);
      if (subsyshistos != SubsysHisto.end())
      {
        for (auto &hiter : subsyshistos->second)
        {
          hiter.second->ServerHost("UNKNOWN");
          hiter.second->ServerPort(0);
          if (hiter.second->Histo() != nullptr)
          {
            hiter.second->Histo()->Delete();
          }
          hiter.second->Histo(nullptr);
        }
      }
      return iret;
    }
  }

  if (getall == 0)
  {
    std::map<std::string, std::map<const std::string, ClientHistoList *>>::const_iterator histos = SubsysHisto.find(subsys);
    for (auto &hiter : histos->second)
    {
      if (requestHistoByName(subsys, hiter.first))
      {
        if (Verbosity() > 2)
        {
          std::cout << "Request for " << hiter.first << " on " << subsys << " failed " << std::endl;
        }
        if (hiter.second->Histo())
        {
          hiter.second->Histo()->Delete();
          hiter.second->Histo(nullptr);
          hiter.second->ServerHost("UNKNOWN");
          hiter.second->ServerPort(0);
        }
        iret = -1;
      }
    }
  }
  else if (getall == 1)
  {
    std::map<std::string, std::list<std::string>> transferlist;
    std::ostringstream host_port;

    int failed_to_locate = 0;
    auto subs = SubsysHisto.find(subsys);
    {
      for (auto &histos : subs->second)
      {
        if (Verbosity() > 2)
        {
          std::cout << "checking for subsystem " << subs->first << ", histo " << histos.first << std::endl;
        }
        if (histos.second->SubSystem() == subsys)
        {
          int unknown_histo = 0;
          std::string hname = histos.first;
          if (histos.second->ServerHost() == "UNKNOWN")
          {
            if (!failed_to_locate)
            {
              if (LocateHistogram(hname, subsys) < 1)
              {
                if (Verbosity() > 2)
                {
                  std::cout << "Subsystem " << subsys << " Histogram " << hname << " cannot be located" << std::endl;
                }
                failed_to_locate = 1;
                unknown_histo = 1;
                iret = -1;
              }
            }
            else
            {
              unknown_histo = 1;
            }
          }
          // reset histogram in case they don't exist on the server anymore
          if (histos.second->Histo())
          {
            histos.second->Histo()->Reset();
          }
          if (!unknown_histo)
          {
            std::string fullhname = subsys + std::string(" ") + hname;
            (transferlist[subsys]).push_back(fullhname);
          }
        }
      }
    }
    std::list<std::string>::const_iterator liter;
    for (auto listiter = transferlist.begin(); listiter != transferlist.end(); ++listiter)
    {
      std::list<std::string> hlist = listiter->second;
      auto hostportiter = MonitorHostPorts.find(listiter->first);
      if (hostportiter == MonitorHostPorts.end())
      {
        std::cout << __PRETTY_FUNCTION__ << "Cannot find MonitorHostPorts entry for " << listiter->first << std::endl;
        std::cout << "existing hosts: " << std::endl;
        for (auto &hport : MonitorHostPorts)
        {
          std::cout << "subsystem " << hport.first << " on host " << hport.second.first
                    << " listening on port " << hport.second.second << std::endl;
        }
        continue;
      }
      if (requestHistoList(listiter->first, hostportiter->second.first, hostportiter->second.second, hlist) != 0)
      {
        for (liter = hlist.begin(); liter != hlist.end(); ++liter)
        {
          if (requestHistoByName(*liter))
          {
            std::cout << "Request for " << *liter << " failed " << std::endl;
            histoiter = Histo.find(*liter);
            if (histoiter->second && histoiter->second->Histo())
            {
              histoiter->second->Histo()->Delete();
              histoiter->second->Histo(nullptr);
              histoiter->second->ServerHost("UNKNOWN");
              histoiter->second->ServerPort(0);
            }
          }
        }
      }
    }
  }
  else if (getall == 2)
  {
    const char *hname = "none";
    for (histoiter = Histo.begin(); histoiter != Histo.end(); ++histoiter)
    {
      if (histoiter->second->SubSystem() == subsys)
      {
        hname = histoiter->first.c_str();
        break;
      }
    }
    // histoiter points here to first histogram of subsystem
    if (histoiter != Histo.end())
    {
      std::string hostname = histoiter->second->ServerHost();
      int moniport = histoiter->second->ServerPort();
      {
        if (hostname == "UNKNOWN")
        {
          if (LocateHistogram(hname, "") < 1)
          {
            std::cout << "Histogram " << hname << " cannot be located" << std::endl;
            return -1;
          }
          histonewiter = Histo.find(hname);
          hostname = histonewiter->second->ServerHost();
          moniport = histonewiter->second->ServerPort();
        }

        // Reset the histograms in case they don't exhist on the host.
        for (histonewiter = Histo.begin(); histonewiter != Histo.end(); ++histonewiter)
        {
          if (histonewiter->second->SubSystem() == subsys)
          {
            if (histonewiter->second->Histo())
            {
              histonewiter->second->Histo()->Reset();
            }
          }
        }

        if (requestHisto("ALL", hostname, moniport))
        {
          std::cout << "Request for all histograms from "
                    << hostname << " failed, trying single " << std::endl;
          iret = requestHistoBySubSystem(subsys, 0);
        }
      }
    }
    else
    {
      std::cout << "No Histogram of subsystem "
                << subsys << " registered" << std::endl;
    }
  }
  m_MonitorFetchedSet.insert(subsys);
  return iret;
}

void OnlMonClient::registerDrawer(OnlMonDraw *Drawer)
{
  std::map<const std::string, OnlMonDraw *>::iterator iter = DrawerList.find(Drawer->Name());
  if (iter != DrawerList.end())
  {
    std::cout << "Drawer " << Drawer->Name() << " already registered, I won't overwrite it" << std::endl;
    std::cout << "Use a different name and try again" << std::endl;
  }
  else
  {
    DrawerList[Drawer->Name()] = Drawer;
    Drawer->Init();
    SetStyleToDefault();
  }
  return;
}

int OnlMonClient::Draw(const char *who, const char *what)
{
  int iret = DoSomething(who, what, "DRAW");
  //  gSystem->ProcessEvents();
  return iret;
}

int OnlMonClient::SavePlot(const std::string &who, const std::string &what)
{
  int iret = DoSomething(who, what, "SAVEPLOT");
  TSeqCollection *allCanvases = gROOT->GetListOfCanvases();
  TCanvas *canvas = nullptr;
  while ((canvas = static_cast<TCanvas *>(allCanvases->First())))
  {
    if (verbosity > 0)
    {
      std::cout << "Deleting Canvas " << canvas->GetName() << std::endl;
    }
    delete canvas;
  }
  //  gSystem->ProcessEvents();
  return iret;
}

int OnlMonClient::MakePS(const char *who, const char *what)
{
  int iret = DoSomething(who, what, "PS");
  //  gSystem->ProcessEvents();
  return iret;
}

int OnlMonClient::MakeHtml(const char *who, const char *what)
{
  mode_t old_umask;
  int runno = RunNumber();
  if (runno <= 0)
  {
    std::cout << "Run Number too small: " << runno
              << " not creating html output" << std::endl;
    return 0;
  }
  char *onlmon_real_html = getenv("ONLMON_REAL_HTML");
  if (!onlmon_real_html)
  {
    old_umask = umask(S_IWOTH);
    std::cout << "Making html output group writable so others can run tests as well" << std::endl;
  }
  fHtml->runNumber(runno);  // do not forget this !
  int iret = DoSomething(who, what, "HTML");

  if (!onlmon_real_html)
  {
    umask(old_umask);
  }
  //  gSystem->ProcessEvents();
  return iret;
}

int OnlMonClient::DoSomething(const std::string &who, const std::string &what, const std::string &opt)
{
  std::map<const std::string, OnlMonDraw *>::iterator iter;
  if (who != "ALL")  // is drawer was specified (!All)
  {
    iter = DrawerList.find(who);
    if (iter != DrawerList.end())
    {
      if (opt == "DRAW")
      {
        iter->second->Draw(what);
      }
      else if (opt == "PS")
      {
        iter->second->MakePS(what);
      }
      else if (opt == "SAVEPLOT")
      {
        iter->second->SavePlot(what, "png");
      }
      else if (opt == "HTML")
      {
        if (verbosity > 0)
        {
          std::cout << __PRETTY_FUNCTION__ << " creating html output for "
                    << iter->second->Name() << std::endl;
        }
        if (iter->second->MakeHtml(what))
        {
          std::cout << "subsystem " << iter->second->Name()
                    << " not in root file, skipping" << std::endl;
        }
      }
      else
      {
        std::cout << "option " << opt << " not implemented" << std::endl;
        return 0;
      }
      SetStyleToDefault();
      return 0;
    }
    else
    {
      std::cout << "Drawer " << who << " not in list" << std::endl;
      Print("DRAWER");
      return -1;
    }
  }
  else
  {
    for (iter = DrawerList.begin(); iter != DrawerList.end(); ++iter)
    {
      if (opt == "DRAW")
      {
        iter->second->Draw(what);
      }
      else if (opt == "PS")
      {
        iter->second->MakePS(what);
      }
      else if (opt == "HTML")
      {
        if (verbosity > 0)
        {
          std::cout << __PRETTY_FUNCTION__ << " creating html output for "
                    << iter->second->Name() << std::endl;
        }
        gROOT->Reset();
        int iret = iter->second->MakeHtml(what);
        if (iret)
        {
          std::cout << "subsystem " << iter->second->Name()
                    << " not in root file, skipping" << std::endl;
          // delete all canvases (no more piling up of 50 canvases)
          // if run for a single subsystem this leaves the canvas intact
          // for debugging
          TSeqCollection *allCanvases = gROOT->GetListOfCanvases();
          TCanvas *canvas = nullptr;
          while ((canvas = static_cast<TCanvas *>(allCanvases->First())))
          {
            std::cout << "Deleting Canvas " << canvas->GetName() << std::endl;
            delete canvas;
          }
        }
      }
      SetStyleToDefault();
    }
  }
  return 0;
}

int OnlMonClient::requestHistoByName(const std::string &subsys, const std::string &what)
{
  std::string hostname = "UNKNOWN";
  int moniport = OnlMonDefs::MONIPORT;
  std::map<std::string, std::map<const std::string, ClientHistoList *>>::const_iterator histos = SubsysHisto.find(subsys);
  auto histoiter = histos->second.find(what);
  if (Verbosity() > 2)
  {
    std::cout << __PRETTY_FUNCTION__ << "checking for " << what << " on monitor " << subsys << std::endl;
  }
  if (histoiter != histos->second.end())
  {
    hostname = histoiter->second->ServerHost();
    moniport = histoiter->second->ServerPort();
    if (hostname == "UNKNOWN")
    {
      if (LocateHistogram(what, subsys) < 1)
      {
        std::cout << "Histogram " << what << " cannot be located" << std::endl;
        return -1;
      }
      // search again since LocateHistogram can change the map
      histoiter = histos->second.find(what);
      hostname = histoiter->second->ServerHost();
      moniport = histoiter->second->ServerPort();
      if (hostname == "UNKNOWN")
      {
        std::cout << __PRETTY_FUNCTION__ << "host UNKNOWN for whatever reason" << std::endl;
        return -3;
      }
    }
  }
  else
  {
    if (LocateHistogram(what, subsys) < 1)
    {
      std::cout << "Histogram " << what << " cannot be located" << std::endl;
      return -2;
    }
    histoiter = histos->second.find(what);
    if (histoiter != histos->second.end())
    {
      hostname = histoiter->second->ServerHost();
      moniport = histoiter->second->ServerPort();
      if (hostname == "UNKNOWN")
      {
        std::cout << __PRETTY_FUNCTION__ << "host UNKNOWN for whatever reason" << std::endl;
        return -3;
      }
    }
    else
    {
      std::cout << __PRETTY_FUNCTION__ << "Problem determining host" << std::endl;
    }
  }
  // Open connection to server
  TSocket sock(hostname.c_str(), moniport);
  TMessage *mess;
  std::string fullhistoname = subsys + std::string(" ") + what;
  if (Verbosity() > 2)
  {
    std::cout << __PRETTY_FUNCTION__ << " sending " << fullhistoname << " to " << hostname << " port " << moniport << std::endl;
  }
  sock.Send(fullhistoname.c_str());
  while (true)
  {
    if (verbosity > 1)
    {
      std::cout << __PRETTY_FUNCTION__ << "Waiting for Message from : " << hostname
                << " on port " << moniport << std::endl;
    }
    sock.Recv(mess);
    if (!mess)  // if server is not up mess is NULL
    {
      std::cout << __PRETTY_FUNCTION__ << "Server not running on " << hostname << std::endl;
      sock.Close();
      return 1;
    }
    if (mess->What() == kMESS_STRING)
    {
      char str[OnlMonDefs::MSGLEN];
      mess->ReadString(str, OnlMonDefs::MSGLEN);
      delete mess;
      if (verbosity > 1)
      {
        std::cout << __PRETTY_FUNCTION__ << "Message: " << str << std::endl;
      }
      if (!strcmp(str, "Finished"))
      {
        break;
      }
      else if (!strcmp(str, "UnknownHisto"))
      {
        break;
      }
      else
      {
        std::cout << __PRETTY_FUNCTION__ << "Unknown Text Message: " << str << std::endl;
        sock.Send("Ack");
      }
    }
    else if (mess->What() == kMESS_OBJECT)
    {
      TH1 *histo = static_cast<TH1 *>(mess->ReadObjectAny(mess->GetClass()));
      delete mess;
      TH1 *maphist = static_cast<TH1 *>(histo->Clone(histo->GetName()));
      if (verbosity > 1)
      {
        std::cout << __PRETTY_FUNCTION__ << "histoname: " << histo->GetName() << " at "
                  << histo << std::endl;
      }

      updateHistoMap(subsys, histo->GetName(), maphist);
      delete histo;
      sock.Send("Ack");
    }
  }
  sock.Send("Finished");  // tell server we are finished
  // Close the socket
  sock.Close();
  return 0;
}

int OnlMonClient::requestHisto(const std::string &what, const std::string &hostname, const int moniport)
{
  // Open connection to server
  TSocket sock(hostname.c_str(), moniport);
  TMessage *mess;
  sock.Send(what.c_str());
  while (true)
  {
    sock.Recv(mess);
    if (!mess)  // if server is not up mess is NULL
    {
      std::cout << __PRETTY_FUNCTION__ << "Server not running on " << hostname << std::endl;
      sock.Close();
      return 1;
    }
    if (mess->What() == kMESS_STRING)
    {
      char str[OnlMonDefs::MSGLEN];
      mess->ReadString(str, OnlMonDefs::MSGLEN);
      delete mess;
      if (verbosity > 1)
      {
        std::cout << __PRETTY_FUNCTION__ << "Message: " << str << std::endl;
      }

      if (!strcmp(str, "Finished"))
      {
        break;
      }
      else if (!strcmp(str, "UnknownHisto"))
      {
        break;
      }
      else
      {
        std::cout << __PRETTY_FUNCTION__ << "Unknown Text Message: " << str << std::endl;
        sock.Send("Ack");
      }
    }
    else if (mess->What() == kMESS_OBJECT)
    {
      // this reads the message and allocate space for new histogram
      TH1 *histo = static_cast<TH1 *>(mess->ReadObjectAny(mess->GetClass()));
      delete mess;
      if (verbosity > 1)
      {
        std::cout << __PRETTY_FUNCTION__ << "histoname: " << histo->GetName() << " at "
                  << histo << std::endl;
      }

      updateHistoMap(what, histo->GetName(), histo);
      sock.Send("Ack");
    }
  }
  sock.Send("Finished");  // tell server we are finished

  // Close the socket
  sock.Close();
  return 0;
}

int OnlMonClient::requestMonitorList(const std::string &hostname, const int moniport)
{
  TSocket sock(hostname.c_str(), moniport);
  TMessage *mess;
  sock.Send("LISTMONITORS");
  sock.Recv(mess);
  if (!mess)  // if server is not up mess is NULL
  {
    std::cout << __PRETTY_FUNCTION__ << "Server not running on " << hostname << std::endl;
    sock.Close();
    return 1;
  }
  delete mess;
  while (true)
  {
    sock.Recv(mess);
    if (mess->What() == kMESS_STRING)
    {
      char strmess[OnlMonDefs::MSGLEN];
      mess->ReadString(strmess, OnlMonDefs::MSGLEN);
      delete mess;
      if (Verbosity() > 1)
      {
        std::cout << "received " << strmess << std::endl;
      }
      std::string str(strmess);
      if (str == "Finished")
      {
        break;
      }
      if (Verbosity() > 2)
      {
        std::cout << "inserting " << str << " on host " << hostname
                  << " listening to " << moniport << std::endl;
      }
      MonitorHostPorts.insert(std::make_pair(str, std::make_pair(hostname, moniport)));
    }
    else
    {
      std::cout << "requestMonitorList: received unexpected message type: " << mess->What() << std::endl;
      break;
    }
  }
  sock.Send("Finished");  // tell server we are finished

  // Close the socket
  sock.Close();
  return 0;
}

int OnlMonClient::requestHistoList(const std::string &subsys, const std::string &hostname, const int moniport, std::list<std::string> &histolist)
{
  // Open connection to server
  TSocket sock(hostname.c_str(), moniport);
  TMessage *mess;
  sock.Send("LIST");
  std::list<std::string>::const_iterator listiter;
  sock.Recv(mess);
  if (!mess)  // if server is not up mess is NULL
  {
    std::cout << __PRETTY_FUNCTION__ << "Server not running on " << hostname << std::endl;
    sock.Close();
    return 1;
  }

  delete mess;
  for (listiter = histolist.begin(); listiter != histolist.end(); ++listiter)
  {
    if (Verbosity() > 2)
    {
      std::cout << __PRETTY_FUNCTION__ << "asking for " << *listiter << std::endl;
    }
    sock.Send((*listiter).c_str());
    sock.Recv(mess);
    if (!mess)
    {
      std::cout << __PRETTY_FUNCTION__ << "Server shut down during getting histo list" << std::endl;
      sock.Close();
      return 1;
    }
    if (mess->What() == kMESS_STRING)
    {
      char str[OnlMonDefs::MSGLEN];
      mess->ReadString(str, OnlMonDefs::MSGLEN);
      delete mess;
      if (verbosity > 1)
      {
        std::cout << __PRETTY_FUNCTION__ << "Message: " << str << std::endl;
      }

      if (!strcmp(str, "Ack"))
      {
        break;
      }
      else if (!strcmp(str, "UnknownHisto"))
      {
        break;
      }
      else
      {
        std::cout << __PRETTY_FUNCTION__ << "Unknown Text Message: " << str << std::endl;
      }
    }
    else if (mess->What() == kMESS_OBJECT)
    {
      // this reads the message and allocate space for new histogram
      TH1 *histo = static_cast<TH1 *>(mess->ReadObjectAny(mess->GetClass()));
      delete mess;
      if (verbosity > 1)
      {
        std::cout << __PRETTY_FUNCTION__ << "histoname: " << histo->GetName() << " at "
                  << histo << std::endl;
      }

      updateHistoMap(subsys, histo->GetName(), histo);
    }
  }
  sock.Send("alldone");
  sock.Recv(mess);
  delete mess;
  sock.Send("Finished");  // tell server we are finished

  // Close the socket
  sock.Close();
  return 0;
}

void OnlMonClient::updateHistoMap(const std::string &subsys, const std::string &hname, TH1 *h1d)
{
  auto subsysiter = SubsysHisto.find(subsys);
  if (subsysiter == SubsysHisto.end())
  {
    std::map<const std::string, ClientHistoList *> newmap;
    ClientHistoList *entry = new ClientHistoList(subsys);
    entry->Histo(h1d);
    newmap.insert(std::make_pair(hname, entry));
    SubsysHisto.insert(std::make_pair(subsys, newmap));
    return;
  }
  auto histoiter = subsysiter->second.find(hname);
  //  std::map<const std::string, ClientHistoList *>::const_iterator histoiter = Histo.find(hname);
  if (histoiter != subsysiter->second.end())
  {
    if (Verbosity() > 2)
    {
      std::cout << "deleting histogram " << hname << " at " << Histo[hname] << std::endl;
    }
    delete histoiter->second->Histo();  // delete old histogram
    histoiter->second->Histo(h1d);
  }
  else
  {
    ClientHistoList *newhisto = new ClientHistoList(subsys);
    newhisto->Histo(h1d);
    subsysiter->second.insert(std::make_pair(hname, newhisto));
    if (Verbosity() > 2)
    {
      std::cout << "new histogram " << hname << " at " << newhisto->Histo() << std::endl;
    }
  }
  return;
}

OnlMonDraw *
OnlMonClient::getDrawer(const std::string &name)
{
  std::map<const std::string, OnlMonDraw *>::iterator iter = DrawerList.find(name);
  if (iter != DrawerList.end())
  {
    return iter->second;
  }
  std::cout << "Could not locate drawer" << name << std::endl;
  return nullptr;
}

TH1 *OnlMonClient::getHisto(const std::string &monitor, const std::string &hname)
{
  auto subsysiter = SubsysHisto.find(monitor);
  if (subsysiter == SubsysHisto.end())
  {
    return nullptr;
  }
  auto hiter = subsysiter->second.find(hname);
  if (hiter == subsysiter->second.end())
  {
    return nullptr;
  }
  return hiter->second->Histo();
}

void OnlMonClient::Print(const char *what)
{
  if (!strcmp(what, "ALL") || !strcmp(what, "DRAWER"))
  {
    // loop over the map and print out the content (name and location in memory)
    std::cout << "--------------------------------------" << std::endl
              << std::endl;
    std::cout << "List of Drawers in OnlMonClient:" << std::endl;

    std::map<const std::string, OnlMonDraw *>::const_iterator hiter;
    for (hiter = DrawerList.begin(); hiter != DrawerList.end(); ++hiter)
    {
      std::cout << hiter->first << " is at " << hiter->second << std::endl;
    }
    std::cout << std::endl;
  }
  if (!strcmp(what, "ALL") || !strcmp(what, "SERVERS"))
  {
    // loop over the map and print out the content (name and location in memory)
    std::cout << "--------------------------------------" << std::endl
              << std::endl;
    std::cout << "List of Servers in OnlMonClient:" << std::endl;

    std::vector<std::string>::iterator hostiter;
    for (hostiter = MonitorHosts.begin(); hostiter != MonitorHosts.end(); ++hostiter)
    {
      std::cout << "ServerHost: " << *hostiter << std::endl;
    }
  }
  if (!strcmp(what, "ALL") || !strcmp(what, "MONITORS"))
  {
    // loop over the map and print out the content (name and location in memory)
    std::cout << "--------------------------------------" << std::endl
              << std::endl;
    std::cout << "List of Monitors in OnlMonClient:" << std::endl;

    for (auto &moniiter : MonitorHostPorts)
    {
      std::cout << "Monitor " << moniiter.first << " runs on " << moniiter.second.first
                << " listening to port " << moniiter.second.second << std::endl;
    }
  }
  if (!strcmp(what, "ALL") || !strcmp(what, "HISTOS"))
  {
    // loop over the map and print out the content (name and location in memory)
    std::cout << "--------------------------------------" << std::endl
              << std::endl;
    std::cout << "List of Histograms in OnlMonClient:" << std::endl;

    std::map<const std::string, ClientHistoList *>::const_iterator hiter;
    for (hiter = Histo.begin(); hiter != Histo.end(); ++hiter)
    {
      std::cout << hiter->first << " Address " << hiter->second->Histo()
                << " on host " << hiter->second->ServerHost()
                << " port " << hiter->second->ServerPort()
                << ", subsystem " << hiter->second->SubSystem() << std::endl;
    }
    std::cout << std::endl;
    for (auto &subs : SubsysHisto)
    {
      auto subiter = MonitorHostPorts.find(subs.first);
      if (subiter != MonitorHostPorts.end())
      {
        std::cout << "Subsystem " << subs.first << " runs on " << subiter->second.first;
        std::cout << " on port " << subiter->second.second << std::endl;
      }
      for (auto &histos : subs.second)
        std::cout << histos.first << " @ " << subs.first
                  << " Address " << histos.second->Histo()
                  << " on host " << histos.second->ServerHost()
                  << " port " << histos.second->ServerPort()
                  << ", subsystem " << histos.second->SubSystem() << std::endl;
    }
    std::cout << std::endl;
  }
  if (!strcmp(what, "ALL") || !strcmp(what, "UNKNOWN"))
  {
    // loop over the map and print out the content (name and location in memory)
    std::cout << "--------------------------------------" << std::endl
              << std::endl;
    std::cout << "List of Unknown Histograms in OnlMonClient:" << std::endl;

    std::map<const std::string, ClientHistoList *>::const_iterator hiter;
    for (hiter = Histo.begin(); hiter != Histo.end(); ++hiter)
    {
      if (hiter->second->ServerHost() == "UNKNOWN" ||
          hiter->second->SubSystem() == "UNKNOWN")
      {
        std::cout << hiter->first << " Address " << hiter->second->Histo()
                  << " on host " << hiter->second->ServerHost()
                  << " port " << hiter->second->ServerPort()
                  << ", subsystem " << hiter->second->SubSystem() << std::endl;
      }
    }
    std::cout << std::endl;
  }
  return;
}

int OnlMonClient::UpdateServerHistoMap(const std::string &hname, const std::string &subsys, const std::string &hostname)
{
  // Open connection to server
  std::string searchstring = subsys + ' ' + hname;
  int MoniPort = OnlMonDefs::MONIPORT;
  int foundit = 0;
  do
  {
    std::cout << "Connecting to " << hostname << ", if it is frozen here - this is the one you need to restart" << std::endl;
    TSocket sock(hostname.c_str(), MoniPort);
    TMessage *mess;
    if (verbosity > 0)
    {
      std::cout << "UpdateServerHistoMap: sending cmd HistoList to "
                << hostname << " on port "
                << MoniPort
                << std::endl;
    }
    if (!sock.Send("HistoList"))
    {
      std::cout << "Server not running on " << hostname
                << " port " << MoniPort << std::endl;
      goto noserver;
    }
    while (true)
    {
      if (verbosity > 0)
      {
        std::cout << "UpdateServerHistoMap: waiting for response on "
                  << hostname << " on port "
                  << MoniPort
                  << std::endl;
      }
      sock.Recv(mess);
      if (!mess)  // if server is not up mess is NULL
      {
        std::cout << "UpdateServerHistoMap: No Recv, Server not running on "
                  << hostname
                  << " on port " << MoniPort << std::endl;
        goto noserver;
      }
      if (mess->What() == kMESS_STRING)
      {
        char strchr[OnlMonDefs::MSGLEN];
        mess->ReadString(strchr, OnlMonDefs::MSGLEN);
        delete mess;
        std::string str = strchr;
        if (verbosity > 0)
        {
          std::cout << "UpdateServerHistoMap Response: " << str << std::endl;
        }
        if (str == "Finished")
        {
          break;
        }

        if (str == searchstring)
        {
          std::cout << "found subsystem " << subsys << " histo " << hname << std::endl;
          foundit = 1;
        }
        unsigned int pos_space = str.find(' ');
        PutHistoInMap(str.substr(pos_space + 1, str.size()), str.substr(0, pos_space), hostname, MoniPort);
        sock.Send("Ack");
      }
    }
    sock.Send("Finished");  // tell server we are finished

    // Close the socket
  noserver:
    sock.Close();
    if (foundit)
    {
      return foundit;
    }
    MoniPort++;
  } while ((MoniPort - OnlMonDefs::MONIPORT) < OnlMonDefs::NUMMONIPORT);  // no more than NUMMONIPORT parallel servers
  return foundit;
}

void OnlMonClient::PutHistoInMap(const std::string &hname, const std::string &subsys, const std::string &hostname, const int port)
{
  auto hiter = SubsysHisto.find(subsys);
  if (hiter == SubsysHisto.end())
  {
    std::map<const std::string, ClientHistoList *> entry;
    // c++11 map.insert returns pair<iterator,bool>
    hiter = SubsysHisto.insert(std::make_pair(subsys, entry)).first;
  }
  auto histoiter = hiter->second.find(hname);
  if (histoiter != hiter->second.end())
  {
    histoiter->second->ServerHost(hostname);
    histoiter->second->ServerPort(port);
  }
  else
  {
    ClientHistoList *newhisto = new ClientHistoList();
    newhisto->ServerHost(hostname);
    newhisto->ServerPort(port);
    hiter->second.insert(std::make_pair(hname, newhisto));
  }
  return;
}

void OnlMonClient::AddServerHost(const std::string &hostname)
{
  if (find(MonitorHosts.begin(), MonitorHosts.end(), hostname) != MonitorHosts.end())
  {
    std::cout << "Host " << hostname << " already in list" << std::endl;
  }
  else
  {
    MonitorHosts.emplace_back(hostname);
  }
  return;
}

int OnlMonClient::LocateHistogram(const std::string &hname, const std::string &subsys)
{
  for (auto &hostiter : MonitorHosts)
  {
    if (UpdateServerHistoMap(hname, subsys, hostiter) > 0)
    {
      return 1;
    }
  }
  return 0;
}

int OnlMonClient::RunNumber()
{
  int runno = -9999;
  for (auto frwrkiter : m_MonitorFetchedSet)
  {
    TH1 *frameworkvars = getHisto(frwrkiter, "FrameWorkVars");
    if (frameworkvars)
    {
      runno = std::max(runno, (int) frameworkvars->GetBinContent(RUNNUMBERBIN));
    }
  }
  return (runno);
}

time_t OnlMonClient::EventTime(const std::string &which)
{
  time_t tret = 0;
  for (auto frwrkiter : m_MonitorFetchedSet)
  {
    tret = std::max(tret, EventTime(frwrkiter,which));
  }

  if (verbosity > 0)
  {
    std::cout << "Time is " << ctime(&tret) << std::endl;
  }
  return (tret);
}

time_t OnlMonClient::EventTime(const std::string &servername, const std::string &which)
{
  time_t tret = 0;
  int ibin = 0;
  if (which == "BOR")
  {
    ibin = BORTIMEBIN;
  }
  else if (which == "CURRENT")
  {
    ibin = CURRENTTIMEBIN;
  }
  else if (which == "EOR")
  {
    ibin = EORTIMEBIN;
  }
  else
  {
    std::cout << "Bad Option for Time: " << which
              << ", implemented are BOR EOR CURRENT" << std::endl;
    ibin = CURRENTTIMEBIN;
  }
  TH1 *frameworkvars = getHisto(servername, "FrameWorkVars");
  if (frameworkvars == nullptr)
  {
    tret = 0;
  }
  else
  {
    tret = (time_t) frameworkvars->GetBinContent(ibin);
  }
  if (verbosity > 0)
  {
    std::cout << "Time is " << ctime(&tret) << std::endl;
  }
  return (tret);
}

int OnlMonClient::ReadHistogramsFromFile(const std::string &filename)
{
  std::string subsys = ExtractSubsystem(filename);
  TDirectory *save = gDirectory;  // save current dir (which will be overwritten by the following fileopen)
  TFile *histofile = new TFile(filename.c_str(), "READ");
  if (!histofile)
  {
    std::cout << "Can't open " << filename << std::endl;
    return -1;
  }
  save->cd();
  TIterator *titer = histofile->GetListOfKeys()->MakeIterator();
  TObject *obj;
  TH1 *histo, *histoptr;
  while ((obj = titer->Next()))
  {
    if (verbosity > 0)
    {
      std::cout << "TObject at " << obj << std::endl;
      std::cout << obj->GetName() << std::endl;
      std::cout << obj->ClassName() << std::endl;
    }
    histofile->GetObject(obj->GetName(), histoptr);
    if (histoptr)
    {
      histo = static_cast<TH1 *>(histoptr->Clone());
      updateHistoMap(subsys, histo->GetName(), histo);
      if (verbosity > 0)
      {
        std::cout << "HistoName: " << histo->GetName() << std::endl;
        std::cout << "HistoClass: " << histo->ClassName() << std::endl;
      }
    }
  }
  delete histofile;
  delete titer;
  return 0;
}

int OnlMonClient::SendCommand(const char *hostname, const int port, const char *cmd)
{
  // Open connection to server
  TSocket sock(hostname, port);
  TMessage *mess;
  if (!sock.Send(cmd))
  {
    std::cout << "Server not running on " << hostname
              << " port " << port << std::endl;
    sock.Close();
    return -1;
  }
  while (true)
  {
    sock.Recv(mess);
    if (!mess)  // if server is not up mess is NULL
    {
      std::cout << __PRETTY_FUNCTION__ << "No Recv, Server not running on " << hostname
                << " on port " << port << std::endl;
      sock.Close();
      return -1;
    }
    if (mess->What() == kMESS_STRING)
    {
      char str[OnlMonDefs::MSGLEN];
      mess->ReadString(str, OnlMonDefs::MSGLEN);
      delete mess;
      if (verbosity > 1)
      {
        std::cout << __PRETTY_FUNCTION__ << " Message: " << str << std::endl;
      }

      if (!strcmp(str, "Finished"))
      {
        break;
      }
    }
  }
  sock.Send("Finished");  // tell server we are finished

  // Close the socket
  sock.Close();
  return 0;
}

//_____________________________________________________________________________
void OnlMonClient::Verbosity(const int v)
{
  verbosity = v;
  if (fHtml)
  {
    fHtml->verbosity(v);
  }
}

//_____________________________________________________________________________
void OnlMonClient::htmlAddMenu(const OnlMonDraw &drawer,
                               const std::string &path,
                               const std::string &relfilename)
{
  fHtml->addMenu(drawer.Name(), path, relfilename);
}

//_____________________________________________________________________________
void OnlMonClient::htmlNamer(const OnlMonDraw &drawer,
                             const std::string &basefilename,
                             const std::string &ext,
                             std::string &fullfilename,
                             std::string &filename)
{
  fHtml->namer(drawer.Name(), basefilename, ext, fullfilename, filename);
}

//_____________________________________________________________________________
std::string
OnlMonClient::htmlRegisterPage(const OnlMonDraw &drawer,
                               const std::string &path,
                               const std::string &basefilename,
                               const std::string &ext)
{
  return fHtml->registerPage(drawer.Name(), path, basefilename, ext);
}

int OnlMonClient::CanvasToPng(TCanvas *canvas, std::string const &pngfilename)
{
  // in order not to collide when running multiple html generators
  // create a unique filename (okay tempnam is not totally safe against
  // multiple procs getting the same name but the local /tmp filesystem should
  // prevent at least multiple machines colliding)
  // returned char array needs to be free'd after use
  if (!canvas)
  {
    std::cout << __PRETTY_FUNCTION__ << " TCanvas is Null Pointer" << std::endl;
    return -2;
  }
  if (pngfilename.empty())
  {
    std::cout << __PRETTY_FUNCTION__ << " emtpy png filename, not saving TCanvas "
              << canvas->GetName() << std::endl;
    return -1;
  }
  uuid_t uu;
  uuid_generate(uu);
  char uuid[50];
  uuid_unparse(uu, uuid);
  std::string tmpname = "/tmp/TC" + std::string(uuid);
  canvas->Print(tmpname.c_str(), "gif");  // write gif format
  TImage *img = TImage::Open(tmpname.c_str());
  img->WriteImage(pngfilename.c_str());
  delete img;
  if (remove(tmpname.c_str()))
  {
    std::cout << "Error removing " << tmpname << std::endl;
  }
  return 0;
}

int OnlMonClient::HistoToPng(TH1 *histo, std::string const &pngfilename, const char *drawopt, const int statopt)
{
  TCanvas *cgiCanv = new TCanvas("cgiCanv", "cgiCanv", 200, 200, 650, 500);
  gStyle->SetOptStat(statopt);
  cgiCanv->SetFillColor(0);
  cgiCanv->SetBorderMode(0);
  cgiCanv->SetBorderSize(2);
  cgiCanv->SetFrameFillColor(0);
  cgiCanv->SetFrameBorderMode(0);
  cgiCanv->SetTickx();
  cgiCanv->SetTicky();
  cgiCanv->cd();
  histo->SetMarkerStyle(8);
  histo->SetMarkerSize(0.15);
  histo->Draw(drawopt);
  uuid_t uu;
  uuid_generate(uu);
  char uuid[50];
  uuid_unparse(uu, uuid);
  std::string tmpname = "/tmp/TC" + std::string(uuid);
  cgiCanv->Print(tmpname.c_str(), "gif");
  TImage *img = TImage::Open(tmpname.c_str());
  img->WriteImage(pngfilename.c_str());
  if (remove(tmpname.c_str()))
  {
    std::cout << "Error removing " << tmpname << std::endl;
  }
  delete cgiCanv;
  return 0;
}

int OnlMonClient::SaveLogFile(const OnlMonDraw &drawer)
{
  // sendfile example shamelessly copied from
  // http://www.linuxgazette.com/issue91/tranter.html
  std::ostringstream logfilename;
  const char *logdir = getenv("ONLMON_LOGDIR");
  if (logdir)
  {
    logfilename << logdir << "/";
  }
  int irun = RunNumber();
  logfilename << drawer.Name() << "_" << irun << ".log.gz";
  std::ifstream infile(logfilename.str().c_str(), std::ios_base::binary);
  if (infile.good())
  {
    std::string outfilename = htmlRegisterPage(drawer, "Logfile", "log", "txt.gz");
    std::ofstream outfile(outfilename.c_str(), std::ios_base::binary);
    outfile << infile.rdbuf();
    infile.close();
    outfile.close();
  }
  return 0;
}

int OnlMonClient::SetStyleToDefault()
{
  defaultStyle->cd();
  defaultStyle->Reset();
  defaultStyle->SetFrameBorderMode(0);
  defaultStyle->SetCanvasColor(0);
  defaultStyle->SetPadBorderMode(0);
  defaultStyle->SetCanvasBorderMode(0);
  defaultStyle->SetPalette(1, nullptr);
  return 0;
}

void OnlMonClient::CacheRunDB(const int runnoinput)
{
  int runno = 221;
  if (runnoinput == 221)
  {
    runno = runnoinput;
  }
  if (runno == cachedrun)
  {
    return;
  }
  std::string TriggerConfig = "UNKNOWN";
  standalone = 0;
  cosmicrun = 0;
  runtype = "UNKNOWN";

  if (runno == 0xFEE2DCB)  // dcm2 standalone runs have this runnumber
  {
    standalone = 1;
    return;
  }
  odbc::Connection *con = nullptr;
  odbc::Statement *query = nullptr;
  std::ostringstream cmd;
  try
  {
    con = odbc::DriverManager::getConnection("daq", "phnxrc", "");
  }
  catch (odbc::SQLException &e)
  {
    printf(" Exception caught during DriverManager::getConnection, Message: %s\n", e.getMessage().c_str());
    return;
  }

  query = con->createStatement();
  cmd << "select runnumber from run where runnumber = " << runno;
  odbc::ResultSet *rs = nullptr;
  int ncount = 10;
  while (ncount > 0)
  {
    try
    {
      rs = query->executeQuery(cmd.str());
    }
    catch (odbc::SQLException &e)
    {
      printf("Exception caught for query %s\nMessage: %s", cmd.str().c_str(), e.getMessage().c_str());
    }
    if (!rs->next())
    {
      printf("run table query did not give any result, run %d not in DB yet countdown %d\n", runno, ncount);
      ncount--;
      sleep(10);
      delete rs;
    }
    else
    {
      delete rs;
      break;
    }
  }
  cmd.str("");
  cmd << "SELECT triggerconfig,brunixtime,runtype FROM RUN  WHERE RUNNUMBER = "
      << runno;
  if (verbosity > 0)
  {
    printf("command: %s\n", cmd.str().c_str());
  }
  try
  {
    rs = query->executeQuery(cmd.str());
  }
  catch (odbc::SQLException &e)
  {
    printf("Exception caught for query %s\nMessage: %s", cmd.str().c_str(), e.getMessage().c_str());
  }
  if (rs->next())
  {
    runtype = rs->getString("runtype");
    TriggerConfig = rs->getString("triggerconfig");
    if (TriggerConfig == "StandAloneMode")
    {
      standalone = 1;
    }
    else
    {
      standalone = 0;
    }
    if (TriggerConfig.find("Cosmic") != std::string::npos)
    {
      cosmicrun = 1;
    }
    else
    {
      cosmicrun = 0;
    }
  }
  delete con;
  //  printf("CacheRunDB: runno: %d\n",runno);
  return;
}

int OnlMonClient::isCosmicRun()
{
  CacheRunDB(RunNumber());
  return cosmicrun;
}

int OnlMonClient::isStandalone()
{
  CacheRunDB(RunNumber());
  return standalone;
}

std::string
OnlMonClient::RunType()
{
  CacheRunDB(RunNumber());
  return runtype;
}

void OnlMonClient::FindAllMonitors()
{
  for (auto &hostiter : MonitorHosts)
  {
    if (Verbosity() > 2)
    {
      std::cout << "checking " << hostiter << std::endl;
    }
    for (unsigned int moniport = OnlMonDefs::MONIPORT; moniport < OnlMonDefs::MONIPORT + OnlMonDefs::NUMMONIPORT; ++moniport)
    {
      requestMonitorList(hostiter, moniport);
    }
  }
  return;
}

int OnlMonClient::FindMonitor(const std::string &name)
{
  // loop over all hosts/ports until we find ours
  int iret = 0;
  for (auto &hostiter : MonitorHosts)
  {
    if (Verbosity() > 2)
    {
      std::cout << "checking " << hostiter << std::endl;
    }
    for (unsigned int moniport = OnlMonDefs::MONIPORT; moniport < OnlMonDefs::MONIPORT + OnlMonDefs::NUMMONIPORT; ++moniport)
    {
      requestMonitorList(hostiter, moniport);
      if (Verbosity() > 2)
      {
        std::cout << "looking for " << name << std::endl;
      }
      auto moniter = MonitorHostPorts.find(name);
      if (moniter != MonitorHostPorts.end())
      {
        if (Verbosity() > 2)
        {
          std::cout << "found " << name << " running on " << moniter->second.first
                    << " listening to port " << moniter->second.second << std::endl;
        }
        return 1;
      }
    }
  }
  return iret;
}

int OnlMonClient::IsMonitorRunning(const std::string &name)
{
  int iret = 0;
  std::string command = std::string("ISRUNNING") + ' ' + name;
  auto moniter = MonitorHostPorts.find(name);
  if (moniter == MonitorHostPorts.end())
  {
    return iret;
  }
  TSocket sock(moniter->second.first.c_str(), moniter->second.second);
  TMessage *mess;
  sock.Send(command.c_str());
  sock.Recv(mess);
  if (!mess)  // if server is not up mess is NULL
  {
    std::cout << __PRETTY_FUNCTION__ << "Server not running on " << moniter->second.first << std::endl;
    sock.Close();
    return iret;
  }
  if (mess->What() == kMESS_STRING)
  {
    char str[OnlMonDefs::MSGLEN];
    mess->ReadString(str, OnlMonDefs::MSGLEN);
    delete mess;
    if (verbosity > 1)
    {
      std::cout << __PRETTY_FUNCTION__ << "Message: " << str << std::endl;
    }
    if (!strcmp(str, "Yes"))
    {
      iret = 1;
    }
  }
  sock.Send("Finished");
  sock.Close();
  return iret;
}

std::string OnlMonClient::ExtractSubsystem(const std::string &fullfilename)
{
  std::string subsys = std::filesystem::path(fullfilename).filename();
  subsys = subsys.substr(subsys.find('-') + 1);
  subsys = subsys.substr(0, subsys.find(".root"));
  m_MonitorFetchedSet.insert(subsys);
  return subsys;
}
