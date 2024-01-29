#include "OnlMonServer.h"

#include "OnlMon.h"
#include "OnlMonStatusDB.h"

#include "MessageSystem.h"

#include <Event/msg_profile.h>  // for MSG_SEV_ERROR, MSG_SEV...

#include <TFile.h>
#include <TH1.h>
#include <TROOT.h>

#include <odbc++/connection.h>
#include <odbc++/drivermanager.h>
#include <odbc++/resultset.h>
#include <odbc++/statement.h>  // for Statement
#include <odbc++/types.h>      // for SQLException

#include <zlib.h>

#include <sys/utsname.h>
#include <unistd.h>   // for sleep
#include <algorithm>  // for max
#include <cstdio>     // for printf
#include <cstdlib>
#include <cstring>  // for strcmp
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <utility>  // for pair
#include <vector>

OnlMonServer *OnlMonServer::__instance = nullptr;

// root TFile compression level
static const int compression_level = 3;

OnlMonServer *OnlMonServer::instance()
{
  if (__instance)
  {
    return __instance;
  }
  __instance = new OnlMonServer("OnlMonServer");
  return __instance;
}

OnlMonServer::OnlMonServer(const std::string &name)
  : OnlMonBase(name)
{
#ifdef USE_MUTEX
  pthread_mutex_init(&mutex, nullptr);
#endif
  MsgSystem[ThisName] = new MessageSystem(ThisName);
  statusDB = new OnlMonStatusDB();
  RunStatusDB = new OnlMonStatusDB("onlmonrunstatus");
  InitAll();
  return;
}

OnlMonServer::~OnlMonServer()
{
#ifdef USE_MUTEX
  pthread_mutex_lock(&mutex);
#endif
  if (int tret = pthread_cancel(serverthreadid))
  {
    std::cout << __PRETTY_FUNCTION__ << "pthread cancel returned error: " << tret << std::endl;
  }
  delete serverrunning;

#ifdef USE_MUTEX
  pthread_mutex_destroy(&mutex);
#endif
  while (MonitorList.begin() != MonitorList.end())
  {
    delete MonitorList.back();
    MonitorList.pop_back();
  }
  delete statusDB;
  delete RunStatusDB;
  while(MonitorHistoSet.begin() !=  MonitorHistoSet.end())
  {
    while(MonitorHistoSet.begin()->second.begin() != MonitorHistoSet.begin()->second.end())
    {
      if (CommonHistoMap.find(MonitorHistoSet.begin()->second.begin()->second->GetName()) == CommonHistoMap.end())
      {
      delete MonitorHistoSet.begin()->second.begin()->second;
      }
      else
      {
	if (Verbosity() > 2)
	{
	std::cout << "not deleting " << MonitorHistoSet.begin()->second.begin()->second->GetName() << std::endl;
	}
      }
      MonitorHistoSet.begin()->second.erase(MonitorHistoSet.begin()->second.begin());
    }
    MonitorHistoSet.erase(MonitorHistoSet.begin());
  }
  while(CommonHistoMap.begin() != CommonHistoMap.end())
  {
    delete CommonHistoMap.begin()->second;
    CommonHistoMap.erase(CommonHistoMap.begin());
  }
  while (MsgSystem.begin() != MsgSystem.end())
  {
    delete MsgSystem.begin()->second;
    MsgSystem.erase(MsgSystem.begin());
  }
  return;
}

void OnlMonServer::InitAll()
{
  if (gROOT->FindObject("ClientRunning"))
  {
    std::ostringstream msg;
    msg << "Don't run Server and Client in same session, exiting";
    send_message(MSG_SEV_FATAL, msg.str(), 1);
    exit(1);
  }
  serverrunning = new TH1F("ServerRunning", "ServerRunning", 1, 0, 1);
  return;
}

void OnlMonServer::dumpHistos(const std::string &filename)
{
  TFile *hfile = TFile::Open(filename.c_str(), "RECREATE", "Created by Online Monitor");
  for (auto &moniiter : MonitorHistoSet)
  {
    std::cout << "saving " << moniiter.first << std::endl;
    for (auto &histiter : moniiter.second)
    {
      std::cout << "saving " << histiter.first << std::endl;
      histiter.second->Write();
    }
  }
  hfile->Close();
  delete hfile;
  return;
}

void OnlMonServer::registerCommonHisto(TH1 *h1d)
{
  CommonHistoMap.insert(std::make_pair(h1d->GetName(), h1d));
  return;
}

void OnlMonServer::registerHisto(const OnlMon *monitor, TH1 *h1d)
{
  registerHisto(monitor->Name(), h1d->GetName(), h1d, 0);
  return;
}

void OnlMonServer::registerHisto(const std::string &monitorname, const std::string &hname, TH1 *h1d, const int replace)
{
  if (hname.find(' ') != std::string::npos)
  {
    std::cout << "No empty spaces in registered histogram names : " << hname << std::endl;
    exit(1);
  }
  auto moniiter = MonitorHistoSet.find(monitorname);
  if (moniiter == MonitorHistoSet.end())
  {
    std::map<std::string, TH1 *> histo;
    histo[hname] = h1d;
    std::cout << __PRETTY_FUNCTION__ << " inserting " << monitorname << " hname " << hname << std::endl;
    MonitorHistoSet.insert(std::make_pair(monitorname, histo));
    return;
  }
  auto histoiter = moniiter->second.find(hname);
  if (histoiter == moniiter->second.end())
  {
    moniiter->second.insert(std::make_pair(hname, h1d));
  }
  else
  {
    if (replace)
    {
      delete histoiter->second;
      histoiter->second = h1d;
    }
    else
    {
      std::cout << "Histogram " << hname << " already registered with " << monitorname
                << ", it will not be overwritten" << std::endl;
    }
  }
  return;
}

void OnlMonServer::registerHisto(const std::string &hname, TH1 *h1d, const int replace)
{
  if (hname.find(' ') != std::string::npos)
  {
    std::cout << "No empty spaces in registered histogram names : " << hname << std::endl;
    exit(1);
  }
  const std::string &tmpstr = hname;
  std::map<const std::string, TH1 *>::const_iterator histoiter = CommonHistoMap.find(tmpstr);
  std::ostringstream msg;
  int histoexist;
  TH1 *delhis;
  if (histoiter != CommonHistoMap.end())
  {
    delhis = histoiter->second;
    histoexist = 1;
  }
  else
  {
    delhis = nullptr;
    histoexist = 0;
  }
  if (histoexist && replace == 0)
  {
    msg << "Histogram " << hname << " already registered, I won't overwrite it";
    send_message(MSG_SEV_WARNING, msg.str(), 2);
    msg.str("");
    msg << "Use a different name and try again" << std::endl;
    send_message(MSG_SEV_WARNING, msg.str(), 2);
  }
  else
  {
    if (verbosity > 1)
    {
      if (strcmp(h1d->GetName(), tmpstr.c_str()))
      {
        msg.str("");
        msg << __PRETTY_FUNCTION__ << "Histogram " << h1d->GetName()
            << " at " << h1d << " renamed to " << tmpstr;
        send_message(MSG_SEV_INFORMATIONAL, msg.str(), 3);
      }
    }
    CommonHistoMap[tmpstr] = h1d;
    if (delhis)
    {
      delete delhis;
    }
    h1d->SetName(tmpstr.c_str());
  }
  return;
}

OnlMon *
OnlMonServer::getMonitor(const std::string &name)
{
  for (OnlMon *mon : MonitorList)
  {
    if (name == mon->Name())
    {
      return mon;
    }
  }
  std::cout << "Could not locate monitor" << name << std::endl;
  return nullptr;
}

void OnlMonServer::registerMonitor(OnlMon *Monitor)
{
  for (OnlMon *mon : MonitorList)
  {
    if (Monitor->Name() == mon->Name())
    {
      std::ostringstream msg;
      msg << "Monitor " << Monitor->Name() << " already registered, I won't overwrite it";
      send_message(MSG_SEV_SEVEREERROR, msg.str(), 4);
      msg.str("");
      msg << "Use a different name and try again";
      send_message(MSG_SEV_SEVEREERROR, msg.str(), 4);
      return;
    }
  }
  MonitorList.push_back(Monitor);
  MsgSystem[Monitor->Name()] = new MessageSystem(Monitor->Name());
  Monitor->InitCommon(this);
  Monitor->Init();
  return;
}

TH1 *OnlMonServer::getHisto(const unsigned int ihisto) const
{
  std::map<const std::string, TH1 *>::const_iterator histoiter = CommonHistoMap.begin();
  unsigned int size = CommonHistoMap.size();
  if (Verbosity() > 3)
  {
    std::ostringstream msg;

    msg << "Map contains " << size << " Elements";
    send_message(MSG_SEV_INFORMATIONAL, msg.str(), 5);
  }
  if (ihisto < size)
  {
    for (unsigned int i = 0; i < ihisto; i++)
    {
      ++histoiter;
    }
    return histoiter->second;
  }
  else
  {
    std::ostringstream msg;
    msg << "OnlMonServer::getHisto: ERROR Invalid histogram number: "
        << ihisto << ", maximum number is " << size;
    send_message(MSG_SEV_ERROR, msg.str(), 6);
  }
  return nullptr;
}

const std::string
OnlMonServer::getHistoName(const unsigned int ihisto) const
{
  std::map<const std::string, TH1 *>::const_iterator histoiter = CommonHistoMap.begin();
  unsigned int size = CommonHistoMap.size();
  if (verbosity > 3)
  {
    std::ostringstream msg;
    msg << "Map contains " << size << " Elements";
    send_message(MSG_SEV_INFORMATIONAL, msg.str(), 5);
  }
  if (ihisto < size)
  {
    for (unsigned int i = 0; i < ihisto; i++)
    {
      ++histoiter;
    }
    return histoiter->first.c_str();
  }
  else
  {
    std::ostringstream msg;
    msg << "OnlMonServer::getHisto: ERROR Invalid histogram number: "
        << ihisto << ", maximum number is " << size;
    send_message(MSG_SEV_ERROR, msg.str(), 6);
  }
  return nullptr;
}

TH1 *OnlMonServer::getHisto(const std::string &subsys, const std::string &hname) const
{
  if (Verbosity() > 2)
  {
    std::cout << __PRETTY_FUNCTION__ << " checking for subsys " << subsys << ", hname " << hname << std::endl;
  }
  auto moniiter = MonitorHistoSet.find(subsys);
  if (moniiter != MonitorHistoSet.end())
  {
    auto histoiter = moniiter->second.find(hname);
    if (histoiter != moniiter->second.end())
    {
      return histoiter->second;
    }
  }
  std::ostringstream msg;

  msg << "OnlMonServer::getHisto: ERROR Unknown Histogram " << hname
      << ", The following are implemented: ";
  send_message(MSG_SEV_ERROR, msg.str(), 7);
  Print("HISTOS");
  return nullptr;
}

TH1 *OnlMonServer::getCommonHisto(const std::string &hname) const
{
  std::map<const std::string, TH1 *>::const_iterator histoiter = CommonHistoMap.find(hname);
  if (histoiter != CommonHistoMap.end())
  {
    return histoiter->second;
  }
  std::ostringstream msg;

  msg << "OnlMonServer::getHisto: ERROR Unknown Histogram " << hname
      << ", The following are implemented: ";
  send_message(MSG_SEV_ERROR, msg.str(), 7);
  Print("HISTOS");
  return nullptr;
}

int OnlMonServer::run_empty(const int nevents)
{
  int iret = 0;
  for (int i = 0; i < nevents; i++)
  {
    iret = process_event(nullptr);
    if (iret)
    {
      break;
    }
  }
  return iret;
}

int OnlMonServer::process_event(Event *evt)
{
  int i = 0;
  std::vector<OnlMon *>::iterator iter;

  for (iter = MonitorList.begin(); iter != MonitorList.end(); ++iter)
  {
    i += (*iter)->process_event_common(evt);
  }
  for (iter = MonitorList.begin(); iter != MonitorList.end(); ++iter)
  {
    i += (*iter)->ResetEvent();
  }

  return i;
}

int OnlMonServer::Reset()
{
  int i = 0;
  std::vector<OnlMon *>::iterator iter;
  for (iter = MonitorList.begin(); iter != MonitorList.end(); ++iter)
  {
    i += (*iter)->Reset();
  }
  std::map<const std::string, TH1 *>::const_iterator hiter;
  for (auto &moniiter : MonitorHistoSet)
  {
     for (auto &histiter : moniiter.second)
    {
      histiter.second->Reset();
    }
  }

  for (hiter = CommonHistoMap.begin(); hiter != CommonHistoMap.end(); ++hiter)
  {
    hiter->second->Reset();
  }
  eventnumber = 0;
  std::map<std::string, MessageSystem *>::const_iterator miter;
  for (miter = MsgSystem.begin(); miter != MsgSystem.end(); ++miter)
  {
    miter->second->Reset();
  }
  return i;
}

int OnlMonServer::BeginRun(const int runno)
{
  int i = 0;
  i = CacheRunDB(runno);
  if (i)
  {
    printf("bad return code from CacheRunDB(%d): %d\n", runno, i);
    i = 0;
  }

  std::vector<OnlMon *>::iterator iter;
  activepacketsinit = 0;
  for (iter = MonitorList.begin(); iter != MonitorList.end(); ++iter)
  {
    (*iter)->BeginRunCommon(runno, this);
    i += (*iter)->BeginRun(runno);
  }
  return i;
}

int OnlMonServer::EndRun(const int runno)
{
  int i = 0;
  std::vector<OnlMon *>::iterator iter;
  for (iter = MonitorList.begin(); iter != MonitorList.end(); ++iter)
  {
    i += (*iter)->EndRun(runno);
  }
  return i;
}

void OnlMonServer::Print(const std::string &what, std::ostream& os) const
{
  if (what == "ALL" || what == "PORT")
  {
    utsname ThisNode;
    uname(&ThisNode);
    os << "--------------------------------------" << std::endl << std::endl;
    os << "Server running on " << ThisNode.nodename
       << " and is listening on port " << PortNumber() << std::endl
       << std::endl;
  }
  if (what == "ALL" || what == "HISTOS")
  {
    os << "--------------------------------------" << std::endl << std::endl;
    os << "List of Assigned histograms in OnlMonServer:" << std::endl << std::endl;
    for (auto &moniiter : MonitorHistoSet)
    {
      os << "Monitor " << moniiter.first << std::endl;
      for (auto &histiter : moniiter.second)
      {
        os <<  moniiter.first << " " << histiter.first
	   << " at " << histiter.second << std::endl;
      }
    }
    // loop over the map and print out the content (name and location in memory)
    os << std::endl << "--------------------------------------" << std::endl << std::endl;
    os << "List of Common Histograms in OnlMonServer"  << std::endl;
    for (auto &hiter : CommonHistoMap)
    {
      os << hiter.first << std::endl;
    }
    os << std::endl;
  }
  if (what == "ALL" || what == "MONITOR")
  {
    // loop over the map and print out the content (name and location in memory)
    os << "--------------------------------------" << std::endl << std::endl;
    os << "List of Monitors with registered histos in OnlMonServer:"  << std::endl;

    for (auto &miter : MonitorList)
    {
      os << miter->Name() << std::endl;
    }
    os << std::endl;
  }
  if (what == "ALL" || what == "ACTIVE")
  {
    os << "--------------------------------------" << std::endl << std::endl;
    os << "List of active packets:" << std::endl;
    std::set<unsigned int>::const_iterator iter;
    for (iter = activepackets.begin(); iter != activepackets.end(); ++iter)
    {
      os << *iter << std::endl;
    }
  }
return;
}

void OnlMonServer::PrintFile(const std::string &fname) const
{
  std::ofstream fout(fname);
  Print("ALL",fout);
  fout.close();
}

void OnlMonServer::RunNumber(const int irun)
{
  runnumber = irun;
  // recoConsts *rc = recoConsts::instance();
  // rc->set_IntFlag("RUNNUMBER", irun);
  return;
}

int OnlMonServer::WriteHistoFile()
{
  for (auto &moniiter : MonitorHistoSet)
  {
    std::string dirname = "./";
    if (getenv("ONLMON_SAVEDIR"))
    {
      dirname = std::string(getenv("ONLMON_SAVEDIR")) + "/";
    }
    std::string filename = dirname + "Run_" + std::to_string(RunNumber()) + "-" + moniiter.first + ".root";
    if (Verbosity() > 2)
    {
      std::cout << "saving histos for " << moniiter.first << " in " << filename << std::endl;
    }
    TFile *hfile = TFile::Open(filename.c_str(), "RECREATE", "Created by Online Monitor");
    for (auto &histiter : moniiter.second)
    {
      histiter.second->Write();
    }
    hfile->Close();
    delete hfile;
  }
  return 0;
}

int OnlMonServer::send_message(const OnlMon *Monitor, const int msgsource, const int severity, const std::string &err_message, const int msgtype) const
{
  int iret = -1;
  std::map<std::string, MessageSystem *>::const_iterator iter = MsgSystem.find(Monitor->Name());
  std::ostringstream msg("Run ");
  int irun = RunNumber();
  msg << "Run " << irun << " Event# " << eventnumber << ": " << err_message;
  if (iter != MsgSystem.end())
  {
    iret = iter->second->send_message(msgsource, severity, msg.str(), msgtype);
  }
  else
  {
    // in case there are messages from the ctor before the subsystem is registered
    // print them out from the server
    iter = MsgSystem.find(ThisName);
    iret = iter->second->send_message(msgsource, severity, msg.str(), 0);
  }
  if (severity > MSG_SEV_WARNING)
  {
    WriteLogFile(Monitor->Name(), err_message);
  }
  return iret;
}

int OnlMonServer::send_message(const int severity, const std::string &err_message, const int msgtype) const
{
  std::map<std::string, MessageSystem *>::const_iterator iter = MsgSystem.find(ThisName);
  int iret = iter->second->send_message(MSG_SOURCE_UNSPECIFIED, severity, err_message, msgtype);
  return iret;
}

int OnlMonServer::WriteLogFile(const std::string &name, const std::string &message) const
{
  int irun = RunNumber();
  // if no run is opened yet, I don't want to produce a -1 file which will
  // increase in size indefinitely
  if (irun <= 0)
  {
    return 0;
  }
  // no logifles for junk and calibration runs
  if (RunType == "JUNK" || RunType == "CALIBRATION")
  {
    return 0;
  }
  std::ostringstream logfilename, msg;
  const char *logdir = getenv("ONLMON_LOGDIR");
  if (logdir)
  {
    logfilename << logdir << "/";
  }
  logfilename << name << "_" << irun << ".log.gz";
  gzFile fout = gzopen(logfilename.str().c_str(), "a9");
  time_t curticks = CurrentTicks();
  std::string timestr = ctime(&curticks);
  std::string::size_type backslpos;
  // get rid of this damn end line of ctime
  if ((backslpos = timestr.find('\n')) != std::string::npos)
  {
    timestr.erase(backslpos);
  }
  msg << timestr
      << ", EventNo " << eventnumber
      << ": " << message;
  gzprintf(fout, "%s\n", msg.str().c_str());
  gzclose(fout);
  return 0;
}

int OnlMonServer::CacheRunDB(const int runnoinput)
{
  int runno = -1;
  if (runnoinput == 221)
  {
    runno = runnoinput;
  }
  else
  {
    runno = 221;
  }
  RunType = "PHYSICS";
  standalone = 0;
  cosmicrun = 0;
  borticks = 0;
  return 0;
  odbc::Connection *con = nullptr;
  odbc::Statement *query = nullptr;
  std::ostringstream cmd;
  int iret = 0;
  try
  {
    con = odbc::DriverManager::getConnection("daq", "phnxrc", "");
  }
  catch (odbc::SQLException &e)
  {
    printf(" Exception caught during DriverManager::getConnection, Message: %s\n", e.getMessage().c_str());
    return -1;
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
  cmd << "SELECT runtype,triggerconfig,brunixtime FROM RUN  WHERE RUNNUMBER = "
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
    RunType = rs->getString("runtype");
    TriggerConfig = rs->getString("triggerconfig");
    borticks = rs->getInt("brunixtime");
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
  else
  {
    iret = -1;
  }
  delete con;
  //  printf("CacheRunDB: runno: %d, RunType %s\n",runno,RunType.c_str());
  return iret;
}

int OnlMonServer::SetSubsystemStatus(OnlMon *Monitor, const int status)
{
  if (GetRunType() != "JUNK")
  {
    statusDB->UpdateStatus(Monitor->Name(), runnumber, status);
  }
  return 0;
}

int OnlMonServer::SetSubsystemRunStatus(OnlMon *Monitor, const int status)
{
  if (GetRunType() == "PHYSICS")
  {
    RunStatusDB->UpdateStatus(Monitor->Name(), runnumber, status);
  }
  return 0;
}

int OnlMonServer::LookAtMe(OnlMon *Monitor, const int level, const std::string &message)
{
  std::cout << "got a LookAtMe from " << Monitor->Name()
            << ", level: " << level
            << ", message: " << message
            << std::endl;
  return 0;
}
