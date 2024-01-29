#include "OnlMon.h"
#include "OnlMonServer.h"

#include <Event/msg_profile.h>

#include <TH1.h>

#include <cstdio>  // for printf
#include <iostream>
#include <sstream>

class Event;

OnlMon::OnlMon(const std::string &name)
  : OnlMonBase(name)
{
  if (name.find('_') != std::string::npos)
  {
    std::cout << "No underscore (_) in online monitoring server name " << name << " allowed" << std::endl;
    exit(1);
  }
  status = OnlMon::ACTIVE;
  return;
}

int OnlMon::process_event_common(Event *evt)
{
  int iret = process_event(evt);
  return iret;
}

int OnlMon::process_event(Event * /*evt*/)
{
  std::cout << "process_event(Event *evt) not implemented by daughter class" << std::endl;
  return -1;
}

int OnlMon::Reset()
{
  //  cout << "Reset() not implemented by daughter class" << endl;
  return -1;
}

void OnlMon::identify(std::ostream &out) const
{
  out << "identify() not implemented by daughter class" << std::endl;
  return;
}

int OnlMon::InitCommon(OnlMonServer *  se)
{
//  m_LocalFrameWorkVars = static_cast<TH1 *>(se->getCommonHisto("FrameWorkVars")->Clone());
  se->registerHisto(this,se->getCommonHisto("FrameWorkVars"));
  return 0;
}

int OnlMon::BeginRunCommon(const int /* runno */, OnlMonServer * /*se*/)
{
  SetStatus(OnlMon::OK);
  return 0;
}

void OnlMon::SetStatus(const int newstatus)
{
  if (newstatus != status)
  {
    OnlMonServer *se = OnlMonServer::instance();
    status = newstatus;
    se->SetSubsystemRunStatus(this, status);
  }
  return;
}

void OnlMon::SetMonitorServerId(unsigned int i)
{
  if (Name().find('_') != std::string::npos)
  {
    std::cout << "Monitor Server Id was already set " << Name() << std::endl;
    return;
  }
  m_MonitorServerId = i;
  Name(Name() + '_' + std::to_string(i));
  return;
}

