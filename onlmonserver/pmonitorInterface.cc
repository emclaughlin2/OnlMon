/*
  This is the interface from the monitoring to pmonitor
*/

#include "pmonitorInterface.h"
#include "HistoBinDefs.h"
#include "OnlMon.h"
#include "OnlMonDefs.h"
#include "OnlMonServer.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <Event/Event.h>
#include <pmonitor.h>
#pragma GCC diagnostic pop

#include <Event/EventTypes.h>
#include <Event/msg_control.h>
#include <Event/msg_profile.h>
#include <Event/packet.h>

#include <MessageTypes.h>  // for kMESS_OBJECT, kMESS_STRING
#include <TClass.h>
#include <TH1.h>
#include <TInetAddress.h>  // for TInetAddress
#include <TMessage.h>
#include <TROOT.h>
#include <TServerSocket.h>
#include <TSocket.h>
#include <TSystem.h>
#include <TThread.h>

#include <pthread.h>
#include <sys/types.h>  // for time_t
#include <unistd.h>     // for sleep
#include <cstdio>       // for printf, NULL
#include <cstdlib>      // for exit
#include <cstring>      // for strcmp
#include <iostream>     // for operator<<, basic_ostream, endl, basic_o...
#include <sstream>
#include <string>

//#define ROOTTHREAD

#ifndef ROOTTHREAD
#define SERVER
#endif

#ifdef SERVER
static void *server(void *);
int ServerThread = 0;
#endif

#ifdef ROOTTHREAD
static void *server(void *);
static TThread *ServerThread = nullptr;
#endif

#ifdef USE_MUTEX
pthread_mutex_t mutex;
#endif

TH1 *FrameWorkVars = nullptr;

//*********************************************************************

int pinit()
{
  OnlMonServer *Onlmonserver = OnlMonServer::instance();
#ifdef USE_MUTEX
  Onlmonserver->GetMutex(mutex);
#endif
  for (int i = 0; i < kMAXSIGNALS; i++)
  {
    gSystem->IgnoreSignal((ESignals) i);
  }
#ifdef USE_MUTEX
  pthread_mutex_lock(&mutex);
#endif

#if defined(SERVER) || defined(ROOTTHREAD)

  pthread_t ThreadId = 0;
  if (!ServerThread)
  {
    // std::cout << "creating server thread" << std::endl;
#ifdef SERVER

    ServerThread = pthread_create(&ThreadId, nullptr, server, (void *) nullptr);
    Onlmonserver->SetThreadId(ThreadId);
#endif
#ifdef ROOTTHREAD

    ServerThread = new TThread(server, (void *) 0);
    ServerThread->Run();
#endif
  }
#endif
  // for the timestamp we need doubles
  FrameWorkVars = new TH1D("FrameWorkVars", "FrameWorkVars", NFRAMEWORKBINS, 0., NFRAMEWORKBINS);
  Onlmonserver->registerCommonHisto(FrameWorkVars);
#ifdef USE_MUTEX
  pthread_mutex_unlock(&mutex);
#endif

  return 0;
}

//*********************************************************************
int process_event(Event *evt)
{
  static uint64_t savetmpticks = 0x7FFFFFFF;
  static uint64_t borticks = 0;
  static uint64_t eorticks = 0;
  static int eventcnt = 0;

  OnlMonServer *se = OnlMonServer::instance();
  uint64_t tmpticks = evt->getTime();

  // first test if a new run has started and call BOR/EOR methods of monitors
  if (se->RunNumber() == -1)
  {
    savetmpticks = 0x7FFFFFFF;
    eventcnt = 0;
    int newrun = evt->getRunNumber();
#ifdef USE_MUTEX
    pthread_mutex_lock(&mutex);
#endif
    FrameWorkVars->SetBinContent(RUNNUMBERBIN, (Stat_t) newrun);
    se->BadEvents(0);
    se->EventNumber(evt->getEvtSequence());
    se->RunNumber(newrun);
    // set ticks to the current event, so the begin run has the
    // up to date time stamp, while the end run has the time stamp from
    // the last event of the previous run
    se->CurrentTicks(tmpticks);
    se->BeginRun(newrun);
    // set trigger mask in et pool frontend
    borticks = se->BorTicks();
    FrameWorkVars->SetBinContent(BORTIMEBIN, (Stat_t) borticks);
#ifdef USE_MUTEX
    pthread_mutex_unlock(&mutex);
#endif
    eorticks = borticks;
  }
  if (evt->getEvtLength() <= 0 || evt->getEvtLength() > 2500000)
  {
    std::ostringstream msg;
    msg << __PRETTY_FUNCTION__ << "Discarding event with length "
        << evt->getEvtLength();
    send_message(MSG_SEV_WARNING, msg.str());
    se->AddBadEvent();
    return 0;
  }

  int oldrun;
  if ((oldrun = se->RunNumber()) != evt->getRunNumber())
  {
    // ROOT crashes when one thread updates histos while they are
    // being saved, need mutex protection here
#ifdef USE_MUTEX
    pthread_mutex_lock(&mutex);
#endif
    FrameWorkVars->SetBinContent(EORTIMEBIN, (Stat_t) eorticks);  // set EOR time
    se->EndRun(oldrun);
    se->WriteHistoFile();
    se->Reset();  // reset all monitors
    int newrun = evt->getRunNumber();
    FrameWorkVars->SetBinContent(RUNNUMBERBIN, (Stat_t) newrun);
    eorticks = tmpticks;  // initialize eorticks
    se->BadEvents(0);
    se->RunNumber(newrun);
    se->EventNumber(evt->getEvtSequence());
    // set ticks to the current event, so the begin run has the
    // up to date time stamp, while the end run has the time stamp from
    // the last event of the previous run
    se->CurrentTicks(tmpticks);
    se->BeginRun(newrun);
    borticks = se->BorTicks();
    FrameWorkVars->SetBinContent(BORTIMEBIN, (Stat_t) borticks);
    eventcnt = 0;
#ifdef USE_MUTEX
    pthread_mutex_unlock(&mutex);
#endif
    savetmpticks = 0x7FFFFFFF;
  }

  se->CurrentTicks(tmpticks);
  // check if we get an event which was earlier than the BOR timestamp
  // save earliest time stamp and number of events with earlier timestamps
  if (tmpticks < borticks)
  {
#ifdef USE_MUTEX
    pthread_mutex_lock(&mutex);
#endif
    FrameWorkVars->AddBinContent(EARLYEVENTNUMBIN);
    if (tmpticks < savetmpticks)
    {
      savetmpticks = tmpticks;
      FrameWorkVars->SetBinContent(EARLYEVENTTIMEBIN, (Stat_t) tmpticks);
    }
#ifdef USE_MUTEX
    pthread_mutex_unlock(&mutex);
#endif
  }
  if (eorticks < se->CurrentTicks())
  {
    eorticks = se->CurrentTicks();
  }
  if (evt->getErrorCode())
  {
    std::ostringstream msg;
    msg << __PRETTY_FUNCTION__ << " Event with error code: "
        << evt->getErrorCode()
        << " discarding event " << evt->getEvtSequence();
    send_message(MSG_SEV_WARNING, msg.str());
    se->AddBadEvent();
    return 0;
  }

  eventcnt++;
#ifdef USE_MUTEX
  pthread_mutex_lock(&mutex);
#endif
  FrameWorkVars->SetBinContent(CURRENTTIMEBIN, (Stat_t) se->CurrentTicks());
  se->EventNumber(evt->getEvtSequence());
  se->process_event(evt);
#ifdef USE_MUTEX
  pthread_mutex_unlock(&mutex);
#endif
  return 0;
}

int setup_server()
{
  return 0;
}

static void *server(void * /* arg */)
{
  OnlMonServer *Onlmonserver = OnlMonServer::instance();
  int MoniPort = OnlMonDefs::MONIPORT;
  //  int thread_arg[5];
#ifdef USE_MUTEX
  pthread_mutex_lock(&mutex);
#endif
  TServerSocket *ss = nullptr;
  sleep(5);
  do
  {
    if (ss)
    {
      delete ss;
    }
    ss = new TServerSocket(MoniPort, kTRUE);
    // Accept a connection and return a full-duplex communication socket.
    Onlmonserver->PortNumber(MoniPort);
    if ((MoniPort - OnlMonDefs::MONIPORT) >= OnlMonDefs::NUMMONIPORT)
    {
      std::ostringstream msg;
      msg << "Too many Online Monitors running on this machine, bailing out";
      send_message(MSG_SEV_FATAL, msg.str());

      exit(1);
    }
    MoniPort++;
    if (!ss->IsValid())
    {
      printf("Ignore ROOT error about socket in use, I try another one\n");
    }
  } while (!ss->IsValid());  // from do{}while

  // root keeps a list of sockets and tries to close them when quitting.
  // this interferes with my own threading and makes valgrind crash
  // The solution is to remove the TServerSocket *ss from roots list of
  // sockets. Then it will leave this socket alone.
  int isock = gROOT->GetListOfSockets()->IndexOf(ss);
  gROOT->GetListOfSockets()->RemoveAt(isock);
  sleep(10);
#ifdef USE_MUTEX
  pthread_mutex_unlock(&mutex);
#endif
again:
  TSocket *s0 = ss->Accept();
  if (!s0)
  {
    std::cout << "Server socket " << OnlMonDefs::MONIPORT
              << " in use, either go to a different node or" << std::endl
              << "change MONIPORT in server/OnlMonDefs.h and recompile" << std::endl
              << "server and client" << std::endl;
    exit(1);
  }
  // mutex protected since writing of histo
  // to outgoing buffer and updating by other thread do not
  // go well together
  if (Onlmonserver->Verbosity() > 2)
  {
    TInetAddress adr = s0->GetInetAddress();
    std::cout << "got connection from " << std::endl;
    adr.Print();
  }
  //  std::cout << "try locking mutex" << std::endl;
#ifdef USE_MUTEX
  pthread_mutex_lock(&mutex);
#endif
  // std::cout << "got mutex" << std::endl;
  handleconnection(s0);
  // std::cout << "try releasing mutex" << std::endl;
#ifdef USE_MUTEX
  pthread_mutex_unlock(&mutex);
#endif
  // std::cout << "mutex released" << std::endl;
  delete s0;
  /*
    if (!aargh)
    {
    std::cout << "making thread" << std::endl;
    aargh = new TThread(handletest,(void *)0);
    aargh->Run();
    }
  */
  // std::cout << "closing socket" << std::endl;
  // s0->Close();
  goto again;
}

void handletest(void * /* arg */)
{
  //  std::cout << "threading" << std::endl;
  return;
}

void handleconnection(void *arg)
{
  TSocket *s0 = (TSocket *) arg;

  OnlMonServer *Onlmonserver = OnlMonServer::instance();
  /*
    int val;
    s0->GetOption(kSendBuffer, val);
    printf("sendbuffer size: %d\n", val);
    s0->GetOption(kRecvBuffer, val);
    printf("recvbuffer size: %d\n", val);
  */
  TMessage *mess = nullptr;
  TMessage outgoing(kMESS_OBJECT);
  while (true)
  {
    if (Onlmonserver->Verbosity() > 2)
    {
      std::cout << "Waiting for message" << std::endl;
    }
    s0->Recv(mess);
    if (!mess)
    {
      std::cout << "Broken Connection, closing socket" << std::endl;
      break;
    }
    if (mess->What() == kMESS_STRING)
    {
      char strchr[OnlMonDefs::MSGLEN];
      mess->ReadString(strchr, OnlMonDefs::MSGLEN);
      delete mess;
      mess = nullptr;
      std::string str = strchr;
      if (Onlmonserver->Verbosity() > 2)
      {
        std::cout << "received message: " << str << std::endl;
      }
      if (str == "Finished")
      {
        break;
      }
      else if (str == "WriteRootFile")
      {
        Onlmonserver->WriteHistoFile();
        s0->Send("Finished");
        break;
      }
      else if (str == "Ack")
      {
        continue;
      }
      else if (str == "HistoList")
      {
        if (Onlmonserver->Verbosity() > 2)
        {
          std::cout << "number of histos: " << Onlmonserver->nHistos() << std::endl;
        }
        for (auto monitors = Onlmonserver->monibegin(); monitors != Onlmonserver->moniend(); ++monitors)
        {
          for (auto &histos : monitors->second)
          {
            std::string subsyshisto = monitors->first + ' ' + histos.first;
            if (Onlmonserver->Verbosity() > 2)
            {
              std::cout << "subsystem: " << monitors->first << ", histo: " << histos.first << std::endl;
              std::cout << " sending: \"" << subsyshisto << "\"" << std::endl;
            }
            s0->Send(subsyshisto.c_str());
            int nbytes = s0->Recv(mess);
            delete mess;
            mess = nullptr;
            if (nbytes <= 0)
            {
              std::ostringstream msg;

              msg << "Problem receiving message: return code: " << nbytes;
              send_message(MSG_SEV_ERROR, msg.str());
            }
          }
        }
        s0->Send("Finished");
      }
      else if (str == "ALL")
      {
        if (Onlmonserver->Verbosity() > 2)
        {
          std::cout << "number of histos: " << Onlmonserver->nHistos() << std::endl;
        }
        for (unsigned int i = 0; i < Onlmonserver->nHistos(); i++)
        {
          TH1 *histo = Onlmonserver->getHisto(i);
          if (histo)
          {
            outgoing.Reset();
            outgoing.WriteObject(histo);
            s0->Send(outgoing);
            outgoing.Reset();
            s0->Recv(mess);
            delete mess;
            mess = nullptr;
          }
        }
        s0->Send("Finished");
      }
      else if (str.find("ISRUNNING") != std::string::npos)
      {
        std::string answer = "No";
        unsigned int pos_space = str.find(' ');
        std::string moniname = str.substr(pos_space + 1, str.size());
        for (auto moniter = Onlmonserver->monitor_vec_begin(); moniter != Onlmonserver->monitor_vec_end(); ++moniter)
        {
          if ((*moniter)->Name() == moniname)
          {
            answer = "Yes";
            break;
          }
        }
        if (Onlmonserver->Verbosity() > 2)
        {
          std::cout << "got " << str << ", replied " << answer << std::endl;
        }
        s0->Send(answer.c_str());
      }
      else if (str == "LISTMONITORS")
      {
        s0->Send("go");
        for (auto moniter = Onlmonserver->monitor_vec_begin(); moniter != Onlmonserver->monitor_vec_end(); ++moniter)
        {
          if (Onlmonserver->Verbosity() > 2)
          {
            std::cout << "sending " << (*moniter)->Name().c_str() << std::endl;
          }
          s0->Send((*moniter)->Name().c_str());
        }
        s0->Send("Finished");
        break;
      }
      else if (str == "LIST")
      {
        s0->Send("go");
        while (true)
        {
          char strmess[OnlMonDefs::MSGLEN];
          s0->Recv(mess);
          if (!mess)
          {
            break;
          }
          if (mess->What() == kMESS_STRING)
          {
            mess->ReadString(strmess, OnlMonDefs::MSGLEN);
            delete mess;
            mess = nullptr;
            if (std::string(strmess) == "alldone")
            {
              break;
            }
          }
          std::string str1(strmess);
          unsigned int pos_space = str1.find(' ');
          if (Onlmonserver->Verbosity() > 2)
          {
            std::cout << __PRETTY_FUNCTION__ << " getting subsystem " << str1.substr(0, pos_space) << ", histo " << str1.substr(pos_space + 1, str1.size()) << std::endl;
          }
          TH1 *histo = Onlmonserver->getHisto(str1.substr(0, pos_space), str1.substr(pos_space + 1, str1.size()));
          if (histo)
          {
            outgoing.Reset();
            outgoing.WriteObject(histo);
            s0->Send(outgoing);
            outgoing.Reset();
          }
          else
          {
            s0->Send("UnknownHisto");
          }
          //		  delete mess;
        }
        s0->Send("Finished");
      }
      else
      {
        std::string strstr(str);
        unsigned int pos_space = str.find(' ');
        TH1 *histo = Onlmonserver->getHisto(strstr.substr(0, pos_space), strstr.substr(pos_space + 1, str.size()));
        if (histo)
        {
          //		  const char *hisname = histo->GetName();
          outgoing.Reset();
          outgoing.WriteObject(histo);
          s0->Send(outgoing);
          outgoing.Reset();
          s0->Recv(mess);
          delete mess;
          s0->Send("Finished");
        }
        else
        {
          s0->Send("UnknownHisto");
        }
      }
    }
    else if (mess->What() == kMESS_OBJECT)
    {
      printf("got object of class: %s\n", mess->GetClass()->GetName());
      delete mess;
    }
    else
    {
      printf("*** Unexpected message ***\n");
      delete mess;
    }
  }

  // Close the socket.
  s0->Close();
  return;
}

int send_message(const int severity, const std::string &msg)
{
  // check $ONLINE_MAIN/include/msg_profile.h for MSG defs
  // if you do not find your subsystem, do not initialize it and drop me a line
  msg_control *Message = new msg_control(MSG_TYPE_MONITORING,
                                         MSG_SOURCE_DAQMON,
                                         severity, "pmonitorInterface");
  std::cout << *Message << msg << std::endl;
  delete Message;
  return 0;
}
