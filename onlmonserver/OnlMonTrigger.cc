#include "OnlMonTrigger.h"
#include "OnlMonServer.h"
#include "TrigInfo.h"

#include <TrigRunLvl1.h>
#include <TriggerHelper.h>

#include <msg_control.h>
#include <msg_profile.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

OnlMonTrigger::OnlMonTrigger()
  : MyRunNumber(-2)
  , MyRunNumberBit(-2)
  , MyRunNumberTriggerHelper(-2)
  , runnumber(-2)
  , th(nullptr)
{
  for (int i = 0; i < NGL1TRIG; i++)
  {
    triginfo[i] = new TrigInfo();
  }
  return;
}

OnlMonTrigger::~OnlMonTrigger()
{
  // deleting triggerhelper here gives invalid reads,
  // the dtor of triggerUtilities is called after this dtor
  // which takes care of the trigger helper object
  for (int i = 0; i < NGL1TRIG; i++)
  {
    delete triginfo[i];
  }
  return;
}

int OnlMonTrigger::SetupTriggerHelper()
{
  if (MyRunNumberTriggerHelper != runnumber)
  {
    MyRunNumberTriggerHelper = runnumber;
    th = tu.getTriggerHelper(runnumber);
    if (!th)
    {
      ostringstream msg;
      msg << "Could not find run in DB, loading ASCII file which is used by GL1";
      send_message(MSG_SEV_WARNING, msg.str());
      msg.str("");
      msg << "Check the l2config table in the Phenix database, run " << runnumber << " is missing from it";
      send_message(MSG_SEV_SEVEREERROR, msg.str());
    }
  }
  return 0;
}

void OnlMonTrigger::Print(const char *what) const
{
  if (!strcmp(what, "ALL") || !strcmp(what, "TRIGGERBITS"))
  {
    printf("--------------------------------------\n\n");
    printf("List of Trigger Bits in OnlMonServer:\n");

    map<const string, unsigned int>::const_iterator iter;
    for (iter = trigmask.begin(); iter != trigmask.end(); ++iter)
    {
      printf("%s is %08x\n", (iter->first).c_str(), iter->second);
    }
    printf("\n");
  }
  if (!strcmp(what, "ALL") || !strcmp(what, "TRIGGER"))
  {
    printf("--------------------------------------\n\n");
    printf("List of Triggers in OnlMonServer:\n");
    for (int i = 0; i < NGL1TRIG; i++)
    {
      printf("Index %d", triginfo[i]->Index());
      printf(", Bit %d", triginfo[i]->Bit());
      printf(", Name: %s", triginfo[i]->Name().c_str());
      printf(", BitMask: %08x", triginfo[i]->BitMask());
      printf(", Scale Down: %d\n", triginfo[i]->ScaleDown());
    }
    printf("\n");
  }
  return;
}

// this is for susbsystems who want to implement their own
// triggerbits selection
void OnlMonTrigger::TrigMask(const string &trigname, const unsigned int bitmask)
{
  map<const string, unsigned int>::const_iterator trigiter = trigmask.find(trigname);
  if (trigiter != trigmask.end())
  {
    ostringstream msg;
    msg << "Replacing Trigmask " << trigname << " of 0x" << hex
        << trigiter->second << " by 0x" << bitmask << dec;
    send_message(MSG_SEV_WARNING, msg.str());
  }
  trigmask[trigname] = bitmask;
  return;
}

unsigned int
OnlMonTrigger::getLevel1Bit(const string &name)
{
  int ibit;
  vector<string> trgnames;
  if (MyRunNumberBit != runnumber)
  {
    MyRunNumberBit = runnumber;
    trigmask.clear();
    SetupTriggerHelper();
  }
  map<const std::string, unsigned int>::const_iterator trigiter = trigmask.find(name);
  if (trigiter != trigmask.end())
  {
    return trigiter->second;
  }
  if (name == "ONLMONBBCLL1")
  {
    set<string>::const_iterator trgiter;
    for (trgiter = BbcLL1TrigNames.begin(); trgiter != BbcLL1TrigNames.end(); ++trgiter)
    {
      printf("Adding trigger: %s  to ONLMONBBCLL1 definition\n", (*trgiter).c_str());
      trgnames.push_back(*trgiter);
    }
  }
  else
  {
    trgnames.push_back(name);
  }
  vector<string>::const_iterator iter;
  unsigned int bitmask = 0x0;
  for (iter = trgnames.begin(); iter != trgnames.end(); ++iter)
  {
    unsigned int bittmp;
    if (th)
    {
      ibit = th->getLevel1BitNumber(*iter);
      if (ibit < 0)
      {
        ibit = 33;
      }
    }
    // if triggerhelper did not work (runnumber not in db),
    // read the ascii file used by GL1
    else
    {
      ibit = GetBitFromTriggerFile((*iter).c_str());
    }
    if (ibit < 32)
    {
      bittmp = 0x1 << ibit;
    }
    else
    {
      bittmp = 0;
    }
    bitmask += bittmp;
  }
  printf("using bitmask %08x for trigger %s\n", bitmask, name.c_str());
  trigmask[name] = bitmask;
  return bitmask;
}

unsigned int
OnlMonTrigger::GetBitFromTriggerFile(const char *trigname) const
{
  string filename = "/export/software/oncs/online_configuration/GL1/GL1TEMPDB_DIR/trigger_0.names";
  ifstream infile(filename.c_str());
  string FullLine;
  getline(infile, FullLine);
  int ibit = 0;
  while (!infile.eof())
  {
    if (FullLine == trigname)
    {
      return ibit;
    }
    ibit++;
    getline(infile, FullLine);
  }
  return 33;  // return impossible trigger bit (31 is maximum)
}

int OnlMonTrigger::FillTrigInfo()
{
  if (MyRunNumber != runnumber)
  {
    printf("Reloading Trig info for run %d\n", runnumber);
    ClearTrigInfo();
    SetupTriggerHelper();
    MyRunNumber = runnumber;
    if (th)
    {
      TrigRunLvl1 *trun = th->get_trigRunLvl1();
      if (trun)
      {
        for (int i = 0; i < 32; i++)
        {
          // if trigger not set this comes a NULL pointer
          // which causes exception when filling it into
          // an stl string
          if (trun->get_lvl1_trig_name(i))
          {
            triginfo[i]->SetInfo(trun->get_lvl1_trig_name(i),
                                 trun->get_lvl1_trig_bit(i),
                                 trun->get_lvl1_trig_bitmask(i),
                                 i,
                                 trun->get_lvl1_trig_scale_down(i));
            trigbits.insert(trun->get_lvl1_trig_bit(i));
          }
        }
      }
    }
    else
    {
      FillTrigInfoFromFile();
    }
  }
  return 0;
}

int OnlMonTrigger::FillTrigInfoFromFile()
{
  // trigger_0.names contains the trigger names. The position in that file
  // gives the trigger bit
  string filename = "/export/software/oncs/online_configuration/GL1/GL1TEMPDB_DIR/trigger_0.names";
  ifstream infiletrig(filename.c_str(), ios_base::in);
  // maskoff_0.msk gives the mask, position in the file is trigger bit
  filename = "/export/software/oncs/online_configuration/GL1/GL1TEMPDB_DIR/maskoff_0.msk";
  ifstream infilemask(filename.c_str(), ios_base::in);
  // scdn_0.sdn gives the scaledowns, position in the file is trigger bit
  filename = "/export/software/oncs/online_configuration/GL1/GL1TEMPDB_DIR/scdn_0.sdn";
  ifstream infilescal(filename.c_str(), ios_base::in);
  string FullLine1;
  string FullLine2;
  string FullLine3;
  getline(infiletrig, FullLine1);
  getline(infilemask, FullLine2);
  getline(infilescal, FullLine3);
  unsigned short ibit = 0;
  unsigned short index = 0;
  while (!infiletrig.eof() && !infilemask.eof() && !infilescal.eof())
  {
    istringstream line2(FullLine2.c_str());
    unsigned short imask;
    line2 >> imask;
    istringstream line3(FullLine3.c_str());
    unsigned int iscal;
    line3 >> iscal;
    triginfo[index]->SetInfo(FullLine1, ibit, imask, index, iscal);
    getline(infiletrig, FullLine1);
    getline(infilemask, FullLine2);
    getline(infilescal, FullLine3);
    ibit++;
    index++;
  }
  infiletrig.close();
  infilemask.close();
  infilescal.close();
  return 0;
}

int OnlMonTrigger::ClearTrigInfo()
{
  for (int i = 0; i < NGL1TRIG; i++)
  {
    triginfo[i]->Reset();
  }
  trigbits.clear();
  return 0;
}

int OnlMonTrigger::send_message(const int severity, const string &err_message) const
{
  // check $ONLINE_MAIN/include/msg_profile.h for MSG defs
  // if you do not find your subsystem, do not initialize it and drop me a line
  msg_control *Message = new msg_control(MSG_TYPE_MONITORING,
                                         MSG_SOURCE_DAQMON,
                                         severity, "OnlMonTrigger");
  cout << *Message << err_message << endl;
  delete Message;
  return 0;
}

unsigned int
OnlMonTrigger::get_lvl1_trig_scale_down(const unsigned int i)
{
  FillTrigInfo();
  return triginfo[i]->ScaleDown();
}

unsigned int
OnlMonTrigger::get_lvl1_trig_scale_down_bybit(const unsigned int ibit)
{
  FillTrigInfo();
  // that's the usual case, triggerindex and trigger bit are identical
  if (triginfo[ibit]->Bit() == ibit)
  {
    return triginfo[ibit]->ScaleDown();
  }
  else
  {
    if (trigbits.find(ibit) != trigbits.end())
    {
      for (short index = 0; index < NGL1TRIG; index++)
      {
        if (triginfo[index]->Bit() == ibit)
          return triginfo[index]->ScaleDown();
      }
      screwupmessage("get_lvl1_trig_scale_down_bybit", ibit);
    }
  }
  return 0xFFFFFFFF;
}

unsigned int
OnlMonTrigger::get_lvl1_trig_bitmask(const unsigned int index)
{
  FillTrigInfo();
  return triginfo[index]->BitMask();
}

unsigned int
OnlMonTrigger::get_lvl1_trig_bitmask_bybit(const unsigned int ibit)
{
  FillTrigInfo();
  // that's the usual case, triggerindex and trigger bit are identical
  if (triginfo[ibit]->Bit() == ibit)
  {
    return triginfo[ibit]->BitMask();
  }
  else
  {
    if (trigbits.find(ibit) != trigbits.end())
    {
      for (short index = 0; index < NGL1TRIG; index++)
      {
        if (triginfo[index]->Bit() == ibit)
          return triginfo[index]->BitMask();
      }
      screwupmessage("get_lvl1_trig_bitmask_bybit", ibit);
    }
  }
  return 0;
}

unsigned int
OnlMonTrigger::get_lvl1_trig_bit(const unsigned int i)
{
  FillTrigInfo();
  return triginfo[i]->Bit();
}

unsigned int
OnlMonTrigger::get_lvl1_trig_index(const unsigned int index)
{
  FillTrigInfo();
  return triginfo[index]->Index();
}

unsigned int
OnlMonTrigger::get_lvl1_trig_index_bybit(const unsigned int ibit)
{
  FillTrigInfo();
  // that's the usual case, triggerindex and trigger bit are identical
  if (triginfo[ibit]->Bit() == ibit)
  {
    return triginfo[ibit]->Index();
  }
  else
  {
    if (trigbits.find(ibit) != trigbits.end())
    {
      for (short index = 0; index < NGL1TRIG; index++)
      {
        if (triginfo[index]->Bit() == ibit)
        {
          return triginfo[index]->Index();
        }
      }
      screwupmessage("get_lvl1_trig_index_bybit", ibit);
    }
  }
  return 0xFFFFFFFF;
}

string
OnlMonTrigger::get_lvl1_trig_name(const unsigned int i)
{
  FillTrigInfo();
  return triginfo[i]->Name();
}

string
OnlMonTrigger::get_lvl1_trig_name_bybit(const unsigned int ibit)
{
  FillTrigInfo();
  // that's the usual case, triggerindex and trigger bit are identical
  if (triginfo[ibit]->Bit() == ibit)
  {
    return triginfo[ibit]->Name();
  }
  else
  {
    if (trigbits.find(ibit) != trigbits.end())
    {
      for (short index = 0; index < NGL1TRIG; index++)
      {
        if (triginfo[index]->Bit() == ibit)
        {
          return triginfo[index]->Name();
        }
      }
      screwupmessage("get_lvl1_trig_name_bybit", ibit);
    }
  }
  return "NONE";
}

int OnlMonTrigger::screwupmessage(const string &method, const unsigned int ibit) const
{
  ostringstream msg;
  msg << method
      << ": trigger screwup, bit mismatch for bit: "
      << ibit
      << ", tell Chris immediately and save the following info";
  send_message(MSG_SEV_FATAL, msg.str());
  Print("TRIGGER");
  return 0;
}

bool OnlMonTrigger::didLevel1TriggerFire(const string &name)
{
  unsigned int triggermask = getLevel1Bit(name);
  OnlMonServer *se = OnlMonServer::instance();
  static const unsigned short itrg = 1;  // use live trigger
  if (triggermask & se->Trigger(itrg))
  {
    return true;
  }
  return false;
}

void OnlMonTrigger::AddBbcLL1TrigName(const std::string &trgname)
{
  BbcLL1TrigNames.insert(trgname);
  return;
}
