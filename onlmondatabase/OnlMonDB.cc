#include "OnlMonDB.h"
#include "OnlMonDBVar.h"
#include "OnlMonDBodbc.h"

#include <onlmon/OnlMonBase.h>  // for OnlMonBase
#include <onlmon/OnlMonServer.h>

#include <ctype.h>  // for tolower
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <utility>  // for pair

OnlMonDB::OnlMonDB(const std::string &thisname)
  : OnlMonBase(thisname)
{
  return;
}

OnlMonDB::~OnlMonDB()
{
  delete db;
  while (varmap.begin() != varmap.end())
  {
    delete varmap.begin()->second;
    varmap.erase(varmap.begin());
  }
}

void OnlMonDB::Print() const
{
  std::cout << "OnlMonDB Name: " << ThisName << std::endl;
  std::map<const std::string, OnlMonDBVar *>::const_iterator iter;
  for (iter = varmap.begin(); iter != varmap.end(); ++iter)
  {
    iter->second->Print();
  }
  return;
}

int OnlMonDB::registerVar(const std::string &varname)
{
  std::string cpstring = varname;
  transform(cpstring.begin(), cpstring.end(), cpstring.begin(), (int (*)(int)) tolower);
  std::map<const std::string, OnlMonDBVar *>::const_iterator iter = varmap.find(cpstring);
  if (iter != varmap.end())
  {
    std::cout << "Variable " << varname << " allready registered in DB" << std::endl;
    return -1;
  }
  varmap[cpstring] = new OnlMonDBVar();
  return 0;
}

int OnlMonDB::SetVar(const std::string &varname, const float var, const float varerr, const float varqual)
{
  float vararray[3];
  vararray[0] = var;
  vararray[1] = varerr;
  vararray[2] = varqual;
  return SetVar(varname, vararray);
}

int OnlMonDB::SetVar(const std::string &varname, const float var[3])
{
  std::string cpstring = varname;
  transform(cpstring.begin(), cpstring.end(), cpstring.begin(), (int (*)(int)) tolower);
  std::map<const std::string, OnlMonDBVar *>::iterator iter = varmap.find(cpstring);
  if (iter != varmap.end())
  {
    iter->second->SetVar(var);
    return 0;
  }
  std::cout << __PRETTY_FUNCTION__ << " Could not find Variable " << varname << " in DB list" << std::endl;
  return -1;
}

int OnlMonDB::DBInit()
{
  if (!db)
  {
    db = new OnlMonDBodbc(ThisName);
  }
  db->Verbosity(Verbosity());
  db->CheckAndCreateTable(varmap);
  return 0;
}

int OnlMonDB::DBcommit()
{
  OnlMonServer *se = OnlMonServer::instance();

  int iret = 0;
  if (!db)
  {
    std::cout << "Data Base not initialized, fix your code." << std::endl;
    std::cout << "You need to call DBInit() after you registered your variables" << std::endl;
    return -1;
  }
  iret = db->AddRow(se->CurrentTicks(), se->RunNumber(), varmap);
  if (iret)
  {
    printf("error in dbcommit, ret code %d\n", iret);
    return iret;
  }
  //db->Dump();
  std::map<const std::string, OnlMonDBVar *>::iterator iter;
  for (iter = varmap.begin(); iter != varmap.end(); ++iter)
  {
    iter->second->resetupdated();
  }
  return iret;
}

int OnlMonDB::DBcommitTest()
{
  static int ifirst = 1;
  int iret = 0;
  static time_t ticks = 0;
  static int runnumber = 90000;
  if (!db)
  {
    std::cout << "Data Base not initialized, fix your code." << std::endl;
    std::cout << "You need to call DBInit() after you registered your variables" << std::endl;
    return -1;
  }
  if (ifirst)
  {
    ticks = time(nullptr);
    ticks -= 2 * 30 * 24 * 60 * 60;
    ifirst = 0;
  }
  else
  {
    ticks += 3 * 60;
    runnumber += 1;
  }
  db->AddRow(ticks, runnumber, varmap);
  if (iret)
  {
    std::cout << "error in dbcommit" << std::endl;
    return iret;
  }
  //db->Dump();
  std::map<const std::string, OnlMonDBVar *>::iterator iter;
  for (iter = varmap.begin(); iter != varmap.end(); ++iter)
  {
    iter->second->resetupdated();
  }
  return iret;
}

int OnlMonDB::GetVar(const time_t begin, const time_t end, const std::string &varname, std::vector<time_t> &timestp, std::vector<int> &runnumber, std::vector<float> &var, std::vector<float> &varerr)
{
  if (!db)
  {
    db = new OnlMonDBodbc(ThisName);
  }
  int iret = db->GetVar(begin, end, varname, timestp, runnumber, var, varerr);
  return iret;
}

void OnlMonDB::Reset()
{
  std::map<const std::string, OnlMonDBVar *>::iterator iter;
  for (iter = varmap.begin(); iter != varmap.end(); ++iter)
  {
    iter->second->resetupdated();
  }
  return;
}
