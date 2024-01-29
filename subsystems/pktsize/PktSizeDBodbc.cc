#include "PktSizeDBodbc.h"

#include <phool/phool.h>

#include <odbc++/connection.h>
#include <odbc++/drivermanager.h>
#include <odbc++/errorhandler.h>
#include <odbc++/resultset.h>
#include <odbc++/resultsetmetadata.h>
#include <odbc++/setup.h>
#include <odbc++/types.h>
#include <sql.h>
#include <iostream>

#include <algorithm>
#include <cctype>
#include <ctime>
#include <fstream>
#include <set>
#include <sstream>

static odbc::Connection* con = nullptr;

PktSizeDBodbc::PktSizeDBodbc(const std::string& name)
  : OnlMonBase(name)
  , tableprefix(name)
{
}

PktSizeDBodbc::~PktSizeDBodbc()
{
  delete con;
  con = nullptr;
}

int PktSizeDBodbc::CheckAndCreateTable(const std::string& name, const std::map<unsigned int, unsigned int>& packetsize)
{
  if (GetConnection())
  {
    return -1;
  }

  // Postgres version
  // std::cout << con->getMetaData()-> getDatabaseProductVersion() << std::endl;
  odbc::Statement* stmt = con->createStatement();
  std::ostringstream cmd;
  //  cmd << "SELECT COUNT(*) FROM " << name << " WHERE 1 = 2" ;
  // cmd << "select " << name << " from pg_tables where schemaname='public";
  std::string lowname = name;
  // The bizarre cast here is needed for newer gccs
  transform(lowname.begin(), lowname.end(), lowname.begin(), (int (*)(int)) tolower);
  cmd << "select * from pg_tables where tablename = '" << lowname << "'";
  //  cmd << "SELECT * FROM " << name << " LIMIT 1" ;
  if (verbosity > 0)
  {
    std::cout << "cmd: " << cmd.str() << std::endl;
  }

  odbc::ResultSet* rs = nullptr;
  try
  {
    rs = stmt->executeQuery(cmd.str());
  }
  catch (odbc::SQLException& e)
  {
    std::cout << name << " does not exist, creating it" << std::endl;
  }
  int iret = 0;

  cmd.str("");
  if (!rs->next())
  {
    delete rs;
    rs = nullptr;
    cmd << "CREATE TABLE " << name << "(runnumber int NOT NULL, events int NOT NULL";
    std::map<unsigned int, unsigned int>::const_iterator iter;
    for (iter = packetsize.begin(); iter != packetsize.end(); ++iter)
    {
      cmd << ", p_" << iter->first << " float DEFAULT 0";
    }
    cmd << ", primary key(runnumber))";
    if (verbosity > 0)
    {
      std::cout << "Executing " << cmd.str() << std::endl;
    }
    try
    {
      iret = stmt->executeUpdate(cmd.str());
    }
    catch (odbc::SQLException& e)
    {
      std::cout << "Exception caught, Message: " << e.getMessage() << std::endl;
    }
  }
  delete stmt;
  return iret;
}

int PktSizeDBodbc::AddRow(const std::string& granulename, const int runnumber, const int nevnts, const std::map<unsigned int, unsigned int>& packetsize)
{
  std::string table = tableprefix + granulename;

  CheckAndCreateTable(table, packetsize);

  int iret = 0;
  std::ostringstream cmd;

  if (GetConnection())
  {
    return -1;
  }

  odbc::Statement* stmt = con->createStatement();

  // check if an entry for this run exists already
  cmd << "SELECT events FROM " << table << " where runnumber = "
      << runnumber;
  odbc::ResultSet* rs;
  try
  {
    rs = stmt->executeQuery(cmd.str());
  }
  catch (odbc::SQLException& e)
  {
    std::cout << "Exception caught, Message: " << e.getMessage() << std::endl;
    return -1;
  }
  if (rs->next())
  {
    int events = rs->getInt("events");
    if (nevnts <= events)
    {
      std::cout << "Run " << runnumber << " already in table "
           << table << " extracted from " << events << " Events"
           << std::endl;
      std::cout << "Run more events than " << events
           << " if you want to overwrite this entry" << std::endl;
      delete rs;
      return 0;
    }
    else
    {
      cmd.str("");
      cmd << "DELETE FROM " << table << " WHERE runnumber = " << runnumber;
      odbc::Statement* stmt2 = con->createStatement();
      stmt2->executeUpdate(cmd.str());
    }
  }
  delete rs;
  // update columns in table in case there is a new packet
  CheckAndAddColumns(table, packetsize);

  // now add the content
  std::map<unsigned int, unsigned int>::const_iterator iter;
  cmd.str("");
  std::ostringstream cmd1, cmd2;
  cmd1 << "INSERT INTO " << table
       << " (runnumber, events";
  cmd2 << " VALUES(" << runnumber << ", " << nevnts;
  float size_in_bytes;
  for (iter = packetsize.begin(); iter != packetsize.end(); ++iter)
  {
    cmd1 << ", p_" << iter->first;
    size_in_bytes = (float) (iter->second) / (float) (nevnts);
    size_in_bytes *= 4;  // convert from 32 bit words to bytes
    cmd2 << ", " << size_in_bytes;
  }
  cmd1 << ")";
  cmd2 << ")";
  cmd << cmd1.str() << cmd2.str();
  try
  {
    stmt->executeUpdate(cmd.str());
  }
  catch (odbc::SQLException& e)
  {
    std::cout << "Exception caught, Message: " << e.getMessage() << std::endl;
  }

  return iret;
}

int PktSizeDBodbc::GetPacketContent(std::map<unsigned int, float>& packetsize, const int runnumber, const std::string& granulename)
{
  if (GetConnection())
  {
    return -1;
  }
  int iret = 0;
  std::string table = tableprefix + granulename;

  odbc::Statement* query = con->createStatement();
  std::ostringstream cmd;
  cmd << "SELECT * FROM " << table << " WHERE runnumber = " << runnumber;

  if (verbosity > 0)
  {
    std::cout << "command: " << cmd.str() << std::endl;
  }

  odbc::ResultSet* rs;
  try
  {
    rs = query->executeQuery(cmd.str());
  }
  catch (odbc::SQLException& e)
  {
    const std::string& errmsg = e.getMessage();
    if (errmsg.find("does not exist") == std::string::npos)
    {
      std::cout << "Exception caught, when accessing table "
           << table << std::endl;
      std::cout << "Message: " << e.getMessage() << std::endl;
    }
    return -1;
  }
  if (rs->next())
  {
    for (int i = 1; i <= rs->getMetaData()->getColumnCount(); i++)
    {
      std::string colname = rs->getMetaData()->getColumnName(i);
      if (colname == "runnumber" || colname == "events")
      {
        continue;
      }
      std::string packet = colname.substr(colname.find_last_of('_') + 1);
      std::istringstream istr(packet);
      unsigned int ipkt;
      istr >> ipkt;
      float size = rs->getFloat(colname);
      if (rs->wasNull() || size <= 1.)
      {
        continue;
      }
      packetsize[ipkt] = size / 4.;  // convert from bytes to long words
    }
    delete rs;
  }
  return iret;
}

int PktSizeDBodbc::GetConnection()
{
  if (con)
  {
    return 0;
  }
  try
  {
    con = odbc::DriverManager::getConnection(dbname.c_str(), dbowner.c_str(), dbpasswd.c_str());
  }
  catch (odbc::SQLException& e)
  {
    std::cout << PHWHERE
         << " Exception caught during DriverManager::getConnection" << std::endl;
    std::cout << "Message: " << e.getMessage() << std::endl;
    if (con)
    {
      delete con;
      con = nullptr;
    }
    return -1;
  }
  printf("opened DB connection\n");
  return 0;
}

int PktSizeDBodbc::CheckAndAddColumns(const std::string& table, const std::map<unsigned int, unsigned int>& packetsize)
{
  if (GetConnection())
  {
    return -1;
  }
  odbc::Statement* stmt = con->createStatement();
  std::ostringstream cmd;
  cmd << "SELECT * FROM " << table << " limit 1 ";
  odbc::ResultSet* rs;
  try
  {
    rs = stmt->executeQuery(cmd.str());
  }
  catch (odbc::SQLException& e)
  {
    std::cout << "Exception caught, Message: " << e.getMessage() << std::endl;
    return -1;
  }
  std::map<unsigned int, unsigned int>::const_iterator iter;
  std::set<unsigned int> packetids;
  for (iter = packetsize.begin(); iter != packetsize.end(); ++iter)
  {
    packetids.insert(iter->first);
  }
  for (int i = 1; i <= rs->getMetaData()->getColumnCount(); i++)
  {
    std::string colname = rs->getMetaData()->getColumnName(i);
    if (colname == "runnumber" || colname == "events")
    {
      continue;
    }
    std::string packet = colname.substr(colname.find_last_of('_') + 1);
    std::istringstream istr(packet);
    unsigned int ipkt;
    istr >> ipkt;
    packetids.erase(ipkt);
  }
  delete rs;
  std::set<unsigned int>::const_iterator siter;
  for (siter = packetids.begin(); siter != packetids.end(); ++siter)
  {
    std::ostringstream newcol;
    newcol.str("");
    newcol << "p_" << *siter;
    cmd.str("");
    cmd << "ALTER TABLE " << table << " ADD COLUMN " << newcol.str()
        << " float DEFAULT 0";
    try
    {
      stmt->executeUpdate(cmd.str());
    }
    catch (odbc::SQLException& e)
    {
      std::cout << "Exception caught, Message: " << e.getMessage() << std::endl;
    }
  }
  return 0;
}
