#include "OnlMonDBodbc.h"
#include "OnlMonDBReturnCodes.h"
#include "OnlMonDBVar.h"

#include <onlmon/OnlMonBase.h>  // for OnlMonBase

#include <odbc++/connection.h>
#include <odbc++/databasemetadata.h>
#include <odbc++/drivermanager.h>
#include <odbc++/resultset.h>  // for ResultSet
#include <odbc++/resultsetmetadata.h>
#include <odbc++/statement.h>  // for Statement
#include <odbc++/types.h>      // for SQLException, odbc

#include <ctype.h>  // for tolower
#include <algorithm>
#include <cstdio>  // for printf
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <utility>  // for pair

static const unsigned int DEFAULTCOLUMNS = 2;
//static const unsigned int COLUMNSPARVAR = 3;
static const std::string addvarname[3] = {"", "err", "qual"};
static const unsigned int MINUTESINTERVAL = 4;

static odbc::Connection* con = nullptr;

// #define VERBOSE

OnlMonDBodbc::OnlMonDBodbc(const std::string& tablename)
  : OnlMonBase(tablename)
  , dbname("OnlMonDB")
  , dbowner("phnxrc")
  , dbpasswd("")
  , table(tablename)
{
  // table names are lower case only, so convert string to lowercase
  // The bizarre cast here is needed for newer gccs
  transform(table.begin(), table.end(), table.begin(), (int (*)(int)) tolower);
}

OnlMonDBodbc::~OnlMonDBodbc()
{
  delete con;
  con = nullptr;
}

int OnlMonDBodbc::CreateTable()
{
  if (GetConnection())
  {
    return -1;
  }

  std::ostringstream cmd;
  cmd << "SELECT * FROM pg_tables where tablename = '" << table << "'";
  odbc::ResultSet* rs = nullptr;
  odbc::Statement* stmt = con->createStatement();
  try
  {
    rs = stmt->executeQuery(cmd.str());
  }
  catch (odbc::SQLException& e)
  {
    std::cout << "caught exception Message: " << e.getMessage() << std::endl;
  }
  int iret = 0;
  if (!rs->next())
  {
    delete rs;
    rs = nullptr;
    if (verbosity > 0)
    {
      std::cout << "need to create table" << std::endl;
    }
    cmd.str("");
    cmd << "CREATE TABLE " << table << "(date timestamp(0) with time zone NOT NULL, run int NOT NULL, primary key(date,run))";
    try
    {
      stmt->executeUpdate(cmd.str());
    }
    catch (odbc::SQLException& e)
    {
      std::cout << "caught exception Message: " << e.getMessage() << std::endl;
    }
  }
  else
  {
    if (verbosity > 0)
    {
      std::cout << "table " << table << " exists" << std::endl;
    }
  }
  delete rs;
  delete stmt;
  return iret;
}

int OnlMonDBodbc::CheckAndCreateTable(const std::map<const std::string, OnlMonDBVar*>& varmap)
{
  if (GetConnection())
  {
    return -1;
  }
  if (CreateTable())  // check and create the table
  {
    return -1;
  }
  odbc::Statement* stmt = con->createStatement();
  odbc::ResultSet* rs = nullptr;
  std::ostringstream cmd;
  int iret = 0;
  cmd << "select * from " << table << " limit 1";
  std::map<const std::string, OnlMonDBVar*>::const_iterator iter;
  try
  {
    rs = stmt->executeQuery(cmd.str());
  }
  catch (odbc::SQLException& e)
  {
    std::cout << "caught exception Message: " << e.getMessage() << std::endl;
  }
  odbc::ResultSetMetaData* meta = rs->getMetaData();
  unsigned int nocolumn = rs->getMetaData()->getColumnCount();
  for (iter = varmap.begin(); iter != varmap.end(); ++iter)
  {
    std::string varname = iter->first;
    // column names are lower case only, so convert string to lowercase
    // The bizarre cast here is needed for newer gccs
    transform(varname.begin(), varname.end(), varname.begin(), (int (*)(int)) tolower);
    for (const auto & j : addvarname)
    {
      std::string thisvar = varname + j;
      for (unsigned int i = DEFAULTCOLUMNS + 1; i <= nocolumn; i++)
      {
        if (meta->getColumnName(i) == thisvar)
        {
          if (verbosity > 0)
          {
            std::cout << thisvar << " is in table" << std::endl;
          }
          goto foundvar;
        }
      }
      cmd.str("");
      cmd << "ALTER TABLE " << table << " ADD COLUMN "
          << thisvar << " real";
      if (verbosity > 0)
      {
        std::cout << "executing sql command: " << cmd.str() << std::endl;
      }

      try
      {
        odbc::Statement* chgtable = con->createStatement();
        iret = chgtable->executeUpdate(cmd.str());
        delete chgtable;
      }
      catch (odbc::SQLException& e)
      {
        std::cout << "Table " << table << " update failed" << std::endl;
        std::cout << "Message: " << e.getMessage() << std::endl;
      }

    foundvar:
      continue;
    }
  }
  delete rs;
  delete stmt;
  // check columns
  return iret;
}

int OnlMonDBodbc::DropTable(const std::string& name)
{
  if (GetConnection())
  {
    return -1;
  }
  std::string tablename = name;
  // table names are lower case only, so convert string to lowercase
  // The bizarre cast here is needed for newer gccs
  transform(tablename.begin(), tablename.end(), tablename.begin(), (int (*)(int)) tolower);
  std::cout << con->getMetaData()->getDatabaseProductVersion() << std::endl;
  std::cout << con->getCatalog() << std::endl;
  odbc::Statement* stmt = con->createStatement();
  std::ostringstream cmd;
  cmd << "DROP TABLE " << tablename;

  int iret = stmt->executeUpdate(cmd.str());
  std::cout << "iret: " << iret << std::endl;
  return iret;
}

void OnlMonDBodbc::Dump()
{
  if (GetConnection())
  {
    return;
  }

  //std::cout << con->getMetaData()-> getDatabaseProductVersion() << std::endl;
  odbc::Statement* stmt = con->createStatement();
  std::ostringstream cmd;
  cmd << "SELECT * FROM " << table;
  std::cout << "Executing " << cmd.str() << std::endl;
  odbc::ResultSet* rs = stmt->executeQuery(cmd.str());
  Dump(rs);
  return;
}

void OnlMonDBodbc::Dump(odbc::ResultSet* rs) const
{
  odbc::ResultSetMetaData* meta = rs->getMetaData();
  unsigned int ncolumn = meta->getColumnCount();
  while (rs->next())
  {
    odbc::Timestamp id = rs->getTimestamp(1);
    int Name = rs->getInt(2);
    std::cout << "TimeStamp: " << (id.toString()) << " Run #= " << Name << std::endl;
    for (unsigned int i = DEFAULTCOLUMNS + 1; i <= ncolumn; i++)
    {
      float rval = rs->getFloat(i);
      if (rs->wasNull())
      {
        std::cout << meta->getColumnName(i) << ": NULL" << std::endl;
      }
      else
      {
        std::cout << meta->getColumnName(i) << ": " << rval << std::endl;
      }
    }
  }
  try
  {
    rs->first();
  }
  catch (odbc::SQLException& e)
  {
    std::cout << "Exception caught" << std::endl;
    std::cout << "Message: " << e.getMessage() << std::endl;
  }

  return;
}

int OnlMonDBodbc::Info(const char* name)
{
  if (GetConnection())
  {
    return -1;
  }

  std::cout << con->getMetaData()->getDatabaseProductVersion() << std::endl;
  // if no argument is given print out list of tables
  if (!name)
  {
    std::cout << "Driver: " << con->getMetaData()->getDriverName() << std::endl;
    std::cout << "User: " << con->getMetaData()->getUserName() << std::endl;
    std::string catalog = con->getMetaData()->getCatalogTerm();
    std::cout << "Catalog: " << catalog << std::endl;
    std::string schemapattern = con->getMetaData()->getSchemaTerm();
    std::cout << "Schema: " << schemapattern << std::endl;
    std::string tablenamepattern = con->getMetaData()->getTableTerm();
    std::cout << "Table: " << tablenamepattern << std::endl;
    std::vector<std::string> types;
    odbc::ResultSet* rs = con->getMetaData()->getTableTypes();
    while (rs->next())
    {
      std::cout << "1: " << rs->getString(1) << std::endl;
      std::cout << "2: " << rs->getString(2) << std::endl;
      std::cout << "3: " << rs->getString(3) << std::endl;
      std::cout << "4: " << rs->getString(4) << std::endl;
      std::cout << "5: " << rs->getString(5) << std::endl;
    }
    odbc::ResultSet* rs1 = con->getMetaData()->getTables(catalog, schemapattern, tablenamepattern, types);

    std::cout << "rs1: " << rs1 << std::endl;
    //        std::cout << rs->getMetaData()->getTableName(1) << std::endl;
    while (rs1->next())
    {
      std::cout << "Table: " << rs->getString("TABLE_CAT") << std::endl;
      std::cout << "Table: " << rs->getString(3) << std::endl;
    }
  }
  else
  {
    odbc::Statement* stmt = con->createStatement();
    std::string tablename = name;
    // table names are lower case only, so convert string to lowercase
    // The bizarre cast here is needed for newer gccs
    transform(tablename.begin(), tablename.end(), tablename.begin(), (int (*)(int)) tolower);
    std::ostringstream cmd;
    cmd << "select * from " << tablename;
    odbc::ResultSet* rs = stmt->executeQuery("select * from inttest");
    std::cout << rs->getMetaData()->getColumnCount() << std::endl;
    while (rs->next())
    {
      int id = rs->getInt(1);
      float Name = rs->getFloat(2);
      std::cout << "Row: " << rs->getRow() << ", id=" << id << "   Name= " << Name << std::endl;
    }
  }
  return 0;
}

void OnlMonDBodbc::identify() const
{
  std::cout << "DB Name: " << dbname << std::endl;
  std::cout << "DB Owner: " << dbowner << std::endl;
  std::cout << "DB Pwd: " << dbpasswd << std::endl;
  std::cout << "DB table: " << table << std::endl;
  return;
}

int OnlMonDBodbc::AddRow(const time_t ticks, const int runnumber, const std::map<const std::string, OnlMonDBVar*>& varmap)
{
  // bail out when restarted before a run was taken - run=0, ticks=0
  if (ticks == 0 || runnumber <= 0)
  {
    return -1;
  }
  std::map<const std::string, OnlMonDBVar*>::const_iterator iter;
  int iret = 0;
  int minutesinterval = MINUTESINTERVAL;
  std::ostringstream cmd, cmd1, datestream;
  odbc::Timestamp thistime(ticks);
  odbc::Timestamp mintime;
  odbc::Timestamp maxtime;

  if (GetConnection())
  {
    return -1;
  }

  odbc::Statement* query = con->createStatement();
searchagain:
  mintime.setTime(ticks - minutesinterval * 60);
  maxtime.setTime(ticks + minutesinterval * 60);
#ifdef VERBOSE

  std::cout << "mintime stp: " << mintime.toString() << std::endl;
  std::cout << "maxtime stp: " << maxtime.toString() << std::endl;
#endif

  cmd.str("");
  datestream.str("");
  datestream << "date > '" << mintime.toString()
             << "' and date < '" << maxtime.toString() << "'";
  cmd << "SELECT COUNT(*) FROM " << table << " WHERE run = "
      << runnumber << " and " << datestream.str();
#ifdef VERBOSE

  std::cout << "cmd: " << cmd.str() << std::endl;
#endif

  odbc::ResultSet* rs = nullptr;

  try
  {
    rs = query->executeQuery(cmd.str());
  }
  catch (odbc::SQLException& e)
  {
    std::cout << "Exception caught" << std::endl;
    std::cout << "Message: " << e.getMessage() << std::endl;
    std::cout << "sql cmd: " << cmd.str() << std::endl;
  }
  int haverow = 0;
  if (rs)
  {
    while (rs->next())
    {
      haverow = rs->getInt(1);
#ifdef VERBOSE

      std::cout << "found rows: " << haverow << std::endl;
#endif
    }
    delete rs;
  }
  if (haverow > 1)
  {
    minutesinterval = minutesinterval / 2;
    goto searchagain;
  }
  else if (haverow == 1)
  {
    cmd.str("");
    cmd << "SELECT * FROM " << table << " WHERE " << datestream.str();
#ifdef VERBOSE

    std::cout << "command: " << cmd.str() << std::endl;
#endif

    rs = query->executeQuery(cmd.str());
    // if the following works, the query returned one or more rows
    // in the given time interval
    while (rs->next())
    {
      cmd.str("");
      for (iter = varmap.begin(); iter != varmap.end(); ++iter)
      {
        if (iter->second->wasupdated())
        {
          std::string varqualname = iter->first + addvarname[2];
          float varqual = iter->second->GetVar(2);
          float sqlvarqual = rs->getFloat(varqualname);
          if (varqual > sqlvarqual || rs->wasNull())
          {
            odbc::Statement* upd = con->createStatement();
            // if there is only 1 row we do not need to worry which one to update
            for (unsigned int j = 0; j < 3; j++)
            {
              cmd.str("");
              cmd << "UPDATE " << table << " SET "
                  << iter->first << addvarname[j] << " = "
                  << iter->second->GetVar(j) << " WHERE "
                  << datestream.str();
#ifdef VERBOSE

              std::cout << "Command: " << cmd.str() << std::endl;
#endif

              int iret2 = upd->executeUpdate(cmd.str());
              if (!iret2)
              {
                std::cout << __PRETTY_FUNCTION__ << "Update failed please send mail to pinkenburg@bnl.gov"
                          << std::endl;
                std::cout << "And include the macro and the following info" << std::endl;
                std::cout << "TableName: " << table << std::endl;
                std::cout << "Variable: " << iter->first << addvarname[j] << std::endl;
                std::cout << "Value: " << iter->second->GetVar(j) << std::endl;
                std::cout << "TimeStamp: " << rs->getTimestamp(1).toString() << std::endl;
              }
            }
          }
        }
      }
    }
  }
  else
  {
    cmd.str("");
    cmd << "INSERT INTO " << table << "(date, run";
    cmd1 << "VALUES('" << thistime.toString() << "'," << runnumber;
    int newval = 0;
    for (iter = varmap.begin(); iter != varmap.end(); ++iter)
    {
      if (iter->second->wasupdated())
      {
        for (unsigned int j = 0; j < 3; j++)
        {
          cmd << ", " << iter->first << addvarname[j];
          cmd1 << ", " << iter->second->GetVar(j);
          newval++;
        }
      }
    }
    if (!newval)
    {
      printf("No updated values\n");
      return -1;
    }
    cmd << ") ";
    cmd1 << ")";
    cmd << cmd1.str();

#ifdef VERBOSE

    std::cout << cmd.str() << std::endl;
#endif

    odbc::Statement* stmt = con->createStatement();
    try
    {
      stmt->executeUpdate(cmd.str());
    }
    catch (odbc::SQLException& e)
    {
      const std::string& errmsg = e.getMessage();
      if (errmsg.find("Cannot insert a duplicate key into unique index") != std::string::npos)
      {
#ifdef VERBOSE
        std::cout << "Identical entry already in DB" << std::endl;
#endif
        iret = 0;
      }
      else
      {
        std::cout << __PRETTY_FUNCTION__ << " DB Error in execute stmt: " << e.getMessage() << std::endl;
        std::ofstream savesql("lostupdates.sql", std::ios_base::app);
        savesql << cmd.str() << std::endl;
        savesql.close();
        iret = -1;
      }
    }
  }
  return iret;
}

int OnlMonDBodbc::GetVar(const time_t begin, const time_t end, const std::string& varname, std::vector<time_t>& timestp, std::vector<int>& runnumber, std::vector<float>& var, std::vector<float>& varerr)
{
  if (GetConnection())
  {
    return DBNOCON;
  }
  int iret = 0;

  odbc::Statement* query = con->createStatement();
  odbc::Timestamp mintime(begin);
  odbc::Timestamp maxtime(end);
  std::string varnameerr = varname + addvarname[1];
  std::ostringstream cmd, datestream;
  datestream << "date > '" << mintime.toString()
             << "' and date < '" << maxtime.toString() << "'";
  cmd << "SELECT COUNT(*) FROM " << table << " WHERE " << datestream.str();

#ifdef VERBOSE
  std::cout << "command: " << cmd.str() << std::endl;
#endif

  odbc::ResultSet* rs;
  try
  {
    rs = query->executeQuery(cmd.str());
  }
  catch (odbc::SQLException& e)
  {
    std::cout << "Exception caught, probably your table "
              << table
              << " does not exist" << std::endl;
    std::cout << "Message: " << e.getMessage() << std::endl;
    varerr.resize(0);
    var.resize(0);
    timestp.resize(0);
    runnumber.resize(0);
    return -1;
  }
  int nrows = 0;
  while (rs->next())
  {
    nrows = rs->getInt(1);
  }
  timestp.resize(nrows);
  runnumber.resize(nrows);
  var.resize(nrows);
  varerr.resize(nrows);
  if (!nrows)
  {
    return DBNOENTRIES;
  }
  delete rs;
  cmd.str("");
  // get data in ascending date order (earliest time first)
  cmd << "SELECT date, run, " << varname << ", " << varnameerr
      << " FROM " << table << " WHERE "
      << datestream.str() << " ORDER BY date ASC";
#ifdef VERBOSE

  std::cout << "Command: " << cmd.str() << std::endl;
#endif

  try
  {
    rs = query->executeQuery(cmd.str());
  }
  catch (odbc::SQLException& e)
  {
    std::cout << "Exception caught, probably your variable "
              << varname << " or the table " << table
              << " does not exist" << std::endl;
    std::cout << "Message: " << e.getMessage() << std::endl;
    varerr.resize(0);
    var.resize(0);
    timestp.resize(0);
    runnumber.resize(0);
    return -3;
  }
  unsigned int index = 0;
  int filledrows = nrows;
  while (rs->next())
  {
    float val = rs->getFloat(3);
    if (rs->wasNull())
    {
      filledrows--;
      continue;
    }
    timestp[index] = rs->getTimestamp(1).getTime();
    runnumber[index] = rs->getInt(2);
    var[index] = val;
    varerr[index] = rs->getFloat(4);
    index++;
  }
  if (filledrows != nrows)
  {
    varerr.resize(filledrows);
    var.resize(filledrows);
    timestp.resize(filledrows);
    runnumber.resize(filledrows);
  }
  return iret;
}

int OnlMonDBodbc::GetConnection()
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
    std::cout << __PRETTY_FUNCTION__
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
