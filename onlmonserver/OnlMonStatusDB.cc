#include "OnlMonStatusDB.h"

#include <odbc++/connection.h>
#include <odbc++/drivermanager.h>
#include <odbc++/resultset.h>
#include <odbc++/statement.h>  // for Statement
#include <odbc++/types.h>      // for SQLException

#include <iostream>
#include <sstream>

static odbc::Connection* con = nullptr;

OnlMonStatusDB::OnlMonStatusDB(const std::string& tablename)
  : table(tablename)
{
}

OnlMonStatusDB::~OnlMonStatusDB()
{
  delete con;
  con = nullptr;
}

int OnlMonStatusDB::CheckAndCreateTable()
{
  if (GetConnection())
  {
    return -1;
  }

  //Postgres version
  //std::cout << con->getMetaData()-> getDatabaseProductVersion() << std::endl;
  odbc::Statement* stmt = con->createStatement();
  std::ostringstream cmd;
  //  cmd << "SELECT COUNT(*) FROM " << name << " WHERE 1 = 2" ;
  cmd << "SELECT * FROM pg_tables where tablename = '" << table << "'";
#ifdef VERBOSE

  std::cout << cmd.str() << std::endl;
#endif

  odbc::ResultSet* rs = nullptr;
  try
  {
    rs = stmt->executeQuery(cmd.str());
  }
  catch (odbc::SQLException& e)
  {
    const std::string& message = e.getMessage();
    if (message.find("does not exist") == std::string::npos)
    {
      std::cout << "Exception caught" << std::endl;
      std::cout << "Message: " << e.getMessage() << std::endl;
    }
  }
  cmd.str("");
  if (!rs->next())
  {
    delete rs;
    cmd << "CREATE TABLE " << table << "(runnumber int, primary key(runnumber))";
    try
    {
      stmt->executeUpdate(cmd.str());
    }
    catch (odbc::SQLException& e)
    {
      std::cout << "caught exception Message: " << e.getMessage() << std::endl;
    }
  }
  return 0;
}

int OnlMonStatusDB::CheckAndCreateMonitor(const std::string& name)
{
  std::ostringstream cmd;
  if (CheckAndCreateTable())
  {
    std::cout << "Problem creating " << table << std::endl;
    return -1;
  }
  cmd << "SELECT * FROM " << table << " LIMIT 1";
  odbc::ResultSet* rs = nullptr;
  odbc::Statement* stmt = con->createStatement();
  try
  {
    rs = stmt->executeQuery(cmd.str());
  }
  catch (odbc::SQLException& e)
  {
    std::cout << __PRETTY_FUNCTION__ << "Exception caught" << std::endl;
    std::cout << "Message: " << e.getMessage() << std::endl;
    return -1;
  }
  try
  {
    rs->findColumn(name);
  }
  catch (odbc::SQLException& e)
  {
    const std::string& exceptionmessage = e.getMessage();
    if (exceptionmessage.find("not found in result set") != std::string::npos)
    {
      cmd.str("");
      cmd << "ALTER TABLE "
          << table
          << " ADD "
          << name
          << " int default 0";
      try
      {
        odbc::Statement* stmtup = con->createStatement();
        stmtup->executeUpdate(cmd.str());
      }
      catch (odbc::SQLException& e2)
      {
        std::cout << "Exception caught: " << e2.getMessage() << std::endl;
        return -1;
      }
    }
  }
  return 0;
}

int OnlMonStatusDB::FindAndInsertRunNum(const int runnumber)
{
  if (GetConnection() != 0)
  {
    std::cout << "problem" << std::endl;
    return -1;
  }
  if (findRunNumInDB(runnumber) == 0)
  {
    return 0;
  }
  odbc::Statement* statement = con->createStatement();
  std::ostringstream cmd;
  cmd << "INSERT INTO "
      << table
      << " (runnumber) VALUES ("
      << runnumber << ")";
  try
  {
    statement->executeUpdate(cmd.str());
  }
  catch (odbc::SQLException& e)
  {
    std::cout << e.getMessage() << std::endl;
    return -1;
  }
  return 0;
}

int OnlMonStatusDB::findRunNumInDB(const int runnumber)
{
  odbc::Statement* statement = nullptr;
  odbc::ResultSet* rs = nullptr;
  std::ostringstream cmd;
  cmd << "SELECT runnumber FROM "
      << table
      << " WHERE runnumber = "
      << runnumber;

  statement = con->createStatement();

  try
  {
    rs = statement->executeQuery(cmd.str());
  }
  catch (odbc::SQLException& e)
  {
    std::cout << "exception caught: " << e.getMessage() << std::endl;
    return -1;
  }

  if (rs->next())
  {
    try
    {
      rs->getInt("runnumber");
    }
    catch (odbc::SQLException& e)
    {
      std::cout << "exception caught: " << e.getMessage() << std::endl;
      return -1;
    }
  }
  else
  {
    return -1;
  }
  return 0;
}

int OnlMonStatusDB::UpdateStatus(const std::string& name, const int runnumber, const int status)
{
  if (CheckAndCreateMonitor(name))
  {
    std::cout << __PRETTY_FUNCTION__ << "Problem encountered, cannot do update" << std::endl;
    return -1;
  }
  if (FindAndInsertRunNum(runnumber) != 0)
  {
    std::cout << __PRETTY_FUNCTION__ << "Problem updating runnumber encountered, cannot do update" << std::endl;
    return -1;
  }

  std::ostringstream cmd;
  cmd << "Update "
      << table
      << " set " << name
      << " = " << status
      << " where runnumber = "
      << runnumber;
  odbc::Statement* stmtupd = nullptr;
  try
  {
    stmtupd = con->createStatement();
  }
  catch (odbc::SQLException& e)
  {
    std::cout << "Cannot create statement" << std::endl;
    std::cout << e.getMessage() << std::endl;
    return -1;
  }

  try
  {
    stmtupd->executeUpdate(cmd.str());
  }
  catch (odbc::SQLException& e)
  {
    std::cout << __PRETTY_FUNCTION__ << "Exception caught" << std::endl;
    std::cout << "Message: " << e.getMessage() << std::endl;
    return -1;
  }
  return 0;
}

int OnlMonStatusDB::GetConnection()
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
  return 0;
}
