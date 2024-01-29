#ifndef ONLMONDATABASE_RUNDBODBC_H
#define ONLMONDATABASE_RUNDBODBC_H

#include <set>
#include <string>

class RunDBodbc
{
 public:
  RunDBodbc() {}
  virtual ~RunDBodbc() {}
  //  void Dump(const int nrows = 0) const;
  void identify() const;
  std::string RunType(const int runno = 0) const;
  void Verbosity(const int i) { verbosity = i; }
  int GetRunNumbers(std::set<int> &result, const std::string &type, const int nruns, const int lastrunexclusive) const;

 private:
  int verbosity = 0;
  std::string dbname = "daq";
  std::string dbowner = "phnxrc";
  std::string dbpasswd = "";
};

#endif
