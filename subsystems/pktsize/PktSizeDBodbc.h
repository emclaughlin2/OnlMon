#ifndef __PKTSIZEDBODBC_H__
#define __PKTSIZEDBODBC_H__

#include <onlmon/OnlMonBase.h>

#include <ctime>
#include <map>
#include <string>
#include <vector>

class PktSizeDBodbc : public OnlMonBase
{
 public:
  PktSizeDBodbc(const std::string &name);
  virtual ~PktSizeDBodbc();
  int CheckAndCreateTable(const std::string &name, const std::map<unsigned int, unsigned int> &packetsize);
  int CheckAndAddColumns(const std::string &table, const std::map<unsigned int, unsigned int> &packetsize);
  int AddRow(const std::string &granulename, const int runnumber, const int nevnts, const std::map<unsigned int, unsigned int> &packetsize);
  int GetPacketContent(std::map<unsigned int, float> &packetsize, const int runnumber, const std::string &table);

 private:
  int GetConnection();
  std::string dbname;
  std::string dbowner;
  std::string dbpasswd;
  std::string tableprefix;
};

#endif
