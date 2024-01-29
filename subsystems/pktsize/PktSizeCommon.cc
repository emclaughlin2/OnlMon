#include "PktSizeCommon.h"
#include <onlmon/OnlMonServer.h>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>

int PktSizeCommon::fillgranules(std::map<std::string, std::pair<unsigned int, unsigned int> > &granulepacketlimits)
{
  std::pair<unsigned int, unsigned int> pktlimits;
  pktlimits.first = 17000;
  pktlimits.second = 17999;
  granulepacketlimits["ACC"] = pktlimits;

  pktlimits.first = 1000;
  pktlimits.second = 1999;
  granulepacketlimits["BBC"] = pktlimits;

  pktlimits.first = 3000;
  pktlimits.second = 3999;
  granulepacketlimits["DCH"] = pktlimits;

  pktlimits.first = 8000;
  pktlimits.second = 8999;
  granulepacketlimits["EMC"] = pktlimits;

  pktlimits.first = 14200;
  pktlimits.second = 14201;
  granulepacketlimits["ERT"] = pktlimits;

  //   pktlimits.first = 16000;
  //   pktlimits.second = 16999;
  //   granulepacketlimits["FCAL"] = pktlimits;

  pktlimits.first = 25000;
  pktlimits.second = 25999;
  granulepacketlimits["FVTX"] = pktlimits;

  pktlimits.first = 14000;
  pktlimits.second = 14024;
  granulepacketlimits["GL1"] = pktlimits;

  pktlimits.first = 14050;
  pktlimits.second = 14199;
  granulepacketlimits["LVL2"] = pktlimits;

  pktlimits.first = 14777;
  pktlimits.second = 14777;
  granulepacketlimits["LVL1"] = pktlimits;

  pktlimits.first = 21100;
  pktlimits.second = 21110;
  granulepacketlimits["MPC"] = pktlimits;

  pktlimits.first = 12000;
  pktlimits.second = 12999;
  granulepacketlimits["MUID"] = pktlimits;

  pktlimits.first = 11000;
  pktlimits.second = 11999;
  granulepacketlimits["MUTR"] = pktlimits;

  //   this one is for testing if no packet exists
  //   pktlimits.first = 2000;
  //   pktlimits.second = 2999;
  //   granulepacketlimits["MVD"] = pktlimits;

  pktlimits.first = 4000;
  pktlimits.second = 4999;
  granulepacketlimits["PAD"] = pktlimits;

  pktlimits.first = 6000;
  pktlimits.second = 6999;
  granulepacketlimits["RICH"] = pktlimits;

  pktlimits.first = 19000;
  pktlimits.second = 19999;
  granulepacketlimits["RPC"] = pktlimits;

  pktlimits.first = 7000;
  pktlimits.second = 7099;
  granulepacketlimits["TOFE"] = pktlimits;

  pktlimits.first = 7100;
  pktlimits.second = 7999;
  granulepacketlimits["TOFW"] = pktlimits;

  pktlimits.first = 21300;
  pktlimits.second = 21360;
  granulepacketlimits["MPCEX"] = pktlimits;

  pktlimits.first = 24000;
  pktlimits.second = 24099;
  granulepacketlimits["VTXP"] = pktlimits;

  pktlimits.first = 24100;
  pktlimits.second = 24199;
  granulepacketlimits["VTXS"] = pktlimits;

  pktlimits.first = 13000;
  pktlimits.second = 13999;
  granulepacketlimits["ZDC"] = pktlimits;

  return 0;
}

int PktSizeCommon::filldcmgroups(std::map<unsigned int, std::string> &dcmgroupmap)
{
  //  OnlMonServer *se = OnlMonServer::instance();
  std::set<std::string> pcffiles;
  //  se->parse_granuleDef(pcffiles);
  std::set<std::string>::const_iterator piter;
  for (piter = pcffiles.begin(); piter != pcffiles.end(); ++piter)
  {
    std::ostringstream filenam;
    if (getenv("ONLINE_CONFIGURATION"))
    {
      filenam << getenv("ONLINE_CONFIGURATION") << "/rc/hw/";
    }
    filenam << *piter;
    std::ifstream infile;
    infile.open(filenam.str().c_str(), std::ifstream::in);
    if (!infile)
    {
      if (filenam.str().find("gl1test.pcf") != std::string::npos)
      {
        std::cout << "  Could not open " << filenam.str() << std::endl;
      }
      continue;
    }
    std::string FullLine;  // a complete line in the config file
    getline(infile, FullLine);
    std::string::size_type pos1;
    std::string::size_type pos2;
    std::string dcmgrp = "NONE";
    while (!infile.eof())

    {
      if (FullLine.find("//") == std::string::npos)
      {
        if ((pos1 = FullLine.find("DCMGROUP")) != std::string::npos && FullLine.find("level1dd") != std::string::npos)
        {
          FullLine.erase(0, pos1);  // erase all before DCMGROUP string
          pos2 = FullLine.find(',');
          dcmgrp = FullLine.substr(0, pos2);
          //                  std::cout << dcmgrp << std::endl;
        }
        if ((pos1 = FullLine.find("packetid")) != std::string::npos)
        {
          FullLine.erase(0, pos1);  // erase all before packetid string
          while ((pos1 = FullLine.find(':')) != std::string::npos)
          {
            pos2 = FullLine.find(',');
            // search the int between the ":" and the ","
            std::string packetidstr = FullLine.substr(pos1 + 1, pos2 - (pos1 + 1));
            std::istringstream line;
            line.str(packetidstr);
            unsigned int packetid;
            line >> packetid;
            if (packetid > 0)
            {
              if (dcmgrp != "NONE")
              {
                dcmgroupmap[packetid] = dcmgrp;
              }
              else
              {
                std::cout << "error assigning packet " << packetid << " to dcm group" << std::endl;
              }
            }
            // erase this entry from the line
            FullLine.erase(0, pos2 + 1);
          }
        }
      }
      getline(infile, FullLine);
    }
    infile.close();
  }
  return 0;
}

int PktSizeCommon::fillfibergroups(std::map<unsigned int, std::string> &fibergroupmap)
{
  //  OnlMonServer *se = OnlMonServer::instance();
  std::set<std::string> pcffiles;
  //  se->parse_granuleDef(pcffiles);
  std::set<std::string>::const_iterator piter;
  for (piter = pcffiles.begin(); piter != pcffiles.end(); ++piter)
  {
    std::ostringstream filenam;
    if (getenv("ONLINE_CONFIGURATION"))
    {
      filenam << getenv("ONLINE_CONFIGURATION") << "/rc/hw/";
    }
    filenam << *piter;
    std::ifstream infile;
    infile.open(filenam.str().c_str(), std::ifstream::in);
    if (!infile)
    {
      std::cout << "  Could not open " << filenam.str() << std::endl;
      continue;
    }
    std::string FullLine;  // a complete line in the config file
    getline(infile, FullLine);
    std::string::size_type pos1;
    std::string::size_type pos2;
    std::string dcmgrp = "NONE";
    std::string unitgrp = "NONE";
    while (!infile.eof())
    {
      if (FullLine.find("//") == std::string::npos)
      {
        if ((pos1 = FullLine.find("DCMGROUP")) != std::string::npos && FullLine.find("level1dd") != std::string::npos)
        {
          FullLine.erase(0, pos1);  // erase all before DCMGROUP string
          pos2 = FullLine.find(',');
          dcmgrp = FullLine.substr(0, pos2);
          //                  std::cout << dcmgrp << std::endl;
        }
        if ((pos1 = FullLine.find("unit")) != std::string::npos)
        {
          FullLine.erase(0, pos1);  // erase all before unit string
          pos1 = FullLine.find(':');
          FullLine.erase(0, pos1 + 1);  // erase all before : string
          pos2 = FullLine.find(',');
          unitgrp = FullLine.substr(0, pos2);
          //                  std::cout << unitgrp << std::endl;
        }
        if ((pos1 = FullLine.find("packetid")) != std::string::npos)
        {
          FullLine.erase(0, pos1);  // erase all before packetid string
          while ((pos1 = FullLine.find(':')) != std::string::npos)
          {
            pos2 = FullLine.find(',');
            // search the int between the ":" and the ","
            std::string packetidstr = FullLine.substr(pos1 + 1, pos2 - (pos1 + 1));
            std::istringstream line;
            line.str(packetidstr);
            unsigned int packetid;
            line >> packetid;
            if (packetid > 0)
            {
              if (dcmgrp != "NONE" && unitgrp != "NONE")
              {
                std::string fibergrp = dcmgrp + ":" + unitgrp;
                fibergroupmap[packetid] = fibergrp;
              }
              else
              {
                std::cout << "error assigning packet " << packetid << " to fiber group" << std::endl;
              }
            }
            // erase this entry from the line
            FullLine.erase(0, pos2 + 1);
          }
        }
      }
      getline(infile, FullLine);
    }
    infile.close();
  }
  return 0;
}
