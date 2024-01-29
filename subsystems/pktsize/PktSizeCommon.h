#ifndef __PKTSIZECOMMON_H__
#define __PKTSIZECOMMON_H__

#include <map>
#include <string>

namespace PktSizeCommon
{
  int fillgranules(std::map<std::string, std::pair<unsigned int, unsigned int> > &granulepacketlimits);
  int filldcmgroups(std::map<unsigned int, std::string> &dcmgroupmap);
  int fillfibergroups(std::map<unsigned int, std::string> &dcmgroupmap);
};  // namespace PktSizeCommon

#endif /* __PKTSIZECOMMON_H__ */
