#include "TrigInfo.h"

#include <string>

TrigInfo::TrigInfo()
{
  Reset();
}

void TrigInfo::SetInfo(const std::string &name,
                       const unsigned short ibit,
                       const unsigned short ibitmask,
                       const unsigned short iindex,
                       const unsigned int iscaledown)
{
  trigname = name;
  bit = ibit;
  bitmask = ibitmask;
  index = iindex;
  scaledown = iscaledown;
  return;
}

void TrigInfo::Reset()
{
  trigname = "NONE";
  bit = 0;
  bitmask = 0;
  index = 0xFFFF;
  scaledown = 0xFFFFFFFF;
  return;
}
