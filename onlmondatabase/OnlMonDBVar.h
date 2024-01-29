#ifndef ONLMONDATABASE_ONLMONDBVAR_H
#define ONLMONDATABASE_ONLMONDBVAR_H

#include <array>

class OnlMonDBVar
{
 public:
  OnlMonDBVar();
  virtual ~OnlMonDBVar(){};

  int SetVar(const float rval[3]);
  void Print() const;
  int wasupdated() { return updated; }
  void resetupdated() { updated = 0; }
  float GetVar(const int i) { return val[i]; }

 protected:
  int updated = 0;
  std::array<float, 3> val;
};

#endif
