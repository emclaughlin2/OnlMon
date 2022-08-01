#ifndef MBDLL1_MBDLL1MONDRAW_H
#define MBDLL1_MBDLL1MONDRAW_H

#include <onlmon/OnlMonDraw.h>

#include <string>  // for allocator, string

class OnlMonDB;
class TCanvas;
class TGraphErrors;
class TPad;

class Mbdll1MonDraw : public OnlMonDraw
{
 public:
  Mbdll1MonDraw(const std::string &name = "MBDLL1MON"); // same name as server!

  ~Mbdll1MonDraw() override {}

  int Init() override;
  int Draw(const std::string &what = "ALL") override;
  int MakePS(const std::string &what = "ALL") override;
  int MakeHtml(const std::string &what = "ALL") override;

 protected:
  int MakeCanvas(const std::string &name);
  int DrawFirst(const std::string &what = "ALL");
  int DrawSecond(const std::string &what = "ALL");
  int DrawHistory(const std::string &what = "ALL");
  int TimeOffsetTicks = -1;
  TCanvas *TC[3] = {nullptr};
  TPad *transparent[3] = {nullptr};
  TPad *Pad[6] = {nullptr};
  TGraphErrors *gr[2] = {nullptr};
  OnlMonDB *dbvars = nullptr;
};

#endif /* MBDLL1_MBDLL1MONDRAW_H */
