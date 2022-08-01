#ifndef EXAMPLE_MYMONDRAW_H
#define EXAMPLE_MYMONDRAW_H

#include <onlmon/OnlMonDraw.h>

#include <string>  // for allocator, string

class OnlMonDB;
class TCanvas;
class TGraphErrors;
class TPad;

class MyMonDraw : public OnlMonDraw
{
 public:
  MyMonDraw(const std::string &name = "MYMON"); // same name as server!
  ~MyMonDraw() override {}

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

#endif /*__MYMONDRAW_H__ */
