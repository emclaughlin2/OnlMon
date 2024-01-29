#ifndef SEPD_SEPDMONDRAW_H
#define SEPD_SEPDMONDRAW_H

#include <onlmon/OnlMonDraw.h>

#include <TH2.h>
#include <TStyle.h>

#include <string>  // for allocator, string

class TCanvas;
class TGraphErrors;
class TPad;

class SepdMonDraw : public OnlMonDraw
{
 public:
  SepdMonDraw(const std::string &name);

  ~SepdMonDraw() override {}

  int Init() override;
  int Draw(const std::string &what = "ALL") override;
  int MakeHtml(const std::string &what = "ALL") override;
  int SavePlot(const std::string &what = "ALL", const std::string &type = "png") override;

 protected:
  int MakeCanvas(const std::string &name);
  int DrawFirst(const std::string &what = "ALL");
  int DrawSecond(const std::string &what = "ALL");
  int DrawThird(const std::string &what = "ALL");
  int DrawFourth(const std::string &what = "ALL");
  int DrawFifth(const std::string &what = "ALL");
  time_t getTime();

  TCanvas* TC[9] = {nullptr};
  TPad* transparent[9] = {nullptr};
  TPad* Pad[18] = {nullptr};
  TPad* warning[18] = {nullptr};
  TGraphErrors *gr[2] = {nullptr};
};

#endif /* SEPD_SEPDMONDRAW_H */
