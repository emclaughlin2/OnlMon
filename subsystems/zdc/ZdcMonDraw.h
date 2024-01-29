#ifndef ZDC_ZDCMONDRAW_H
#define ZDC_ZDCMONDRAW_H

#include <onlmon/OnlMonDraw.h>

#include <string>  // for allocator, string

class OnlMonDB;
class TCanvas;
class TGraphErrors;
class TPad;

class ZdcMonDraw : public OnlMonDraw
{
 public:
  ZdcMonDraw(const std::string &name);

  ~ZdcMonDraw() override {}

  int Init() override;
  int Draw(const std::string &what = "ALL") override;
  int MakeHtml(const std::string &what = "ALL") override;
  int SavePlot(const std::string &what = "ALL", const std::string &type = "png") override;

 protected:
  int MakeCanvas(const std::string &name);
  int DrawFirst(const std::string &what = "ALL");
  int DrawSecond(const std::string &what = "ALL");
  TCanvas *TC[2] = {nullptr};
  TPad *transparent[2] = {nullptr};
  TPad *Pad[4] = {nullptr};
};

#endif /* ZDC_ZDCMONDRAW_H */
