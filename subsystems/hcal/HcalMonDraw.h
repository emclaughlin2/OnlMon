#ifndef HCAL_HCALMONDRAW_H
#define HCAL_HCALMONDRAW_H

#include <onlmon/OnlMonDraw.h>

#include <TH2.h>
#include <TStyle.h>

#include <string>  // for allocator, string

class TCanvas;
class TGraphErrors;
class TPad;

class HcalMonDraw : public OnlMonDraw
{
 public:
  explicit HcalMonDraw(const std::string& name);

  ~HcalMonDraw() override {}

  int Init() override;
  int Draw(const std::string& what = "ALL") override;
  int MakeHtml(const std::string& what = "ALL") override;
  int SavePlot(const std::string &what = "ALL", const std::string &type = "png") override;
  void HandleEvent(int, int, int, TObject*);
  

 private:
  int MakeCanvas(const std::string& name);
  int DrawFirst(const std::string& what = "ALL");
  int DrawSecond(const std::string& what = "ALL");
  int DrawThird(const std::string& what = "ALL");
  int DrawFourth(const std::string& what = "ALL");
  int FindHotTower(TPad* warn, TH2*);
  void DrawTowerAvg();
  void DrawHitMap();
  time_t getTime();
  
  // int DrawDeadServer(TPad *transparent);
  std::string prefix = "HCALMON";
  TCanvas* TC[9] = {nullptr};
  TPad* transparent[9] = {nullptr};
  TPad* Pad[18] = {nullptr};
  TPad* warning[18] = {nullptr};
  TGraphErrors* gr[2] = {nullptr};
  TStyle* hcalStyle = nullptr;
};

#endif /* HCAL_HCALMONDRAW_H */
