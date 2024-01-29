#ifndef CEMC_CEMCMONDRAW_H
#define CEMC_CEMCMONDRAW_H

#include <onlmon/OnlMonDraw.h>

#include <string>  // for allocator, string
#include <TH2.h>
#include <TStyle.h>

class TCanvas;
class TGraphErrors;
class TPad;

class CemcMonDraw : public OnlMonDraw
{
 public:
  CemcMonDraw(const std::string &name = "CEMCMON"); // same name as server!

  ~CemcMonDraw() override {}

  int Init() override;
  int Draw(const std::string &what = "ALL") override;
  int MakeHtml(const std::string &what = "ALL") override;
  int SavePlot(const std::string &what = "ALL", const std::string &type = "png") override;
  void setSave(int s) {save = s;}

 private:
  int MakeCanvas(const std::string &name);
  int DrawFirst(const std::string &what = "ALL");
  int DrawSecond(const std::string &what = "ALL");
  int DrawThird(const std::string &what = "ALL");
  int DrawFourth(const std::string &what = "ALL");
  int DrawHistory(const std::string &what = "ALL");
  int FindHotTower(TPad *warn,TH2* );
  time_t getTime();
  std::vector<int>  getBadPackets(TH1 *hist, int what, float cutoff);


  const int nTowersEta = 96;
  const int nTowersPhi = 256;
  const int templateDepth=10000;
  int save = 0;
  TCanvas *TC[9] = {nullptr};
  TPad *transparent[9] = {nullptr};
  TPad *Pad[18] = {nullptr};
  TPad *warning[18] = {nullptr};
  TGraphErrors *gr[2] = {nullptr};
  TStyle* cemcStyle = nullptr;
  const int nSEBs = 8;
};

#endif /* CEMC_CEMCMONDRAW_H */
