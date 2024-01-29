#ifndef TPC_TPCMONDRAW_H
#define TPC_TPCMONDRAW_H

#include <onlmon/OnlMonDraw.h>

#include <string>  // for allocator, string

class TCanvas;
class TGraphErrors;
class TPad;
class TH2;
class TPaveLabel;

class TpcMonDraw : public OnlMonDraw
{
 public:
  TpcMonDraw(const std::string &name);

  ~TpcMonDraw() override {}

  int Init() override;
  int Draw(const std::string &what = "ALL") override;
  int MakeHtml(const std::string &what = "ALL") override;
  int SavePlot(const std::string &what = "ALL", const std::string &type = "png") override;

 protected:
  int MakeCanvas(const std::string &name);
  int DrawFirst(const std::string &what = "ALL");
  int DrawSecond(const std::string &what = "ALL");
  int DrawTPCModules(const std::string &what = "ALL");
  int DrawTPCSampleSize(const std::string &what = "ALL");
  int DrawTPCCheckSum(const std::string &what = "ALL");
  int DrawTPCADCSample(const std::string &what = "ALL");
  int DrawTPCADCSampleLarge(const std::string &what = "ALL");
  int DrawTPCMaxADCModule(const std::string &what = "ALL");
  int DrawTPCRawADC1D(const std::string &what = "ALL");
  int DrawTPCMaxADC1D(const std::string &what = "ALL");
  int DrawTPCXYclusters(const std::string &what = "ALL");
  int DrawTPCXYclusters_unweighted(const std::string &what = "ALL");
  int DrawTPCZYclusters(const std::string &what = "ALL");
  int DrawTPCZYclusters_unweighted(const std::string &what = "ALL");
  time_t getTime();
  
  TCanvas *TC[15] = {nullptr};
  TPad *transparent[14] = {nullptr};
  TPad *Pad[10] = {nullptr};
  TGraphErrors *gr[2] = {nullptr};
  //TPC Module
  TH2 *dummy_his1 = nullptr;
  TH2 *dummy_his2 = nullptr;

  //TPC Module
  TH2 *dummy_his1_XY = nullptr;
  TH2 *dummy_his2_XY = nullptr;

  TH2 *dummy_his1_ZY = nullptr;

  TH2 *dummy_his1_XY_unw = nullptr;
  TH2 *dummy_his2_XY_unw = nullptr;

  TH2 *dummy_his1_ZY_unw = nullptr;

  TPaveLabel* NS18 = nullptr; //North Side labels
  TPaveLabel* NS17 = nullptr;
  TPaveLabel* NS16 = nullptr;
  TPaveLabel* NS15 = nullptr;
  TPaveLabel* NS14 = nullptr;
  TPaveLabel* NS13 = nullptr;
  TPaveLabel* NS12 = nullptr;
  TPaveLabel* NS23 = nullptr;
  TPaveLabel* NS22 = nullptr;
  TPaveLabel* NS21 = nullptr;
  TPaveLabel* NS20 = nullptr;
  TPaveLabel* NS19 = nullptr;

  TPaveLabel* SS00 = nullptr; //South Side labels
  TPaveLabel* SS01 = nullptr;
  TPaveLabel* SS02 = nullptr;
  TPaveLabel* SS03 = nullptr;
  TPaveLabel* SS04 = nullptr;
  TPaveLabel* SS05 = nullptr;
  TPaveLabel* SS06 = nullptr;
  TPaveLabel* SS07 = nullptr;
  TPaveLabel* SS08 = nullptr;
  TPaveLabel* SS09 = nullptr;
  TPaveLabel* SS10 = nullptr;
  TPaveLabel* SS11 = nullptr;
};

#endif /* TPC_TPCMONDRAW_H */
