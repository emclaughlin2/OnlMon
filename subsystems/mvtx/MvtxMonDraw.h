#ifndef MVTX_MVTXMONDRAW_H
#define MVTX_MVTXMONDRAW_H

#include <onlmon/OnlMonDraw.h>

#include <string>  // for allocator, string
#include <TColor.h>

class OnlMonDB;
class TCanvas;
class TGraphErrors;
class TPad;
class OnlMonClient;
class TPaveText;
class TH1;
class TH2Poly;


class MvtxMonDraw : public OnlMonDraw
{
 public:
  MvtxMonDraw(const std::string &name);

  ~MvtxMonDraw() override {}

  int Init() override;
  int Draw(const std::string &what = "ALL") override;
  int MakeHtml(const std::string &what = "ALL") override;
  int SavePlot(const std::string &what = "ALL", const std::string &type = "png") override;

     const static int NSTAVE = 48;
  const static int NCHIP = 9;
      const static int NLAYERS = 3;
  static constexpr int NFlags = 3;
  const int NStaves[NLAYERS] = { 12, 16, 20 };
  const int StaveBoundary[NLAYERS + 1] = { 0, 12, 28, 48};
  std::string mLaneStatusFlag[NFlags] = { "WARNING", "ERROR", "FAULT" };
  int PublishHistogram(TCanvas *c, int pad, TH1 *h, const char* opt = "");
  int PublishHistogram(TPad *p, TH1 *h, const char* opt = "");
  //template <typename T>
  //int PublishHistogram(TPad *p, int pad, T h, const char* opt = "");
  int PublishHistogram(TPad *p, int pad, TH1 *h, const char* opt = "");
  void PublishStatistics(TCanvas *c,OnlMonClient *cl );
  template <typename T>
  int MergeServers(T *h);
  void formatPaveText(TPaveText* aPT, float aTextSize, Color_t aTextColor, short aTextAlign, const char* aText);
  enum Quality {Good, Medium, Bad};
  std::vector<Quality> analyseForError(TH2Poly* over1, TH2Poly* over2,TH2Poly* over3, TH1* decErr);
  void DrawPave(std::vector<MvtxMonDraw::Quality> status, int position, const char* what = "");

 private:
static const int NFlx = 6;
 static constexpr int NCols = 1024;
  static constexpr int NRows = 512;
  static constexpr int NPixels = NRows * NCols;
  const int chipmapoffset[3] = { 0, 12, 28};

  int MakeCanvas(const std::string &name);
  int DrawFirst(const std::string &what = "ALL");
  int DrawSecond(const std::string &what = "ALL");
  int DrawHitMap(const std::string &what = "ALL");
  int DrawGeneral(const std::string &what = "ALL");
  int DrawFEE(const std::string &what = "ALL");
  int DrawOCC(const std::string &what = "ALL");
  int DrawFHR(const std::string &what = "ALL");
  int DrawHistory(const std::string &what = "ALL");
  time_t getTime();
  int TimeOffsetTicks = -1;
  TCanvas *TC[6] = {nullptr};
  TPad *transparent[1] = {nullptr};
  TPad *Pad[6] = {nullptr};
  TGraphErrors *gr[6] = {nullptr};
  OnlMonDB *dbvars[NFlx] = {nullptr};


 
int maxbadchips = 2;



};

#endif /* MVTX_MVTXMONDRAW_H */
