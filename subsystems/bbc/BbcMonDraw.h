#ifndef BBC_BBCMONDRAW_H
#define BBC_BBCMONDRAW_H

#include "BbcMonDefs.h"

#include <onlmon/OnlMonDraw.h>

#include <string>  // for allocator, string

class OnlMonDB;
class TCanvas;
class TGraphErrors;
class TPad;
class TStyle;
class TPaveText;
class TText;
class TArc;
class TLine;
class TBox;
class TH1;
class TH2;
class TGraph;
class TArrow;
class TF1;
class TLatex;
class TSpectrum;

#define MAX_WARNING 16

class BbcMonDraw : public OnlMonDraw
{
 public:
  explicit BbcMonDraw(const std::string &name);

  ~BbcMonDraw() override;

  int Init() override;
  int Draw(const std::string &what = "ALL") override;
  int MakeHtml(const std::string &what = "ALL") override;
  int SavePlot(const std::string &what = "ALL", const std::string &type = "png") override;

 protected:
  TStyle *bbcStyle = nullptr;

  int MakeCanvas(const std::string &name);
  int DrawFirst(const std::string &what = "ALL");
  int DrawSecond(const std::string &what = "ALL");
  int DrawHistory(const std::string &what = "ALL");
  int TimeOffsetTicks = -1;

  int Warning(TPad *pad, const float x, const float y, const int r, const std::string &msg);
  int ClearWarning(void);

  OnlMonDB *dbvars = nullptr;

  TCanvas *TC[nCANVAS] = {nullptr};
  TPad *transparent[nCANVAS] = {nullptr};

  TPad *Pad[6] = {nullptr};

  // for all Pages
  TPad *PadTop[nCANVAS] = {};
  TPaveText *PaveTop = nullptr;
  TText *TextTop = nullptr;

  TPaveText *PaveWarning[MAX_WARNING] = {};
  TArc *ArcWarning[MAX_WARNING] = {};
  TPad *PadWarning[MAX_WARNING] = {};
  int nPadWarning[MAX_WARNING] = {};
  int nWarning = 0;

  // for the 1st Page
  TGraphErrors *TdcOver[nSIDE] = {};
  TPad *PadTdcOver[nSIDE] = {};
  TLine *LineTdcOver[nSIDE][2] = {};
  TBox *BoxTdcOver[nSIDE] = {};
  TH1 *FrameTdcOver[nSIDE] = {};

  TGraph *nHit[nTRIGGER][nSIDE] = {};
  TPad *PadnHit[nSIDE] = {};
  TLine *LinenHit[nTRIGGER][nSIDE][2] = {};
  TBox *BoxnHit[nTRIGGER][nSIDE] = {};
  TPad *PadnHitStatus = nullptr;
  TText *TextnHitStatus = nullptr;
  TH1 *FramenHit[nSIDE] = {};

  // for 2nd Page
  TH1 *HitTime[nSIDE] = {};
  TPad *PadHitTime[nSIDE] = {};
  TLine *LineHitTime[nSIDE][2] = {};
  TArrow *ArrowHitTime[nSIDE] = {};
  TText *TextHitTime[nSIDE] = {};
  TF1 *FitHitTime[nSIDE] = {};

  TH1 *AvrHitTime = nullptr;
  TPad *PadAvrHitTime = nullptr;
  TLine *LineAvrHitTime[2] = {};
  TArrow *ArrowAvrHitTime = nullptr;
  TText *TextAvrHitTime = nullptr;
  TF1 *FitAvrHitTime = nullptr;

  TH1 *SouthHitTime = nullptr;
  TPad *PadSouthHitTime = nullptr;
  TLine *LineSouthHitTime[2] = {};
  TArrow *ArrowSouthHitTime = nullptr;
  TText *TextSouthHitTime = nullptr;
  TF1 *FitSouthHitTime = nullptr;

  TH1 *NorthHitTime = nullptr;
  TPad *PadNorthHitTime = nullptr;
  TLine *LineNorthHitTime[2] = {};
  TArrow *ArrowNorthHitTime = nullptr;
  TText *TextNorthHitTime = nullptr;
  TF1 *FitNorthHitTime = nullptr;

  TH1 *SouthChargeSum = nullptr;
  TH1 *NorthChargeSum = nullptr;
  TPad *PadChargeSum = nullptr;
  TText *TextSouthChargeSum = nullptr;
  TText *TextNorthChargeSum = nullptr;

  TH2 *TimeWave = nullptr;
  TPad *PadTimeWave = nullptr;

  TH2 *ChargeWave = nullptr;
  TPad *PadChargeWave = nullptr;

  TH2 *SouthHitMap = nullptr;
  TPad *PadSouthHitMap = nullptr;

  TH2 *NorthHitMap = nullptr;
  TPad *PadNorthHitMap = nullptr;

  TLine *LineTzeroZvtx[4] = {};
  TText *TextTzeroZvtx = nullptr;

  TPad *PadBbcSummary = nullptr;
  // TText * TextBbcSummaryHitTime = nullptr;
  TText *TextBbcSummaryHitTime[2] = {};
  TText *TextBbcSummaryGlobalOffset[2] = {};
  TText *TextBbcSummaryZvertex = nullptr;
  // TText * TextBbcSummaryTrigRate = nullptr;
  TLatex *TextBbcSummaryTrigRate = nullptr;

  TH1 *Zvtx = nullptr;
  TH1 *Zvtx_bbll1 = nullptr;
  //  TH1 * Zvtx_zdc = nullptr;
  //  TH1 * Zvtx_zdc_scale3 = nullptr;
  TH1 *Zvtx_bbll1_novtx = nullptr;
  TH1 *Zvtx_bbll1_narrowvtx = nullptr;
  //  TH1 * Zvtx_bbll1_zdc = nullptrNor;
  TPad *PadZvtx = nullptr;
  TF1 *FitZvtx = nullptr;
  TLine *LineZvtx[2] = {};
  TArrow *ArrowZvtx = nullptr;
  TText *TextZvtx = nullptr;
  TText *TextZvtxNorth = nullptr;
  TText *TextZvtxSouth = nullptr;

  TH2 *ArmHit = nullptr;
  TPad *PadArmHit = nullptr;
  TText *TextArmHit = nullptr;
  TArc *ArcArmHit = nullptr;

  // for 3rd Page
  // It is for Only Expert
  TH2 *Adc = nullptr;
  TPad *PadAdc = nullptr;

  TPad *PadButton = nullptr;
  TPad *PadMultiView = nullptr;
  TH1 *MultiView1F = nullptr;
  TH2 *MultiView2F = nullptr;
  TPad *PadWarnings = nullptr;
  TPaveText *PaveWarnings = nullptr;

  // for 4th Page
  // Vertex Monitor
  TPad *PadZVertex = nullptr;
  TPad *PadTzeroZVertex = nullptr;
  TPad *PadZVertexSummary = nullptr;
  TText *TextZVertexExpress = nullptr;
  TText *TextZVertexNotice = nullptr;

  //  TText * TextZVertex[3] = {};
  // TText * TextZVertex_scale[3] = {};
  // TText * TextZVertex[4] = {};
  // TText * TextZVertex_scale[4] = {};
  // TText * TextZVertex_mean[4] = {};
  TText *TextZVertex[5] = {};        // RUN11 pp
  TText *TextZVertex_scale[5] = {};  // RUN11 pp
  TText *TextZVertex_mean[5] = {};   // RUN11 pp

  TLatex *TextZvtxStatus[2] = {};

  TH2 *TzeroZvtx = nullptr;

  // Scale down factor
  TH1 *Prescale_hist = nullptr;
  TSpectrum *tspec = nullptr;

private:
  int DrawDeadServer(TPad *transparent_pad) override;

};

#endif /* BBC_BBCMONDRAW_H */
