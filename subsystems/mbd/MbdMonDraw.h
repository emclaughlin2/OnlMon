#ifndef MBD_MBDMONDRAW_H
#define MBD_MBDMONDRAW_H

#include <onlmon/OnlMonDraw.h>
#include "MbdMonDefs.h"

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

using namespace mbd_onlmon;

#define MAX_WARNING 16

class MbdMonDraw : public OnlMonDraw
{
 public:
  MbdMonDraw(const std::string &name = "MBDMON"); // same name as server!

  ~MbdMonDraw() override;

  int Init() override;
  int Draw(const std::string &what = "ALL") override;
  int MakePS(const std::string &what = "ALL") override;
  int MakeHtml(const std::string &what = "ALL") override;

 protected:
  TStyle *mbdStyle;
  
  int MakeCanvas(const std::string &name);
  int DrawFirst(const std::string &what = "ALL");
  int DrawSecond(const std::string &what = "ALL");
  int DrawHistory(const std::string &what = "ALL");
  int TimeOffsetTicks = -1;

  int Warning( TPad * pad, const float x, const float y, const int r, const std::string& msg);
  int ClearWarning( void );

  OnlMonDB *dbvars = nullptr;

  TCanvas *TC[nCANVAS] = {nullptr};
  TPad *transparent[nCANVAS] = {nullptr};

  //TPad *Pad[6] = {nullptr};


  // for all Pages
  TPad * PadTop[nCANVAS];
  TPaveText * PaveTop;
  TText * TextTop;

  TPaveText * PaveWarning[MAX_WARNING];
  TArc * ArcWarning[MAX_WARNING];
  TPad * PadWarning[MAX_WARNING];
  int nPadWarning[MAX_WARNING];
  int nWarning;

  // for the 1st Page
  TGraphErrors * TdcOver[nSIDE];
  TPad * PadTdcOver[nSIDE];
  TLine * LineTdcOver[nSIDE][2];
  TBox * BoxTdcOver[nSIDE];
  TH1 *FrameTdcOver[nSIDE];

  TGraph * nHit[nTRIGGER][nSIDE];
  TPad * PadnHit [nSIDE];
  TLine *LinenHit[nTRIGGER][nSIDE][2];
  TBox * BoxnHit[nTRIGGER][nSIDE];
  TPad * PadnHitStatus;
  TText *TextnHitStatus;
  TH1 *FramenHit[nSIDE];

  // for 2nd Page
  TH1 * HitTime[nSIDE];
  TPad * PadHitTime[nSIDE];
  TLine * LineHitTime[nSIDE][2];
  TArrow *ArrowHitTime[nSIDE];
  TText * TextHitTime[nSIDE];
  TF1 * FitHitTime[nSIDE];

  TH1        * AvrHitTime;
  TPad    * PadAvrHitTime;
  TLine  * LineAvrHitTime[2];
  TArrow *ArrowAvrHitTime;
  TText  * TextAvrHitTime;
  TF1     * FitAvrHitTime;

  TH1        * SouthHitTime;
  TPad    * PadSouthHitTime;
  TLine  * LineSouthHitTime[2];
  TArrow *ArrowSouthHitTime;
  TText  * TextSouthHitTime;
  TF1     * FitSouthHitTime;

  TH1        * NorthHitTime;
  TPad    * PadNorthHitTime;
  TLine  * LineNorthHitTime[2];
  TArrow *ArrowNorthHitTime;
  TText  * TextNorthHitTime;
  TF1     * FitNorthHitTime;

  TH1        * SouthChargeSum;
  TH1        * NorthChargeSum;
  TPad    * PadChargeSum;
  TText  * TextSouthChargeSum;
  TText  * TextNorthChargeSum;

  TLine  * LineTzeroZvtx[4];
  TText  * TextTzeroZvtx;

  TPad * PadBbcSummary;
  //TText * TextBbcSummaryHitTime;
  TText * TextBbcSummaryHitTime[2];
  TText * TextBbcSummaryGlobalOffset[2];
  TText * TextBbcSummaryZvertex;
  //TText * TextBbcSummaryTrigRate;
  TLatex * TextBbcSummaryTrigRate;

  TH1 * Zvtx;
  TH1 * Zvtx_bbll1;
//  TH1 * Zvtx_zdc;
//  TH1 * Zvtx_zdc_scale3;
  TH1 * Zvtx_bbll1_novtx;
  TH1 * Zvtx_bbll1_narrowvtx;
//  TH1 * Zvtx_bbll1_zdc;
  TPad * PadZvtx;
  TF1 * FitZvtx;
  TLine * LineZvtx[2];
  TArrow *ArrowZvtx;
  TText * TextZvtx;
  TText * TextZvtxNorth;
  TText * TextZvtxSouth;

  TH2 * ArmHit;
  TPad * PadArmHit;
  TText * TextArmHit;
  TArc * ArcArmHit;

  // for 3rd Page
  // It is for Only Expert
  TH2 * Adc;
  TPad *PadAdc;

  TPad * PadButton;
  TPad * PadMultiView;
  TH1 * MultiView1F;
  TH2 * MultiView2F;
  TPad * PadWarnings;
  TPaveText * PaveWarnings;

  // for 4th Page
  // Vertex Monitor
  TPad  * PadZVertex;
  TPad  * PadTzeroZVertex;
  TPad  * PadZVertexSummary;
  TText * TextZVertexExpress;
  TText * TextZVertexNotice;

  //  TText * TextZVertex[3];
  //TText * TextZVertex_scale[3];
  //TText * TextZVertex[4];
  //TText * TextZVertex_scale[4];
  //TText * TextZVertex_mean[4];
  TText * TextZVertex[5];//RUN11 pp
  TText * TextZVertex_scale[5];//RUN11 pp
  TText * TextZVertex_mean[5];//RUN11 pp

  TF1 *FitZvtxBBLL1NoVtx;
  TLatex *TextZVtxStatus[2];

  TH2 * TzeroZvtx;

  // Scale down factor
  TH1 * Prescale_hist;
  TSpectrum *tspec;


};

#endif /* MBD_MBDMONDRAW_H */
