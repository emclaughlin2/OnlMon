#ifndef MVTX_MVTXMON_H
#define MVTX_MVTXMON_H

#include <onlmon/OnlMon.h>

#include <TString.h>
#include <TPaveText.h>

#include <map>
#include <cmath>


class Event;
class TH1;
class TH2;
class TH1I;
class TH2I;
class TH3I;
class TH1D;
class TH2D;
class TH2Poly;
class map;
class pair;

class MvtxRawHit;
class Packet;

class MvtxMon : public OnlMon
{
 public:
  MvtxMon(const std::string &name);
  virtual ~MvtxMon();

  int process_event(Event *evt);
  int Init();
  int BeginRun(const int runno);
  int Reset();


 protected:
  int evtcnt = 0;
  int idummy = 0;

  static constexpr int NSTAVE = 48;
  static constexpr int NCHIP = 9;
  static constexpr int NLAYERS = 3;
  static constexpr int NFees = 48 * 3;
  static constexpr int NTrigger = 13;
  static constexpr int NLanesMax = 9;
  static constexpr int NCols = 1024;
  static constexpr int NRows = 512;
  int mHitPerChip[NLAYERS][NSTAVE][NCHIP] = {};
  static constexpr int NFlags = 3;

  int mMaxGeneralAxisRange = -3;  // the range of TH2Poly plots z axis range, pow(10, mMinGeneralAxisRange) ~ pow(10, mMaxGeneralAxisRange)
  int mMinGeneralAxisRange = -12;
  std::string mLaneStatusFlag[NFlags] = { "WARNING", "ERROR", "FAULT" };

  const int StaveBoundary[4] = { 0, 12, 28, 48};
  const int LayerBoundaryFEE[NLAYERS - 1] = { 35, 83};

  TH2	*mvtxmon_ChipStaveOcc= nullptr;
  TH1	*mvtxmon_ChipStave1D= nullptr;
  TH1	*mvtxmon_ChipFiredHis= nullptr;
  TH1	*mvtxmon_EvtHitChip= nullptr;
  TH1 *mvtxmon_EvtHitDis= nullptr;
  TH2 *mvtxmon_HitMap[NSTAVE][NCHIP] = {nullptr}; 
  TH2Poly* mvtxmon_GeneralOccupancy = nullptr;
  TH2Poly* mvtxmon_LaneStatusOverview[NFlags] = {nullptr}; 

  TH1I* mRCDAQevt= nullptr;

  TH1I *hChipStrobes = nullptr;
  TH1I *hChipL1 = nullptr;
 //fee
  TH2I* mTriggerVsFeeId= nullptr;
  TH1I* mTrigger= nullptr;
  TH2I* mLaneInfo= nullptr;
  TH2I* mLaneStatus[NFlags] = {nullptr};
  TH2I* mLaneStatusCumulative[NFlags] = {nullptr};
  TH1I* mLaneStatusSummary[NLAYERS] = {nullptr};
  TH1I* mLaneStatusSummaryIB= nullptr;
  TString mTriggerType[NTrigger] = { "ORBIT", "HB", "HBr", "HC", "PHYSICS", "PP", "CAL", "SOT", "EOT", "SOC", "EOC", "TF", "INT" };

  //raw task
  TH1D* hErrorPlots= nullptr;
  TH2D* hErrorFile= nullptr;
  TH1D* hOccupancyPlot[NLAYERS] = {nullptr};
  TH2I* hEtaPhiHitmap[NLAYERS] = {nullptr};
  TH2D* hChipStaveOccupancy[NLAYERS] = {nullptr};
  TH3I* hChipHitmap = nullptr;
  TH3I* hChipHitmap_evt = nullptr;

  //fhr
  TH2I* mErrorVsFeeid= nullptr;
  TH2Poly* mGeneralOccupancy= nullptr;
  TH2Poly* mGeneralNoisyPixel= nullptr;
  TH2D* mDeadChipPos[NLAYERS] = {nullptr};
  TH2D* mAliveChipPos[NLAYERS] = {nullptr};
  TH2D* mTotalDeadChipPos= nullptr;
  TH2D* mTotalAliveChipPos= nullptr;
  TH2D* mChipStaveOccupancy[NLAYERS] = {nullptr};
  TH1D* mOccupancyPlot[NLAYERS] = {nullptr};
  TH2D* hChipStaveNoisy[NLAYERS] = {nullptr};

  float mOccupancyCutForNoisyPixel = 0.2;
  int mNoisyPixelNumber[3][20][9] = { { 0 } };

  static constexpr int NError = 11;
  static constexpr int NErrorExtended = 19;
 
 /* TString ErrorType[NError] = { "Error ID 1: ErrPageCounterDiscontinuity", "Error ID 2: ErrRDHvsGBTHPageCnt",
                                "Error ID 3: ErrMissingGBTHeader", "Error ID 4: ErrMissingGBTTrailer", "Error ID 5: ErrNonZeroPageAfterStop",
                                "Error ID 6: ErrUnstoppedLanes", "Error ID 7: ErrDataForStoppedLane", "Error ID 8: ErrNoDataForActiveLane",
                                "Error ID 9: ErrIBChipLaneMismatch", "Error ID 10: ErrCableDataHeadWrong",
                                "Error ID 11: Jump in RDH_packetCounter" };*/

  TString ErrorType[NError] = { "Error ID 1", "Error ID 2",
                                "Error ID 3", "Error ID 4", "Error ID 5",
                                "Error ID 6", "Error ID 7", "Error ID 8",
                                "Error ID 9", "Error ID 10",
                                "Error ID 11" };


  const int NFiles = 6;
  int ntriggers = 0;


  const int NStaves[3] = { 12, 16, 20 };
  const float StartAngle[3] = { 16.997 / 360 * (M_PI * 2.), 17.504 / 360 * (M_PI * 2.), 17.337 / 360 * (M_PI * 2.)}; // start angle of first stave in each layer
  const float MidPointRad[3] = { 23.49, 31.586, 39.341}; 

  const int mapstave[3][20] = {{10,11,12,1,2,3,4,5,6,7,8,9,0,0,0,0,0,0,0,0},{25,26,27,28,13,14,15,16,17,18,19,20,21,22,23,24,0,0,0,0},{44,45,46,47,48,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43}};


  typedef struct linkId
{
  uint32_t layer = 0xFF;
  uint32_t stave = 0xFF;
  uint32_t gbtid = 0xFF;
} LinkId_t;

  LinkId_t DecodeFeeid(const uint16_t& feeid)
  {
    LinkId_t ret = {};
    ret.layer = (feeid >> 12) & 0x7;
    ret.stave = feeid & 0x1F;
    ret.gbtid = (feeid >> 8) & 0x3;
    return ret;
  }

  private:
  unsigned short decode_row(int hit){	return hit >> 16;}
  unsigned short decode_col(int hit){	return hit & 0xffff;}
  void getStavePoint(int layer, int stave, double* px, double* py);
  void drawLayerName(TH2* histo2D);
  void createPoly(TH2Poly *h);

  Packet **plist = nullptr;
  /*unsigned int m_NumSpecialEvents = 0;
  std::map<uint64_t, std::set<int>> m_BeamClockFEE;
  std::map<uint64_t, std::vector<MvtxRawHit *>> m_MvtxRawHitMap;
  std::map<int, uint64_t> m_FEEBclkMap;
  std::map<int, uint64_t> m_FeeStrobeMap;
  std::set<uint64_t> m_BclkStack;
  std::set<uint64_t> gtmL1BcoSet;*/


};


#endif /* MVTX_MVTXMON_H */
