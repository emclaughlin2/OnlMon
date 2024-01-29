#ifndef TPC_TPCMON_H
#define TPC_TPCMON_H

#include <onlmon/OnlMon.h>

#include <map>
#include <tpc/TpcMap.h> //this needs to be pointed to coresoftware - not sure how to do that on EBDCXX...
#include <memory>
#include <string>
#include <cmath>
#include <vector>


class Event;
class TH1;
class TH2;
class TTree;

class TpcMon : public OnlMon
{
 public:
  TpcMon(const std::string &name);
  virtual ~TpcMon();

  int process_event(Event *evt);
  int Init();
  int BeginRun(const int runno);
  int Reset();

 protected:
  int evtcnt = 0;
  int idummy = 0;

  static const int N_rBins = 4; //(one inner bin not filled, 2nd bin is R1, 3rd bin is R2, 4th bin is R3)
  const int N_thBins = 12; //(12 theta bins of uniform angle (360/12 = 30 degrees = TMath::Pi()/6 ~= 0.523 rad)
  const double rBin_edges[N_rBins+1] = {0.0, 0.256, 0.504, 0.752, 1.00}; //variable edges for radial dims

  //for X-Y channel plot
  static const int N_rBins_XY = 66; //From Evgeny code - global to all modules
  const double r_bins[N_rBins_XY + 1] = {217.83,
                                 223.83, 229.83, 235.83, 241.83, 247.83, 253.83, 259.83, 265.83, 271.83, 277.83, 283.83, 289.83, 295.83, 301.83, 306.83,
                                 311.05, 317.92, 323.31, 329.27, 334.63, 340.59, 345.95, 351.91, 357.27, 363.23, 368.59, 374.55, 379.91, 385.87, 391.23, 397.19, 402.49,
                                 411.53, 421.70, 431.90, 442.11, 452.32, 462.52, 472.73, 482.94, 493.14, 503.35, 513.56, 523.76, 533.97, 544.18, 554.39, 564.59, 574.76,
                                 583.67, 594.59, 605.57, 616.54, 627.51, 638.48, 649.45, 660.42, 671.39, 682.36, 693.33, 704.30, 715.27, 726.24, 737.21, 748.18, 759.11}; // From Evgeny code - global to all modules

  static const int N_phi_binx_XY_R1 = 1152; // from Evgeny code (0 to 2pi)
  static const int N_phi_binx_XY_R2 = 1536; // from Evgeny code (0 to 2pi)
  static const int N_phi_binx_XY_R3 = 2304; // from Evgeny code (0 to 2pi)

  static const int N_z_bins = 3350; // (50 ns time bins * ~12.3 micron/ns drift time = 615 microns). Drift Length = 1.03 * 10^6 microns. 2 *1.03 *10^6 microns / 615 microns = 3350 bins

  TH2 *NorthSideADC_clusterXY_R1 = nullptr;
  TH2 *NorthSideADC_clusterXY_R2 = nullptr;
  TH2 *NorthSideADC_clusterXY_R3 = nullptr;

  TH2 *SouthSideADC_clusterXY_R1 = nullptr;
  TH2 *SouthSideADC_clusterXY_R2 = nullptr;
  TH2 *SouthSideADC_clusterXY_R3 = nullptr;

  TH2 *NorthSideADC_clusterXY_R1_unw = nullptr;
  TH2 *NorthSideADC_clusterXY_R2_unw = nullptr;
  TH2 *NorthSideADC_clusterXY_R3_unw = nullptr;

  TH2 *SouthSideADC_clusterXY_R1_unw = nullptr;
  TH2 *SouthSideADC_clusterXY_R2_unw = nullptr;
  TH2 *SouthSideADC_clusterXY_R3_unw = nullptr;

  TH2 *NorthSideADC_clusterZY = nullptr;
  TH2 *SouthSideADC_clusterZY = nullptr;

  TH2 *NorthSideADC_clusterZY_unw = nullptr;
  TH2 *SouthSideADC_clusterZY_unw = nullptr;

  TH1 *tpchist1 = nullptr;
  TH2 *tpchist2 = nullptr;
  
  TH2 *NorthSideADC = nullptr;
  TH2 *SouthSideADC = nullptr;

  TH2 *ADC_vs_SAMPLE = nullptr;
  TH2 *ADC_vs_SAMPLE_large = nullptr;

  TH1 *sample_size_hist = nullptr;
  TH1 *Check_Sum_Error = nullptr;
  TH1 *Check_Sums = nullptr;

  TH2 *MAXADC = nullptr;

  TH1 *RAWADC_1D_R1= nullptr;
  TH1 *MAXADC_1D_R1 = nullptr;

  TH1 *RAWADC_1D_R2= nullptr;
  TH1 *MAXADC_1D_R2 = nullptr;

  TH1 *RAWADC_1D_R3= nullptr;
  TH1 *MAXADC_1D_R3 = nullptr;

  TpcMap M; //declare Martin's map

  int starting_BCO;
  int rollover_value;
  int current_BCOBIN;

  int serverid;

  void Locate(int id, float *rbin, float *thbin);
  int Index_from_Module(int sec_id, int fee_id);
  int Module_ID(int fee_id);
  int Max_Nine(int one, int two, int three, int four, int five, int six, int seven, int eight, int nine);
  bool side(int server_id);
  std::pair<float, float> calculateMeanAndStdDev(const std::vector<int>& values);
};

#endif /* TPC_TPCMON_H */
