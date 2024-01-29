#ifndef __BBCMONDEFS_H__
#define __BBCMONDEFS_H__

#define AUAU  // PROTONPROTON //PAL //PAU //DAU //HE3AU //CUCU //AUAU

namespace bbc_onlmon
{

  const int num_pmt = 128;                 // number of PMT
  const int tdc_max = 4096;                // TDC number of max
  const int tdc_min_overflow = 2600;       // TDC number of max
  const int tdc_max_overflow = 3299;       // TDC number of max
  const float max_armhittime = 22;         // Minimum of Tdc-ArmHitTime
  const float min_armhittime = 0;          // Maximum od Tdc-ArmHitTime
  static const double TDC_FIT_MAX = 19.0;  // compatible type with TF::GetParameter(1)
  static const double TDC_FIT_MIN = 1.0;
  //static const int zvtnbin = 180;
  static const int zvtnbin = 600;
  static const float max_zvertex = 300;
  static const float min_zvertex = -300;

#define nPMT_1SIDE_BBC 64
#define nPMT_BBC (nPMT_1SIDE_BBC * 2)

#define nCANVAS 4

#define nSIDE 2
#define nTRIGGER 2

  static const char *SIDE_Str[nSIDE] __attribute__((unused)) = {"South", "North"};
  static const char *TRIGGER_str[nTRIGGER] __attribute__((unused)) = {"minimum_bias", "laser"};

#define nTDC 1
  // static const char *TDC_STR[nTDC] = {"TDC0","TDC1"};
  // static const char *TDC_Str[nTDC] = {"Tdc0","Tdc1"};
  // static const char *TDC_str[nTDC] = {"tdc0","tdc1"};

  static const int nBIN_TDC = 64;

  static const float TDC_CONVERSION_FACTOR = 0.007;
  static const int nBIN_ADC = 128;
  static const int MAX_ADC = 4096;

  static const float VIEW_OVERFLOW_MAX = 20;
  static const float VIEW_OVERFLOW_MIN = -20;
  static const float BBC_TDC_OVERFLOW_REGULAR_MAX = 3;        // Sigma
  static const float BBC_TDC_OVERFLOW_REGULAR_MIN = -3;       // Sigma
  static const float BBC_TDC_OVERFLOW_REGULAR_RMS_MAX = 3.0;  // sigma
  static const float BBC_ZVERTEX_OVERFLOW_REGULAR_MAX = 50;   // cm
  static const float BBC_ZVERTEX_OVERFLOW_REGULAR_MIN = -50;  // cm

  static const float BBC_nHIT_LASER_MAX[nSIDE] = {1.00, 1.00};
  static const float BBC_nHIT_LASER_MIN[nSIDE] = {0.95, 0.95};
  // changed at Jan. 25. 2003
  // Due to the laser timing in ch25 is faster than normal ch,
  // the ich25 tend to have under-flow of the tdc(tdc<0ch).
  // Then, I changed the BBC_nHIT_LASER_MIN[nSIDE] to 0.95%
  // Takashi Hachiya

  /*
    0 : for 6 devided pads
  */
  static const float CANVAS_XAXIS_TITLE_SIZE[] = {0.060};
  static const float CANVAS_XAXIS_TITLE_OFFSET[] = {1.000};
  static const float CANVAS_XAXIS_LABEL_SIZE[] = {0.050};
  static const float CANVAS_XAXIS_LABEL_OFFSET[] = {0.005};
  static const float CANVAS_XAXIS_PAD_OFFSET[] = {1.000};

  static const float CANVAS_YAXIS_TITLE_SIZE[] = {0.070};
  static const float CANVAS_YAXIS_TITLE_OFFSET[] = {0.500};
  static const float CANVAS_YAXIS_LABEL_SIZE[] = {0.050};
  static const float CANVAS_YAXIS_LABEL_OFFSET[] = {0.005};
  static const float CANVAS_YAXIS_PAD_OFFSET[] = {0.150};

  static const float BBC_TYPICAL_REGULAR_TDC0_MEAN = 10.5;
  static const float BBC_MAX_REGULAR_TDC0_MEAN = 13.5;  // for 7.7GeV run
  static const float BBC_MIN_REGULAR_TDC0_MEAN = 7.5;   // for 7.7GeV run
  // static const float BBC_MAX_REGULAR_TDC0_MEAN = 12.5;
  // static const float BBC_MIN_REGULAR_TDC0_MEAN =  8.5;
  static const float BBC_TYPICAL_REGULAR_TDC1_MEAN = 10.5;
  static const float BBC_MAX_REGULAR_TDC1_MEAN = 13.5;  // for 7.7GeV run
  static const float BBC_MIN_REGULAR_TDC1_MEAN = 7.5;   // for 7.7GeV run
  // static const float BBC_MAX_REGULAR_TDC1_MEAN = 12.5;
  // static const float BBC_MIN_REGULAR_TDC1_MEAN =  8.5;
  static const float BBC_MAX_REGULAR_ZVERTEX_MEAN = 20.0;
  static const float BBC_MIN_REGULAR_ZVERTEX_MEAN = -20.0;
  static const float BBC_MIN_WORNING_STATISTICS_FOR_ZVERTEX_MEAN = 5000;

  static const int BBC_PACKET_ID = 1001;
  static const int BBCN_PACKET_ID = 1002;
  static const int BBCS_PACKET_ID = 1003;
  static const int BBC_ZDC_PACKET_ID = 13001;
  static const int BBC_TRIG_PACKET_ID = 14001;

#ifdef HE3AU
  static const unsigned long BBC_ZDC_TRIG = 0x4;     // 0x2;
  static const unsigned long BBC_BBC_TRIG = 0x2;     // 0x1;
  static const unsigned long BBC_BBZDC_TRIG = 0x8;   // 0x8;
  static const unsigned long BBC_BBLL1_TRIG = 0x40;  // 0x1000;
  // for He3Au
  static const float BBC_nHIT_MB_MAX[nSIDE] = {0.90, 0.90};
  static const float BBC_nHIT_MB_MIN[nSIDE] = {0.05, 0.05};
  static const float ZDC_BBC_TRIGSCALE = 1;
  // Charge Sum plot
  static const int MAX_ADC_MIP = 100;     // for He3Au
  static const int MAX_CHARGE_SUM = 500;  // for He3Au
#endif

#ifdef CUCU
  static const unsigned long BBC_ZDC_TRIG = 0x4;     // 0x2;
  static const unsigned long BBC_BBC_TRIG = 0x2;     // 0x1;
  static const unsigned long BBC_BBZDC_TRIG = 0x8;   // 0x8;
  static const unsigned long BBC_BBLL1_TRIG = 0x40;  // 0x1000;
  // for CuCu
  static const float BBC_nHIT_MB_MAX[nSIDE] = {0.90, 0.90};
  static const float BBC_nHIT_MB_MIN[nSIDE] = {0.10, 0.10};
  static const float ZDC_BBC_TRIGSCALE = 10;
  // Charge Sum plot
  static const int MAX_ADC_MIP = 50;       // for AuAu
  static const int MAX_CHARGE_SUM = 1500;  // for AuAu
#endif
                                           // last year
#ifdef AUAU
  static const unsigned long BBC_ZDC_TRIG = 0x4;     // 0x2;
  static const unsigned long BBC_BBC_TRIG = 0x2;     // 0x1;
  static const unsigned long BBC_BBZDC_TRIG = 0x8;   // 0x8;
  static const unsigned long BBC_BBLL1_TRIG = 0x40;  // 0x1000;
  // for AuAu
  static const float BBC_nHIT_MB_MAX[nSIDE] = {0.90, 0.90};
  // static const float BBC_nHIT_MB_MIN[nSIDE] = { 0.10, 0.10 };
  static const float BBC_nHIT_MB_MIN[nSIDE] = {0.05, 0.05};
  static const float ZDC_BBC_TRIGSCALE = 1;
  // Charge Sum plot
  static const int MAX_ADC_MIP = 100;      // for AuAu
  static const int MAX_CHARGE_SUM = 1500;  // for AuAu 200GeV
// static const int MAX_ADC_MIP  = 100;     // for AuAu 15GeV
// static const int MAX_CHARGE_SUM = 150; // for AuAu 15GeV
// static const int MAX_ADC_MIP  = 50;     // for AuAu
// static const int MAX_CHARGE_SUM = 400; // for AuAu
// static const int MAX_CHARGE_SUM = 1500; // for AuAu 200GeV
// static const int MAX_CHARGE_SUM = 750; // for AuAu 62GeV
// static const int MAX_CHARGE_SUM = 150; // for AuAu 7.7GeV
#endif

#ifdef DAU
  static const unsigned long BBC_ZDC_TRIG = 0x4;     // 0x2;
  static const unsigned long BBC_BBC_TRIG = 0x2;     // 0x1;
  static const unsigned long BBC_BBZDC_TRIG = 0x8;   // 0x8;
  static const unsigned long BBC_BBLL1_TRIG = 0x40;  // 0x1000;
  // for dAu
  static const float BBC_nHIT_MB_MAX[nSIDE] = {0.80, 0.40};
  static const float BBC_nHIT_MB_MIN[nSIDE] = {0.02, 0.005};
  static const float ZDC_BBC_TRIGSCALE = 3;
  // Charge Sum plot
  static const int MAX_ADC_MIP = 50;      // for dAu
  static const int MAX_CHARGE_SUM = 150;  // for dAu
#endif

#ifdef PAU
  static const unsigned long BBC_ZDC_TRIG = 0x4;     // 0x2;
  static const unsigned long BBC_BBC_TRIG = 0x2;     // 0x1;
  static const unsigned long BBC_BBZDC_TRIG = 0x8;   // 0x8;
  static const unsigned long BBC_BBLL1_TRIG = 0x40;  // 0x1000;
  // for pAu
  static const float BBC_nHIT_MB_MAX[nSIDE] = {0.80, 0.60};
  static const float BBC_nHIT_MB_MIN[nSIDE] = {0.02, 0.01};
  static const float ZDC_BBC_TRIGSCALE = 3;
  // Charge Sum plot
  static const int MAX_ADC_MIP = 100;     // for pAu
  static const int MAX_CHARGE_SUM = 250;  // for pAu
#endif

#ifdef PAL
  static const unsigned long BBC_ZDC_TRIG = 0x4;     // 0x2;
  static const unsigned long BBC_BBC_TRIG = 0x2;     // 0x1;
  static const unsigned long BBC_BBZDC_TRIG = 0x8;   // 0x8;
  static const unsigned long BBC_BBLL1_TRIG = 0x40;  // 0x1000;
  // for pAl
  static const float BBC_nHIT_MB_MAX[nSIDE] = {0.80, 0.60};
  static const float BBC_nHIT_MB_MIN[nSIDE] = {0.02, 0.01};
  static const float ZDC_BBC_TRIGSCALE = 3;
  // Charge Sum plot
  static const int MAX_ADC_MIP = 100;     // for pAu
  static const int MAX_CHARGE_SUM = 200;  // for pAu
#endif

#ifdef PROTONPROTON

  static const unsigned long BBC_ZDC_TRIG = 0x4000000;           // 27 zdcns;
  static const unsigned long BBC_BBC_TRIG = 0x4;                 // 3bit bbll1>1 only
  static const unsigned long BBC_BBZDC_TRIG = 0x8;               // 0x8;
  static const unsigned long BBC_BBLL1_TRIG = 0x40;              // 0x1000;
  static const unsigned long BBC_BBLL1_NOVTX_TRIG = 0x80000000;  //;
  static const unsigned long BBC_RBIB_NTC_WIDE = 0x100000;       //;

  // for pp
  static const float BBC_nHIT_MB_MAX[nSIDE] = {0.60, 0.60};
  static const float BBC_nHIT_MB_MIN[nSIDE] = {0.01, 0.01};
  static const float ZDC_BBC_TRIGSCALE = 55;
  // Charge Sum plot
  static const int MAX_ADC_MIP = 100;     // for pp
  static const int MAX_CHARGE_SUM = 100;  // for pp
#endif

#ifdef UU
  static const unsigned long BBC_ZDC_TRIG = 0x4;     // 0x2;
  static const unsigned long BBC_BBC_TRIG = 0x2;     // 0x1;
  static const unsigned long BBC_BBZDC_TRIG = 0x8;   // 0x8;
  static const unsigned long BBC_BBLL1_TRIG = 0x40;  // 0x1000;
  static const float BBC_nHIT_MB_MAX[nSIDE] = {0.90, 0.90};
  static const float BBC_nHIT_MB_MIN[nSIDE] = {0.05, 0.05};
  static const float ZDC_BBC_TRIGSCALE = 1;
  // Charge Sum plot
  static const int MAX_ADC_MIP = 50;       // for UU
  static const int MAX_CHARGE_SUM = 1500;  // for UU 200GeV
#endif

  static const unsigned long BBC_BBLL1_COPY_TRIG = 0x1000000;  // new trigger
  // static const int BBC_BBLL1_TRIG    = 0x10;       //0x1000;
  //static const unsigned long BBC_LASER_TRIG = 0x40000000;  // 0x40000000;

  // the mean of TDC between south and north should be this value.
  // static const int BBC_DEFAULT_OFFSET = 1500;
  static const float BBC_DEFAULT_OFFSET = 10.5;

  // for trigger rate between BBCLL1 and ZDC within certain vertex cut
  // the number should be within 3-5
  static const float BBC_ZVERTEX_CUT_FOR_TRIG_RATE = 30.0;  // [cm]

  // color definition
  // static const int BBC_COLOR_ZVTX[4] = {4,2,7,1};
  static const int BBC_COLOR_ZVTX[5] = {4, 2, 7, 6, 1};

  static const int BBC_NSAMPLES = 31;  // number of digitizer samples

}  // namespace bbc_onlmon

#endif /*__BBCMONDEFS_H__ */
