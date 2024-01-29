#ifndef LL1_LL1MON_H
#define LL1_LL1MON_H

#include <onlmon/OnlMon.h>

class Event;
class TH1;
class TH2;

class LL1Mon : public OnlMon
{
 public:
  LL1Mon(const std::string &name = "LL1MON");
  virtual ~LL1Mon();

  int process_event(Event *evt);
  int Init();
  int BeginRun(const int runno);
  int Reset();

 protected:
  int evtcnt = 0;
  int idummy = 0;
  int thresh=2;
  TH2* h_line_up=nullptr;
  TH2* h_nhit_corr=nullptr;
  TH1 *h_hit_n= nullptr;
  TH1 *h_hit_s= nullptr;
  TH1* h_nhit_n=nullptr;
  TH1* h_nhit_s=nullptr;
  TH1* h_nhit_n1=nullptr;
  TH1* h_nhit_n2=nullptr;
  TH1* h_nhit_s1=nullptr;
  TH1* h_nhit_s2=nullptr;
  TH2* h_nhit_n_corr=nullptr;
  TH2* h_nhit_s_corr=nullptr;
  TH1* h_time_diff=nullptr;
  TH2* h_hit_check=nullptr;
  
  float binstart = -0.5;
  float binend = 64.5;
  float binend2 = 32.5;
  const int nSamples = 20;
  const int nChannels =52;
  const int nChargeChannels = 8;
  const int nHitSampleIdx = 9;

  const int nhitbins = 65;
};

#endif /* LL1_LL1MON_H */
