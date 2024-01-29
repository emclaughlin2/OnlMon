#ifndef CEMC_CEMCMON_H
#define CEMC_CEMCMON_H


#include <onlmon/OnlMon.h>

#include <vector>

class CaloWaveformFitting;
class TowerInfoContainer;
class Event;
class TH1;
class TH2;
class Packet;
class runningMean;

class CemcMon : public OnlMon
{
 public:
  CemcMon(const std::string &name);
  virtual ~CemcMon();

  int process_event(Event *evt);
  int Init();
  int BeginRun(const int runno);
  int Reset();

 protected:
  std::vector<float> getSignal(Packet *p, const int channel);
  std::vector<float> anaWaveformFast(Packet *p, const int channel);
  std::vector<float> anaWaveformTemp(Packet *p, const int channel);

  int idummy = 0;
  TH2 *cemc_occupancy = nullptr;
  TH2 *cemc_runningmean = nullptr;
  TH1 *cemc_signal = nullptr;
  TH1 *h1_cemc_adc = nullptr;
  
  const int Nsector = 64;
  const int Ntower = 64*2*192;
  const int packetlow = 6001;
  const int packethigh = 6128;
  const int m_nSamples = 16;
  const int m_nChannels = 192;
  const int templateDepth = 10000;
  int eventCounter = 0;
  
  TH1* h1_packet_chans = nullptr;
  TH1* h1_packet_length = nullptr;
  TH1* h1_packet_number = nullptr;
  TH1* h1_cemc_fitting_sigDiff = nullptr;
  TH1* h1_cemc_fitting_pedDiff = nullptr;
  TH1* h1_cemc_fitting_timeDiff = nullptr;
  TH2* h2_cemc_hits = nullptr;
  TH2* h2_waveform_twrAvg = nullptr;
  TH1* h1_waveform_time = nullptr;
  TH1* h1_waveform_pedestal = nullptr;
  TH2* h2_cemc_rm = nullptr;
  TH2* h2_cemc_mean = nullptr;
  TH1* h1_sectorAvg_total = nullptr;
  TH1* h1_event = nullptr;
  TH1* h1_rm_sectorAvg[100] = {nullptr};

  std::vector<runningMean*> rm_vector_twr;
  std::vector<runningMean*> rm_vector_sectAvg;

  std::string runtypestr = "Unknown";
  std::string id_string;

  CaloWaveformFitting* WaveformProcessingFast = nullptr;
  CaloWaveformFitting* WaveformProcessingTemp = nullptr;


  std::vector<runningMean*> rm_vector; 

};

#endif /* CEMC_CEMCMON_H */
