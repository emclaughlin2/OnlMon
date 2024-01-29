#ifndef HCAL_HCALMON_H
#define HCAL_HCALMON_H

#include <onlmon/OnlMon.h>
#include <vector>

#include <vector>

class CaloWaveformFitting;
class TowerInfoContainer;
class Event;
class TH1;
class TH2;
class Packet;
class runningMean;

class HcalMon : public OnlMon
{
 public:
  HcalMon(const std::string& name);
  virtual ~HcalMon();

  int process_event(Event* evt);
  int Init();
  int BeginRun(const int runno);
  int Reset();
  std::vector<float> getSignal(Packet* p, const int channel);
  std::vector<float> anaWaveform(Packet* p, const int channel);

 protected:
  int evtcnt = 0;
  int idummy = 0;
  const int Nsector = 32;
  const int Ntower = 1536;
  const int m_nChannels = 192;
  int packetlow = 8001;
  int packethigh = 8008;
  // TH1 *hcalhist1 = nullptr;
  // TH2 *hcalhist2 = nullptr;
  TH2* h2_hcal_hits = nullptr;
  TH1* h_waveform_twrAvg = nullptr;
  TH1* h_waveform_time = nullptr;
  TH1* h_waveform_pedestal = nullptr;
  TH2* h2_hcal_rm = nullptr;
  TH2* h2_hcal_mean = nullptr;
  TH2* h2_hcal_waveform = nullptr;
  TH2* h2_hcal_correlation = nullptr;
  TH1* h_sectorAvg_total = nullptr;
  TH1* h_event = nullptr;
  TH1* h_ntower = nullptr;
  TH1* h1_packet_chans = nullptr;
  TH1* h1_packet_length = nullptr;
  TH1* h1_packet_number = nullptr;
  TH1* h1_packet_event = nullptr;
  TH1* h_rm_sectorAvg[32] = {nullptr};
  TH1* h_rm_tower[24][64] = {nullptr};

  CaloWaveformFitting* WaveformProcessing = nullptr;

  std::vector<runningMean*> rm_vector_sectAvg;
  std::vector<runningMean*> rm_vector_twr;
  std::vector<runningMean*> rm_packet_number;
  std::vector<runningMean*> rm_packet_length;
  std::vector<runningMean*> rm_packet_chans;
};

#endif /* HCAL_HCALMON_H */
