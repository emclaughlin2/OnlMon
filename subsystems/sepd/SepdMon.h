#ifndef SEPD_SEPDMON_H
#define SEPD_SEPDMON_H

#include <onlmon/OnlMon.h>
#include <cmath>
#include <vector>

class CaloWaveformFitting;
class TowerInfoContainer;
class Event;
class TH1;
class TH2;
class Packet;
class runningMean;

class SepdMon : public OnlMon
{
 public:
  SepdMon(const std::string &name);
  virtual ~SepdMon();

  int process_event(Event *evt);
  int Init();
  int BeginRun(const int runno);
  int SepdMapChannel(int ch);
  int Reset();

 protected:
  std::vector<float> getSignal(Packet *p, const int channel);
  std::vector<float> anaWaveformFast(Packet *p, const int channel);
  std::vector<float> anaWaveformTemp(Packet *p, const int channel);
  int evtcnt = 0;
  int idummy = 0;

  const int Nsector = 24;
  const int Nchannel = 192 * 4;
  const int packetlow = 9001;
  const int packethigh = 9006;
  const int m_nChannels = 192;

  const int packet_depth = 1000;
  const float hit_threshold = 10;
  const int n_samples_show = 31;

  const int nChannels = 744;
  int nPhi0 = 12;
  int nPhi = 24;
  int nRad = 24;
  double axislimit = M_PI;

  TH2 *h_ADC0_s = nullptr;
  TH2 *h_ADC0_n = nullptr;
  TH2 *h_ADC_s = nullptr;
  TH2 *h_ADC_n = nullptr;

  TH2 *h_hits0_s = nullptr;
  TH2 *h_hits0_n = nullptr;
  TH2 *h_hits_s = nullptr;
  TH2 *h_hits_n = nullptr;

  TH2 *h_ADC_corr = nullptr;
  TH2 *h_hits_corr = nullptr;

  TH1 *h1_waveform_twrAvg = nullptr;
  TH1 *h1_waveform_time = nullptr;
  TH1 *h1_waveform_pedestal = nullptr;
  TH1 *h_event = nullptr;
  TH2 *h2_sepd_waveform = nullptr;

  TH1 *h1_sepd_fitting_sigDiff = nullptr;
  TH1 *h1_sepd_fitting_pedDiff = nullptr;
  TH1 *h1_sepd_fitting_timeDiff = nullptr;

  TH1 *h1_packet_chans = nullptr;
  TH1 *h1_packet_length = nullptr;
  TH1 *h1_packet_number = nullptr;
  TH1 *h1_packet_event = nullptr;

  TH1 *h_ADC_channel[744] = {nullptr};

  std::string runtypestr = "Unknown";
  std::string id_string;

  CaloWaveformFitting *WaveformProcessingFast = nullptr;
  CaloWaveformFitting *WaveformProcessingTemp = nullptr;

  std::vector<runningMean *> rm_packet_number;
  std::vector<runningMean *> rm_packet_length;
  std::vector<runningMean *> rm_packet_chans;
};

#endif /* SEPD_SEPDMON_H */
