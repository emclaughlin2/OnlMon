#ifndef TPOT_TPOTMON_H
#define TPOT_TPOTMON_H

#include "MicromegasGeometry.h"

#include <micromegas/MicromegasDefs.h>
#include <micromegas/MicromegasMapping.h>
#include <micromegas/MicromegasCalibrationData.h>
#include <onlmon/OnlMon.h>

#include <array>
#include <memory>

class Event;
class TH1;
class TH2;
class TH2Poly;

class TpotMon : public OnlMon
{
  public:
  TpotMon(const std::string &name);
  ~TpotMon() override = default;

  int process_event(Event *evt) override;
  int Init() override;
  int BeginRun(const int runno) override;
  int Reset() override;

  //! calibration file
  void set_calibration_file( const std::string& value ) 
  { m_calibration_filename = value; }

  //! max sample
  void set_max_sample( int value ) 
  { m_max_sample = value; }
  
  /// set number of RMS sigma used to defined static threshold on a given channel
  void set_n_sigma( double value ) { m_n_sigma = value; }

  /// set minimum ADC value, disregarding pedestal and RMS. This removes channels for which calibration has failed
  void set_min_adc( double value ) { m_min_adc = value; }
  
  // define signal sample window
  using sample_window_t = std::pair<int, int>;
  void set_sample_window_signal( const sample_window_t& value ) 
  { m_sample_window_signal = value; }

  private:
    
  //! setup bins in a TH2Poly. One bin per detector
  void setup_detector_bins( TH2Poly* );
  
  //! setup bins in TH2Poly, one bin per Resist sector
  void setup_resist_bins( TH2Poly*, MicromegasDefs::SegmentationType );
  
  //! keep track of number of events
  int m_evtcnt = 0;
  
  //! keep track of number of 'full' triggers
  /* 
   * it is estimated by summing the number of recorded waveforms/max_waveform/trigger
   * this will break when zero suppression is implemented
   */
  double m_fullevtcnt = 0;
  
  //! mapping
  MicromegasMapping m_mapping;
    
  //! geometry
  MicromegasGeometry m_geometry;

  //! calibration filename
  std::string m_calibration_filename;

  //! calibration data
  MicromegasCalibrationData m_calibration_data;

  //! max sample
  int m_max_sample = 100;

  // sample window
  sample_window_t m_sample_window_signal = {20, 40};

  //! number of RMS sigma used to define threshold
  double m_n_sigma = 5;

  //! minimum ADC value, disregarding pedestal and RMS. This removes channels for which calibration has failed
  double m_min_adc = 50;
  
  //! counter
  TH1* m_counters = nullptr;
   
  //! TPOT per/detector multiplicity
  TH2Poly* m_detector_multiplicity_z = nullptr;
  TH2Poly* m_detector_multiplicity_phi = nullptr;
  
  //! TPOT per/detector occupancy
  TH2Poly* m_detector_occupancy_z = nullptr;
  TH2Poly* m_detector_occupancy_phi = nullptr;

  //! TPOT per/detector multiplicity
  TH2Poly* m_resist_multiplicity_z = nullptr;
  TH2Poly* m_resist_multiplicity_phi = nullptr;
  
  //! TPOT per/detector occupancy
  TH2Poly* m_resist_occupancy_z = nullptr;
  TH2Poly* m_resist_occupancy_phi = nullptr;

  //@name per detector structure
  //@{
  class detector_histograms_t 
  {
    public:

    /// counts avoce threshold vs sample id in each detector
    TH1* m_counts_sample = nullptr;

    /// adc counts vs sample id in each detector
    TH2* m_adc_sample = nullptr;
    
    /// total charge per hit in each detector
    TH1* m_hit_charge = nullptr;
    
    /// number of hits per event in each detector
    TH1* m_hit_multiplicity = nullptr;
    
    /// m_hit_vs_channel
    TH1* m_hit_vs_channel = nullptr;    
  };
  //@}
  
  //@name map tile centers (from MicromegasGeometry) to fee_id
  std::map<int, MicromegasGeometry::point_t> m_tile_centers;
  
  //@name map detector histograms to fee id
  std::map<int, detector_histograms_t> m_detector_histograms;
  
};

#endif /* TPOT_TPOTMON_H */
