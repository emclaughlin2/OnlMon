#ifndef TPOT_TPOTMONGEOMETRY_H
#define TPOT_TPOTMONGEOMETRY_H

#include <micromegas/MicromegasDefs.h>

#include <cstddef>
#include <string>
#include <vector>
#include <map>

/// handle simplified TPOT 2D geometry, for histogram representation
/**
 * positioning follows the same tile_id order as in coresoftware/simulation/g4simulation/g4micromegas/PHG4MicromegasDetector::setup_tiles
 * it is ordered by tile number as defined in MicromegasMapping.cc
 */
class MicromegasGeometry
{
  public:

  using point_t = std::pair<double, double>;
  using point_list_t = std::vector<point_t>;

  /// constructor
  MicromegasGeometry();
  
  //! get number of defined tiles
  size_t get_ntiles() const 
  { return m_tile_centers.size(); }
  
  //! get tile center for given tile id
  point_t get_tile_center( size_t tile_id )  const 
  { return m_tile_centers[tile_id]; }

  //! get detector name for a given tile id
  std::string get_detname_sphenix( size_t tile_id ) const
  { return m_detnames_sphenix[tile_id]; }
  
  //! get tile boundaries for a given tile id
  point_list_t get_tile_boundaries( size_t /* tile_id */ ) const;

  //! get resist region boundaries for a given tile id
  /** there are 4 resistive layer regions for a given tile, each oriented along the strips */
  point_list_t get_resist_boundaries( size_t /* tile_id */, size_t /* resist_id */, MicromegasDefs::SegmentationType ) const;

  //! tile definitions
  static constexpr double m_tile_length = 54.2; // cm
  static constexpr double m_tile_width = 31.6;  // cm

  /// pitch for z views (cm)
  static constexpr double m_pitch_z = 0.2;
  
  /// pitch for phi views (cm)
  static constexpr double m_pitch_phi = 0.1;

  //! number of resist sectors
  static constexpr int m_nresist = 4;
      
  private:
  
  //! tile centers
  /** tile centers are ordered by tile_ids as defined in MicromegasMapping.cc */
  point_list_t m_tile_centers;

  //! matching detector name
  std::vector<std::string> m_detnames_sphenix;
  
};

#endif
