#include "MicromegasGeometry.h"

//__________________________________________________
MicromegasGeometry::MicromegasGeometry():
  m_detnames_sphenix( { "SCO", "SCI", "NCI", "NCO", "SEI", "NEI", "SWI", "NWI" } )
{
  /*
   * to convert sphenix coordinates into a x,y 2D histogram, 
   * we transform z(3D) = x(2D)
   * and x (3D) = y (2D)
   */
  {
    const double tile_x = 0;
    for( const double& tile_z:{ -84.6, -28.2, 28.2, 84.6 } )
    { m_tile_centers.push_back( {tile_z, tile_x} ); }
  }
    
  {
    // neighbor sectors have two modules, separated by 10cm
    for( const double& tile_x: { -m_tile_width - 2, m_tile_width+2 } )
      for( const double& tile_z:{ -37.1, 37.1 } )
    { m_tile_centers.push_back( {tile_z, tile_x} ); }
  }  
}

//__________________________________________________
MicromegasGeometry::point_list_t MicromegasGeometry::get_tile_boundaries( size_t index ) const
{
  const auto center = get_tile_center( index );
  return 
  {    
    { center.first-m_tile_length/2,center.second-m_tile_width/2 },
    { center.first-m_tile_length/2,center.second+m_tile_width/2 },     
    { center.first+m_tile_length/2,center.second+m_tile_width/2 },     
    { center.first+m_tile_length/2,center.second-m_tile_width/2 }     
  };                                                                                       
}

//__________________________________________________
MicromegasGeometry::point_list_t MicromegasGeometry::get_resist_boundaries( size_t tile_index, size_t resist_index, MicromegasDefs::SegmentationType segmentation ) const
{
  const auto tile_center = get_tile_center( tile_index );
  switch( segmentation )
  {
    case MicromegasDefs::SegmentationType::SEGMENTATION_Z:
    {
      const double& resist_width = m_tile_width;
      const double resist_length = m_tile_length/m_nresist;
      const point_t resist_center({ tile_center.first - (m_tile_length - (2*resist_index+1)*resist_length)/2, tile_center.second });
      return 
      {    
        { resist_center.first-resist_length/2,resist_center.second-resist_width/2 },
        { resist_center.first-resist_length/2,resist_center.second+resist_width/2 },     
        { resist_center.first+resist_length/2,resist_center.second+resist_width/2 },     
        { resist_center.first+resist_length/2,resist_center.second-resist_width/2 }     
      };                                                                                       
    }
    case MicromegasDefs::SegmentationType::SEGMENTATION_PHI:
    {
      const double resist_width = m_tile_width/m_nresist;
      const double& resist_length = m_tile_length;
      const point_t resist_center({tile_center.first, tile_center.second - (m_tile_width - (2*resist_index+1)*resist_width)/2 });
      return 
      {    
        { resist_center.first-resist_length/2,resist_center.second-resist_width/2 },
        { resist_center.first-resist_length/2,resist_center.second+resist_width/2 },     
        { resist_center.first+resist_length/2,resist_center.second+resist_width/2 },     
        { resist_center.first+resist_length/2,resist_center.second-resist_width/2 }     
      };                                                                                       
    }
  }

  // unreached
  return point_list_t();
}
