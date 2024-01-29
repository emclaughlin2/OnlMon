#include "TpotMonDraw.h"
#include "TpotMonDefs.h"

#include <onlmon/OnlMonClient.h>
#include <onlmon/OnlMonDB.h>

#include <TAxis.h>  // for TAxis
#include <TCanvas.h>
#include <TDatime.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TLine.h>
#include <TPad.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TText.h>

#include <cstring>  // for memset
#include <ctime>
#include <fstream>
#include <iostream>  // for operator<<, basic_ostream, basic_os...
#include <sstream>
#include <vector>  // for vector

namespace
{
  //! make canvas editable in creator, and non-editable in destructor
  class CanvasEditor
  {
    public:
    CanvasEditor( TCanvas* cv ):m_cv(cv)
    { if( m_cv ) m_cv->SetEditable(true); }

    ~CanvasEditor()
    // {}
    { if( m_cv ) m_cv->SetEditable(false); }

    private:
    TCanvas* m_cv = nullptr;
  };

  TPad* create_transparent_pad( const std::string& name )
  {
    auto transparent = new TPad( (name+"_transparent").c_str(), "", 0, 0, 1, 1);
    transparent->SetFillStyle(4000);
    transparent->Draw();
    return transparent;
  };

  TPad* get_transparent_pad( TPad* parent, const std::string& name, bool clear = true)
  {
    if( !parent ) return nullptr;
    const std::string transparent_name = name+"_transparent";
    auto out = dynamic_cast<TPad*>( parent->FindObject( transparent_name.c_str() ) );

    if( !out ) std::cout << "get_transparent_pad - " << transparent_name << " not found" << std::endl;
    if( out && clear ) out->Clear("D");
    return out;

  }
 
  // draw an vertical line that extends automatically from both sides of a canvas
  [[maybe_unused]] TLine* vertical_line( TVirtualPad* pad, Double_t x )
  {
    Double_t yMin = pad->GetUymin();
    Double_t yMax = pad->GetUymax();
    
    if( pad->GetLogy() )
    {
      yMin = std::pow( 10, yMin );
      yMax = std::pow( 10, yMax );
    }
    
    return new TLine( x, yMin, x, yMax );
  }
  
  // draw an horizontal line that extends automatically from both sides of a canvas
  [[maybe_unused]] TLine* horizontal_line( TVirtualPad* pad, Double_t y )
  {
    Double_t xMin = pad->GetUxmin();
    Double_t xMax = pad->GetUxmax();
    
    if( pad->GetLogx() )
    {
      xMin = std::pow( 10, xMin );
      xMax = std::pow( 10, xMax );
    }
    
    return new TLine( xMin, y, xMax, y );    
  }
  
  // divide canvas, adjusting canvas positions to leave room for a banner at the top
  void divide_canvas( TCanvas* cv, int ncol, int nrow )
  {
    static constexpr double max_height = 0.94;

    cv->Divide( ncol, nrow );
    for( int i = 0; i < ncol*nrow; ++i )
    {
      auto pad = cv->GetPad( i+1 );
      int col = i%ncol;
      int row = i/ncol;
      const double xmin = double(col)/ncol;
      const double xmax = double(col+1)/ncol;

      const double ymin = max_height*(1. - double(row+1)/nrow);
      const double ymax = max_height*(1. - double(row)/nrow);
      pad->SetPad( xmin, ymin, xmax, ymax );
    }

  }
}

//__________________________________________________________________________________
TpotMonDraw::TpotMonDraw(const std::string &name)
  : OnlMonDraw(name)
{
  // setup default filename for reference histograms
  const auto tpotcalibref = getenv("TPOTCALIBREF");
  if( tpotcalibref )
  {
    m_ref_histograms_filename = std::string(tpotcalibref) + "/" + "Run_00000-TPOTMON_0.root";
    std::cout << "TpotMon::TpotMon - reading reference histograms from: " << m_ref_histograms_filename << std::endl;
    m_ref_histograms_tfile.reset( TFile::Open( m_ref_histograms_filename.c_str(), "READ" ) );
  } else {
    m_ref_histograms_filename = "Run_00000-TPOTMON_0.root";
    std::cout << "TpotMon::TpotMon - TPOTCALIBREF environment variable not set. Reading reference histograms from: " << m_ref_histograms_filename << std::endl;
    m_ref_histograms_tfile.reset( TFile::Open( m_ref_histograms_filename.c_str(), "READ" ) );
  }
  
  // this TimeOffsetTicks is neccessary to get the time axis right
  TDatime T0(2003, 01, 01, 00, 00, 00);
  TimeOffsetTicks = T0.Convert();
  dbvars.reset( new OnlMonDB(ThisName) );

  // initialize local list of detector names
  for( const auto& fee_id:m_mapping.get_fee_id_list() )
  { m_detnames_sphenix.push_back( m_mapping.get_detname_sphenix( fee_id ) ); }

}

//__________________________________________________________________________________
int TpotMonDraw::Init()
{ return 0; }

//__________________________________________________________________________________
int TpotMonDraw::DrawDeadServer( TPad* pad )
{
  if( !pad ) 
  {
    if( Verbosity() ) std::cout << "TpotMonDraw::DrawDeadServer - invalid pad" << std::endl;
    return 0;
  } else {  
    pad->SetPad(0,0,1,1);
    return OnlMonDraw::DrawDeadServer( pad );
  }
}

//__________________________________________________________________________________
void TpotMonDraw::draw_time( TPad* pad )
{
  if( !pad ) 
  {
    if( Verbosity() ) std::cout << "TpotMonDraw::draw_time - invalid pad" << std::endl;
    return;
  }
  
  pad->SetPad( 0, 0.94, 1, 1 );
  pad->Clear();
  TText PrintRun;
  PrintRun.SetTextFont(62);
  PrintRun.SetTextSize(0.6);
  PrintRun.SetNDC();          // set to normalized coordinates
  PrintRun.SetTextAlign(23);  // center/top alignment

  std::ostringstream runnostream;
  auto cl = OnlMonClient::instance();
  time_t evttime = cl->EventTime("CURRENT");

  runnostream
    << ThisName << " Run " << cl->RunNumber()
    << ", Time: " << ctime(&evttime);

  pad->cd();
  PrintRun.DrawText(0.5, 0.5, runnostream.str().c_str());
}

//__________________________________________________________________________________
TCanvas* TpotMonDraw::get_canvas(const std::string& name, bool clear )
{
  auto cv = dynamic_cast<TCanvas*>( gROOT->FindObject( name.c_str() ) );
  if( !cv ) cv = create_canvas( name );
  if( cv && clear ) cv->Clear("D");
  return cv;
}

//__________________________________________________________________________________
TCanvas* TpotMonDraw::create_canvas(const std::string &name)
{

  if( Verbosity() )
  { std::cout << "TpotMonDraw::create_canvas - name: " << name << std::endl; }

  OnlMonClient *cl = OnlMonClient::instance();
  int xsize = cl->GetDisplaySizeX();
  int ysize = cl->GetDisplaySizeY();

  if (name == "TPOT_counters")
  {

    // xpos (-1) negative: do not draw menu bar
    auto cv = new TCanvas(name.c_str(), "TPOT event counters", -1, 0, xsize / 2, ysize);
    gSystem->ProcessEvents();
    divide_canvas( cv, 1, 1 );
    create_transparent_pad(name);
    cv->SetEditable(false);
    m_canvas.push_back( cv );
    return cv;
  } else if (name == "TPOT_detector_occupancy") {

    // xpos (-1) negative: do not draw menu bar
    auto cv = new TCanvas(name.c_str(), "TPOT detector occupancy", -1, 0, xsize / 2, ysize);
    gSystem->ProcessEvents();
    divide_canvas( cv, 1, 2 );
    create_transparent_pad(name);
    cv->SetEditable(false);
    m_canvas.push_back( cv );
    return cv;

  } else if (name == "TPOT_resist_occupancy") {

    // xpos (-1) negative: do not draw menu bar
    auto cv = new TCanvas(name.c_str(), "TPOT resist occupancy", -1, 0, xsize / 2, ysize);
    gSystem->ProcessEvents();
    divide_canvas( cv, 1, 2 );
    create_transparent_pad(name);
    cv->SetEditable(false);
    m_canvas.push_back( cv );
    return cv;

  } else if (name == "TPOT_adc_vs_sample") {

    auto cv = new TCanvas(name.c_str(), "TpotMon adc vs sample", -1, 0, xsize / 2, ysize);
    gSystem->ProcessEvents();
    divide_canvas( cv, 4, 4 );
    create_transparent_pad(name);
    for( int i = 0; i < 16; ++i )
    {
      cv->GetPad(i+1)->SetLeftMargin(0.15);
      cv->GetPad(i+1)->SetRightMargin(0.02);
    }
    cv->SetEditable(false);
    m_canvas.push_back( cv );
    return cv;

  } else if (name == "TPOT_counts_vs_sample") {

    auto cv = new TCanvas(name.c_str(), "TpotMon counts vs sample", -1, 0, xsize / 2, ysize);
    gSystem->ProcessEvents();
    divide_canvas( cv, 4, 4 );
    create_transparent_pad(name);
    for( int i = 0; i < 16; ++i )
    {
      auto&& pad = cv->GetPad(i+1);
      pad->SetLeftMargin(0.15);
      pad->SetRightMargin(0.02);    
    }
        
    cv->SetEditable(false);
    m_canvas.push_back( cv );
    return cv;

  } else if (name == "TPOT_hit_charge") {

    auto cv = new TCanvas(name.c_str(), "TpotMon hit charge", -1, 0, xsize / 2, ysize);
    gSystem->ProcessEvents();
    divide_canvas( cv, 4, 4 );
    create_transparent_pad(name);
    cv->SetEditable(false);
    m_canvas.push_back( cv );
    return cv;

  } else if (name == "TPOT_hit_multiplicity") {

    auto cv = new TCanvas(name.c_str(), "TpotMon hit multiplicity", -1, 0, xsize / 2, ysize);
    gSystem->ProcessEvents();
    divide_canvas( cv, 4, 4 );
    create_transparent_pad(name);
    cv->SetEditable(false);
    m_canvas.push_back( cv );
    return cv;

  } else if (name == "TPOT_hit_vs_channel") {

    auto cv = new TCanvas(name.c_str(), "TpotMon hit vs channel", -1, 0, xsize / 2, ysize);
    gSystem->ProcessEvents();
    divide_canvas( cv, 4, 4 );
    create_transparent_pad(name);
    cv->SetEditable(false);
    m_canvas.push_back( cv );
    return cv;

  }
  return nullptr;
}

//_______________________________________________________________________________
int TpotMonDraw::Draw(const std::string &what)
{
  if( Verbosity() )
  { std::cout << "TpotMonDraw::Draw - what: " << what << std::endl; }

  int iret = 0;
  int idraw = 0;

  {
    // get counters
    const auto m_counters = get_histogram( "m_counters");
    if( m_counters && Verbosity() )
    {
      const int events = m_counters->GetBinContent( TpotMonDefs::kEventCounter );
      const int valid_events = m_counters->GetBinContent( TpotMonDefs::kValidEventCounter );
      std::cout << "TpotMonDraw::Draw - events: " << events << " valid events: " << valid_events << std::endl;
    }
  }

  if( what == "ALL" || what == "TPOT_counters" )
  {
    iret += draw_counters();
    ++idraw;
  }

  if( what == "ALL" || what == "TPOT_detector_occupancy" )
  {
    iret += draw_detector_occupancy();
    ++idraw;
  }

  if( what == "ALL" || what == "TPOT_resist_occupancy" )
  {
    iret += draw_resist_occupancy();
    ++idraw;
  }

  if (what == "ALL" || what == "TPOT_adc_vs_sample")
  {
    iret += draw_array("TPOT_adc_vs_sample", get_histograms( "m_adc_sample" ), DrawOptions::Colz );
    auto cv = get_canvas("TPOT_adc_vs_sample");
    if( cv )
    {
      CanvasEditor cv_edit(cv);
      cv->Update();
      for( int i = 0; i < 16; ++i )
      {
        // draw vertical lines that match sample window
        auto&& pad = cv->GetPad(i+1);
        pad->cd();
        for( const auto line:{vertical_line( pad, m_sample_window_signal.first ), vertical_line( pad, m_sample_window_signal.second ) } )
        {
          line->SetLineStyle(2);
          line->SetLineColor(2);
          line->SetLineWidth(2);
          line->Draw();
        }
      }
    }
    ++idraw;
  }

  if (what == "ALL" || what == "TPOT_counts_vs_sample")
  {
    iret += draw_array("TPOT_counts_vs_sample", get_histograms( "m_counts_sample" ), get_ref_histograms_scaled( "m_counts_sample" ) );
    auto cv = get_canvas("TPOT_counts_vs_sample");
    if( cv )
    {
      std::cout << "TpotMonDraw::Draw - draw vertical lines" << std::endl;
      CanvasEditor cv_edit(cv);
      cv->Update();
      for( int i = 0; i < 16; ++i )
      {
        // draw vertical lines that match sample window
        auto&& pad = cv->GetPad(i+1);
        pad->cd();
        for( const auto line:{vertical_line( pad, m_sample_window_signal.first ), vertical_line( pad, m_sample_window_signal.second ) } )
        {
          line->SetLineStyle(2);
          line->SetLineColor(2);
          line->SetLineWidth(2);
          line->Draw();
        }
      }
    }
    ++idraw;
  }

  if (what == "ALL" || what == "TPOT_hit_charge")
  {
    iret += draw_array("TPOT_hit_charge", get_histograms( "m_hit_charge" ), get_ref_histograms_scaled( "m_hit_charge" ), DrawOptions::Logy );
    ++idraw;
  }

  if (what == "ALL" || what == "TPOT_hit_multiplicity")
  {
    iret += draw_array("TPOT_hit_multiplicity", get_histograms( "m_hit_multiplicity" ), get_ref_histograms_scaled( "m_hit_multiplicity" ), DrawOptions::Logy );
    ++idraw;
  }

  if (what == "ALL" || what == "TPOT_hit_vs_channel")
  {
    iret += draw_array("TPOT_hit_vs_channel", get_histograms( "m_hit_vs_channel" ), get_ref_histograms_scaled( "m_hit_vs_channel" ) );
    auto cv = get_canvas("TPOT_hit_vs_channel");
    if( cv )
    {
      CanvasEditor cv_edit(cv);
      cv->Update();
      for( int i = 0; i < 16; ++i )
      {
        // draw vertical lines that match sample window
        auto&& pad = cv->GetPad(i+1);
        pad->cd();
        for( const int& channel:{64, 128, 196} ) 
        {
          const auto line = vertical_line( pad, channel );
          line->SetLineStyle(2);
          line->SetLineColor(1);
          line->SetLineWidth(1);
          line->Draw();
        }
      }
    }
    ++idraw;
  }

  if (!idraw)
  {
    std::cout << "TpotMonDraw::Draw - Unimplemented Drawing option: " << what << std::endl;
    iret = -1;
  }
  return iret;
}

//_______________________________________________________________________________
int TpotMonDraw::SavePlot(const std::string &what, const std::string &type)
{
  auto cl = OnlMonClient::instance();
  const int iret = Draw(what);
  if (iret) return iret;

  for( std::size_t i =0; i < m_canvas.size(); ++i )
  {
    const auto& cv = m_canvas[i];
    if( cv )
    {
    std::string filename = ThisName + "_" + std::to_string(i+1) + "_" +
      std::to_string(cl->RunNumber()) + "." + type;
     cl->CanvasToPng(cv, filename);
    }
  }

  return 0;
}

//__________________________________________________________________________________
int TpotMonDraw::MakeHtml(const std::string &what)
{
  const int iret = Draw(what);
  if (iret) return iret;

  auto cl = OnlMonClient::instance();

  // Register the 1st canvas png file to the menu and produces the png file.
  for( std::size_t i =0; i < m_canvas.size(); ++i )
  {
    const auto& cv = m_canvas[i];
    if( cv )
    {
      const auto pngfile = cl->htmlRegisterPage(*this, cv->GetName(), Form("%lu", i+1), "png");
      cl->CanvasToPng(cv, pngfile);
    }
  }

  // log
  {
    const std::string logfile = cl->htmlRegisterPage(*this, "EXPERTS/Log", "log", "html");
    std::ofstream out(logfile.c_str());
    out
      << "<HTML><HEAD><TITLE>Log file for run " << cl->RunNumber()
      << "</TITLE></HEAD>"
      << std::endl
      << "<P>Some log file output would go here."
      << std::endl;
    out.close();
  }

  // status
  {
    const std::string status = cl->htmlRegisterPage(*this, "EXPERTS/Status", "status", "html");
    std::ofstream out(status.c_str());
    out
      << "<HTML><HEAD><TITLE>Status file for run " << cl->RunNumber()
      << "</TITLE></HEAD>"
      << std::endl
      << "<P>Some status output would go here."
      << std::endl;
    out.close();
    cl->SaveLogFile(*this);
  }

  return 0;
}

//__________________________________________________________________________________
int TpotMonDraw::draw_counters()
{

  if( Verbosity() ) std::cout << "TpotMonDraw::draw_counters" << std::endl;

  // get histograms
  auto m_counters =  get_histogram( "m_counters");  
  std::unique_ptr<TH1> m_counters_ref( normalize( get_ref_histogram( "m_counters" ), get_ref_scale_factor() ) );

  auto cv = get_canvas("TPOT_counters");
  auto transparent = get_transparent_pad( cv, "TPOT_counters");
  if( !cv )
  {
    if( Verbosity() ) std::cout << "TpotMonDraw::draw_counters - no canvas" << std::endl;
    return -1;
  }

  CanvasEditor cv_edit(cv);

  if( m_counters )
  {
    cv->cd(1);
    gPad->SetLeftMargin( 0.07 );
    gPad->SetRightMargin( 0.15 );
    m_counters->SetFillStyle(1001);
    m_counters->SetFillColor(kYellow );
    m_counters->DrawCopy( "h" );
    
    if( m_counters_ref ) 
    {
      m_counters_ref->SetLineColor(2);
      m_counters_ref->DrawCopy( "hist same" );
    }
    
    draw_time(transparent);
    return 0;
  } else {
    
    DrawDeadServer(transparent);
    return -1;
  }
}

//__________________________________________________________________________________
int TpotMonDraw::draw_detector_occupancy()
{

  if( Verbosity() ) std::cout << "TpotMonDraw::draw_detector_occupancy" << std::endl;

  // get histograms
  auto m_detector_occupancy_phi =  get_histogram( "m_detector_occupancy_phi");
  auto m_detector_occupancy_z =  get_histogram( "m_detector_occupancy_z");

  // turn off stat panel
  for( const auto& h:{m_detector_occupancy_phi,m_detector_occupancy_z} )
  { if(h) h->SetStats(0); }
  
  auto cv = get_canvas("TPOT_detector_occupancy");
  auto transparent = get_transparent_pad( cv, "TPOT_detector_occupancy");
  if( !cv )
  {
    if( Verbosity() ) std::cout << "TpotMonDraw::draw_detector_occupancy - no canvas" << std::endl;
    return -1;
  }

  CanvasEditor cv_edit(cv);

  if( m_detector_occupancy_phi && m_detector_occupancy_z )
  {
    cv->cd(1);
    gPad->SetLeftMargin( 0.07 );
    gPad->SetRightMargin( 0.15 );
    m_detector_occupancy_z->DrawCopy( "colz" );
    draw_detnames_sphenix( "Z" );

    cv->cd(2);
    gPad->SetLeftMargin( 0.07 );
    gPad->SetRightMargin( 0.15 );
    m_detector_occupancy_phi->DrawCopy( "colz" );
    draw_detnames_sphenix( "P" );

    draw_time(transparent);
    return 0;

  } else {

    DrawDeadServer(transparent);
    return -1;

  }
}

//__________________________________________________________________________________
int TpotMonDraw::draw_resist_occupancy()
{

  if( Verbosity() ) std::cout << "TpotMonDraw::draw_resist_occupancy" << std::endl;

  // get histograms
  auto m_resist_occupancy_phi =  get_histogram( "m_resist_occupancy_phi");
  auto m_resist_occupancy_z =  get_histogram( "m_resist_occupancy_z");

  // turn off stat panel
  for( const auto& h:{m_resist_occupancy_phi,m_resist_occupancy_z} )
  { if(h) h->SetStats(0); }
    
  auto cv = get_canvas("TPOT_resist_occupancy");
  auto transparent = get_transparent_pad( cv, "TPOT_resist_occupancy");
  if( !cv )
  {
    if( Verbosity() ) std::cout << "TpotMonDraw::draw_resist_occupancy - no canvas" << std::endl;
    return -1;
  }

  CanvasEditor cv_edit(cv);

  if( m_resist_occupancy_phi && m_resist_occupancy_z )
  {
    cv->cd(1);
    gPad->SetLeftMargin( 0.07 );
    gPad->SetRightMargin( 0.15 );
    m_resist_occupancy_z->DrawCopy( "colz" );
    draw_detnames_sphenix( "Z" );

    cv->cd(2);
    gPad->SetLeftMargin( 0.07 );
    gPad->SetRightMargin( 0.15 );
    m_resist_occupancy_phi->DrawCopy( "colz" );
    draw_detnames_sphenix( "P" );

    draw_time(transparent);
    return 0;

  } else {
    DrawDeadServer(transparent);
    return -1;
  }
}

//__________________________________________________________________________________
void TpotMonDraw::draw_detnames_sphenix( const std::string& suffix)
{
  gPad->Update();
  for( size_t i = 0; i < m_geometry.get_ntiles(); ++i )
  {
    const auto name = m_geometry.get_detname_sphenix(i)+suffix;
    const auto [x,y] = m_geometry.get_tile_center(i);
    auto text = new TText();
    text->DrawText( x-0.8*m_geometry.m_tile_length/2, y-0.8*m_geometry.m_tile_width/2, name.c_str() );
    text->Draw();
  }

}

//__________________________________________________________________________________
TH1* TpotMonDraw::get_histogram( const std::string& name ) const
{
  auto cl = OnlMonClient::instance();
  return cl->getHisto("TPOTMON_0", name );
}

//__________________________________________________________________________________
TpotMonDraw::histogram_array_t TpotMonDraw::get_histograms( const std::string& name ) const
{
  histogram_array_t out{{nullptr}};
  for( size_t i=0; i<m_detnames_sphenix.size(); ++i)
  {
    const auto& detector_name=m_detnames_sphenix[i];
    const auto hname = name + "_" + detector_name;
    out[i] =  get_histogram(  hname );
    if( Verbosity() )
    { std::cout << "TpotMonDraw::get_histograms - " << hname << (out[i]?" found":" not found" ) << std::endl; }
  }

  return out;
}

//__________________________________________________________________________________
TH1* TpotMonDraw::get_ref_histogram( const std::string& name ) const
{ return m_ref_histograms_tfile ? static_cast<TH1*>( m_ref_histograms_tfile->Get( name.c_str() ) ):nullptr; }

//__________________________________________________________________________________
TpotMonDraw::histogram_array_t TpotMonDraw::get_ref_histograms( const std::string& name ) const
{
  histogram_array_t out{{nullptr}};
  for( size_t i=0; i<m_detnames_sphenix.size(); ++i)
  {
    const auto& detector_name=m_detnames_sphenix[i];
    const auto hname = name + "_" + detector_name;
    out[i] =  get_ref_histogram(  hname );
    if( Verbosity() )
    { std::cout << "TpotMonDraw::get_ref_histograms - " << hname << (out[i]?" found":" not found" ) << std::endl; }
  }

  return out;
}

//__________________________________________________________________________________
TpotMonDraw::histogram_array_t TpotMonDraw::get_ref_histograms_scaled( const std::string& name ) const
{
  histogram_array_t source( get_ref_histograms( name ) );
  histogram_array_t out{{nullptr}};
  
  const double scale = get_ref_scale_factor();
  for( size_t i=0; i<source.size(); ++i)
  { if( source[i] ) out[i]=normalize( source[i], scale ); } 

  return out;
}

//__________________________________________________________________________________
double TpotMonDraw::get_ref_scale_factor() const
{
  if( !m_ref_histograms_tfile ) return 0;
  const auto m_counters = get_histogram( "m_counters");
  const auto m_counters_ref = get_ref_histogram( "m_counters");
  if( !( m_counters && m_counters_ref ) ) return 0;
  
  const double full_events = m_counters->GetBinContent( TpotMonDefs::kFullEventCounter );
  const double full_events_ref = m_counters_ref->GetBinContent( TpotMonDefs::kFullEventCounter );
  return full_events_ref > 0 ? full_events/full_events_ref : 0;
}

//__________________________________________________________________________________
TH1* TpotMonDraw::normalize( TH1* source, double scale ) const
{
  if( !source ) return nullptr;
  auto destination = static_cast<TH1*>( source->Clone() );
  destination->SetName( TString( source->GetName() )+"_scaled" );
  destination->Scale( scale );
  return destination;
}

//__________________________________________________________________________________
int TpotMonDraw::draw_array( const std::string& name, const TpotMonDraw::histogram_array_t& histograms, const TpotMonDraw::histogram_array_t& ref_histograms, unsigned int options )
{
  if( Verbosity() ) std::cout << "TpotMonDraw::draw_array - name: " << name << std::endl;

  auto cv = get_canvas(name);
  auto transparent = get_transparent_pad( cv, name);
  if( !cv ) return -1;

  bool drawn = false;
  CanvasEditor cv_edit(cv);
  for( size_t i = 0; i < histograms.size(); ++i )
  {
    if( histograms[i] )
    {
      cv->cd(i+1);
      if( options&DrawOptions::Colz ) histograms[i]->DrawCopy( "col" );
      else {
        
        histograms[i]->SetFillStyle(1001);
        histograms[i]->SetFillColor(kYellow );
        histograms[i]->DrawCopy( "h" );
        histograms[i]->DrawCopy();
      
      }
      
      // also draw reference
      if( ref_histograms[i] )
      {
        ref_histograms[i]->SetLineColor(2);
        ref_histograms[i]->Draw("hist same" );
      }
      
      gPad->SetBottomMargin(0.12);
      if( options&DrawOptions::Logx ) gPad->SetLogx( true );
      if( options&DrawOptions::Logy && histograms[i]->GetEntries() > 0 ) gPad->SetLogy( true );
      if( options&DrawOptions::Logz ) gPad->SetLogz( true );
      drawn = true;
    }
  }

  if( drawn )
  {
    draw_time(transparent);
    return 0;
  } else {
    DrawDeadServer(transparent);
    return -1;
  }
  
  // need to delete reference histograms to avoid leak
  for( auto h:ref_histograms ) { delete h; }
  
}
