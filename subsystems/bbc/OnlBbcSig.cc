#include "OnlBbcSig.h"
//#include "RunningStats.h"

#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TH2.h>
//#include <THnSparse.h>
#include <TSpline.h>
#include <TPad.h>
#include <TMath.h>
#include <TGraphErrors.h>
#include <TSpectrum.h>

#include <iostream>
#include <fstream>
#include <limits>

using namespace std;


OnlBbcSig::OnlBbcSig(const int chnum, const int nsamp) :
  ch{chnum},
  nsamples{nsamp},
  f_ampl{0},
  f_time{0},
  f_time_offset{4.0}, // time shift from fit
  f_integral{0.},     // time shift from fit
  hRawPulse{nullptr},
  hSubPulse{nullptr},
  hpulse{nullptr},
  gRawPulse{nullptr},
  gSubPulse{nullptr},
  gpulse{nullptr},
  hPed0{nullptr},
  ped0{0},            // ped average
  ped0rms{0},
  use_ped0{0},
  minped0samp{-9999},
  maxped0samp{-9999},
  minped0x{0.},
  maxped0x{0.},
  time_calib{0},
  h2Template{nullptr},
  h2Residuals{nullptr},
  // range of good amplitudes for templates
  // units are usually in ADC counts
  hAmpl{nullptr},
  hTime{nullptr},
  template_npointsx{0},
  template_npointsy{0},
  template_begintime{0},
  template_endtime{0},
  template_min_good_amplitude{20.},
  template_max_good_amplitude{4080},
  template_min_xrange{0},
  template_max_xrange{0},
  template_fcn{nullptr},
  verbose{0}
{
  //cout << "In OnlBbcSig::OnlBbcSig(" << ch << "," << nsamples << ")" << endl;

}


void OnlBbcSig::Init()
{
  TString name;

  name = "hrawpulse"; name += ch;
  hRawPulse = new TH1F(name,name,nsamples,-0.5,nsamples-0.5);
  name = "hsubpulse"; name += ch;
  hSubPulse = new TH1F(name,name,nsamples,-0.5,nsamples-0.5);

  //gRawPulse = new TGraphErrors(nsamples);
  gRawPulse = new TGraphErrors();
  name = "grawpulse"; name += ch;
  gRawPulse->SetName(name);
  //gSubPulse = new TGraphErrors(nsamples);
  gSubPulse = new TGraphErrors();
  name = "gsubpulse"; name += ch;
  gSubPulse->SetName(name);

  hpulse = hRawPulse;   // hpulse,gpulse point to raw by default
  gpulse = gRawPulse;   // we switch to sub for default if ped is applied

  //ped0stats = new RunningStats();
  name = "hPed0_"; name += ch;
  hPed0 = new TH1F(name,name,16384,-0.5,16383.5);
  //hPed0 = new TH1F(name,name,10000,1,0); // automatically determine the range

  SetTemplateSize(120,2048,-2,9.9);
}

void  OnlBbcSig::SetTemplateSize(const Int_t nptsx, const Int_t nptsy, const Double_t begt, const Double_t endt)
{
  template_npointsx = nptsx;
  template_npointsy = nptsy;
  template_begintime = begt;
  template_endtime = endt;

  template_y.resize(template_npointsx);
  template_yrms.resize(template_npointsx);

  Double_t xbinwid = (template_endtime - template_begintime)/(template_npointsx-1);
  Double_t ybinwid = (1.1+0.1)/template_npointsy;  // yscale... should we vary this?
  if ( h2Template ) delete h2Template;
  if ( h2Residuals ) delete h2Residuals;

  TString name = "h2Template"; name += ch;
  h2Template = new TH2F(name,name,template_npointsx,template_begintime-xbinwid/2.,template_endtime+xbinwid/2,
      template_npointsy,-0.1+ybinwid/2.0,1.1+ybinwid/2.0);
 
  name = "h2Residuals"; name += ch;
  h2Residuals = new TH2F(name,name,template_npointsx,template_begintime-xbinwid/2.,template_endtime+xbinwid/2,
      80,-20,20);
 
  /*
  int nbins[] = { template_npointsx, nbinsy };
  Double_t lowrange[] = { template_begintime-xbinwid/2.0, -0.1+ybinwid/2.0 };
  Double_t highrange[] = { template_endtime+xbinwid/2.0, 1.1+ybinwid/2.0 };
  h2Template = new THnSparseF(name,name,2,nbins,lowrange,highrange);
  */
  //h2Template->cd( gDirectory );

}

OnlBbcSig::~OnlBbcSig()
{
  delete hRawPulse;
  delete hSubPulse;
  delete gRawPulse;
  delete gSubPulse;
  //delete ped0stats;
  delete hPed0;
  //h2Template->Write();
  delete h2Template;
  delete h2Residuals;
  delete hAmpl;
  delete hTime;
  delete template_fcn;

}

void  OnlBbcSig::SetTemplateMinMaxGoodADC(const Double_t min, const Double_t max)
{
  template_min_good_amplitude = min;
  template_max_good_amplitude = max;
}

void  OnlBbcSig::SetTemplateMinMaxFitRange(const Double_t min, const Double_t max)
{
  template_min_xrange = min;
  template_max_xrange = max;
}

// This sets y, and x to sample number (starts at 0)
void OnlBbcSig::SetY(const Float_t *y, const int invert)
{
  if ( hRawPulse == nullptr )
  {
    Init();
  }

  hpulse->Reset();
  f_ampl = -9999.;
  f_time = -9999.;

  for (int isamp=0; isamp<nsamples; isamp++)
  {
    hRawPulse->SetBinContent( isamp+1, y[isamp] );
    gRawPulse->SetPoint( isamp, Double_t(isamp), y[isamp] );
  }

  // Apply pedestal
  if ( use_ped0 != 0 || minped0samp >= 0 || minped0x != maxped0x )
  {
    //cout << "sub" << endl;

    if ( minped0samp >= 0 )
    {
      CalcEventPed0(minped0samp,maxped0samp);
    }
    else if ( minped0x != maxped0x )
    {
      CalcEventPed0(minped0x,maxped0x);
    }

    for (int isamp=0; isamp<nsamples; isamp++)
    {
      hSubPulse->SetBinContent( isamp+1, invert*(y[isamp]-ped0) );
      hSubPulse->SetBinError( isamp+1, ped0rms );
      gSubPulse->SetPoint( isamp, (Double_t)isamp, invert*(y[isamp]-ped0) );
      gSubPulse->SetPointError( isamp, 0., ped0rms );
    }
  }
}

void OnlBbcSig::SetXY(const Float_t *x, const Float_t *y, const int invert)
{
  if ( hRawPulse == nullptr )
  {
    Init();
  }

  hRawPulse->Reset();
  hSubPulse->Reset();

  f_ampl = -9999.;
  f_time = -9999.;

  //cout << "nsamples " << nsamples << endl;
  //cout << "use_ped0 " << use_ped0 << "\t" << ped0 << endl;

  for( int isamp=0; isamp<nsamples; isamp++ )
  {
    //cout << "aaa\t" << isamp << "\t" << x[isamp] << "\t" << y[isamp] << endl;
    hRawPulse->SetBinContent( isamp+1, y[isamp] );
    gRawPulse->SetPoint( isamp, x[isamp], y[isamp] );
  }

  if ( use_ped0 != 0 || minped0samp >= 0 || minped0x != maxped0x )
  {
    if ( minped0samp >= 0 )
    {
      CalcEventPed0(minped0samp,maxped0samp);
    }
    else if ( minped0x != maxped0x )
    {
      CalcEventPed0(minped0x,maxped0x);
    }

    for (int isamp=0; isamp<nsamples; isamp++)
    {
      {
        // How do we handle data which is not in samples, but is in time,
        // such as DRS4 data
        //cout << "bbb\t" << isamp << "\t" << x[isamp] << "\t" << invert*(y[isamp]-ped0) << endl;
        hSubPulse->SetBinContent( isamp+1, invert*(y[isamp]-ped0) );
        hSubPulse->SetBinError( isamp+1, ped0rms );
        gSubPulse->SetPoint( isamp, x[isamp], invert*(y[isamp]-ped0) );
        gSubPulse->SetPointError( isamp, 0., ped0rms );
      }
    }
  }
}

Double_t OnlBbcSig::GetSplineAmpl()
{
  TSpline3 s3("s3",gSubPulse);

  // First find maximum, to rescale
  f_ampl = -999999.;
  double step_size = 0.01;
  //cout << "step size " << step_size << endl;
  for (double ix=0; ix<nsamples; ix += step_size)
  {
    Double_t val = s3.Eval(ix);
    if ( val > f_ampl )
    {
      f_ampl = val;
    }
  }

  return f_ampl;
}

// This does a straight line fit for now...
/*
Double_t OnlBbcSig::FitPulse()
{
  const Double_t pedcut[] = {1650,1560};
  const Double_t maxcut[] = {12000,12600};

  TF1 f("f","pol1",0,31);
  f.SetParameter(0,0);
  f.SetParameter(1,3000.);

  Double_t start = 0;
  Double_t stop = 31;
  Double_t x, y;
  for (int isamp=0; isamp<nsamples; isamp++)
  {
    gpulse->GetPoint(isamp,x,y);
    if ( y>pedcut[ch] )
    {
      start = x;
      break;
    }
  }
  
  for (int isamp=(int)start; isamp<nsamples; isamp++)
  {
    gpulse->GetPoint(isamp,x,y);
    if ( y>maxcut[ch] )
    {
      stop = x-1;
      break;
    }
  }

  f.SetRange(start,stop);
  gpulse->Fit(&f,"R");

  Double_t slope = f.GetParameter(1);

  //cout << "xxx " << slope << endl;
  return slope;
}
*/

void OnlBbcSig::FillPed0(const Int_t sampmin, const Int_t sampmax)
{
  Double_t x, y;
  for (int isamp=sampmin; isamp<=sampmax; isamp++)
  {
    gRawPulse->GetPoint(isamp,x,y);
    //gRawPulse->Print("all");
    hPed0->Fill( y );

    /*
// chiu taken out
    ped0stats->Push( y );
    ped0 = ped0stats->Mean();
    ped0rms = ped0stats->RMS();
    */

    //cout << "ped0 " << ch << " " << n << "\t" << ped0 << endl;
    //cout << "ped0 " << ch << "\t" << ped0 << endl;
  }

}


void OnlBbcSig::FillPed0(const Double_t begin, const Double_t end)
{
  Double_t x, y;
  Int_t n = gRawPulse->GetN();
  for (int isamp=0; isamp<n; isamp++)
  {
    gRawPulse->GetPoint(isamp,x,y);
    if ( x>=begin && x<=end )
    {
      hPed0->Fill( y );

      /*
      ped0stats->Push( y );
      ped0 = ped0stats->Mean();
      ped0rms = ped0stats->RMS();
      */

      //cout << "ped0 " << ch << " " << n << "\t" << x << "\t" << y << endl;
    }

    // quit if we are past the ped region
    if ( x>end ) break;
  }

}


void OnlBbcSig::SetPed0(const Double_t mean, const Double_t rms)
{
  ped0 = mean;
  ped0rms = rms;
  use_ped0 = 1;
  hpulse = hSubPulse;
  gpulse = gSubPulse;
  //if ( ch==8 ) cout << "ch " << ch << " Ped = " << ped0 << endl;
}

// Get Event by Event Ped0 if requested
void OnlBbcSig::CalcEventPed0(const Int_t minpedsamp, const Int_t maxpedsamp)
{
  //if (ch==8) cout << "In OnlBbcSig::CalcEventPed0(int,int)" << endl;
  hPed0->Reset();
  //ped0stats->Clear();

  Double_t x, y;
  for (int isamp=minpedsamp; isamp<=maxpedsamp; isamp++)
  {
    gRawPulse->GetPoint(isamp,x,y);

    hPed0->Fill(y);
    //ped0stats->Push( y );
    //if ( ch==8 ) cout << "ped0stats " << isamp << "\t" << y << endl;
  }


  // use straight mean for pedestal
  // Could consider using fit to hPed0 to remove outliers
  float mean = hPed0->GetMean();
  float rms = hPed0->GetRMS();

  //SetPed0( ped0stats->Mean(), ped0stats->RMS() );

  SetPed0( mean, rms );
  //if (ch==8) cout << "ped0stats mean, rms " << mean << "\t" << rms << endl;
}

// Get Event by Event Ped0 if requested
void OnlBbcSig::CalcEventPed0(const Double_t minpedx, const Double_t maxpedx)
{
  hPed0->Reset();
  //ped0stats->Clear();

  Double_t x, y;
  Int_t n = gRawPulse->GetN();

  for (int isamp=0; isamp<n; isamp++)
  {
    gRawPulse->GetPoint(isamp,x,y);

    if ( x>= minpedx && x<= maxpedx)
    {
      hPed0->Fill(y);
      //ped0stats->Push( y );
    }
  }

  // use straight mean for pedestal
  // Could consider using fit to hPed0 to remove outliers
  SetPed0( hPed0->GetMean(), hPed0->GetRMS() );
  /*
  Double_t mean = ped0stats->Mean();
  Double_t rms = ped0stats->RMS();
  SetPed0( mean, rms );
  */
  //cout << "ped0stats " << mean << "\t" << rms << endl;
}

Double_t OnlBbcSig::LeadingEdge(const Double_t threshold)
{
  // Find first point above threshold
  // We also make sure the next point is above threshold
  // to get rid of a high fluctuation
  int n = gSubPulse->GetN();
  Double_t *x = gSubPulse->GetX();
  Double_t *y = gSubPulse->GetY();

  int sample = -1;
  for (int isamp=0; isamp<n; isamp++)
  {
    if ( y[isamp] > threshold )
    {
      if ( isamp==(n-1) || y[isamp+1] > threshold )
      {
        sample = isamp;
        break;
      }
    }
  }
  if ( sample == -1 ) return -9999.;  // no signal above threshold

  // Linear Interpolation of start time
  Double_t dx = x[sample] - x[sample-1];
  Double_t dy = y[sample] - y[sample-1];
  Double_t dt1 = y[sample] - threshold;

  Double_t t0 = x[sample] - dt1*(dx/dy);

  return t0;
}

Double_t OnlBbcSig::dCFD(const Double_t fraction_threshold)
{
  // Find first point above threshold
  // We also make sure the next point is above threshold
  // to get rid of a high fluctuation
  int n = gSubPulse->GetN();
  Double_t *x = gSubPulse->GetX();
  Double_t *y = gSubPulse->GetY();

  // Get max amplitude
  Double_t ymax = TMath::MaxElement(n,y);
  if ( f_ampl == -9999. ) f_ampl = ymax;

  Double_t threshold = fraction_threshold * ymax; // get fraction of amplitude
  //cout << "threshold = " << threshold << "\tymax = " << ymax <<endl;

  int sample = -1;
  for (int isamp=0; isamp<n; isamp++)
  {
    if ( y[isamp] > threshold )
    {
      if ( isamp==(n-1) || y[isamp+1] > threshold )
      {
        sample = isamp;
        break;
      }
    }
  }
  if ( sample == -1 ) return -9999.;  // no signal above threshold

  // Linear Interpolation of start time
  Double_t dx = x[sample] - x[sample-1];
  Double_t dy = y[sample] - y[sample-1];
  Double_t dt1 = y[sample] - threshold;

  Double_t t0 = x[sample] - dt1*(dx/dy);

  return t0;
}

Double_t OnlBbcSig::MBD(const Int_t max_samp)
{
  // Get the amplitude of the sample number to get time
  Double_t *y = gSubPulse->GetY();

  if ( y==0 ) { 
      cout << "ERROR y == 0" << endl; 
      return 0;
  }

  // SHOULD INCLUDE TIME CALIBRATION HERE
  Double_t t0 = y[max_samp];

  // Get max amplitude, and set it if it hasn't already been set
  int n = gSubPulse->GetN();
  Double_t ymax = TMath::MaxElement(n,y);
  if ( f_ampl == -9999. ) f_ampl = ymax;

  return t0;
}

Double_t OnlBbcSig::Integral(const Double_t xmin, const Double_t xmax)
{
  Int_t n = gSubPulse->GetN();
  Double_t* x = gSubPulse->GetX();
  Double_t* y = gSubPulse->GetY();

  f_integral = 0.;
  for (int ix=0; ix<n; ix++)
  {
    if (x[ix]>=xmin && x[ix]<=xmax)
    {
      // Get dx
      Double_t dx = (x[ix+1]-x[ix-1])/2.0;
      f_integral += (y[ix]*dx);
    }
  }

  return f_integral;
}

void OnlBbcSig::LocMax(Double_t& x_at_max, Double_t& ymax, Double_t xminrange, Double_t xmaxrange)
{
  // Find index of maximum peak
  Int_t n = gSubPulse->GetN();
  Double_t* x = gSubPulse->GetX();
  Double_t* y = gSubPulse->GetY();

  // if flipped or equal, we search the whole range
  if ( xmaxrange <= xminrange )
  {
    xminrange = -DBL_MAX;
    xmaxrange = DBL_MAX;
  }

  ymax = -DBL_MAX;

  for (int i=0; i<n; i++)
  {
    // Skip if out of range
    if ( x[i] < xminrange ) continue;
    if ( x[i] > xmaxrange ) break;

    if ( y[i] > ymax )
    {
      ymax = y[i];
      x_at_max = x[i];
    }
  }

}

void OnlBbcSig::LocMin(Double_t& x_at_max, Double_t& ymin, Double_t xminrange, Double_t xmaxrange)
{
  // Find index of minimum peak (for neg signals)
  Int_t n = gSubPulse->GetN();
  Double_t* x = gSubPulse->GetX();
  Double_t* y = gSubPulse->GetY();

  // if flipped or equal, we search the whole range
  if ( xmaxrange <= xminrange )
  {
    xminrange = -DBL_MAX;
    xmaxrange = DBL_MAX;
  }

  ymin = DBL_MAX;

  for (int i=0; i<n; i++)
  {
    // Skip if out of range
    if ( x[i] < xminrange ) continue;
    if ( x[i] > xmaxrange ) break;

    if ( y[i] < ymin )
    {
      ymin = y[i];
      x_at_max = x[i];
    }
  }

  // old way of getting locmax
  //int locmax = TMath::LocMin(n,y);
}

void OnlBbcSig::Print()
{
  Double_t x, y;
  cout << "CH " << ch << endl;
  for (int isamp=0; isamp<nsamples; isamp++)
  {
    gpulse->GetPoint(isamp,x,y);
    cout << isamp << "\t" << x << "\t" << y << endl;
  }
}

void OnlBbcSig::PadUpdate()
{
  // Make sure TCanvas is created externally!
  gPad->Modified();
  gPad->Update();
  cout << ch << " ? ";
  TString junk;
  cin >> junk;

  if (junk[0] == 'w' || junk[0] == 's')
  {
    TString name = "ch"; name += ch; name += ".png";
    gPad->SaveAs( name );
  }
}

Double_t OnlBbcSig::TemplateFcn(const Double_t *x, const Double_t *par)
{
  // par[0] is the amplitude (relative to the spline amplitude)
  // par[1] is the start time (in sample number)
  // x[0] units are in sample number
  Double_t xx = x[0]-par[1];
  Double_t f = 0.;

  //verbose = 100;

  // When fit is out of limits of good part of spline, ignore fit
  if ( xx<template_begintime || xx>template_endtime )
    {
      TF1::RejectPoint();
      if ( xx < template_begintime )
        {
          //Double_t x0,y0;
          Double_t y0 = template_y[0];
          return par[0]*y0;
        }
      else if ( xx > template_endtime )
        {
          //Double_t x0,y0;
          Double_t y0 = template_y[template_npointsx-1];
          return par[0]*y0;
        }
    }

  // Linear Interpolation of template
  Double_t x0 = 0.;
  Double_t y0 = 0.;
  Double_t x1 = 0.;
  Double_t y1 = 0.;

  // find the index in the vector which is closest to xx
  Double_t step = (template_endtime - template_begintime) / (template_npointsx-1);
  Double_t index = (xx - template_begintime)/step;

  int ilow = TMath::FloorNint( index );
  int ihigh = TMath::CeilNint( index );
  if ( ilow < 0 || ihigh >= template_npointsx )
    {
      if ( verbose>0 )
      {
        cout << "ERROR, ilow ihigh " << ilow << "\t" << ihigh << endl;
        cout << " " << xx << " " << x[0] << " " << par[1] << endl;
      }

      if ( ilow<0 )
      {
        ilow = 0;
      }
      else if ( ihigh >= template_npointsx )
      {
        ihigh = template_npointsx - 1;
      }
    }

  if ( ilow==ihigh )
    {
      f = par[0]*template_y[ilow];
    }
  else
    {
      x0 = template_begintime + ilow*step;
      y0 = template_y[ilow];
      x1 = template_begintime + ihigh*step;
      y1 = template_y[ihigh];
      f = par[0]*(y0+((y1-y0)/(x1-x0))*(xx-x0));  // linear interpolation
    }

  return f;
}

int OnlBbcSig::FitTemplate()
{
  //verbose = 100;	// uncomment to see fits
  if ( verbose>0 ) cout << "Fitting ch " << ch << endl;
 
  /*
  Int_t maxadc = -1;
  Int_t peak_samp = -1;
  for (int isamp=0; isamp<NSAMPLES; isamp++)
    {
      x[isamp] = (Float_t)isamp;
      adcsub[isamp] = adc[isamp] - ped;
      adcerr[isamp] = pedrms;
      if ( adcsub[isamp]>maxadc )
        {
          maxadc = adcsub[isamp];
          peak_samp = isamp;
        }
    }

  TGraphErrors gSubpulse(NSAMPLES,x,adcsub,0,adcerr);
  if ( verbose>10 )
    {
      gSubPulse->SetMarkerStyle(20);
      gSubPulse->SetMarkerColor(2);
      gSubPulse->SetLineColor(2);
    }

  fit_shape = fit_pshape[ch];
  fit_sherr = fit_psherr[ch];
  //if ( verbose>10 ) cout << "CH is " << ch << "\t" << fit_shape << endl;
  */
  
  // Check if channel is empty
  if ( gSubPulse->GetN() == 0 )
  {
    f_ampl = -9999.;
    f_time = -9999.;
    //cout << "gSubPulse empty" << endl;
    return 1;
  }

  // Get x-position of maximum
  Double_t x_at_max, ymax;
  LocMax(x_at_max, ymax);

  template_fcn->SetParameters(ymax, x_at_max);

  //template_fcn->SetParLimits(1,-5.,4.);
  template_fcn->SetRange(template_min_xrange,template_max_xrange);
  if ( verbose==0 ) gSubPulse->Fit(template_fcn,"RNQ");
  else              gSubPulse->Fit(template_fcn,"R");

  // Get fit parameters
  f_ampl = template_fcn->GetParameter(0);
  f_time = template_fcn->GetParameter(1);

  if ( verbose>0 && fabs(f_ampl) > 0. )
  {
    cout << "FitTemplate " << f_ampl << "\t" << f_time << endl;
    gSubPulse->Draw("ap");
    template_fcn->SetLineColor(4);
    template_fcn->Draw("same");
    PadUpdate();
  }

  /*
  Double_t chi2ndf = template_fcn->GetChisquare()/template_fcn->GetNDF();
  // Store fit values
  amp = static_cast<Float_t>( template_fcn->GetParameter(0) );
  fquality = (Short_t)chi2ndf;	// Need to define this still

  // For the tdc, we choose 360 tdc ticks per sample
  Float_t samp_number = 7.0+template_fcn->GetParameter(1);
  if ( samp_number<0. ) 
    {
      samp_number = 0.;
      amp = 0;    // for now, if the time is bad, we zero out the channel
    }
  else if ( samp_number>18.0 )
    {
      samp_number = 18.0;
      amp = 0;
    }
  tdc = static_cast<Short_t>( samp_number*360. );
  */

  return 1;
}

int OnlBbcSig::FillSplineTemplate()
{
  //verbose = 100;
  
  if ( ! hAmpl )
  {
    TString name = "hAmpl"; name += ch;
    hAmpl = new TH1F(name,name,17100,-100,17000);
  }
  if ( ! hTime )
  {
    TString name = "hTime"; name += ch;
    hTime = new TH1F(name,name,3100,0,31);
  }


  Double_t max = TMath::MaxElement(gSubPulse->GetN(),gSubPulse->GetY());

  if ( max < template_min_good_amplitude ) return 0;

  if ( verbose ) gSubPulse->Draw("ap");
  TSpline3 s3("s3",gSubPulse);

  // First find maximum, to rescale
  f_ampl = -999999.;
  double step_size = h2Template->GetXaxis()->GetBinWidth(1);
  //cout << "step size " << step_size << endl;
  for (double ix=0; ix<nsamples; ix += step_size)
  {
    Double_t val = s3.Eval(ix);
    if ( val > f_ampl )
    {
      f_ampl = val;
    }
  }

  //cout << f_ampl << endl;
  if ( f_ampl<template_min_good_amplitude || f_ampl>template_max_good_amplitude ) return 0; 

  if ( verbose>0 )
  {
    //cout << ch << endl;
    s3.SetLineColor(2);
    s3.Draw("same");
    gSubPulse->Draw("p");
    PadUpdate();
  }


  // Now go back to find the time by finding x at midpoint of rise
  for (double ix=0; ix<nsamples; ix += step_size)
  {
    Double_t val = s3.Eval(ix);
    if ( val > 0.5*f_ampl )
    {
      // interpolate midpoint
      Double_t dy_mid = val - 0.5*f_ampl;
      Double_t dy_prev = val - s3.Eval(ix-step_size);
      Double_t dx_prev = step_size;

      f_time = ix - (dy_mid/dy_prev)*dx_prev;
      break;
    }
  }

  // correct the pulse back
  f_time -= f_time_offset;

  // Get Time and Max of spline to rescale
  //cout << f_ampl << "\t" << time << endl;
  hAmpl->Fill( f_ampl );
  hTime->Fill( f_time );

  //cout << "nsamples " << nsamples << endl;
  for (int isamp=0; isamp<nsamples; isamp++)
  {
    Double_t x, y;
    gSubPulse->GetPoint(isamp,x,y);

    Double_t fillvalues[2] = {0.};
    fillvalues[0] = x - f_time; //corr_time
    if ( f_ampl != 0. )
    {
      fillvalues[1] = y/f_ampl; //scaled_ampl
    }

    h2Template->Fill(fillvalues[0],fillvalues[1]);
    //h2Template->Fill( fillvalues );

    // For splines, by defn they go through the data pts so residuals = 0
    //h2Residuals->Fill( fillvalues[0], y - s3.Eval(x) );
  }

  return 1;
}


void OnlBbcSig::MakeAndWriteTemplate(ostream& out, ostream& oerr)
{
  //verbose = 100;
  if ( verbose ) cout << "In  OnlBbcSig::MakeAndWriteTemplate" << endl;

  //Int_t nbinsy = h2Template->GetNbinsX();
 
  TString name = h2Template->GetName(); name += "_py";
  TH1 *hprojy = h2Template->ProjectionY(name,1,1,"e");
  TF1 *gaus = new TF1("gaussian","gaus",template_begintime,template_endtime);
  gaus->SetLineColor(2);

  if ( verbose )
  {
    h2Template->Draw("colz");
    //TH2 *h = h2Template->Projection(1,0);
    //h->Draw("colz");
    PadUpdate();
    //delete h;
  }

  TString fitarg = "R";
  if ( verbose==0 ) fitarg += "NQ";

  for (int ibin=0; ibin<template_npointsx; ibin++)
  {
    hprojy = h2Template->ProjectionY(name,ibin+1,ibin+1,"e");
   
    //h2Template->GetAxis(0)->SetRange(ibin+1,ibin+1);
    //TH1 *hprojy = h2Template->Projection(1,"e");
    //hprojy->Sumw2();

    Double_t ncounts = hprojy->Integral();
    if ( ncounts < 10. )
    {
      //cout << "ERROR, " << ch << "\t" << ibin << "\tToo few entries " << hprojy->Integral() << endl;
      fitarg += "WLL";
    }
    else if ( ncounts<100. )
    {
      fitarg += "WLL";
    }

    Double_t peak = hprojy->GetBinContent( hprojy->GetMaximumBin() );
    Double_t mean = hprojy->GetBinCenter( hprojy->GetMaximumBin() );
    Double_t rms = hprojy->GetRMS();
    //cout << ch << "\t" << ibin << "\t" << peak << "\t" << mean << endl;
    //gaus->SetParameters(peak,mean,0.01);
    gaus->SetParameters(peak,mean,rms);

    if ( ncounts>=10 ) hprojy->Fit(gaus,fitarg);

    // See the fits
    if ( verbose )
    {
      hprojy->Draw();
      PadUpdate();
    }

    /*
    template_y[ibin] = gaus->GetParameter(1);
    template_yrms[ibin] = gaus->GetParameter(2);
    */

    template_y[ibin] = hprojy->GetMean();
    template_yrms[ibin] = hprojy->GetRMS();

    if ( verbose ) cout << "ibin " << ibin << "\t" << template_y[ibin] << "\t" << template_yrms[ibin] << endl;
    //delete hprojy;

  }
  delete hprojy;
  delete gaus;

  /* //thnsparse
  h2Template->GetAxis(0)->SetRange(1,template_npointsx);
  h2Template->Write();
  */

  out << ch << "\t" << template_npointsx << "\t" << template_begintime << "\t" << template_endtime << endl;
  oerr << ch << "\t" << template_npointsx << "\t" << template_begintime << "\t" << template_endtime << endl;
  for (int ibin=0; ibin<template_npointsx; ibin++)
  {
    // Write out the template value for bin i
    out << template_y[ibin] << "\t";
    if (ibin%10==9) out << endl;

    // Now write out the rms
    oerr << template_yrms[ibin] << "\t";
    if (ibin%10==9) oerr << endl;
  }
  out << endl;
  oerr << endl;

}

void OnlBbcSig::FillFcnTemplate()
{
  //verbose = 100;
  if ( ! hAmpl )
  {
    TString name = "hAmpl"; name += ch;
    hAmpl = new TH1F(name,name,17100,-100,17000);
  }
  if ( ! hTime )
  {
    TString name = "hTime"; name += ch;
    hTime = new TH1F(name,name,3100,0,31);
  }


  FitTemplate();

  if ( f_ampl<template_min_good_amplitude || f_ampl>template_max_good_amplitude ) return; 

  // Get Time and Max of spline to rescale
  if ( verbose ) cout << "ampl time " << ch << "\t" << f_ampl << "\t" << f_time << endl;
  hAmpl->Fill( f_ampl );
  hTime->Fill( f_time );

  // Rescale and shift time using template fit
  for (int isamp=0; isamp<nsamples; isamp++)
  {
    Double_t x, y;
    gSubPulse->GetPoint(isamp,x,y);

    Double_t fillvalues[2] = {0};
    fillvalues[0] = x - f_time; //corr_time
    //Double_t scaled_ampl = 0.;
    if ( f_ampl != 0 )
    {
      fillvalues[1] = y/f_ampl; //scaled_ampl
    }
    h2Template->Fill( fillvalues[0], fillvalues[1] );
    //h2Template->Fill( fillvalues );
  }
}


int OnlBbcSig::ReadTemplate(ifstream& shapefile, ifstream& sherrfile)
{
  //verbose = 100;
  Int_t temp_ch = -9999;
  Int_t temp_nsamples;
  Double_t temp_begintime;
  Double_t temp_endtime;

  template_y.clear();
  template_yrms.clear();

  // Template 
  while ( shapefile >> temp_ch >> temp_nsamples >> temp_begintime >> temp_endtime )
  {
    if ( verbose ) cout << "shape " << temp_ch << "\t" <<  temp_nsamples << "\t" <<  temp_begintime << "\t" <<  temp_endtime << endl;
    if ( temp_ch != ch )
    {
      cerr << "ERROR in shape: ch is " << temp_ch << "but should be " << ch << endl;
      return -1;
    }

    Double_t temp_val;
    for (int isamp=0; isamp<temp_nsamples; isamp++)
    {
      shapefile >> temp_val;
      template_y.push_back( temp_val );
      if ( verbose )
      {
        cout << template_y[isamp] << " ";
        if ( isamp%10==9 ) cout << endl;
      }
    }
    if ( verbose ) cout << endl;
    break;
  }

  // Now get the errors
  while ( sherrfile >> temp_ch >> temp_nsamples >> temp_begintime >> temp_endtime )
  {
    if ( verbose ) cout << "sherr " << temp_ch << "\t" <<  temp_nsamples << "\t" <<  temp_begintime << "\t" <<  temp_endtime << endl;
    if ( temp_ch != ch )
    {
      cerr << "ERROR in sherr: ch is " << temp_ch << " but should be " << ch << endl;
      return -1;
    }

    Double_t temp_val;
    for (int isamp=0; isamp<temp_nsamples; isamp++)
    {
      sherrfile >> temp_val;
      template_yrms.push_back( temp_val );
      if ( verbose )
      {
        cout << template_yrms[isamp] << " ";
        if ( isamp%10==9 ) cout << endl;
      }
    }
    if ( verbose ) cout << endl;
    break;
  }

  TString name = "template_fcn"; name += ch;
  template_fcn = new TF1(name,this,&OnlBbcSig::TemplateFcn,0,nsamples,2,"OnlBbcSig","TemplateFcn");
  template_fcn->SetParameters(1,8);

  return 1;
}

