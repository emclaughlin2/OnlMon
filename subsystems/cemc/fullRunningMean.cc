#include "fullRunningMean.h"

fullRunningMean::fullRunningMean( const int n, const int d)
{
  int i,j;
  depth = d;
  NumberofChannels  = n;
  current_depth =0;
  array = new double *[NumberofChannels];
  for (i=0; i< NumberofChannels; i++) 
    {
      array[i] = new double[depth];
      for (j=0; j< depth; j++) array[i][j] = 0.;
    }
}

fullRunningMean::~fullRunningMean()
{
  int i;
  for (i=0; i< NumberofChannels; i++) 
    {
      delete [] array[i];
    }
  delete [] array;
}

int fullRunningMean::Add (const int iarr[])
{
  int i;

  for (i = 0; i<NumberofChannels; i++) addChannel(i, double(iarr[i]) );
  if ( current_depth < depth) current_depth++;
  return 0;
						  
}
int fullRunningMean::Add (const float farr[])
{
  int i;

  for (i = 0; i<NumberofChannels; i++) addChannel(i, farr[i]);
  if ( current_depth < depth) current_depth++;
  return 0;
						  
}
int fullRunningMean::Add (const double darr[])
{
  int i;

  for (i = 0; i<NumberofChannels; i++) addChannel(i, darr[i]);
  if ( current_depth < depth) current_depth++;
  return 0;
}


int fullRunningMean::Reset()
{
  int i,j;
  for (i=0; i< NumberofChannels; i++) 
    {
      for (j=0; j< depth; j++) array[i][j] = 0.;
    }
  return 0;
}

double fullRunningMean::getMean(const int ich) const
{
  int j;
  double x = 0.;

  if (current_depth == 0) return 0.;

  for (j=0; j<current_depth; j++) x += array[ich][j];
  return  x / double(current_depth);
}


int fullRunningMean::addChannel(const int channel, const double k)
{
  int j;
  if ( current_depth < depth)
    {
       array[channel][current_depth] = k;
    }
  else
    {
      for (j=1; j<depth; j++) array[channel][j-1] = array[channel][j];
      array[channel][depth-1] = k;
    }
  return 0;
}

