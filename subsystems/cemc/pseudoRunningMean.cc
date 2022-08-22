#include "pseudoRunningMean.h"

pseudoRunningMean::pseudoRunningMean( const int n, const int d)
{
  int i;
  depth = d;
  NumberofChannels  = n;
  current_depth =0;
  array = new double[NumberofChannels];
  for (i=0; i< NumberofChannels; i++) array[i] = 0.;
}

pseudoRunningMean::~pseudoRunningMean()
{
  delete [] array;
}


int pseudoRunningMean::Add (const int iarr[])
{
  int i;

  for (i = 0; i<NumberofChannels; i++) addChannel(i, iarr[i]);
  if ( current_depth < depth) current_depth++;
  return 0;
						  
}
int pseudoRunningMean::Add (const float farr[])
{
  int i;

  for (i = 0; i<NumberofChannels; i++) addChannel(i, farr[i]);
  if ( current_depth < depth) current_depth++;
  return 0;
						  
}
int pseudoRunningMean::Add (const double darr[])
{
  int i;

  for (i = 0; i<NumberofChannels; i++) addChannel(i, darr[i]);
  if ( current_depth < depth) current_depth++;
  return 0;
}


int pseudoRunningMean::Reset()
{
  int i;
  for (i=0; i<NumberofChannels; i++) array[i] = 0.;
  return 0;
}

double pseudoRunningMean::getMean(const int ich) const
{
  if (current_depth ==0)  return 0;
  return array[ich] / double(current_depth);
}


int pseudoRunningMean::addChannel(const int channel, const double k)
{

  if ( current_depth < depth)
    {
       array[channel] += k;
    }
  else
    {
      double x,y,ratio;
      x = current_depth -1;
      y = current_depth;
      ratio = x/y;
      array[channel] = array[channel] * ratio + k;
    }
  return 0;
}

