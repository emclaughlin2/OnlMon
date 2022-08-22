#include <runningMean.h>

runningMean::runningMean()
{
  refArray = 0;
}

runningMean::~runningMean()
{
  if (refArray) delete [] refArray;
}

int runningMean::setRefArray(const double darr[])
{
  int i;
  if (refArray ==0) 
    {
      refArray = new double[getNumberofChannels()];
    }

  for (i=0; i< getNumberofChannels(); i++) refArray[i] = darr[i];

  return 0;
}

int runningMean::setRefChannel(const int channel, const double refvalue)
{
  int i;
  if (refArray ==0) 
    {
      refArray = new double[getNumberofChannels()];
      for (i=0; i< getNumberofChannels(); i++) refArray[i] = 0;
    }
  
  if (channel < 0 || channel >= getNumberofChannels() ) return -1;
  refArray[channel] = refvalue;
  
  return 0;
}

int runningMean::setAsReference()
{
  int i;
  if (refArray ==0) 
    {
      refArray = new double[getNumberofChannels()];
    }
  for (i=0; i< getNumberofChannels(); i++) refArray[i] = getMean(i);
  return 0;
}

double runningMean::getReference(const int channel) const
{
  if (refArray ==0) return 0;

  if (channel < 0 || channel >= getNumberofChannels() ) return 0;
  return refArray[channel];
}

double runningMean::getPercentDeviation(const int channel) const
{
  if (refArray ==0) return 0;

   if (channel < 0 || channel >= getNumberofChannels() ) return 0;

   if (getReference(channel) == 0.) return 0.;
   return ( 100.* ( getMean(channel) - getReference(channel) ) / getReference(channel) );
}

