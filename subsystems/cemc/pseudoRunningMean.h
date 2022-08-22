

#ifndef __PSEUDORUNNINGMEAN_H__
#define __PSEUDORUNNINGMEAN_H__

/** 
This is the pseudo running mean class. 

In order to calculate the "real" (mathematically correct) running 
mean value of a certain depth d, you will need to store the d most recent 
entries to the running mean. For large values of d, this can lead to excessive
amount of memory allocated.

This pseudo running mean is an approximation of the running mean. It will
work very well if you want to monitor values which don't vary too much, and 
makes it very well-suited for gain monitoring and similar applications, where
there are no dramatic changes in the values under normal conditions. It will
calculate the new running mean value rn from the existing value ro as

rn = ( ro * (d-1)/d + x ) /d

where x is the new value in the running mean calculation. If there
have been less than d readings so far, it will return the actual mean
value.


This class is mean to monitor lots of  values (such as all channels of a given
detector) simultaneously; In the cnstructor you specify the "width" (how many channels) 
and the depth of the running mean. 

This class is meant to be lightweight, so there is not much in the way of bounds checking 
of the input data going on. 

As a simple example, let's assume that you pre-select "laser trigger" events for the EmCal, and you 
want the running mean of the (144) channels calculated. You could construct a pseudoRunningMean object
with 144 channels and depth 50:

\begin{verbatim}
 pseudoRunningMean *pm = new pseudoRunningMean(144,50);
\end{verbatim}

and the pass your "laser" events and the pm object on to a routine that add the values:

\begin{verbatim}
int caculate_running_mean ( event * evt, RunningMean * rm)
{
   int array[144];
   Packet *p = evt->getPacket(8002);
   if (p) 
      {
          // yes, we got the packet, and we ask it now to get channel 57 for us
	  p->fillIntArray(array, 144, &nw);
          rm->Add(array);
	  delete p;
          return 0;
      }
   return 1;
}
\end{verbatim}

Then, outside of your routine, you could ask the pm object for the (pseudo-) running mean value 
of channel i with

\begin{verbatim}
pm->getmean(i);
\end{verbatim}

 

*/


#include "runningMean.h"


class pseudoRunningMean : public runningMean {

public:
  pseudoRunningMean( const int /*NumberofChannels*/, const int /*depth*/);
  virtual ~pseudoRunningMean();


  /// the getMean(i) funtion returns the current mean value of channel i 
  virtual double getMean(const int /*ich*/) const;

  /// Reset will reset th whole class
  virtual int Reset();

  /**Add will add a new list of readings. It is your responsibility 
     to provide an approriate array of readings. (Typically you can get 
     the array of int's from the Packet object's fillIntArray function).
   */ 
  int Add (const int /*iarr*/[]);
  int Add (const float /*farr*/[]);
  int Add (const double /*darr*/[]);

protected:

  int addChannel(const int /*channel*/, const double /*x*/);
  int depth;
  int current_depth;
  double * array;

};
#endif


