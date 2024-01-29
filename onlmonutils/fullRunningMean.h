#ifndef __FULLRUNNINGMEAN_H__
#define __FULLRUNNINGMEAN_H__

/** 
This is the full running mean class. 

It will calculate the "real" (mathematically correct) running 
mean value of a certain depth d, which comes at a price. For many online monitoring 
applicatiomns, you will be better off with the "pseudoRunningMean" class, which is 
a good approximation of the true running mean value. This class is provided to allow
you to check that your pseudo running mean is a good enough approximation of the true value.
Also, if the series of input values has large variations, the pseudo running mean will 
be off by a few percent.   
 
Since you will need to store the d most recent entries to the running
mean this can lead to excessive amount of memory allocated.


This class is mean to monitor lots of  values (such as all channels of a given
detector) simultaneously; In the constructor you specify the "width" (how many channels) 
and the depth of the running mean. 

Since both this fullRunningMean and pseudoRunningMean inherit from the abstract runningMean class,
you can change your choice by just instantiating a different object, as in the following
example (the calculate_running_mean function will accept either class) 

\begin{verbatim}
 runningMean *pm = new fullRunningMean(144,50);
\end{verbatim}

and the pass your "laser" events and the pm object on to a routine that add the values:

\begin{verbatim}
int calculate_running_mean ( Event * evt, runningMean * rm)
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

Then, outside of your routine, you could ask the pm object for the running mean value 
of channel i with

\begin{verbatim}
pm->getmean(i);
\end{verbatim}

 

*/




#include "runningMean.h"


class fullRunningMean : public runningMean {

public:
  fullRunningMean( const int /*NumberofChannels*/, const int /*depth*/);
  ~fullRunningMean() override;

 // delete copy ctor and assignment operator (cppcheck)
  explicit fullRunningMean(const fullRunningMean&) = delete;
  fullRunningMean& operator=(const fullRunningMean&) = delete;

  /// the getMean(i) funtion returns the current mean value of channel i 
  double getMean(const int /*ich*/) const override;

  /// Reset will reset th whole class
  int Reset() override;

  /**Add will add a new list of readings. It is your responsibility 
     to provide an approriate array of readings. (Typically you can get 
     the array of int's from the Packet object's fillIntArray function).
   */ 
  int Add (const int /*iarr*/[]) override;
  int Add (const float /*farr*/[]) override;
  int Add (const double /*darr*/[]) override;

protected:

  int addChannel(const int /*channel*/, const double /*x*/);
  int depth;
  int current_depth;
  double ** array;

};
#endif


