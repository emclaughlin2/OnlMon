#ifndef __RUNNINGMEAN_H__
#define __RUNNINGMEAN_H__



/** 
This is the abstract running mean parent class. 

We have two different running mean classes derive from it, a "real"
(mathematically correct) running mean value, and a "pseudo" value of the running
mean, which is a lot more efficient and fully adequate for most monitoring purposes.


These classes are meant to monitor lots of  values (such as all channels 
of a given detector) simultaneously; In the constructor you specify the 
"width" (how many channels) and the depth of the running mean. 

This class is meant to be lightweight, so there is not much in the way of 
bounds checking of the input data going on. 

*/



class runningMean {

public:
  runningMean();
  virtual ~runningMean();


  /// the getMean(i) funtion returns the current mean value of channel i 
  virtual double getMean(const int /*ich*/) const = 0;
  virtual double getReference(const int /*ich*/) const;
  virtual int getNumberofChannels() const {return NumberofChannels;};

  /// Reset will reset the mean values (not the references)
  virtual int Reset() = 0;

  /**Add will add a new list of readings. It is your responsibility 
     to provide an approriate array of readings. (Typically you can get 
     the array of int's from the Packet object's fillIntArray function).
   */ 
  virtual int Add (const int /*iarr*/[]) = 0;
  virtual int Add (const float /*farr*/[]) = 0;
  virtual int Add (const double /*darr*/[]) = 0;

  // provide an externAl, ready-made array of the means (e.g, from the
  // database) to be used as a reference
  int setRefArray(const double /*darr*/[]);

  // set the reference for individual channels
  int setRefChannel(const int /*channel*/, const double /*refvalue*/);

  // establish thecurrent mean values as the reference
  int setAsReference();

  double getPercentDeviation(const int /*channel*/) const;

 protected:

  int NumberofChannels = 0;
  double * refArray;


};
#endif


