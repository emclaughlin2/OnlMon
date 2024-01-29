#ifndef INTT_MON_H
#define INTT_MON_H

#include "InttMonConstants.h"
#include "InttFelixMap.h" 

#include <onlmon/OnlMon.h>
#include <onlmon/OnlMonDB.h>
#include <onlmon/OnlMonServer.h>

//#include <pmonitor/pmonitor.h>
#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/msg_profile.h>

#include <TH1D.h>
#include <TH2D.h>
#include <TRandom.h> //for rng; remove later

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

class InttMon : public OnlMon
{
public:
	InttMon(const std::string &name);
	virtual ~InttMon();
	
	int Init();
	int BeginRun(const int);
	int process_event(Event*);
	int Reset();

	//for testing/debugging without unpacker, remove later
	int MiscDebug();
	void RandomEvent(int);
private:
	//// for testing/debugging without unpacker, remove later
	TRandom* rng = nullptr;
	//int InitExpectationHists();
	////~for testing/debugging without unpacker, remove later

	int DBVarInit();
	int DBVarUpdate();

	OnlMonDB* dbvars = nullptr;
	int evtcnt = 0;

	TH1D* NumEvents = nullptr;
	TH1D* HitMap = nullptr;
	//...
};

#endif
