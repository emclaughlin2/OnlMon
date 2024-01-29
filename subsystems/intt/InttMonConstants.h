#ifndef INTT_MON_CONSTANTS_H
#define INTT_MON_CONSTANTS_H

#include <iostream>

namespace INTT
{
	constexpr int LAYER_OFFSET = 0;
	constexpr int CHIP_OFFSET = 1;

	const int FELIX = 8;
	const int FELIX_CHANNEL = 14;

	constexpr int LAYER = 4;
	constexpr int LADDER[LAYER] = {12, 12, 16, 16};
	constexpr int ARM = 2;

	constexpr int CHIP = 26;
	constexpr int CHANNEL = 128;
	constexpr int ADC = 8;

	constexpr int LADDERS_()
	{
		int L = 0;
	
		for(int layer = 0; layer < LAYER; layer++)
		{
			L += LADDER[layer];
		}
	
		return L;
	}
	
	constexpr int LADDERS = LADDERS_();

	constexpr int FELIX_CHANNELS = FELIX * FELIX_CHANNEL;
	constexpr int CHIPS =  CHIP * FELIX_CHANNELS;
	constexpr int CHANNELS = CHANNEL * CHIPS;
	constexpr int ADCS = CHANNELS * ADC;

	struct Indexes_s
	{
		int lyr = 0;
		int ldr = 0;
		int arm = 0;
		int chp = 0;
		int chn = 0;
		int adc = 0;
	};

	void GetFelixBinFromIndexes(int&, int, struct Indexes_s const&);
	void GetFelixIndexesFromBin(int, int&, struct Indexes_s&);

	void GetGlobalChipBinXYFromIndexes(int&, int&, struct Indexes_s const&);
	void GetIndexesFromGlobalChipBinXY(int const&, int const&, struct Indexes_s&);

	void GetLocalChannelBinXYFromIndexes(int&, int&, struct Indexes_s const&);
	void GetIndexesFromLocalChannelBinXY(int const&, int const&, struct Indexes_s&);

	void GetGlobalLadderBinXYFromIndexes(int&, int&, struct Indexes_s const&);
	void GetIndexesFromGlobalLadderBinXY(int const&, int const&, struct Indexes_s&);

	void GetLocalChipBinXYFromIndexes(int&, int&, struct Indexes_s const&);
	void GetIndexesFromLocalChipBinXY(int const&, int const&, struct Indexes_s&);
}

//      Ladder Structure                //
//======================================//
//      U14     U1      Type B  North   //
//      U15     U2        .       .     //
//      U16     U3        .       .     //
//      U17     U4        .       .     //
//      U18     U5      Type B    .     //
//------------------------------  .     //
//      U19     U6      Type A    .     //
//      U20     U7        .       .     //
//      U21     U8        .       .     //
//      U22     U9        .       .     //
//      U23     U10       .       .     //
//      U24     U11       .       .     //
//      U25     U12       .       .     //
//      U26     U13     Type A  North   //
//--------------------------------------//
//      U13     U26     Type A  South   //
//      U12     U25       .       .     //
//      U11     U24       .       .     //
//      U10     U23       .       .     //
//      U9      U22       .       .     //
//      U8      U21       .       .     //
//      U7      U20       .       .     //
//      U6      U19     Type A    .     //
//------------------------------  .     //
//      U5      U18     Type B    .     //
//      U4      U17       .       .     //
//      U3      U16       .       .     //
//      U2      U15       .       .     //
//      U1      U14     Type B  South   //
//======================================//

#endif
