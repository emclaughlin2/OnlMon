#ifndef INTT_MON_DRAW_H
#define INTT_MON_DRAW_H

#include "InttMonConstants.h"
#include "InttFelixMap.h"
#include "InttMon.h"

#include <onlmon/OnlMonDraw.h>
#include <onlmon/OnlMonClient.h>
#include <onlmon/OnlMonDB.h>

#include <TStyle.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TSystem.h>
#include <TROOT.h>

#include <TLine.h>
#include <TText.h>

#include <TH2D.h>
#include <TH1D.h>

#include <cctype>
#include <cmath>
#include <ctime>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

class InttMonDraw : public OnlMonDraw
{
public:
	InttMonDraw(std::string const&);
	~InttMonDraw() override;

	int Init() override;
	int Draw(std::string const& = "ALL") override;
	int MakeHtml(std::string const& = "ALL") override;
	int SavePlot(std::string const& = "ALL", std::string const& = "png") override;

	typedef void(*HeadFunc_t)(std::string const&);
	typedef void(*GlobalHistsPrepFunc_t)(std::string const&, TH2D**);
	typedef void(*LocalHistPrepFunc_t)(std::string const&, TH2D**, struct INTT::Indexes_s&);
	struct Option_s
	{
		HeadFunc_t head = nullptr;
		GlobalHistsPrepFunc_t global = nullptr;
		LocalHistPrepFunc_t local = nullptr;
	};
	typedef std::map<std::string, struct Option_s> Options_t;
	static Options_t OPTIONS;

	//GlobalChip-Channel idiom
	static void GlobalChipLocalChannelHead(std::string const&);
	static void DrawGlobalChipMap(std::string const&);
	static void InttGlobalChipExec(const std::string&, int);
	static void DrawLocalChannelMap(std::string const&, struct INTT::Indexes_s);
	static void InttLocalChannelExec(const std::string&, int, int, int, int);

	//GlobalLadder-Chip idiom
	static void GlobalLadderLocalChipHead(std::string const&);
	static void DrawGlobalLadderMap(std::string const&);
	static void InttGlobalLadderExec(const std::string&, int);
	static void DrawLocalChipMap(std::string const&, struct INTT::Indexes_s);
	static void InttLocalChipExec(const std::string&, int, int, int);

	//GlobalChip-Channel methods
	static void PrepGlobalChipHists_Hitmap(std::string const&, TH2D**);
	static void PrepLocalChannelHists_Hitmap(std::string const&, TH2D**, struct INTT::Indexes_s&);
	static void PrepGlobalChipHists_NLL(std::string const&, TH2D**);

	//GlobalLadder-Chip methods
	static void PrepGlobalLadderHists_Interface(std::string const&, TH2D**);
	static void PrepLocalChipHists_Hitmap(std::string const&, TH2D**, struct INTT::Indexes_s&);

private:
	OnlMonDB *dbvars = nullptr;

	//===	Constants for Drawing	===//
	static constexpr double NUM_SIG = 2.0;

	static constexpr int CNVS_WIDTH = 1280;
	static constexpr int CNVS_HEIGHT = 720;

	static constexpr double T_MARGIN = 0.01;
	static constexpr double B_MARGIN = 0.01;
	static constexpr double L_MARGIN = 0.01;
	static constexpr double R_MARGIN = 0.08;

	static constexpr double TOP_FRAC = 0.05;
	static constexpr double DISP_FRAC = 0.1;
	static constexpr double DISP_TEXT_SIZE = 0.20;

	static constexpr double KEY_FRAC = 0.0;
	static constexpr double Y_LABEL_FRAC = 0.05;
	static constexpr double X_LABEL_FRAC = 0.15;
	static constexpr double Y_LABEL_TEXT_SIZE = 0.5;
	static constexpr double X_LABEL_TEXT_SIZE = 0.5;
	//===	~Constants for Drawing	===//

	//===	Drawing Methods		===//
	static void DrawPad(TPad*, TPad*);
	//===	~Drawing Methods	===//
};

#endif
