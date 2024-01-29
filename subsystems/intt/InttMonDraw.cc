#include "InttMonDraw.h"

InttMonDraw::Options_t InttMonDraw::OPTIONS =
{
	//Chip-Channel
	{"chip_hitmap",		(struct Option_s)
		{
			.head =		&InttMonDraw::GlobalChipLocalChannelHead,
			.global =	&InttMonDraw::PrepGlobalChipHists_Hitmap,
			.local =	&InttMonDraw::PrepLocalChannelHists_Hitmap
		}
	},

	{"chip_nll",		(struct Option_s)
		{
			.head =		&InttMonDraw::GlobalChipLocalChannelHead,
			.global =	&InttMonDraw::PrepGlobalChipHists_NLL,
			.local =	&InttMonDraw::PrepLocalChannelHists_Hitmap
		}
	},


	//Ladder-Chip	
	{"ladder_hitmap",		(struct Option_s)
		{
			.head =		&InttMonDraw::GlobalLadderLocalChipHead,
			.global =	&InttMonDraw::PrepGlobalLadderHists_Interface,
			.local =	&InttMonDraw::PrepLocalChipHists_Hitmap
		}
	},
};


InttMonDraw::InttMonDraw(const std::string &name) : OnlMonDraw(name)
{
	return;
}


InttMonDraw::~InttMonDraw()
{
	return;
}

//===		Inherited Functions		===//
int InttMonDraw::Init()
{
	return 0;
}

int InttMonDraw::Draw(const std::string &what)
{
	bool b = false;
	bool found = false;

	std::string temp = "";
	for(std::size_t s = 0; s < what.length(); ++s)
	{
		temp += (char)std::tolower(what[s]);
	}

	for(Options_t::iterator itr = OPTIONS.begin(); itr != OPTIONS.end(); ++itr)
	{
		b = false;

		if(temp == "all")b = true;
		if(temp == itr->first)b = true;

		if(!b)continue;

		found = true;
		(*(itr->second.head))(itr->first);
	}

	if(!found)
	{
		std::cout << "Option \"" << what << "\" not found" << std::endl;
		std::cout << "Try \"all\" or one of the following (case insensitive)" << std::endl;
		for(Options_t::iterator itr = OPTIONS.begin(); itr != OPTIONS.end(); ++itr)
		{
			std::cout << "\t" << itr->first << std::endl;

			return 1;
		}
	}

	return 0;
}

int InttMonDraw::MakeHtml(const std::string &what)
{
	OnlMonClient* cl = OnlMonClient::instance();
	if(Draw(what))return 1;

	TSeqCollection* canvases = gROOT->GetListOfCanvases();
	TCanvas* canvas = nullptr;
	std::string pngfile;

	bool b = false;
	bool found = false;

	std::string temp = "";
	std::string name = "";
	for(std::size_t s = 0; s < what.length(); ++s)
	{
		temp += (char)std::tolower(what[s]);
	}
	for(Options_t::iterator itr = OPTIONS.begin(); itr != OPTIONS.end(); ++itr)
	{
		b = false;

		if(temp == "all")b = true;
		if(temp == itr->first)b = true;

		if(!b)continue;

		found = true;

		canvas = nullptr;
		name = Form("Intt_%s_Global_Canvas", (itr->first).c_str());
		for(TIter t_itr = canvases->begin(); t_itr != canvases->end(); ++t_itr)
		{
			if(std::string((*t_itr)->GetName()).find(name) == std::string::npos)continue;

			canvas = (TCanvas*)(*t_itr);
			break;
		}

		if(!canvas)continue;

		pngfile = cl->htmlRegisterPage(*this, canvas->GetTitle(), itr->first, "png");
		cl->CanvasToPng(canvas, pngfile);
	}

	if(!found)
	{
		std::cout << "Option \"" << what << "\" not found" << std::endl;
		std::cout << "Try \"all\" or one of the following (case insensitive)" << std::endl;
		for(Options_t::iterator itr = OPTIONS.begin(); itr != OPTIONS.end(); ++itr)
		{
			std::cout << "\t" << itr->first << std::endl;
		}
	}

	std::string logfile = cl->htmlRegisterPage(*this, "EXPERTS/Log", "log", "html");
	std::ofstream log_out(logfile.c_str());
	log_out << "<HTML><HEAD><TITLE>Log file for run " << cl->RunNumber();
	log_out << "</TITLE></HEAD>" << std::endl;
	log_out << "<P>[Implement me]" << std::endl;
	log_out.close();

	std::string stsfile = cl->htmlRegisterPage(*this, "EXPERTS/Status", "status", "html");
	std::ofstream sts_out(stsfile.c_str());
	sts_out << "<HTML><HEAD><TITLE>Status file for run " << cl->RunNumber();
	sts_out << "</TITLE></HEAD>" << std::endl;
	sts_out << "<P>[Implement me]" << std::endl;
	sts_out.close();

	cl->SaveLogFile(*this);

	return 0;
}

int InttMonDraw::SavePlot(std::string const& what, std::string const& type)
{
	OnlMonClient* cl = OnlMonClient::instance();
	if(Draw(what))return 1;

	TSeqCollection* canvases = gROOT->GetListOfCanvases();
	TCanvas* canvas = nullptr;
	std::string filename;

	bool b = false;
	bool found = false;

	std::string temp = "";
	std::string name = "";
	for(std::size_t s = 0; s < what.length(); ++s)
	{
		temp += (char)std::tolower(what[s]);
	}
	for(Options_t::iterator itr = OPTIONS.begin(); itr != OPTIONS.end(); ++itr)
	{
		b = false;

		if(temp == "all")b = true;
		if(temp == itr->first)b = true;

		if(!b)continue;

		found = true;

		canvas = nullptr;
		name = Form("Intt_%s_Global_Canvas", (itr->first).c_str());
		for(TIter t_itr = canvases->begin(); t_itr != canvases->end(); ++t_itr)
		{
			if(std::string((*t_itr)->GetName()).find(name) == std::string::npos)continue;

			canvas = (TCanvas*)(*t_itr);
			break;
		}

		if(!canvas)continue;

		filename = ThisName + "_" + itr->first + "_" + cl->RunNumber() + "." + type;
		cl->CanvasToPng(canvas, filename);
	}

	if(!found)
	{
		std::cout << "Option \"" << what << "\" not found" << std::endl;
		std::cout << "Try \"all\" or one of the following (case insensitive)" << std::endl;
		for(Options_t::iterator itr = OPTIONS.begin(); itr != OPTIONS.end(); ++itr)
		{
			std::cout << "\t" << itr->first << std::endl;
		}
	}

	return 0;
}

void InttMonDraw::DrawPad(TPad* base, TPad* pad)
{
	pad->SetFillStyle(4000); //transparent
	pad->Range(0.0, 0.0, 1.0, 1.0);
	pad->SetTopMargin(T_MARGIN);
	pad->SetBottomMargin(B_MARGIN);
	pad->SetLeftMargin(L_MARGIN);
	pad->SetRightMargin(R_MARGIN);

	base->cd();
	pad->Draw();
	pad->cd();
}

//GlobalChip-Channel idiom
void InttMonDraw::GlobalChipLocalChannelHead(std::string const& option)
{
	DrawGlobalChipMap(option);
}

void InttMonDraw::DrawGlobalChipMap(std::string const& option)
{
	double x_lower = KEY_FRAC + Y_LABEL_FRAC;
	double y_lower = 0.0;
	double x_upper = 1.0;
	double y_upper = 1.0;

	std::string name;
	struct INTT::Indexes_s indexes;

	TStyle* style;
	TCanvas* cnvs;
	TPad* disp_pad;
	TPad* hist_pad;
	TPad* grid_pad;
	TPad* ylbl_pad;
	TPad* xlbl_pad;
	TPad* exec_pad;

	TH2D* client_hists[4] = {nullptr};

	name = Form("Intt_%s_Global_Style", option.c_str());
	style = (TStyle*)gROOT->FindObject(name.c_str());
	if(!style)
	{
		style = new TStyle
		(
			name.c_str(),
			name.c_str()
		);
		style->SetOptStat(0);
		style->SetOptTitle(0);
		//...
	}
	style->cd();

	name = Form("Intt_%s_Global_Canvas", option.c_str());
	cnvs = (TCanvas*)gROOT->FindObject(name.c_str());
	if(!cnvs)
	{
		cnvs = new TCanvas
		(
			name.c_str(),
			name.c_str(),
			0,
			0,
			CNVS_WIDTH,
			CNVS_HEIGHT
		);
		//cnvs->...
		//...
	}

	name = Form("Intt_%s_Global_DispPad", option.c_str());
	disp_pad = (TPad*)gROOT->FindObject(name.c_str());
	if(!disp_pad)
	{
		disp_pad = new TPad
		(
			name.c_str(),
			name.c_str(),
			x_lower,
			0.0,
			x_upper,
			DISP_FRAC
		);
		DrawPad(cnvs, disp_pad);

		TText* disp_text = new TText
		(
			0.5,
			0.80,
			Form("Layer: %2d Ladder: %2d (%s) Chip: %2d", INTT::LAYER_OFFSET, 0, "South", INTT::CHIP_OFFSET)
		);
		disp_text->SetName(Form("Intt_%s_Global_DispText", option.c_str()));
		disp_text->SetTextAlign(22);
		disp_text->SetTextSize(DISP_TEXT_SIZE);
		disp_text->Draw();

		OnlMonClient* cl = OnlMonClient::instance();
		TText* run_text = new TText
		(
			0.5,
			0.50,
			Form("Run: %08d", cl->RunNumber())
		);
		run_text->SetName(Form("Intt_%s_Global_RunText", option.c_str()));
		run_text->SetTextAlign(22);
		run_text->SetTextSize(DISP_TEXT_SIZE);
		run_text->Draw();

		std::time_t now = std::time(nullptr);
		struct std::tm* time_s = std::localtime(&now);
		TText* time_text = new TText
		(
			0.5,
			0.20,
			Form("%02d/%02d/%4d (mm/dd/yyy)", time_s->tm_mon + 1, time_s->tm_mday, time_s->tm_year + 1900)
		);
		time_text->SetName(Form("Intt_%s_Global_TimeText", option.c_str()));
		time_text->SetTextAlign(22);
		time_text->SetTextSize(DISP_TEXT_SIZE);
		time_text->Draw();
	}

	Options_t::const_iterator itr = OPTIONS.find(option);
	if(itr == OPTIONS.end())return;
	(*(itr->second.global))(option, client_hists);

	for(indexes.lyr = 0; indexes.lyr < INTT::LAYER; ++indexes.lyr)
	{
		y_lower = (INTT::LAYER - indexes.lyr - (1.0 - X_LABEL_FRAC)) / INTT::LAYER * (1.0 - TOP_FRAC - DISP_FRAC) + DISP_FRAC;
		y_upper = (INTT::LAYER - indexes.lyr - 0.0) / INTT::LAYER * (1.0 - TOP_FRAC - DISP_FRAC) + DISP_FRAC;

		name = Form("Intt_%s_Global_HistPad_%d", option.c_str(), indexes.lyr);
		hist_pad = (TPad*)gROOT->FindObject(name.c_str());
		if(!hist_pad)
		{
			hist_pad = new TPad
			(
				name.c_str(),
				name.c_str(),
				x_lower,
				y_lower,
				x_upper,
				y_upper
			);
			DrawPad(cnvs, hist_pad);
		}

		hist_pad->cd();
		client_hists[indexes.lyr]->DrawCopy("COLZ");

		name = Form("Intt_%s_Global_GridPad_%d", option.c_str(), indexes.lyr);
		grid_pad = (TPad*)gROOT->FindObject(name.c_str());
		if(!grid_pad)
		{
			grid_pad = new TPad
			(
				name.c_str(),
				name.c_str(),
				x_lower,
				y_lower,
				x_upper,
				y_upper
			);
			DrawPad(cnvs, grid_pad);

			int i;
			double temp;

			for(i = 0; i < 2 * INTT::LADDER[indexes.lyr] + 1; ++i)
			{
				temp = L_MARGIN + (i / 2.0) * (1.0 - L_MARGIN - R_MARGIN) / INTT::LADDER[indexes.lyr];

				TLine* line = new TLine(temp, B_MARGIN, temp, 1.0 - T_MARGIN);			
				line->SetLineStyle((i % 2) ? 3 : 1);
				line->SetLineWidth((i % 2) ? 1 : 2);
				line->Draw();
			}
	
			for(i = 0; i < INTT::CHIP + 1; ++i)
			{
				temp = B_MARGIN + i * (1.0 - T_MARGIN - B_MARGIN) / INTT::CHIP;
		
				TLine* line = new TLine(L_MARGIN, temp, 1.0 - R_MARGIN, temp);	
				line->SetLineStyle((i % (INTT::CHIP / 2)) ? 3 : 1);
				line->SetLineWidth((i % (INTT::CHIP / 2)) ? 1 : 2);
				line->Draw();
			}
		}

		name = Form("Intt_%s_Global_YLblPad_%d", option.c_str(), indexes.lyr);
		ylbl_pad = (TPad*)gROOT->FindObject(name.c_str());
		if(!ylbl_pad)
		{
			ylbl_pad = new TPad
			(
				name.c_str(),
				name.c_str(),
				x_lower - Y_LABEL_FRAC,
				y_lower,
				x_lower,
				y_upper
			);
			DrawPad(cnvs, ylbl_pad);

			TText* label = new TText
			(
				0.25,
				0.5,
				Form("B%01dL%01d", indexes.lyr / 2, indexes.lyr % 2)
			);
			label->SetTextSize(Y_LABEL_TEXT_SIZE);
			label->SetTextAlign(22);
			label->SetTextAngle(90);
			label->Draw();

			int i = 0;
			for(i = 0; i < INTT::ARM; ++i)
			{
				label = new TText
				(
					0.75,
					0.25 + 0.5 * i,
					Form("%s", i ? "North" : "South")
				);
				label->SetTextSize(Y_LABEL_TEXT_SIZE / 2.0);
				label->SetTextAlign(22);
				label->SetTextAngle(90);
				label->Draw();
			}
		}

		name = Form("Intt_%s_Global_XLblPad_%d", option.c_str(), indexes.lyr);
		xlbl_pad = (TPad*)gROOT->FindObject(name.c_str());
		if(!xlbl_pad)
		{
			xlbl_pad = new TPad
			(
				name.c_str(),
				name.c_str(),
				x_lower,
				y_lower - X_LABEL_FRAC / INTT::LAYER * (1.0 - TOP_FRAC - DISP_FRAC),
				x_upper,
				y_lower
			);
			DrawPad(cnvs, xlbl_pad);

			int i;
			double temp;

			for(i = 0; i < INTT::LADDER[indexes.lyr]; ++i)
			{
				temp = L_MARGIN + (i + 0.5) * (1.0 - L_MARGIN - R_MARGIN) / INTT::LADDER[indexes.lyr];
				TText* label = new TText
				(
					temp,
					0.75,
					Form("L%02d", i)
				);
				label->SetTextSize(X_LABEL_TEXT_SIZE);
				label->SetTextAlign(22);
				label->Draw();
			}
		}

		name = Form("Intt_%s_Global_ExecPad_%d", option.c_str(), indexes.lyr);
		exec_pad = (TPad*)gROOT->FindObject(name.c_str());
		if(!exec_pad)
		{
			exec_pad = new TPad
			(
				name.c_str(),
				name.c_str(),
				x_lower + (x_upper - x_lower) * L_MARGIN,
				y_lower + (y_upper - y_lower) * B_MARGIN,
				x_upper - (x_upper - x_lower) * R_MARGIN,
				y_upper - (y_upper - y_lower) * T_MARGIN
			);
			DrawPad(cnvs, exec_pad);
			exec_pad->AddExec
			(
				Form("Intt_%s_Global_Exec_%d",option.c_str(), indexes.lyr),
				Form("InttMonDraw::InttGlobalChipExec(\"%s\", %d)", option.c_str(), indexes.lyr)
			);
		}
	}

	cnvs->Show();
	cnvs->SetEditable(0);
}

void InttMonDraw::InttGlobalChipExec(std::string const& option, int layer)
{
	int bin_x = gPad->AbsPixeltoX(gPad->GetEventX()) * (2 * INTT::LADDER[layer]);
	int bin_y = gPad->AbsPixeltoY(gPad->GetEventY()) * (INTT::CHIP);
	std::string name;

	if(bin_x < 0)bin_x = 0;
	if(bin_y < 0)bin_y = 0;

	if(bin_x >= 2 * INTT::LADDER[layer])bin_x = 2 * INTT::LADDER[layer] - 1;
	if(bin_y >= INTT::CHIP)bin_y = INTT::CHIP - 1;

	INTT::Indexes_s indexes;
	indexes.lyr = layer;
	INTT::GetIndexesFromGlobalChipBinXY(bin_x, bin_y, indexes);

	name = Form("Intt_%s_Global_DispText", option.c_str());
	TText* disp_text = (TText*)gROOT->FindObject(name.c_str());
	if(disp_text)disp_text->SetTitle(Form("Layer: %2d Ladder: %2d (%s) Chip: %2d", indexes.lyr + INTT::LAYER_OFFSET, indexes.ldr, indexes.arm ? "North" : "South", indexes.chp + INTT::CHIP_OFFSET));

	name = Form("Intt_%s_Global_DispPad", option.c_str());
	TPad* disp_pad = (TPad*)gROOT->FindObject(name.c_str());
	if(disp_pad)disp_pad->Update();

	if(gPad->GetEvent() != 11)return; //left click

	DrawLocalChannelMap(option, indexes);
}

void InttMonDraw::DrawLocalChannelMap(std::string const& option, struct INTT::Indexes_s indexes)
{
	double x_lower = 0.0;
	double y_lower = DISP_FRAC;
	double x_upper = 1.0;
	double y_upper = 1.0;

	std::string name;

	TStyle* style;
	TCanvas* cnvs;
	TPad* disp_pad;
	TPad* hist_pad;
	TPad* grid_pad;
	TPad* exec_pad;

	TH2D* client_hists[1] = {nullptr};
	Options_t::const_iterator itr = OPTIONS.find(option);
	if(itr == OPTIONS.end())return;
	(*(itr->second.local))(option, client_hists, indexes);

	name = Form("Intt_%s_Local_Style_Lyr%02d_Ldr%02d_Arm%02d_Chp%02d", option.c_str(), indexes.lyr, indexes.ldr, indexes.arm, indexes.chp);
	style = (TStyle*)gROOT->FindObject(name.c_str());
	if(!style)
	{
		style = new TStyle
		(
			name.c_str(),
			name.c_str()
		);
		style->SetOptStat(0);
		style->SetOptTitle(0);
		//...
	}
	style->cd();

	name = Form("Intt_%s_Local_Canvas_Lyr%02d_Ldr%02d_Arm%02d_Chp%02d", option.c_str(), indexes.lyr, indexes.ldr, indexes.arm, indexes.chp);
	cnvs = (TCanvas*)gROOT->FindObject(name.c_str());
	if(!cnvs)
	{
		cnvs = new TCanvas
		(
			name.c_str(),
			name.c_str(),
			0,
			0,
			CNVS_WIDTH,
			CNVS_HEIGHT
		);
	}

	name = Form("Intt_%s_Local_DispPad_Lyr%02d_Ldr%02d_Arm%02d_Chp%02d", option.c_str(), indexes.lyr, indexes.ldr, indexes.arm, indexes.chp);
	disp_pad = (TPad*)gROOT->FindObject(name.c_str());
	if(!disp_pad)
	{
		disp_pad = new TPad
		(
			name.c_str(),
			name.c_str(),
			x_lower,
			0.0,
			x_upper,
			DISP_FRAC
		);
		DrawPad(cnvs, disp_pad);
		TText* disp_text = new TText
		(
			0.5,
			0.5,
			Form("Channel: %3d ADC: %d", 0, 0)
		);
		disp_text->SetName(Form("Intt_%s_Local_DispText_Lyr%02d_Ldr%02d_Arm%02d_Chp%02d", option.c_str(), indexes.lyr, indexes.ldr, indexes.arm, indexes.chp));
		disp_text->SetTextAlign(22);
		disp_text->SetTextSize(DISP_TEXT_SIZE);
		disp_text->Draw();
	}

	name = Form("Intt_%s_Local_HistPad_Lyr%02d_Ldr%02d_Arm%02d_Chp%02d", option.c_str(), indexes.lyr, indexes.ldr, indexes.arm, indexes.chp);
	hist_pad = (TPad*)gROOT->FindObject(name.c_str());
	if(!hist_pad)
	{
		hist_pad = new TPad
		(
			name.c_str(),
			name.c_str(),
			x_lower,
			y_lower,
			x_upper,
			y_upper
		);
		DrawPad(cnvs, hist_pad);
	}

	hist_pad->cd();
	if(client_hists[0])
	{
		client_hists[0]->DrawCopy("COLZ");
		delete client_hists[0];
	}

	name = Form("Intt_%s_Local_GridPad_%d", option.c_str(), indexes.lyr);
	grid_pad = (TPad*)gROOT->FindObject(name.c_str());
	if(!grid_pad)
	{
		grid_pad = new TPad
		(
			name.c_str(),
			name.c_str(),
			x_lower,
			y_lower,
			x_upper,
			y_upper
		);
		DrawPad(cnvs, grid_pad);

		int i;
		double temp;

		for(i = 0; i < INTT::CHANNEL + 1; ++i)
		{
			temp = L_MARGIN + i * (1.0 - L_MARGIN - R_MARGIN) / INTT::CHANNEL;

			TLine* line = new TLine(temp, B_MARGIN, temp, 1.0 - T_MARGIN);			
			line->SetLineStyle(3);
			line->SetLineWidth(1);
			line->Draw();
		}

		for(i = 0; i < INTT::ADC + 1; ++i)
		{
			temp = B_MARGIN + i * (1.0 - T_MARGIN - B_MARGIN) / INTT::ADC;
	
			TLine* line = new TLine(L_MARGIN, temp, 1.0 - R_MARGIN, temp);	
			line->SetLineStyle(3);
			line->SetLineWidth(1);
			line->Draw();
		}
	}

	name = Form("Intt_%s_Local_ExecPad_Lyr%02d_Ldr%02d_Arm%02d_Chp%02d", option.c_str(), indexes.lyr, indexes.ldr, indexes.arm, indexes.chp);
	exec_pad = (TPad*)gROOT->FindObject(name.c_str());
	if(!exec_pad)
	{
		exec_pad = new TPad
		(
			name.c_str(),
			name.c_str(),
			x_lower + (x_upper - x_lower) * L_MARGIN,
			y_lower + (y_upper - y_lower) * B_MARGIN,
			x_upper - (x_upper - x_lower) * R_MARGIN,
			y_upper - (y_upper - y_lower) * T_MARGIN
		);
		DrawPad(cnvs, exec_pad);
		exec_pad->AddExec
		(
			Form("Intt_%s_Local_Exec_Lyr%02d_Ldr%02d_Arm%02d_Chp%02d", option.c_str(), indexes.lyr, indexes.ldr, indexes.arm, indexes.chp),
			Form("InttMonDraw::InttLocalChannelExec(\"%s\", %d, %d, %d, %d)", option.c_str(), indexes.lyr, indexes.ldr, indexes.arm, indexes.chp)
		);
	}
}

void InttMonDraw::InttLocalChannelExec(const std::string& option, int layer, int ladder, int northsouth, int chip)
{
	int bin_x = gPad->AbsPixeltoX(gPad->GetEventX()) * INTT::CHANNEL;
	int bin_y = gPad->AbsPixeltoY(gPad->GetEventY()) * INTT::ADC;
	std::string name;

	if(bin_x < 0)bin_x = 0;
	if(bin_y < 0)bin_y = 0;

	if(bin_x >= INTT::CHANNEL - 1)bin_x = INTT::CHANNEL - 1;
	if(bin_y >= INTT::ADC - 1)bin_y = INTT::ADC - 1;

	INTT::Indexes_s indexes;
	indexes.lyr = layer;
	indexes.ldr = ladder;
	indexes.arm = northsouth;
	indexes.chp = chip;
	INTT::GetIndexesFromLocalChannelBinXY(bin_x, bin_y, indexes);

	name = Form("Intt_%s_Local_DispText_Lyr%02d_Ldr%02d_Arm%02d_Chp%02d", option.c_str(), indexes.lyr, indexes.ldr, indexes.arm, indexes.chp);
	TText* disp_text = (TText*)gROOT->FindObject(name.c_str());
	if(disp_text)disp_text->SetTitle(Form("Channel: %3d ADC: %d", indexes.chn, indexes.adc));

	name = Form("Intt_%s_Local_DispPad_Lyr%02d_Ldr%02d_Arm%02d_Chp%02d", option.c_str(), indexes.lyr, indexes.ldr, indexes.arm, indexes.chp);
	TPad* disp_pad = (TPad*)gROOT->FindObject(name.c_str());
	if(disp_pad)disp_pad->Update();
}

//GlobalLadder-Chip idiom
void InttMonDraw::GlobalLadderLocalChipHead(std::string const& option)
{
	DrawGlobalLadderMap(option);
}

void InttMonDraw::DrawGlobalLadderMap(std::string const& option)
{
	double x_lower = KEY_FRAC + Y_LABEL_FRAC;
	double y_lower = 0.0;
	double x_upper = 1.0;
	double y_upper = 1.0;

	std::string name;
	struct INTT::Indexes_s indexes;

	TStyle* style;
	TCanvas* cnvs;
	TPad* disp_pad;
	TPad* hist_pad;
	TPad* grid_pad;
	TPad* ylbl_pad;
	TPad* xlbl_pad;
	TPad* exec_pad;

	TH2D* client_hists[4] = {nullptr};

	name = Form("Intt_%s_Global_Style", option.c_str());
	style = (TStyle*)gROOT->FindObject(name.c_str());
	if(!style)
	{
		style = new TStyle
		(
			name.c_str(),
			name.c_str()
		);
		style->SetOptStat(0);
		style->SetOptTitle(0);
		//...
	}
	style->cd();

	name = Form("Intt_%s_Global_Canvas", option.c_str());
	cnvs = (TCanvas*)gROOT->FindObject(name.c_str());
	if(!cnvs)
	{
		cnvs = new TCanvas
		(
			name.c_str(),
			name.c_str(),
			0,
			0,
			CNVS_WIDTH,
			CNVS_HEIGHT
		);
		//cnvs->...
		//...
	}

	name = Form("Intt_%s_Global_DispPad", option.c_str());
	disp_pad = (TPad*)gROOT->FindObject(name.c_str());
	if(!disp_pad)
	{
		disp_pad = new TPad
		(
			name.c_str(),
			name.c_str(),
			x_lower,
			0.0,
			x_upper,
			DISP_FRAC
		);
		DrawPad(cnvs, disp_pad);

		TText* disp_text = new TText
		(
			0.5,
			0.80,
			Form("Layer: %2d Ladder: %2d (%s)", INTT::LAYER_OFFSET, 0, "South")
		);
		disp_text->SetName(Form("Intt_%s_Global_DispText", option.c_str()));
		disp_text->SetTextAlign(22);
		disp_text->SetTextSize(DISP_TEXT_SIZE);
		disp_text->Draw();

		OnlMonClient* cl = OnlMonClient::instance();
		TText* run_text = new TText
		(
			0.5,
			0.50,
			Form("Run: %08d", cl->RunNumber())
		);
		run_text->SetName(Form("Intt_%s_Global_RunText", option.c_str()));
		run_text->SetTextAlign(22);
		run_text->SetTextSize(DISP_TEXT_SIZE);
		run_text->Draw();

		std::time_t now = std::time(nullptr);
		struct std::tm* time_s = std::localtime(&now);
		TText* time_text = new TText
		(
			0.5,
			0.20,
			Form("%02d/%02d/%4d (mm/dd/yyy)", time_s->tm_mon + 1, time_s->tm_mday, time_s->tm_year + 1900)
		);
		time_text->SetName(Form("Intt_%s_Global_TimeText", option.c_str()));
		time_text->SetTextAlign(22);
		time_text->SetTextSize(DISP_TEXT_SIZE);
		time_text->Draw();
	}

	Options_t::const_iterator itr = OPTIONS.find(option);
	if(itr == OPTIONS.end())return;
	(*(itr->second.global))(option, client_hists);

	for(indexes.lyr = 0; indexes.lyr < INTT::LAYER; ++indexes.lyr)
	{
		y_lower = (INTT::LAYER - indexes.lyr - (1.0 - X_LABEL_FRAC)) / INTT::LAYER * (1.0 - TOP_FRAC - DISP_FRAC) + DISP_FRAC;
		y_upper = (INTT::LAYER - indexes.lyr - 0.0) / INTT::LAYER * (1.0 - TOP_FRAC - DISP_FRAC) + DISP_FRAC;

		name = Form("Intt_%s_Global_HistPad_%d", option.c_str(), indexes.lyr);
		hist_pad = (TPad*)gROOT->FindObject(name.c_str());
		if(!hist_pad)
		{
			hist_pad = new TPad
			(
				name.c_str(),
				name.c_str(),
				x_lower,
				y_lower,
				x_upper,
				y_upper
			);
			DrawPad(cnvs, hist_pad);
		}

		hist_pad->cd();
		client_hists[indexes.lyr]->DrawCopy("COL");

		name = Form("Intt_%s_Global_GridPad_%d", option.c_str(), indexes.lyr);
		grid_pad = (TPad*)gROOT->FindObject(name.c_str());
		if(!grid_pad)
		{
			grid_pad = new TPad
			(
				name.c_str(),
				name.c_str(),
				x_lower,
				y_lower,
				x_upper,
				y_upper
			);
			DrawPad(cnvs, grid_pad);

			int i;
			double temp;

			for(i = 0; i < INTT::LADDER[indexes.lyr] + 1; ++i)
			{
				temp = L_MARGIN + i * (1.0 - L_MARGIN - R_MARGIN) / INTT::LADDER[indexes.lyr];
				TLine* line = new TLine(temp, B_MARGIN, temp, 1.0 - T_MARGIN);
				line->SetLineStyle(1);
				line->SetLineWidth(2);
				line->Draw();
			}

			for(i = 0; i < 2; ++i)
			{
				temp = B_MARGIN + i * (1.0 - T_MARGIN - B_MARGIN) / 2;
				TLine* line = new TLine(L_MARGIN, temp, 1.0 - R_MARGIN, temp);
				line->SetLineStyle(1);
				line->SetLineWidth(2);
				line->Draw();
			}
		}

		name = Form("Intt_%s_Global_YLblPad_%d", option.c_str(), indexes.lyr);
		ylbl_pad = (TPad*)gROOT->FindObject(name.c_str());
		if(!ylbl_pad)
		{
			ylbl_pad = new TPad
			(
				name.c_str(),
				name.c_str(),
				x_lower - Y_LABEL_FRAC,
				y_lower,
				x_lower,
				y_upper
			);
			DrawPad(cnvs, ylbl_pad);

			TText* label = new TText
			(
				0.25,
				0.5,
				Form("B%01dL%01d", indexes.lyr / 2, indexes.lyr % 2)
			);
			label->SetTextSize(Y_LABEL_TEXT_SIZE);
			label->SetTextAlign(22);
			label->SetTextAngle(90);
			label->Draw();

			int i = 0;
			for(i = 0; i < INTT::ARM; ++i)
			{
				label = new TText
				(
					0.75,
					0.25 + 0.5 * i,
					Form("%s", i ? "North" : "South")
				);
				label->SetTextSize(Y_LABEL_TEXT_SIZE / 2.0);
				label->SetTextAlign(22);
				label->SetTextAngle(90);
				label->Draw();
			}
		}

		name = Form("Intt_%s_Global_XLblPad_%d", option.c_str(), indexes.lyr);
		xlbl_pad = (TPad*)gROOT->FindObject(name.c_str());
		if(!xlbl_pad)
		{
			xlbl_pad = new TPad
			(
				name.c_str(),
				name.c_str(),
				x_lower,
				y_lower - X_LABEL_FRAC / INTT::LAYER * (1.0 - TOP_FRAC - DISP_FRAC),
				x_upper,
				y_lower
			);
			DrawPad(cnvs, xlbl_pad);

			int i;
			double temp;

			for(i = 0; i < INTT::LADDER[indexes.lyr]; ++i)
			{
				temp = L_MARGIN + (i + 0.5) * (1.0 - L_MARGIN - R_MARGIN) / INTT::LADDER[indexes.lyr];
				TText* label = new TText
				(
					temp,
					0.75,
					Form("L%02d", i)
				);
				label->SetTextSize(X_LABEL_TEXT_SIZE);
				label->SetTextAlign(22);
				label->Draw();
			}
		}

		name = Form("Intt_%s_Global_ExecPad_%d", option.c_str(), indexes.lyr);
		exec_pad = (TPad*)gROOT->FindObject(name.c_str());
		if(!exec_pad)
		{
			exec_pad = new TPad
			(
				name.c_str(),
				name.c_str(),
				x_lower + (x_upper - x_lower) * L_MARGIN,
				y_lower + (y_upper - y_lower) * B_MARGIN,
				x_upper - (x_upper - x_lower) * R_MARGIN,
				y_upper - (y_upper - y_lower) * T_MARGIN
			);
			DrawPad(cnvs, exec_pad);
			exec_pad->AddExec
			(
				Form("Intt_%s_Global_Exec_%d",option.c_str(), indexes.lyr),
				Form("InttMonDraw::InttGlobalLadderExec(\"%s\", %d)", option.c_str(), indexes.lyr)
			);
		}
	}

	cnvs->Show();
	cnvs->SetEditable(0);
}

void InttMonDraw::InttGlobalLadderExec(std::string const& option, int layer)
{
	int bin_x = gPad->AbsPixeltoX(gPad->GetEventX()) * INTT::LADDER[layer];
	int bin_y = gPad->AbsPixeltoY(gPad->GetEventY()) * INTT::ARM;
	std::string name;

	if(bin_x < 0)bin_x = 0;
	if(bin_y < 0)bin_y = 0;

	if(bin_x >= INTT::LADDER[layer])bin_x = INTT::LADDER[layer] - 1;
	if(bin_y >= INTT::ARM)bin_y = INTT::ARM - 1;

	INTT::Indexes_s indexes;
	indexes.lyr = layer;
	INTT::GetIndexesFromGlobalLadderBinXY(bin_x, bin_y, indexes);

	name = Form("Intt_%s_Global_DispText", option.c_str());
	TText* disp_text = (TText*)gROOT->FindObject(name.c_str());
	if(disp_text)disp_text->SetTitle(Form("Layer: %2d Ladder: %2d (%s)", indexes.lyr + INTT::LAYER_OFFSET, indexes.ldr, indexes.arm ? "North" : "South"));

	name = Form("Intt_%s_Global_DispPad", option.c_str());
	TPad* disp_pad = (TPad*)gROOT->FindObject(name.c_str());
	if(disp_pad)disp_pad->Update();

	if(gPad->GetEvent() != 11)return; //left click

	DrawLocalChipMap(option, indexes);
}

void InttMonDraw::DrawLocalChipMap(std::string const& option, struct INTT::Indexes_s indexes)
{
	double x_lower = 0.0;
	double y_lower = DISP_FRAC;
	double x_upper = 1.0;
	double y_upper = 1.0;

	std::string name;

	TStyle* style;
	TCanvas* cnvs;
	TPad* disp_pad;
	TPad* hist_pad;
	TPad* grid_pad;
	TPad* exec_pad;
	
	TH2D* client_hists[1] = {nullptr};
	Options_t::const_iterator itr = OPTIONS.find(option);
	if(itr == OPTIONS.end())return;
	(*(itr->second.local))(option, client_hists, indexes);

	name = Form("Intt_%s_Local_Style_Lyr%02d_Ldr%02d_Arm%02d_Chp%02d", option.c_str(), indexes.lyr, indexes.ldr, indexes.arm, indexes.chp);
	style = (TStyle*)gROOT->FindObject(name.c_str());
	if(!style)
	{
		style = new TStyle
		(
			name.c_str(),
			name.c_str()
		);
		style->SetOptStat(0);
		style->SetOptTitle(0);
		//...
	}
	style->cd();

	name = Form("Intt_%s_Local_Canvas_Lyr%02d_Ldr%02d_Arm%02d", option.c_str(), indexes.lyr, indexes.ldr, indexes.arm);
	cnvs = (TCanvas*)gROOT->FindObject(name.c_str());
	if(!cnvs)
	{
		cnvs = new TCanvas
		(
			name.c_str(),
			name.c_str(),
			0,
			0,
			CNVS_WIDTH,
			CNVS_HEIGHT
		);
	}

	name = Form("Intt_%s_Local_DispPad_Lyr%02d_Ldr%02d_Arm%02d", option.c_str(), indexes.lyr, indexes.ldr, indexes.arm);
	disp_pad = (TPad*)gROOT->FindObject(name.c_str());
	if(!disp_pad)
	{
		disp_pad = new TPad
		(
			name.c_str(),
			name.c_str(),
			x_lower,
			0.0,
			x_upper,
			DISP_FRAC
		);
		DrawPad(cnvs, disp_pad);
		TText* disp_text = new TText
		(
			0.5,
			0.5,
			Form("(%s) Chip: %2d Channel: %3d", "South", 0, INTT::CHIP_OFFSET)
		);
		disp_text->SetName(Form("Intt_%s_Local_DispText_Lyr%02d_Ldr%02d_Arm%02d", option.c_str(), indexes.lyr, indexes.ldr, indexes.arm));
		disp_text->SetTextAlign(22);
		disp_text->SetTextSize(DISP_TEXT_SIZE);
		disp_text->Draw();
	}

	name = Form("Intt_%s_Local_HistPad_Lyr%02d_Ldr%02d_Arm%02d", option.c_str(), indexes.lyr, indexes.ldr, indexes.arm);
	hist_pad = (TPad*)gROOT->FindObject(name.c_str());
	if(!hist_pad)
	{
		hist_pad = new TPad
		(
			name.c_str(),
			name.c_str(),
			x_lower,
			y_lower,
			x_upper,
			y_upper
		);
		DrawPad(cnvs, hist_pad);
	}

	hist_pad->cd();
	if(client_hists[0])
	{
		client_hists[0]->DrawCopy("COLZ");
		delete client_hists[0];
	}

	name = Form("Intt_%s_Local_GridPad_Lyr%d_Ldr%02d_Arm%02d", option.c_str(), indexes.lyr, indexes.ldr, indexes.arm);
	grid_pad = (TPad*)gROOT->FindObject(name.c_str());
	if(!grid_pad)
	{
		grid_pad = new TPad
		(
			name.c_str(),
			name.c_str(),
			x_lower,
			y_lower,
			x_upper,
			y_upper
		);
		DrawPad(cnvs, grid_pad);

		int i;
		double temp;

		for(i = 0; i < 2 * INTT::CHANNEL + 1; ++i)
		{
			temp = L_MARGIN + (i / 2.0) * (1.0 - R_MARGIN - L_MARGIN) / INTT::CHANNEL;
		
			TLine* line = new TLine(temp, B_MARGIN, temp, 1.0 - T_MARGIN);
			line->SetLineStyle(i == INTT::CHANNEL ? 1 : 3);
			line->SetLineWidth(i == INTT::CHANNEL ? 2 : 1);
			line->Draw();
		}

		for(i = 0; i < INTT::CHIP / 2 + 1; ++i)
		{
			temp = B_MARGIN + 2 * i * (1.0 - T_MARGIN - B_MARGIN) / INTT::CHIP;
	
			TLine* line = new TLine(L_MARGIN, temp, 1.0 - R_MARGIN, temp);	
			line->SetLineStyle(i == INTT::CHIP / 2 ? 1 : 3);
			line->SetLineWidth(i == INTT::CHIP / 2 ? 2 : 1);
			line->Draw();
		}
	}

	name = Form("Intt_%s_Local_ExecPad_Lyr%02d_Ldr%02d_Arm%02d", option.c_str(), indexes.lyr, indexes.ldr, indexes.arm);
	exec_pad = (TPad*)gROOT->FindObject(name.c_str());
	if(!exec_pad)
	{
		exec_pad = new TPad
		(
			name.c_str(),
			name.c_str(),
			x_lower + (x_upper - x_lower) * L_MARGIN,
			y_lower + (y_upper - y_lower) * B_MARGIN,
			x_upper - (x_upper - x_lower) * R_MARGIN,
			y_upper - (y_upper - y_lower) * T_MARGIN
		);
		DrawPad(cnvs, exec_pad);
		exec_pad->AddExec
		(
			Form("Intt_%s_Local_Exec_Lyr%02d_Ldr%02d_Arm%02d", option.c_str(), indexes.lyr, indexes.ldr, indexes.arm),
			Form("InttMonDraw::InttLocalChipExec(\"%s\", %d, %d, %d)", option.c_str(), indexes.lyr, indexes.ldr, indexes.arm)
		);
	}
}

void InttMonDraw::InttLocalChipExec(const std::string& option, int layer, int ladder, int arm)
{
	int bin_x = gPad->AbsPixeltoX(gPad->GetEventX()) * 2 * INTT::CHANNEL;
	int bin_y = gPad->AbsPixeltoY(gPad->GetEventY()) * INTT::CHIP / 2;
	std::string name;

	if(bin_x < 0)bin_x = 0;
	if(bin_y < 0)bin_y = 0;

	if(bin_x >= 2 * INTT::CHANNEL - 1)bin_x = 2 * INTT::CHANNEL - 1;
	if(bin_y >= INTT::CHIP / 2 - 1)bin_y = INTT::CHIP / 2 - 1;

	INTT::Indexes_s indexes;
	indexes.lyr = layer;
	indexes.ldr = ladder;
	indexes.arm = arm;
	INTT::GetIndexesFromLocalChipBinXY(bin_x, bin_y, indexes);
	//indexes.chp = indexes.arm * (INTT::CHIP / 2 - 1) - (2 * indexes.arm - 1) * bin_y + ((indexes.arm + bin_x / INTT::CHANNEL) % 2) * (INTT::CHIP / 2);
	//indexes.chn = bin_x % INTT::CHANNEL;

	name = Form("Intt_%s_Local_DispText_Lyr%02d_Ldr%02d_Arm%02d", option.c_str(), indexes.lyr, indexes.ldr, indexes.arm);
	TText* disp_text = (TText*)gROOT->FindObject(name.c_str());
	if(disp_text)disp_text->SetTitle(Form("Chip: %2d Channel: %3d", indexes.chp + INTT::CHIP_OFFSET, indexes.chn));

	name = Form("Intt_%s_Local_DispPad_Lyr%02d_Ldr%02d_Arm%02d", option.c_str(), indexes.lyr, indexes.ldr, indexes.arm);
	TPad* disp_pad = (TPad*)gROOT->FindObject(name.c_str());
	if(disp_pad)disp_pad->Update();
}

//idiom specific methods/implementations
//GlobalChip/LocalChannel methods/implementations
void InttMonDraw::PrepGlobalChipHists_Hitmap(std::string const& option, TH2D** client_hists)
{
	int bin;
	double temp;
	double min = DBL_MAX;
	double max = -DBL_MAX;
	int prev_felix = -1;
	int felix = 0;
	int felix_channel = 0;
	struct INTT_Felix::Ladder_s ldr_struct;

	int bin_x;
	int bin_y;
	std::string name;

	struct INTT::Indexes_s indexes;

	OnlMonClient* cl = OnlMonClient::instance();
	TH1D* server_hist = nullptr;

	for(indexes.lyr = 0; indexes.lyr < INTT::LAYER; ++indexes.lyr)
	{
		name = Form("Intt_%s_Global_Hist_%d", option.c_str(), indexes.lyr);
		client_hists[indexes.lyr] = (TH2D*)gROOT->FindObject(name.c_str());
		if(!client_hists[indexes.lyr])
		{
			client_hists[indexes.lyr] = new TH2D
			(
				name.c_str(),
				name.c_str(),
				2 * INTT::LADDER[indexes.lyr],
				-0.5,
				2 * INTT::LADDER[indexes.lyr] - 0.5,
				INTT::CHIP,
				-0.5,
				INTT::CHIP - 0.5
			);
			client_hists[indexes.lyr]->GetXaxis()->SetNdivisions(INTT::LADDER[indexes.lyr], true);
			client_hists[indexes.lyr]->GetYaxis()->SetNdivisions(INTT::CHIP, true);

			client_hists[indexes.lyr]->GetXaxis()->SetLabelSize(0.0);
			client_hists[indexes.lyr]->GetYaxis()->SetLabelSize(0.0);

			client_hists[indexes.lyr]->GetXaxis()->SetTickLength(0.0);
			client_hists[indexes.lyr]->GetYaxis()->SetTickLength(0.0);
		}

		client_hists[indexes.lyr]->Reset();
	}

	indexes.lyr = 0;
	indexes.ldr = 0;
	indexes.arm = 0;
	indexes.chp = 0;
	indexes.chn = 0;
	indexes.adc = 0;

	prev_felix = -1;
	felix = 0;
	felix_channel = 0;
	while(felix < INTT::FELIX)
	{
		if(prev_felix != felix)
		{
			server_hist = (TH1D*)cl->getHisto(Form("INTTMON_%d", felix), "InttMap");
			prev_felix = felix;
		}

		INTT_Felix::FelixMap(felix, felix_channel, ldr_struct);
		indexes.lyr = ldr_struct.barrel * 2 + ldr_struct.layer;
		indexes.ldr = ldr_struct.ladder;
		indexes.arm = felix / 4;
		INTT::GetGlobalChipBinXYFromIndexes(bin_x, bin_y, indexes);

		//for debugging
		//if(server_hist)temp = 0;
		//temp = indexes.chp + 1;
		//bin = client_hists[indexes.lyr]->GetBin(bin_x, bin_y);
		//client_hists[indexes.lyr]->SetBinContent(bin, temp);

		//actual implementation
		INTT::GetFelixBinFromIndexes(bin, felix_channel, indexes);
		temp = server_hist ? server_hist->GetBinContent(bin) : 0;
		bin = client_hists[indexes.lyr]->GetBin(bin_x, bin_y);
		temp += client_hists[indexes.lyr]->GetBinContent(bin);
		client_hists[indexes.lyr]->SetBinContent(bin, temp);

		++indexes.adc;
		if(indexes.adc < INTT::ADC)continue;
		indexes.adc = 0;

		++indexes.chn;
		if(indexes.chn < INTT::CHANNEL)continue;
		indexes.chn = 0;

		++indexes.chp;
		if(indexes.chp < INTT::CHIP)continue;
		indexes.chp = 0;

		++felix_channel;
		if(felix_channel < INTT::FELIX_CHANNEL)continue;
		felix_channel = 0;

		++felix;
	}

	for(indexes.lyr = 0; indexes.lyr < INTT::LAYER; ++indexes.lyr)
	{
		if(client_hists[indexes.lyr]->GetMinimum() < min)min = client_hists[indexes.lyr]->GetMinimum();
		if(client_hists[indexes.lyr]->GetMaximum() > max)max = client_hists[indexes.lyr]->GetMaximum();
	}

	for(indexes.lyr = 0; indexes.lyr < INTT::LAYER; ++indexes.lyr)
	{
		client_hists[indexes.lyr]->SetMinimum(min);
		client_hists[indexes.lyr]->SetMaximum(max);
	}
}

void InttMonDraw::PrepLocalChannelHists_Hitmap(std::string const& option, TH2D** client_hists, struct INTT::Indexes_s& indexes)
{
	int bin;
	double temp;
	int prev_felix = -1;
	int felix = 0;
	int felix_channel = 0;
	struct INTT_Felix::Ladder_s ldr_struct;

	int bin_x;
	int bin_y;
	std::string name;

	OnlMonClient* cl = OnlMonClient::instance();
	TH1D* server_hist = nullptr;

	name = Form("Intt_%s_Local_Hist_Lyr%02d_Ldr%02d_Arm%02d_Chp%02d", option.c_str(), indexes.lyr, indexes.ldr, indexes.arm, indexes.chp);
	client_hists[0] = (TH2D*)gROOT->FindObject(name.c_str());
	if(!client_hists[0])
	{
		client_hists[0] = new TH2D
		(
			name.c_str(),
			name.c_str(),
			INTT::CHANNEL,
			-0.5,
			INTT::CHANNEL - 0.5,
			INTT::ADC,
			-0.5,
			INTT::ADC - 0.5
		);

		client_hists[0]->GetXaxis()->SetNdivisions(INTT::CHANNEL, true);
		client_hists[0]->GetYaxis()->SetNdivisions(INTT::ADC, true);

		client_hists[0]->GetXaxis()->SetLabelSize(0.0);
		client_hists[0]->GetYaxis()->SetLabelSize(0.0);

		client_hists[0]->GetXaxis()->SetTickLength(0.0);
		client_hists[0]->GetYaxis()->SetTickLength(0.0);
	}
	client_hists[0]->Reset();

	indexes.chn = 0;
	indexes.adc = 0;

	prev_felix = -1;
	felix = 0;
	felix_channel = 0;
	while(felix < INTT::FELIX)
	{
		if(prev_felix != felix)
		{
			server_hist = (TH1D*)cl->getHisto(Form("INTTMON_%d", felix), "InttMap");
			prev_felix = felix;
		}

		INTT_Felix::FelixMap(felix, felix_channel, ldr_struct);
		if(2 * ldr_struct.barrel + ldr_struct.layer == indexes.lyr && ldr_struct.ladder == indexes.ldr && felix / 4 == indexes.arm)
		{
			INTT::GetLocalChannelBinXYFromIndexes(bin_x, bin_y, indexes);

			//for debugging
			//client_hist->SetBinContent(client_hist->GetBin(bin_x, bin_y), indexes.chn);
			//client_hist->SetBinContent(client_hist->GetBin(bin_x, bin_y), indexes.adc);

			INTT::GetFelixBinFromIndexes(bin, felix_channel, indexes);
			temp = server_hist ? server_hist->GetBinContent(bin) : 0;
			bin = client_hists[0]->GetBin(bin_x, bin_y);
			temp += client_hists[0]->GetBinContent(bin);
			client_hists[0]->SetBinContent(bin, temp);
		}

		++indexes.adc;
		if(indexes.adc < INTT::ADC)continue;
		indexes.adc = 0;

		++indexes.chn;
		if(indexes.chn < INTT::CHANNEL)continue;
		indexes.chn = 0;

		++felix_channel;
		if(felix_channel < INTT::FELIX_CHANNEL)continue;
		felix_channel = 0;

		++felix;
	}
}

void InttMonDraw::PrepGlobalChipHists_NLL(std::string const& option, TH2D** client_hists)
{
	int bin;
	double temp;
	double min = DBL_MAX;
	double max = -DBL_MAX;
	int prev_felix = -1;
	int felix = 0;
	int felix_channel = 0;
	struct INTT_Felix::Ladder_s ldr_struct;

	int bin_x;
	int bin_y;
	std::string name;

	struct INTT::Indexes_s indexes;
	double adc_counts[INTT::ADC] = {0};

	OnlMonClient* cl = OnlMonClient::instance();
	TH1D* server_hist = nullptr;

	for(indexes.lyr = 0; indexes.lyr < INTT::LAYER; ++indexes.lyr)
	{
		name = Form("Intt_%s_Global_Hist_%d", option.c_str(), indexes.lyr);
		client_hists[indexes.lyr] = (TH2D*)gROOT->FindObject(name.c_str());
		if(!client_hists[indexes.lyr])
		{
			client_hists[indexes.lyr] = new TH2D
			(
				name.c_str(),
				name.c_str(),
				2 * INTT::LADDER[indexes.lyr],
				-0.5,
				2 * INTT::LADDER[indexes.lyr] - 0.5,
				INTT::CHIP,
				-0.5,
				INTT::CHIP - 0.5
			);
			client_hists[indexes.lyr]->GetXaxis()->SetNdivisions(INTT::LADDER[indexes.lyr], true);
			client_hists[indexes.lyr]->GetYaxis()->SetNdivisions(INTT::CHIP, true);

			client_hists[indexes.lyr]->GetXaxis()->SetLabelSize(0.0);
			client_hists[indexes.lyr]->GetYaxis()->SetLabelSize(0.0);

			client_hists[indexes.lyr]->GetXaxis()->SetTickLength(0.0);
			client_hists[indexes.lyr]->GetYaxis()->SetTickLength(0.0);
		}

		client_hists[indexes.lyr]->Reset();
	}

	indexes.lyr = 0;
	indexes.ldr = 0;
	indexes.arm = 0;
	indexes.chp = 0;
	indexes.chn = 0;
	indexes.adc = 0;

	prev_felix = -1;
	felix = 0;
	felix_channel = 0;
	while(felix < INTT::FELIX)
	{
		if(prev_felix != felix)
		{
			server_hist = (TH1D*)cl->getHisto(Form("INTTMON_%d", felix), "InttMap");
			prev_felix = felix;
		}

		INTT::GetFelixBinFromIndexes(bin, felix_channel, indexes);
		temp = server_hist ? server_hist->GetBinContent(bin) : 0;
		adc_counts[indexes.adc] += temp;

		++indexes.adc;
		if(indexes.adc < INTT::ADC)continue;
		indexes.adc = 0;

		++indexes.chn;
		if(indexes.chn < INTT::CHANNEL)continue;
		indexes.chn = 0;

		++indexes.chp;
		if(indexes.chp < INTT::CHIP)continue;
		indexes.chp = 0;

		++felix_channel;
		if(felix_channel < INTT::FELIX_CHANNEL)continue;
		felix_channel = 0;

		++felix;
	}

	for(indexes.adc = 0; indexes.adc < INTT::ADC; ++indexes.adc)
	{
		adc_counts[indexes.adc] /= INTT::CHANNELS;
	}

	indexes.lyr = 0;
	indexes.ldr = 0;
	indexes.arm = 0;
	indexes.chp = 0;
	indexes.chn = 0;
	indexes.adc = 0;

	prev_felix = -1;
	felix = 0;
	felix_channel = 0;
	while(felix < INTT::FELIX)
	{
		if(prev_felix != felix)
		{
			server_hist = (TH1D*)cl->getHisto(Form("INTTMON_%d", felix), "InttMap");
			prev_felix = felix;
		}

		INTT_Felix::FelixMap(felix, felix_channel, ldr_struct);
		indexes.lyr = ldr_struct.barrel * 2 + ldr_struct.layer;
		indexes.ldr = ldr_struct.ladder;
		indexes.arm = felix / 4;
		INTT::GetGlobalChipBinXYFromIndexes(bin_x, bin_y, indexes);

		INTT::GetFelixBinFromIndexes(bin, felix_channel, indexes);
		temp = server_hist ? server_hist->GetBinContent(bin) : 0;

		//NLL of having temp counts from a poisson distribution with expectation adc_counts[indexes.adc]
		temp = adc_counts[indexes.adc] ? lgamma(temp + 1) - temp * log(adc_counts[indexes.adc]) + adc_counts[indexes.adc] : 0;

		bin = client_hists[indexes.lyr]->GetBin(bin_x, bin_y);
		temp += client_hists[indexes.lyr]->GetBinContent(bin);
		client_hists[indexes.lyr]->SetBinContent(bin, temp);

		++indexes.adc;
		if(indexes.adc < INTT::ADC)continue;
		indexes.adc = 0;

		++indexes.chn;
		if(indexes.chn < INTT::CHANNEL)continue;
		indexes.chn = 0;

		++indexes.chp;
		if(indexes.chp < INTT::CHIP)continue;
		indexes.chp = 0;

		++felix_channel;
		if(felix_channel < INTT::FELIX_CHANNEL)continue;
		felix_channel = 0;

		++felix;
	}

	for(indexes.lyr = 0; indexes.lyr < INTT::LAYER; ++indexes.lyr)
	{
		if(client_hists[indexes.lyr]->GetMinimum() < min)min = client_hists[indexes.lyr]->GetMinimum();
		if(client_hists[indexes.lyr]->GetMaximum() > max)max = client_hists[indexes.lyr]->GetMaximum();
	}

	for(indexes.lyr = 0; indexes.lyr < INTT::LAYER; ++indexes.lyr)
	{
		client_hists[indexes.lyr]->SetMinimum(min);
		client_hists[indexes.lyr]->SetMaximum(max);
	}
}

//GlobalLadder/LocalChip methods/implementations
void InttMonDraw::PrepGlobalLadderHists_Interface(std::string const& option, TH2D** client_hists)
{
	struct INTT::Indexes_s indexes;
	std::string name;

	for(indexes.lyr = 0; indexes.lyr < INTT::LAYER; ++indexes.lyr)
	{
		name = Form("Intt_%s_Global_Hist_%d", option.c_str(), indexes.lyr);
		client_hists[indexes.lyr] = (TH2D*)gROOT->FindObject(name.c_str());
		if(!client_hists[indexes.lyr])
		{
			client_hists[indexes.lyr] = new TH2D
			(
				name.c_str(),
				name.c_str(),
				INTT::LADDER[indexes.lyr],
				-0.5,
				INTT::LADDER[indexes.lyr] - 0.5,
				INTT::ARM,
				-0.5,
				INTT::ARM - 0.5
			);
			client_hists[indexes.lyr]->GetXaxis()->SetNdivisions(INTT::LADDER[indexes.lyr], true);
			client_hists[indexes.lyr]->GetYaxis()->SetNdivisions(INTT::ARM, true);

			client_hists[indexes.lyr]->GetXaxis()->SetLabelSize(0.0);
			client_hists[indexes.lyr]->GetYaxis()->SetLabelSize(0.0);

			client_hists[indexes.lyr]->GetXaxis()->SetTickLength(0.0);
			client_hists[indexes.lyr]->GetYaxis()->SetTickLength(0.0);

			client_hists[indexes.lyr]->SetMinimum(-1.0);
			client_hists[indexes.lyr]->SetMaximum(20.0);
		}

		client_hists[indexes.lyr]->Reset();

		for(int n = 0; n < client_hists[indexes.lyr]->GetNcells(); ++n)
		{
			client_hists[indexes.lyr]->SetBinContent(n, 0);
		}
	}
}

void InttMonDraw::PrepLocalChipHists_Hitmap(std::string const& option, TH2D** client_hists, struct INTT::Indexes_s& indexes)
{
	int bin;
	double temp;
	int prev_felix = -1;
	int felix = 0;
	int felix_channel = 0;
	struct INTT_Felix::Ladder_s ldr_struct;

	int bin_x;
	int bin_y;
	std::string name;

	OnlMonClient* cl = OnlMonClient::instance();
	TH1D* server_hist = nullptr;

	name = Form("Intt_%s_Local_ClientHist_Lyr%02d_Ldr%02d_Arm%02d", option.c_str(), indexes.lyr, indexes.ldr, indexes.arm);
	client_hists[0] = (TH2D*)gROOT->FindObject(name.c_str());
	if(!client_hists[0])
	{
		client_hists[0] = new TH2D
		(
			name.c_str(),
			name.c_str(),
			2 * INTT::CHANNEL,
			-0.5,
			2 * INTT::CHANNEL - 0.5,
			INTT::CHIP / 2,
			-0.5,
			INTT::CHIP / 2 - 0.5
		);

		client_hists[0]->GetXaxis()->SetNdivisions(2 * INTT::CHANNEL, true);
		client_hists[0]->GetYaxis()->SetNdivisions(INTT::CHIP / 2, true);

		client_hists[0]->GetXaxis()->SetLabelSize(0.0);
		client_hists[0]->GetYaxis()->SetLabelSize(0.0);

		client_hists[0]->GetXaxis()->SetTickLength(0.0);
		client_hists[0]->GetYaxis()->SetTickLength(0.0);
	}

	client_hists[0]->Reset();

	indexes.chp = 0;
	indexes.chn = 0;
	indexes.adc = 0;

	prev_felix = -1;
	felix = 0;
	felix_channel = 0;
	while(felix < INTT::FELIX)
	{
		if(prev_felix != felix)
		{
			server_hist = (TH1D*)cl->getHisto(Form("INTTMON_%d", felix), "InttMap");
			prev_felix = felix;
		}

		INTT_Felix::FelixMap(felix, felix_channel, ldr_struct);
		if(2 * ldr_struct.barrel + ldr_struct.layer == indexes.lyr && ldr_struct.ladder == indexes.ldr && felix / 4 == indexes.arm)
		{
			INTT::GetLocalChipBinXYFromIndexes(bin_x, bin_y, indexes);

			//for debugging
			//temp = server_hist ? server_hist->GetBinContent(bin) : 0;
			//client_hists[0]->SetBinContent(client_hists[0]->GetBin(bin_x, bin_y), indexes.chp);

			INTT::GetFelixBinFromIndexes(bin, felix_channel, indexes);
			temp = server_hist ? server_hist->GetBinContent(bin) : 0;
			INTT::GetLocalChipBinXYFromIndexes(bin_x, bin_y, indexes);
			bin = client_hists[0]->GetBin(bin_x, bin_y);
			temp += client_hists[0]->GetBinContent(bin);
			client_hists[0]->SetBinContent(bin, temp);
		}

		++indexes.adc;
		if(indexes.adc < INTT::ADC)continue;
		indexes.adc = 0;

		++indexes.chn;
		if(indexes.chn < INTT::CHANNEL)continue;
		indexes.chn = 0;

		++indexes.chp;
		if(indexes.chp < INTT::CHIP)continue;
		indexes.chp = 0;

		++felix_channel;
		if(felix_channel < INTT::FELIX_CHANNEL)continue;
		felix_channel = 0;

		++felix;
	}
}
