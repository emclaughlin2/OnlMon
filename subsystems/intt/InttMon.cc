#include "InttMon.h"

//===	Public Methods		===//
InttMon::InttMon(const std::string &name) : OnlMon(name)
{
	//leaving fairly empty
	return;
}

InttMon::~InttMon()
{
	delete dbvars;
}

int InttMon::Init()
{
	OnlMonServer *omc = OnlMonServer::instance();

	//dbvars
	dbvars = new OnlMonDB(ThisName);
	DBVarInit();

	//histograms
	NumEvents = new TH1D(Form("InttNumEvents"), Form("InttNumEvents"), 1, 0, 1);
	HitMap = new TH1D(Form("InttMap"), Form("InttMap"), INTT::ADCS, 0, INTT::ADCS);
	//...

	omc->registerHisto(this, NumEvents);
	omc->registerHisto(this, HitMap);
	//...

	//Read in calibrartion data from InttMonData.dat
	const char *inttcalib = getenv("INTTCALIB");
	if (!inttcalib)
	{
	  std::cout << "INTTCALIB environment variable not set" << std::endl;
	  exit(1);
	}
	std::string fullfile = std::string(inttcalib) + "/" + "InttMonData.dat";
	std::ifstream calib(fullfile);
	//probably need to do stuff here (maybe write to expectation maps)
	//or reimplment in BeginRun()
	calib.close();

	// for testing/debugging without unpacker, remove later
	rng = new TRandom(1234);
	//~for testing/debugging without unpacker, remove later

	Reset();

	return 0;
}

int InttMon::BeginRun(const int /* run_num */)
{
	//per-run calibrations; don't think we need to do anything here yet

	return 0;
}

int InttMon::process_event(Event* evt)
{
	int bin;
	int N;
	int n;

	int pid = 3001;
	int felix;
	int felix_channel;
	struct INTT_Felix::Ladder_s lddr_s;

	struct INTT::Indexes_s indexes;

	for(pid = 3001; pid < 3009; ++pid)
	{
		felix = pid - 3001;

		Packet* p = evt->getPacket(pid);
		if(!p)continue;

		N = p->iValue(0, "NR_HITS");

		//p->identify();
		//if(N)std::cout << N << std::endl;

		for(n = 0; n < N; ++n)
		{
			felix_channel = p->iValue(n, "FEE");

			INTT_Felix::FelixMap(felix, felix_channel, lddr_s);

			indexes.lyr = lddr_s.barrel * 2 + lddr_s.layer;
			indexes.ldr = lddr_s.ladder;

			indexes.arm = (felix / 4) % 2;
			indexes.chp = p->iValue(n, "CHIP_ID") % 26;
			indexes.chn = p->iValue(n, "CHANNEL_ID");
			indexes.adc = p->iValue(n, "ADC");

			//std::cout << "\t" << indexes.lyr << std::endl;
			//std::cout << "\t" << indexes.ldr << std::endl;
			//std::cout << "\t" << indexes.arm << std::endl;
			//std::cout << "\t" << indexes.chp << std::endl;
			//std::cout << "\t" << indexes.chn << std::endl;
			//std::cout << std::endl;

			INTT::GetFelixBinFromIndexes(bin, felix_channel, indexes);
			if(bin < 0)
			{
				std::cout << "n: " << n << std::endl;
				std::cout << "bin: " << bin << std::endl;
				std::cout << "lyr: " << indexes.lyr << std::endl;
				std::cout << "ldr: " << indexes.ldr << std::endl;
				std::cout << "arm: " << indexes.arm << std::endl;
				std::cout << "chp: " << indexes.chp << std::endl;
				std::cout << "chn: " << indexes.chn << std::endl;
				std::cout << "adc: " << indexes.adc << std::endl;

				break;
			}
			HitMap->AddBinContent(bin);
		}

		delete p;
	}

	NumEvents->AddBinContent(1);

	DBVarUpdate();

	return 0;
}

int InttMon::Reset()
{
	//reset our DBVars
	evtcnt = 0;
	
	//clear our histogram entries
	NumEvents->Reset();
	HitMap->Reset();

	return 0;
}
//===	~Public Methods		===//

//===	Private Methods		===//
int InttMon::DBVarInit()
{
	std::string var_name;

	var_name = "intt_evtcnt";
	dbvars->registerVar(var_name);

	dbvars->DBInit();

	return 0;
}

int InttMon::DBVarUpdate()
{
	dbvars->SetVar("intt_evtcnt", (float)evtcnt, 0.1 * evtcnt, (float)evtcnt);

	return 0;
}
//===	~Private Methods		===//

// for testing/debugging
void InttMon::RandomEvent(int felix)
{
	int bin;

	int felix_channel;
	struct INTT::Indexes_s indexes;
	struct INTT_Felix::Ladder_s ldr_struct;

	int hits = rng->Poisson(16);
	for(int hit = 0; hit < hits; ++hit)
	{
		felix_channel = rng->Uniform(INTT::FELIX_CHANNEL);
		if(felix_channel == INTT::FELIX_CHANNEL)felix_channel -= 1;

		INTT_Felix::FelixMap(felix, felix_channel, ldr_struct);
		indexes.lyr = 2 * ldr_struct.barrel + ldr_struct.layer;
		indexes.ldr = ldr_struct.ladder;
		indexes.arm = felix / 4;

		indexes.chp = rng->Uniform(INTT::CHIP);
		if(indexes.chp == INTT::CHIP)indexes.chp -= 1;

		indexes.chn = rng->Uniform(INTT::CHANNEL);
		if(indexes.chn == INTT::CHANNEL)indexes.chn -= 1;

		indexes.adc = rng->Uniform(INTT::ADC);
		if(indexes.adc == INTT::ADC)indexes.adc -= 1;

		INTT::GetFelixBinFromIndexes(bin, felix_channel, indexes);
		HitMap->SetBinContent(bin, HitMap->GetBinContent(bin) + 1);
		NumEvents->AddBinContent(1);

		printf("Layer:%2d\tLadder:%3d (%s)\tChip:%3d\tChannel:%4d\n", indexes.lyr, indexes.ldr, indexes.arm ? "North" : "South", indexes.chp, indexes.chn);
	}
}

int InttMon::MiscDebug()
{
	int b = 0;
	int c = -1;

	int felix_channel = 0;
	int gelix_channel = -1;

	struct INTT::Indexes_s indexes = (struct INTT::Indexes_s){.lyr = 0, .ldr = 0, .arm = 0, .chp = 0, .chn = 0, .adc = 0};
	struct INTT::Indexes_s jndexes = (struct INTT::Indexes_s){.lyr = -1, .ldr = -1, .arm = -1, .chp = -1, .chn = -1, .adc = -1};

	while(true)
	{
		INTT::GetFelixBinFromIndexes(b, felix_channel, indexes);
		INTT::GetFelixIndexesFromBin(b, gelix_channel, jndexes);
		INTT::GetFelixBinFromIndexes(c, gelix_channel, jndexes);

		if(b != c)
		{
			std::cout << "Round trip failed" << std::endl;
			std::cout << "bin: " << b << " -> " << c << std::endl;
			std::cout << "felix_channel: " << felix_channel << " -> " << gelix_channel << std::endl;
			std::cout << "chp: " << indexes.chp << " -> " << jndexes.chp << std::endl;
			std::cout << "chn: " << indexes.chn << " -> " << jndexes.chn << std::endl;
			std::cout << "adc: " << indexes.adc << " -> " << jndexes.adc << std::endl;

			return 0;
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

		break;
	}

	std::cout << "Felix Round trip worked" << std::endl;

	return 0;
}


