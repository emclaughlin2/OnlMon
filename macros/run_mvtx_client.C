#include "CommonFuncs.C"

#include <onlmon/mvtx/MvtxMonDraw.h>

#include <onlmon/OnlMonClient.h>

// cppcheck-suppress unknownMacro
R__LOAD_LIBRARY(libonlmvtxmon_client.so)

void mvtxDrawInit(const int online = 0)
{
  OnlMonClient *cl = OnlMonClient::instance();
  // register histos we want with monitor name
  //cl->registerHisto("mvtxmon_hist1", "MVTXMON");
  //cl->registerHisto("mvtxmon_hist2", "MVTXMON");
  std::string mLaneStatusFlag[3] = { "WARNING", "ERROR", "FAULT" };
  const int NStaves[3] = { 12, 16, 20 };

  for(int iflx = 0; iflx < 6; iflx++){
    std::string instanceName = "MVTXMON_" + std::to_string(iflx);
    std::cout<<instanceName<<std::endl;
    cl->registerHisto("MVTXMON_RawData_ChipStaveOcc", instanceName);
    cl->registerHisto("MVTXMON_RawData_ChipStave1D", instanceName);
    cl->registerHisto("MVTXMON_RawData_ChipFiredHis", instanceName);
    cl->registerHisto("MVTXMON_RawData_EvtHitChip", instanceName);
    cl->registerHisto("MVTXMON_RawData_EvtHitDis", instanceName);
    cl->registerHisto("MVTXMON_General_Occupancy", instanceName);
    cl->registerHisto("MVTXMON_LaneStatus_laneStatusOverviewFlagWARNING", instanceName);
    cl->registerHisto("MVTXMON_LaneStatus_laneStatusOverviewFlagERROR", instanceName);
    cl->registerHisto("MVTXMON_LaneStatus_laneStatusOverviewFlagFAULT", instanceName);

    cl->registerHisto("hChipStrobes", instanceName);
    cl->registerHisto("hChipL1", instanceName);

    //FEE
    cl->registerHisto("MVTXMON_FEE_TriggerVsFeeid", instanceName);
    cl->registerHisto("MVTXMON_FEE_TriggerFlag", instanceName);
    //cl->registerHisto("MVTXMON/FEE/LaneInfo", instanceName);

    for (int i = 0; i < MvtxMonDraw::NFlags; i++) {
      cl->registerHisto(Form("MVTXMON_LaneStatus_laneStatusFlag%s", mLaneStatusFlag[i].c_str()), instanceName);
      cl->registerHisto(Form("MVTXMON_LaneStatus_laneStatusFlagCumulative%s", mLaneStatusFlag[i].c_str()), instanceName);
    }

    for(int i = 0; i < 3; i++){
      cl->registerHisto(Form("MVTXMON_LaneStatusSummary_LaneStatusSummaryL%i", i), instanceName);
    }

    cl->registerHisto("MVTXMON_LaneStatusSummary_LaneStatusSummary", instanceName);

    //raw task
    cl->registerHisto("MVTXMON_General_ErrorPlots", instanceName);
    cl->registerHisto("MVTXMON_General_ErrorFile", instanceName);

    for (int aLayer = 0; aLayer < 3; aLayer++) {
      cl->registerHisto(Form("MVTXMON_Occupancy_Layer%dOccupancy", aLayer), instanceName);
      cl->registerHisto(Form("MVTXMON_Occupancy_Layer%d_Layer%dChipStave", aLayer, aLayer), instanceName);
      cl->registerHisto(Form("MVTXMON_Noisy_Layer%d_ChipStave", aLayer), instanceName);
    }


    //fhr
    cl->registerHisto("MVTXMON_General_ErrorVsFeeid", instanceName);
    //cl->registerHisto("MVTXMON_General_Occupancy", instanceName);
    cl->registerHisto("MVTXMON_General_Noisy_Pixel", instanceName);
    cl->registerHisto("RCDAQ_evt", instanceName);

    for (int mLayer = 0; mLayer < 3; mLayer++) {
      cl->registerHisto(Form("MVTXMON_Occupancy_Layer%d_Layer%dDeadChipPos", mLayer, mLayer), instanceName);
      cl->registerHisto(Form("MVTXMON_Occupancy_Layer%d_Layer%dAliveChipPos", mLayer, mLayer), instanceName);
      //cl->registerHisto(Form("MVTXMON/Occupancy/Layer%d/Layer%dChipStaveC", mLayer, mLayer), instanceName);
      cl->registerHisto(Form("MVTXMON_Occupancy_Layer%dOccupancy_LOG", mLayer), instanceName);
    }

    cl->registerHisto("MVTXMON_Occupancy_TotalDeadChipPos", instanceName);
    cl->registerHisto("MVTXMON_Occupancy_TotalAliveChipPos", instanceName);

     cl->registerHisto(Form("MVTXMON_chipHitmapFLX%d", iflx), instanceName);
  }

   


  //cl->AddServerHost("localhost");  // check local host first
  CreateHostList(online);
  // get my histos from server, the second parameter = 1
  // says I know they are all on the same node

  cl->requestHistoBySubSystem("MVTXMON_0", 1);
  cl->requestHistoBySubSystem("MVTXMON_1", 1);
  cl->requestHistoBySubSystem("MVTXMON_2", 1);
  cl->requestHistoBySubSystem("MVTXMON_3", 1);
  cl->requestHistoBySubSystem("MVTXMON_4", 1);
  cl->requestHistoBySubSystem("MVTXMON_5", 1);
  OnlMonDraw *mvtxmon = new MvtxMonDraw("MVTXMONDRAW");  // create Drawing Object
  cl->registerDrawer(mvtxmon);              // register with client framework
}

void mvtxDraw(const char *what = "ALL")
{
  OnlMonClient *cl = OnlMonClient::instance();  // get pointer to framewrk
  cl->requestHistoBySubSystem("MVTXMON_0", 1);
  cl->requestHistoBySubSystem("MVTXMON_1", 1);
  cl->requestHistoBySubSystem("MVTXMON_2", 1);
  cl->requestHistoBySubSystem("MVTXMON_3", 1);
  cl->requestHistoBySubSystem("MVTXMON_4", 1);
  cl->requestHistoBySubSystem("MVTXMON_5", 1);
  cl->Draw("MVTXMONDRAW", what);                    // Draw Histos of registered Drawers
}

void mvtxSavePlot()
{
  OnlMonClient *cl = OnlMonClient::instance();  // get pointer to framewrk
  cl->SavePlot("MVTXMONDRAW");                        // Save Plots
  return;
}

void mvtxHtml()
{
  OnlMonClient *cl = OnlMonClient::instance();  // get pointer to framewrk
  cl->MakeHtml("MVTXMONDRAW");                      // Create html output
  return;
}
