#include <onlmon/OnlMonServer.h>
#include <pmonitor/pmonitor.h>

// cppcheck-suppress unknownMacro
R__LOAD_LIBRARY(libonlmonserver.so)
// cppcheck-suppress unknownMacro
R__LOAD_LIBRARY(libonlmonserver_funcs.so)
void CleanUpServer();

void start_server(const std::string &prdffile = "")
{
  OnlMonServer *se = OnlMonServer::instance();  // get pointer to Server Framework
  if (prdffile.empty())
  {
    cout << "No Input file given" << endl;
    return;
  }
  if (prdffile.find("seb") == 0 || prdffile.find("ebdc") == 0 || prdffile.find("intt") == 0 || prdffile.find("mvtx") == 0)
  {
    pidentify(0);
    rcdaqopen(prdffile.c_str());
    prun();
    //  if the rcdaq server is terminated we execute the EndRun and then
    // save the histograms
    se->EndRun(se->RunNumber());  //
    se->WriteHistoFile();
    //      delete enablecorbabuf;
    CleanUpServer();
  }
  else
  {
    pfileopen(prdffile.c_str());
  }
  return;
}

void CleanUpServer()
{
  pclose();
  OnlMonServer *se = OnlMonServer::instance();
  delete se;
  gSystem->Exit(0);
}
