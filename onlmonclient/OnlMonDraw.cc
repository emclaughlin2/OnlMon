#include "OnlMonDraw.h"

#include <TPad.h>
#include <TText.h>

#include <iostream>

OnlMonDraw::OnlMonDraw(const std::string &name)
  : ThisName(name)
{
}

int OnlMonDraw::Draw(const std::string & /* what */)
{
  std::cout << "Draw not implemented by daughter class" << std::endl;
  return -1;
}

int OnlMonDraw::SavePlot(const std::string & /* what */, const std::string & /* type */)
{
  std::cout << "SavePlot not implemented by daughter class" << std::endl;
  return -1;
}

int OnlMonDraw::MakePS(const std::string & /* what */)
{
  std::cout << "MakePS not implemented by daughter class" << std::endl;
  return -1;
}

int OnlMonDraw::MakeHtml(const std::string & /* what */)
{
  std::cout << "MakeHtml not implemented by daughter class" << std::endl;
  return -1;
}

int OnlMonDraw::DrawDeadServer(TPad *transparentpad)
{
  transparentpad->cd();
  TText FatalMsg;
  FatalMsg.SetTextFont(62);
  FatalMsg.SetTextSize(0.1);
  FatalMsg.SetTextColor(4);
  FatalMsg.SetNDC();          // set to normalized coordinates
  FatalMsg.SetTextAlign(23);  // center/top alignment
  FatalMsg.DrawText(0.5, 0.9,  Name().c_str());
  FatalMsg.SetTextAlign(22);  // center/center alignment
  FatalMsg.DrawText(0.5, 0.5, "SERVER");
  FatalMsg.SetTextAlign(21);  // center/bottom alignment
  FatalMsg.DrawText(0.5, 0.1, "DEAD");
  transparentpad->Update();
  return 0;
}
