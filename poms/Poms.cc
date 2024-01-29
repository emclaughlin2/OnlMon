/****************************************************************************\
 *
 *    PomsClient.C -- A program to graphically launch and control all online
 *        monitoring for PHENIX subsystems.
 *    This program requires the ROOT runtime environment and each subsystems
 *        individual online monitoring client macros.
 *    
 *    Program developed at Brookhaven National Laboratory for the PHENIX
 *        collaboration to be used by offline shifters.
 *
 *
 *
 *    Copyright 2003 Phenix Collaboration
 *
 *
 *    M. McCain (mcm99c@acu.edu)    January 2003
 *
 *    Credit to Ryan Roth (LANL) for coding style and example of ROOT GUI
 *    in his implementation of HVDP (High Voltage Display Program)
 ****************************************************************************
 */
#include "Poms.h"

#include <onlmon/OnlMonClient.h>

#include <GuiTypes.h>  // for kHorizontalFrame, kRaisedFrame
#include <TCanvas.h>
#include <TGButton.h>
#include <TGClient.h>  // for TGClient, gClient
#include <TGLayout.h>  // for TGLayoutHints, kLHintsTop, kLHintsE...
#include <TGMenu.h>
#include <TGMsgBox.h>
#include <TGString.h>  // for TGHotString
#include <TList.h>     // for TList
#include <TROOT.h>
#include <TSeqCollection.h>      // for TSeqCollection
#include <TString.h>             // for TString
#include <WidgetMessageTypes.h>  // for GET_MSG, GET_SUBMSG, kCM_BUTTON

#include <cmath>
#include <cstdlib>
#include <cstring>  // for strlen
#include <iostream>

const int SUBSYSTEM_ACTION_ID_BEGIN = 10001;
//const ULong_t SHUTTER_ITEM_BG_COLOR = 0xffaacc;  // color is rrggbb (red, green, blue in hex)

const char* pomsFileTypes[] =
    {
        "PHENIX raw data files", "*.prdf",
        "All files", "*",
        nullptr, nullptr};

/////////////////////////////////////////////////////////////////////////////
// PomsMainFrame Implementation                                            //
/////////////////////////////////////////////////////////////////////////////

// Init Singleton instance
PomsMainFrame* PomsMainFrame::_instance = nullptr;

PomsMainFrame* PomsMainFrame::Instance()
{
  if (!_instance)
    _instance = new PomsMainFrame(gClient->GetRoot(), 1, 1);

  return _instance;
}

PomsMainFrame::PomsMainFrame(const TGWindow* p, UInt_t w, UInt_t h)
  : TGMainFrame(p, w, h)
  , _rootHorizPad(10)
  , _rootVertPad(0)
  , _windowPad(20)
{
  // Constructor sets up window widgets

  std::cout << POMS_VER << "PomsMainFrame constructor called..." << std::endl;

  // Find macro directory
  if (getenv("ONLMON_MACROS"))
  {
    _macroPath = getenv("ONLMON_MACROS");
  }
  else
  {
    std::cout << "Environment variable ONLMON_MACROS not set, using current dir" << std::endl;
    _macroPath = "./";
  }
  // Discover root window properties
  TGFrame* rootWin = (TGFrame*) gClient->GetRoot();
  _rootWidth = rootWin->GetDefaultWidth();
  if (_rootWidth > 2000)
  {
    _rootWidth = _rootWidth / 2;
  }
  _rootHeight = rootWin->GetDefaultHeight();

  std::cout << POMS_VER << "Screen Width: " << _rootWidth << std::endl;
  std::cout << POMS_VER << "Screen Height: " << _rootHeight << std::endl;

  TGLayoutHints* menuLayout = new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0);

  // TODO:
  //   Add close button

  // Setup MenuBar
  _menuFile = new TGPopupMenu(gClient->GetRoot());
  _menuFile->AddEntry("&Exit", M_FILE_EXIT);
  _menuFile->Associate(this);
  _menuBar = new TGMenuBar(this, 1, 1, kHorizontalFrame | kRaisedFrame);
  _menuBar->AddPopup("&File", _menuFile, menuLayout);

  _menuWindow = new TGPopupMenu(gClient->GetRoot());
  _menuWindow->AddEntry("&Align Control Bar to Right", M_WINDOW_ALIGNRIGHT);
  _menuWindow->AddEntry("&Tile All", M_WINDOW_TILEALL);
  _menuWindow->Associate(this);
  _menuBar->AddPopup("&Window", _menuWindow, menuLayout);
  _menuBar->AddPopup("&DontpushThis", _menuFile, menuLayout);
  AddFrame(_menuBar, new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX,
                                       0, 0, 1, 1));

  // Standard Window Setup functions
  SetWindowName("POMS: PHENIX Online Monitoring System");
  _shutter = nullptr;
}

void PomsMainFrame::SetMacroPath(const char* path)
{
  _macroPath = path;
}

void PomsMainFrame::Draw()
{
  OnlMonClient* cl = OnlMonClient::instance();
  // Let macro know what's happening
  std::cout << POMS_VER << "Attempting Build of PomsMainFrame..." << std::endl;

  // Build Subsystem Shutter
  if (_shutter)
  {
    delete _shutter;
  }
  _shutter = BuildShutter();
  AddFrame(_shutter, new TGLayoutHints(kLHintsTop | kLHintsExpandX | kLHintsExpandY));

  MapSubwindows();
  MapWindow();
  AlignRight();
  int xsize = cl->GetDisplaySizeX();
  if (xsize > 2000)
  {
    xsize = xsize / 2;
  }
  int ysize = cl->GetDisplaySizeY();
  xsize -= (_shutter->GetDefaultWidth()) + 30;
  ysize -= 100;
  cl->SetDisplaySizeX(xsize);
  cl->SetDisplaySizeY(ysize);
}

PomsMainFrame::~PomsMainFrame()
{
  // Destructor gets rid of widgets

  delete _closeButton;
  delete _menuBar;
  delete _menuFile;
  delete _shutter;
  delete _menuWindow;

  //TODO:  Add code to clean up threads

  _instance = nullptr;
}

void PomsMainFrame::CloseWindow()
{
  // Standard close window routine
  TGMainFrame::CloseWindow();
  std::cout << "\n\n"
       << std::endl;  // Prevent prompt from displaying at end of output
  gROOT->ProcessLine(".q");
}

Bool_t PomsMainFrame::ProcessMessage(Long_t msg, Long_t parm1, Long_t /* parm2 */)
{
  // Function that defines what to do when window widgets are activated

  int status = 1;

#ifdef DEBUG

  // Debug widgets
  std::cout << "Msg = " << msg << std::endl
       << "GET_MSG = " << GET_MSG(msg) << std::endl
       << "SUB_MSG = " << GET_SUBMSG(msg) << std::endl
       << "parm1 = " << parm1 << std::endl
       << "parm2 = " << parm2 << std::endl;

#endif

  switch (GET_MSG(msg))
  {
  case kC_COMMAND:  //command type event

    switch (GET_SUBMSG(msg))
    {
    case kCM_MENU:
      switch (parm1)
      {
      case M_FILE_EXIT:
        CloseWindow();
        break;

      case M_WINDOW_ALIGNRIGHT:
        AlignRight();
        break;
      case M_WINDOW_TILEALL:
        TileAllCanvases();
        break;

      default:
        break;
      }
      break;

    case kCM_BUTTON:
      switch (parm1)
      {
      default:
        HandleButtonPoms(parm1);
        break;
      }
    }
    break;
  }

  int retval;
  if (status != 1)  // ERROR -- Create Message Window
  {
    new TGMsgBox(gClient->GetRoot(), this,
                 POMS_VER, "General error in executing widget handler!",
                 kMBIconStop, kMBOk, &retval);
  }

  return kTRUE;
}

int PomsMainFrame::HandleButtonPoms(Long_t parm1)
{
  // Check to see if button belongs to SubSystemAction
  if (parm1 > (SUBSYSTEM_ACTION_ID_BEGIN - 1))
  {
    SubSystemAction* action = SubSystemAction::FindById(parm1);

    if (!action)
      return -1;

    return action->Execute();
  }

  return -1;
}

SubSystem* PomsMainFrame::RegisterSubSystem(const char* name, const char* prefix,
                                            int addDefaultActions, int loadLibrary)
{
  SubSystem* sub = nullptr;

  try
  {
    sub = new SubSystem(name, prefix, loadLibrary);

    if (addDefaultActions != 0)
      sub->AddDefaultActions();

    _subSystemList.push_back(sub);
    std::cout << POMS_VER << "SubSystem " << name << " added..." << std::endl;
  }
  catch (char* str)
  {
    std::cout << POMS_VER << "Unable to add subsystem " << name << "!" << std::endl;
    std::cout << "\t" << str << std::endl;
    delete str;
    delete sub;
  }
  return sub;
}

SubSystem* PomsMainFrame::RegisterSubSystem(SubSystem* subSystem)
{
  _subSystemList.push_back(subSystem);
  std::cout << POMS_VER << "SubSystem " << subSystem->GetName() << " added..." << std::endl;
  return subSystem;
}

TGShutter* PomsMainFrame::BuildShutter()
{
  TGShutter* shutter = new TGShutter(this);
  TGShutterItem* shutterItem;
  TGCompositeFrame* container;
  TGTextButton* button;
  TGLayoutHints* layout;

  SubSystemList::iterator subSystem;
  SubSystemActionList* actionList;
  SubSystemActionList::iterator action;
  int shutterItemId = 101;

  std::cout << POMS_VER << "Building SubSystem Shutter..." << std::endl;

  // Describe layout for buttons in shutter
  layout = new TGLayoutHints(kLHintsExpandX | kLHintsTop, 5, 5, 0, 0);

  for (subSystem = _subSystemList.begin(); subSystem != _subSystemList.end(); ++subSystem)
  {
    actionList = (*subSystem)->GetActions();
    action = actionList->begin();

    // Only add SubSystem if it has actions associated with it.
    if (action != actionList->end())
    {
      std::cout << POMS_VER << "\tAdding " << (*subSystem)->GetName() << " to shutter" << std::endl;

      shutterItem = new TGShutterItem(shutter,
                                      new TGHotString((*subSystem)->GetName().c_str()),
                                      shutterItemId++);
      container = (TGCompositeFrame*) shutterItem->GetContainer();

      for (; action != actionList->end(); ++action)
      {
        std::cout << POMS_VER << "\t\tAdding \"" << (*action)->GetDescription() << "\" button..." << std::endl;
        button = new TGTextButton(container,
                                  (*action)->GetDescription().c_str(),
                                  (*action)->GetId());
        button->SetTextColor(0xCC00FF);
        container->AddFrame(button, layout);
        button->Associate(this);
      }
      shutter->AddItem(shutterItem);
    }
  }
  return shutter;
}

void PomsMainFrame::AlignRight()
{
  Int_t x = (_rootWidth - GetDefaultWidth() - _rootHorizPad);
  Int_t y = _rootVertPad;

  UInt_t width = GetDefaultWidth();
  UInt_t height = _rootHeight - (2 * _rootVertPad) - 90;

  Resize(width, height);
  Move(x, y);
}

void PomsMainFrame::TileCanvases(TList* canvasList)
{
  if (!canvasList)
  {
    std::cout << POMS_VER << "cannot tile canvases, canvas list empty!" << std::endl;
    return;
  }
  AlignRight();

  TCanvas* canvas = (TCanvas*) canvasList->First();
  int windowCount = canvasList->LastIndex() + 1;

  int windowCountHoriz = (int) ceil(sqrt((double) windowCount));
  int windowCountVert = (int) floor(sqrt((double) windowCount));

  int i;
  int j;

  int currX = _rootHorizPad;

  int width = (int) ((_rootWidth - (2 * _rootHorizPad) -
                      (windowCountHoriz * _windowPad) - GetDefaultWidth()) /
                     windowCountHoriz);
  int height = (int) ((_rootHeight - (2 * _rootVertPad) - (windowCountVert * _windowPad)) / windowCountVert);

  for (i = 0; i < windowCountHoriz; i++)
  {
    // Reset Vertical Height
    int currY = _rootVertPad;
    currX += (width * i) + _windowPad;

    for (j = 0; j < windowCountVert; j++)
    {
      currY += (height * j) + _windowPad;

      if (canvas)
      {
        canvas->SetWindowSize(width, height);
        canvas->SetWindowPosition(currX, currY);

        canvas = (TCanvas*) canvasList->After(canvas);
      }
    }
  }
}

void PomsMainFrame::CascadeCanvases(TList* /* canvasList */)
{
}

void PomsMainFrame::TileAllCanvases()
{
  // Scan all subsystems, grab their canvases, tile all together
  SubSystemList::iterator subSystem;
  TList* canvasList = new TList();

  for (subSystem = _subSystemList.begin(); subSystem != _subSystemList.end(); ++subSystem)
  {
    canvasList->AddAll((*subSystem)->GetCanvases());
    (*subSystem)->ShowCanvases();
  }
  TileCanvases(canvasList);

  delete canvasList;
}

/////////////////////////////////////////////////////////////////////////////
// SubSystem Implementation                                                //
/////////////////////////////////////////////////////////////////////////////

SubSystem::SubSystem(const char* name, const char* prefix, int loadLibrary)
  : _name(name)
  , _prefix(prefix)
  , _canvasList(nullptr)
  , _initialized(0)
{
  std::string macroPath;

  if ((strlen(name) < 1) || (strlen(prefix) < 1))
  {
    const char* error = "ERROR: name and prefix must not be null!";
    throw error;
  }

  if (loadLibrary)
  {
    macroPath = PomsMainFrame::Instance()->GetMacroPath();
    if (macroPath.size() == 0)
      macroPath = ".";
    gROOT->LoadMacro((macroPath + "/run_" + _prefix + "_client.C").c_str());
  }
}

SubSystem::~SubSystem()
{
  delete _canvasList;
}

TList* SubSystem::GetCanvases(int forceReQuery)
{
  // If canvasList already exists, only requery if forced
  if (_canvasList && !forceReQuery)
    return _canvasList;

  delete _canvasList;

  TSeqCollection* allCanvases = gROOT->GetListOfCanvases();
  TCanvas* canvas = (TCanvas*) allCanvases->First();
  TString* prefix = new TString(_prefix.c_str());

  prefix->ToLower();  // Avoid case sensitive searching

  while (canvas)
  {
    TString* canvasName = new TString(canvas->GetName());
    canvasName->ToLower();

    if (canvasName->Contains(*prefix))
    {
      // Canvas belongs to this sub system
      if (!_canvasList)
      {
        // only create canvasList if we have canvases
        _canvasList = new TList();
      }
      _canvasList->Add(canvas);
    }

    delete canvasName;
    canvas = (TCanvas*) allCanvases->After(canvas);
  }

  delete prefix;

  if (!_canvasList)
    std::cout << POMS_VER << "Canvas list empty for subsystem " << _name << "!" << std::endl;

  return _canvasList;
}

void SubSystem::PrintCanvasList()
{
  TList* canvasList;
  TCanvas* canvas;

  if (!(canvasList = GetCanvases()))
    return;

  std::cout << POMS_VER << "Querying subsystem " << _name << " for canvases:" << std::endl;
  canvas = (TCanvas*) canvasList->First();
  while (canvas)
  {
    std::cout << "\t" << canvas->GetName() << std::endl;
    canvas = (TCanvas*) canvasList->After(canvas);
  }
  std::cout << "End of Canvas List" << std::endl;
  return;
}

void SubSystem::ShowCanvases()
{
  TList* canvasList;
  TCanvas* canvas;

  if (!(canvasList = GetCanvases()))
    return;

  canvas = (TCanvas*) canvasList->First();
  while (canvas)
  {
    canvas->Show();
    canvas = (TCanvas*) canvasList->After(canvas);
  }
  return;
}

SubSystemAction* SubSystem::AddAction(const char* cmd, const char* description)
{
  SubSystemAction* action = nullptr;

  try
  {
    action = new SubSystemAction(this, cmd, description);
    _actions.push_back(action);

    std::cout << POMS_VER << "Action " << cmd << " added to " << _name << "..." << std::endl;
  }
  catch (char* str)
  {
    std::cout << POMS_VER << "Unable to add action " << cmd << "!" << std::endl;
    std::cout << "\t" << str << std::endl;
    delete str;
    delete action;
  }
  return action;
}

SubSystemAction* SubSystem::AddAction(const std::string& cmd, const std::string& description)
{
  return AddAction(cmd.c_str(), description.c_str());
}

SubSystemAction* SubSystem::AddAction(SubSystemAction* action)
{
  _actions.push_back(action);

  std::cout << POMS_VER << "Action \"" << action->GetDescription() << "\" added to " << _name << "..." << std::endl;
  return action;
}

void SubSystem::AddDefaultActions()
{
  std::cout << POMS_VER << "Add default actions to subsystem " << _name << "..." << std::endl;

  // Add Draw Actions
  AddAction(new SubSystemActionDraw(this));
  AddAction(new SubSystemActionDrawPS(this));
  //  AddAction(new SubSystemActionDrawHtml(this));

  // Add Show Canvases Action
  //AddAction(new SubSystemActionShowCanvases(this));

  // Add Tile Canvases Action
  // AddAction(new SubSystemActionTileCanvases(this));
}

void SubSystem::TileCanvases()
{
  PomsMainFrame* pmf = PomsMainFrame::Instance();
  pmf->TileCanvases(GetCanvases());
  ShowCanvases();
}

void SubSystem::CascadeCanvases()
{
  PomsMainFrame* pmf = PomsMainFrame::Instance();
  pmf->CascadeCanvases(GetCanvases());
  ShowCanvases();
}

/////////////////////////////////////////////////////////////////////////////
//   SubSystemAction Implementation                                        //
/////////////////////////////////////////////////////////////////////////////

int SubSystemAction::_nextId = SUBSYSTEM_ACTION_ID_BEGIN;
SubSystemActionMap SubSystemAction::_map;

SubSystemAction::SubSystemAction(SubSystem* parent)
  : _running(false)
  , _parent(parent)
{
  if (!parent)
  {
    const char* error = "ERROR: Action must have parent!";
    throw error;
  }

  _id = NextId();
  _map[_id];
}

SubSystemAction::SubSystemAction(SubSystem* parent, const char* description)
  : _running(false)
  , _parent(parent)
  , _description(description)
{
  if (!parent)
  {
    const char* error = "ERROR: Action must have parent!";
    throw error;
  }

  _id = NextId();    // Id used to assign to buttons, menu items, etc
  _map[_id] = this;  // Add action to lookup table, see PomsMainFrame::HandleButtonPoms()
}

SubSystemAction::SubSystemAction(SubSystem* parent, const char* cmd, const char* description)
  : _running(false)
  , _parent(parent)
  , _cmd(cmd)
  , _description(description)
{
  if (!parent || (strlen(cmd) < 1))
  {
    const char* error = "ERROR: Action must have parent and command string!";
    throw error;
  }

  _id = NextId();
  _map[_id] = this;
}

SubSystemAction::~SubSystemAction()
{
  _map[_id] = nullptr;
}

int SubSystemAction::Execute()
{
  if (_running)
    return 0;

  _running = true;
  if (!_parent->isInitialized())  // Check to see if DrawInit() has been executed
  {
    gROOT->ProcessLine((_parent->GetPrefix() + "DrawInit(1)").c_str());
    _parent->setInitialized(1);
  }

  TSeqCollection* allCanvases = gROOT->GetListOfCanvases();
  TCanvas* canvas = nullptr;
  while ((canvas = (TCanvas*) allCanvases->First()))
  {
    std::cout << "Deleting Canvas " << canvas->GetName() << std::endl;
    delete canvas;
  }
  gROOT->ProcessLine(_cmd.c_str());
  _running = false;
  return 0;
}

/////////////////////////////////////////////////////////////////////////////
//    SubSystemActionDraw  Implementation                                  //
/////////////////////////////////////////////////////////////////////////////

SubSystemActionDraw::SubSystemActionDraw(SubSystem* parent)
  : SubSystemAction(parent, "Draw")
{
}

int SubSystemActionDraw::Execute()
{
  if (_running)
    return 0;

  _running = true;

  if (!_parent->isInitialized())  // Check to see if DrawInit() has been executed
  {
    gROOT->ProcessLine((_parent->GetPrefix() + "DrawInit(1)").c_str());
    _parent->setInitialized(1);
  }

  TSeqCollection* allCanvases = gROOT->GetListOfCanvases();
  TCanvas* canvas = nullptr;
  while ((canvas = (TCanvas*) allCanvases->First()))
  {
    std::cout << "Deleting Canvas " << canvas->GetName() << std::endl;
    delete canvas;
  }
  gROOT->ProcessLine((_parent->GetPrefix() + "Draw()").c_str());
  _running = false;
  return 0;
}

/////////////////////////////////////////////////////////////////////////////
//    SubSystemActionSavePlot  Implementation                                //
/////////////////////////////////////////////////////////////////////////////

SubSystemActionSavePlot::SubSystemActionSavePlot(SubSystem* parent)
  : SubSystemAction(parent, "Save Plots")
{
}

int SubSystemActionSavePlot::Execute()
{
  if (_running)
    return 0;

  _running = true;
  gROOT->ProcessLine((_parent->GetPrefix() + "SavePlot()").c_str());
  _running = false;
  return 0;
}

/////////////////////////////////////////////////////////////////////////////
//    SubSystemActionDrawPS  Implementation                                //
/////////////////////////////////////////////////////////////////////////////

SubSystemActionDrawPS::SubSystemActionDrawPS(SubSystem* parent)
  : SubSystemAction(parent, "Save Postscript")
{
}

int SubSystemActionDrawPS::Execute()
{
  if (_running)
    return 0;

  _running = true;
  gROOT->ProcessLine((_parent->GetPrefix() + "PS()").c_str());
  _running = false;
  return 0;
}

/////////////////////////////////////////////////////////////////////////////
//    SubSystemActionDrawHtml  Implementation                              //
/////////////////////////////////////////////////////////////////////////////

SubSystemActionDrawHtml::SubSystemActionDrawHtml(SubSystem* parent)
  : SubSystemAction(parent, "Save to HTML")
{
}

int SubSystemActionDrawHtml::Execute()
{
  if (_running)
    return 0;

  _running = true;
  gROOT->ProcessLine((_parent->GetPrefix() + "Html()").c_str());
  _running = false;
  return 0;
}

/////////////////////////////////////////////////////////////////////////////
//    SubSystemActionShowCanvases  Implementation                          //
/////////////////////////////////////////////////////////////////////////////

SubSystemActionShowCanvases::SubSystemActionShowCanvases(SubSystem* parent)
  : SubSystemAction(parent, "Show Canvases")
{
}

int SubSystemActionShowCanvases::Execute()
{
  if (_running)
    return 0;

  _running = true;
  _parent->ShowCanvases();
  _running = false;
  return 0;
}

/////////////////////////////////////////////////////////////////////////////
//    SubSystemActionTileCanvases  Implementation                          //
/////////////////////////////////////////////////////////////////////////////

SubSystemActionTileCanvases::SubSystemActionTileCanvases(SubSystem* parent)
  : SubSystemAction(parent, "Tile Canvases")
{
}

int SubSystemActionTileCanvases::Execute()
{
  if (_running)
    return 0;

  _running = true;
  _parent->TileCanvases();
  _running = false;
  return 0;
}

// -------------------------------------------------------------
// ColorShutterItem Implementation
// -------------------------------------------------------------

ColorShutterItem::ColorShutterItem(const ULong_t bgColor, const TGWindow* p, TGHotString* s,
                                   Int_t id, UInt_t options)
  : TGShutterItem(p, s, id, options)
{
  fButton->ChangeBackground(bgColor);
}
