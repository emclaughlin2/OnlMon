/****************************************************************************\
 *    Poms.h, header for:
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

/* ROOT headers */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"
#include <TGFrame.h>
#include <TGShutter.h>
#pragma GCC diagnostic pop

/* Standard C++ headers */

#include <list>
#include <map>
#include <string>

/* forward declarations to speed up compilation */
class TList;
class TGMenuBar;
class TGButton;
class TGPopupMenu;
class TGHotString;
class TGWindow;

#define POMS_VER "POMS Ver 1.0: "

enum EMessageID
{
  M_FILE_TEST,
  M_FILE_EXIT,
  M_WINDOW_ALIGNRIGHT,
  M_WINDOW_TILEALL,
  B_QUIT
};

// Declare all classes to preven compiler errors
class SubSystem;
class SubSystemAction;

// Declare list types
typedef std::list<SubSystem*> SubSystemList;
typedef std::list<SubSystemAction*> SubSystemActionList;
typedef std::less<int> lessp;
typedef std::map<int, SubSystemAction*, lessp> SubSystemActionMap;

/////////////////////////////////////////////////////////////////////////////
//    Class that defines the main Window                                   //
//                                                                         //
//    (actual main window obtained from gClient->GetRoot() since we        //
//    are running in CINT environment)                                     //
//                                                                         //
//    Implements Singleton Pattern                                         //
/////////////////////////////////////////////////////////////////////////////

class PomsMainFrame : public TGMainFrame
{
 private:
  static PomsMainFrame* _instance;

  //Paths
  std::string _macroPath;

  int looping;
  //Root Screen Properties
  /*** See constructor to change default values ***/
  UInt_t _rootWidth;
  UInt_t _rootHeight;
  UInt_t _rootHorizPad;
  UInt_t _rootVertPad;
  UInt_t _windowPad;

  //GUI OBJECTS -- Main Window
  TGButton* _closeButton;
  TGMenuBar* _menuBar;
  TGPopupMenu* _menuFile;
  TGPopupMenu* _menuWindow;
  TGShutter* _shutter;

  TGButton* startloop;
  TGButton* stoploop;
  //Collections
  SubSystemList _subSystemList;

  //MEMBER FUNCTIONS
  TGShutter* BuildShutter();
  int HandleButtonPoms(Long_t parm1);
  void LoopDiLoop();

  // Singleton, private constructure
  PomsMainFrame(const TGWindow* p, UInt_t w, UInt_t h);

 public:
  static PomsMainFrame* Instance();
  virtual ~PomsMainFrame();
  void StopLoop();

  void SetMacroPath(const char* path);
  virtual void CloseWindow();
  virtual Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t /* parm2 */);
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"
  void Draw();
#pragma GCC diagnostic pop
  void AlignRight();
  void TileCanvases(TList* canvasList);
  void CascadeCanvases(TList* canvasList);
  void TileAllCanvases();

  //Accessors
  UInt_t GetRootWidth() { return _rootWidth; };
  UInt_t GetRootHeight() { return _rootHeight; };
  const std::string& GetMacroPath() { return _macroPath; };

  //CONFIGURATION FUNCTIONS
  SubSystem* RegisterSubSystem(const char* name, const char* prefix,
                               int addDefaultActions = 1, int loadLibrary = 1);
  SubSystem* RegisterSubSystem(SubSystem* subSystem);
};

/////////////////////////////////////////////////////////////////////////////
//    Class to store information on each subsystem that is registered      //
/////////////////////////////////////////////////////////////////////////////

class SubSystem
{
 private:
  std::string _name;
  std::string _prefix;
  TList* _canvasList;
  SubSystemActionList _actions;
  int _initialized;

 public:
  SubSystem(const char* name, const char* prefix, int loadLibrary = 1);
  virtual ~SubSystem();

  // Public Functions
  TList* GetCanvases(int forceReQuery = 0);
  void PrintCanvasList();
  void ShowCanvases();

  SubSystemAction* AddAction(const char* cmd, const char* description);
  SubSystemAction* AddAction(const std::string& cmd, const std::string& description);
  SubSystemAction* AddAction(SubSystemAction* action);
  void AddDefaultActions();

  void TileCanvases();
  void CascadeCanvases();

  // Accessor Methods
  const std::string& GetName() { return _name; };
  const std::string& GetPrefix() { return _prefix; };
  SubSystemActionList* GetActions() { return &_actions; };
  int isInitialized() { return _initialized; }
  void setInitialized(const int i) { _initialized = i; }
};

/////////////////////////////////////////////////////////////////////////////
//    Class to store subsystem actions, commands to be added to shutter    //
/////////////////////////////////////////////////////////////////////////////

class SubSystemAction
{
 protected:
  bool _running;
  int _id;
  static int _nextId;
  static SubSystemActionMap _map;

  // Member Functions
  int NextId() { return _nextId++; };

  SubSystem* _parent;
  std::string _cmd;
  std::string _description;

 public:
  SubSystemAction(SubSystem* parent);
  SubSystemAction(SubSystem* parent, const char* description);
  SubSystemAction(SubSystem* parent, const char* cmd, const char* description);
  virtual ~SubSystemAction();

  virtual int Execute();

  // Accessor Methods
  const std::string& GetCmd() { return _cmd; };
  const std::string& GetDescription() { return _description; };
  int GetId() const { return _id; };
  static SubSystemAction* FindById(int id) { return _map[id]; };
};

class SubSystemActionDraw : public SubSystemAction
{
 public:
  SubSystemActionDraw(SubSystem* parent);
  virtual ~SubSystemActionDraw(){};
  int Execute();
};

class SubSystemActionSavePlot : public SubSystemAction
{
 public:
  SubSystemActionSavePlot(SubSystem* parent);
  virtual ~SubSystemActionSavePlot(){};
  int Execute();
};

class SubSystemActionDrawPS : public SubSystemAction
{
 public:
  SubSystemActionDrawPS(SubSystem* parent);
  virtual ~SubSystemActionDrawPS(){};
  int Execute();
};

class SubSystemActionDrawHtml : public SubSystemAction
{
 public:
  SubSystemActionDrawHtml(SubSystem* parent);
  virtual ~SubSystemActionDrawHtml(){};
  int Execute();
};

class SubSystemActionShowCanvases : public SubSystemAction
{
 public:
  SubSystemActionShowCanvases(SubSystem* parent);
  virtual ~SubSystemActionShowCanvases(){};
  int Execute();
};

class SubSystemActionTileCanvases : public SubSystemAction
{
 public:
  SubSystemActionTileCanvases(SubSystem* parent);
  virtual ~SubSystemActionTileCanvases(){};
  int Execute();
};

class ColorShutterItem : public TGShutterItem
{
 public:
  ColorShutterItem(const ULong_t bgColor, const TGWindow* p, TGHotString* s,
                   Int_t id = -1, UInt_t options = 0);
  virtual ~ColorShutterItem(){};
};
