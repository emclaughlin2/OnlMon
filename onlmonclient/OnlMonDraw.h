#ifndef ONLMONCLIENT_ONLMONDRAW_H
#define ONLMONCLIENT_ONLMONDRAW_H

#include <string>

class TPad;

class OnlMonDraw
{
 public:
  OnlMonDraw(const std::string &name = "NONE");
  virtual ~OnlMonDraw() {}

  virtual int Init() { return 0; }
  virtual int Draw(const std::string &what = "ALL");
  virtual int MakePS(const std::string &what = "ALL");
  virtual int MakeHtml(const std::string &what = "ALL");
  const std::string Name() const { return ThisName; }
  void Verbosity(const int i) { verbosity = i; }
  int Verbosity() const { return verbosity; }

 protected:
  virtual int DrawDeadServer(TPad *transparent);
  std::string ThisName;
  int verbosity = 0;
};

#endif /* ONLMONCLIENT_ONLMONDRAW_H */
