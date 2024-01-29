#ifndef ONLMONHTML_H
#define ONLMONHTML_H

#include <set>
#include <string>

class RunDBodbc;

/**@name Helper class to centralize HTML generation code.
 */

class OnlMonHtml
{
 public:
  explicit OnlMonHtml(const std::string& topdir);
  //! delete copy ctor and assignment opertor (cppcheck)
  explicit OnlMonHtml(const OnlMonHtml&) = delete;
  OnlMonHtml& operator=(const OnlMonHtml&) = delete;

  virtual ~OnlMonHtml();

  /** Generate a bit of the navigation menu for a given file (link).
   *  @param path the path as it will appear in the menu
   *  @param relfilename the filename that will be served when using
   *  path in the menu (must not be a fullpathname, but a plain filename).
   */
  void addMenu(const std::string& header, const std::string& path,
               const std::string& relfilename);

  /** Generate filenames, to be used to produce e.g. gif or html files.
   *  @param drawer the OnlMonDraw child class for which filename must be built
   *  @param basefilename the beginning of the filename 
   *  (will be completed e.g. with the run number by this method)
   *  @param ext the extension of the file (e.g. gif or html)
   *  @return fullfilename = path+filename of the generated filename.
   *  @return filename = only filename of the generated filename
   */
  void namer(const std::string& header, const std::string& basefilename,
             const std::string& ext,
             std::string& fullfilename,
             std::string& filename);

  /** Generate a full filename from specified pieces, *and* register
    * the file to the navigation menu
    * @param drawer the OnlMonDraw child class for which filename must be built
    * @param path the path that will be used in the *navigation menu*
    * @param basefilename the beginning of the filename 
    * (will be completed e.g. with the run number by this method)
    * @param ext the extension of the file (e.g. gif or html)
    * @return the full filename (i.e. full path + filename)
    */
  std::string registerPage(const std::string& header,
                           const std::string& path,
                           const std::string& basefilename,
                           const std::string& ext);

  void runNumber(const int runnumber);
  int runNumber() const { return fRunNumber; }

  void verbosity(const int v) { fVerbosity = v; }
  int verbosity() const { return fVerbosity; }

 protected:
  void plainHtmlMenu(const std::set<std::string>&);
  void runInit();
  std::string runRange();

  RunDBodbc* rundb = nullptr;

  int fVerbosity = 0;
  int fRunNumber = 0;

  std::string fHtmlDir;
  std::string fHtmlRunDir;
};

#endif
