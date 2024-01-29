AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(ONLINE_MAIN)/include \
  -I$(OPT_SPHENIX)/include \
  -isystem$(ROOTSYS)/include 

RINCLUDES = \
  -I$(includedir) \
  -I$(ONLINE_MAIN)/include \
  -I$(OPT_SPHENIX)/include \
  -I$(ROOTSYS)/include 

BUILT_SOURCES = \
  testexternals.cc

lib_LTLIBRARIES = libonlmonclient.la   

noinst_HEADERS = \
  ClientHistoList.h \
  OnlMonHtml.h

pkginclude_HEADERS = \
  OnlMonClient.h \
  OnlMonDraw.h

libonlmonclient_la_SOURCES = \
  OnlMonClient.cc \
  OnlMonDraw.cc \
  OnlMonHtml.cc \
  ClientHistoList.cc

libonlmonclient_la_LDFLAGS = \
  -L$(libdir) \
  -lstdc++fs \
  -lonlmondb \
  -luuid \
  `root-config --glibs`

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = \
  testexternals.cc

testexternals_LDADD = \
  libonlmonclient.la


testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
