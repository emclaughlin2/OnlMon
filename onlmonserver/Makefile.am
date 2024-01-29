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

lib_LTLIBRARIES = \
  libonlmonserver.la \
  libonlmonserver_funcs.la

libonlmonserver_funcs_la_LIBADD = \
  libonlmonserver.la \
  -L$(libdir) \
  -L$(ONLINE_MAIN)/lib \
  -lpmonitor \
  -lmessage

libonlmonserver_la_LDFLAGS = \
  -L$(libdir) \
  -L$(ONLINE_MAIN)/lib \
  -L$(OPT_SPHENIX)/lib \
  `root-config --libs` \
  -lmessage \
  -lodbc++ \
  -lz


noinst_HEADERS = \
  pmonitorInterface.h

pkginclude_HEADERS = \
  HistoBinDefs.h \
  OnlMon.h \
  OnlMonBase.h \
  OnlMonDefs.h \
  OnlMonServer.h \
  OnlMonStatus.h

libonlmonserver_funcs_la_SOURCES = \
  pmonitorInterface.cc

libonlmonserver_la_SOURCES = \
  MessageSystem.cc \
  OnlMon.cc \
  OnlMonBase.cc \
  OnlMonServer.cc \
  OnlMonStatusDB.cc

BUILT_SOURCES = \
  testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = \
  testexternals.cc

testexternals_LDADD = \
  libonlmonserver.la \
  libonlmonserver_funcs.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@


clean-local:
	rm -f $(BUILT_SOURCES)
