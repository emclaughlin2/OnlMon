AUTOMAKE_OPTIONS = foreign

BUILT_SOURCES = \
  testexternals.cc

SUBDIRS = calib

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(ONLMON_MAIN)/include \
  -I$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libonlcemcmon_client.la \
  libonlcemcmon_server.la

libonlcemcmon_server_la_LIBADD = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lonlmonserver \
  -lonlmondb

libonlcemcmon_client_la_LIBADD = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lonlmonclient \
  -lonlmondb


cemcincludedir=$(pkgincludedir)/cemc

cemcinclude_HEADERS = \
  CemcMon.h \
  CemcMonDraw.h \
  runningMean.h

libonlcemcmon_server_la_SOURCES = \
  CemcMon.cc \
  runningMean.cc \
  pseudoRunningMean.cc \
  fullRunningMean.cc

libonlcemcmon_client_la_SOURCES = \
  CemcMonDraw.cc

bin_SCRIPTS = \
  CemcMonSetup.csh \
  CemcMonSetup.sh

noinst_PROGRAMS = \
  testexternals_server \
  testexternals_client


testexternals_server_SOURCES = \
  testexternals.cc

testexternals_client_SOURCES = \
  testexternals.cc

testexternals_server_LDADD = \
  libonlcemcmon_server.la

testexternals_client_LDADD = \
  libonlcemcmon_client.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
