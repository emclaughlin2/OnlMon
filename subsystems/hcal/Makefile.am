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
  libonlhcalmon_client.la \
  libonlhcalmon_server.la

libonlhcalmon_server_la_LIBADD = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lonlmonserver \
  -lonlmondb

libonlhcalmon_client_la_LIBADD = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lonlmonclient \
  -lonlmondb


hcalincludedir=$(pkgincludedir)/hcal

hcalinclude_HEADERS = \
  HcalMon.h \
  HcalMonDraw.h \
  runningMean.h

libonlhcalmon_server_la_SOURCES = \
  HcalMon.cc \
  runningMean.cc \
  pseudoRunningMean.cc \
  fullRunningMean.cc

libonlhcalmon_client_la_SOURCES = \
  HcalMonDraw.cc

bin_SCRIPTS = \
  HcalMonSetup.csh \
  HcalMonSetup.sh

noinst_PROGRAMS = \
  testexternals_server \
  testexternals_client


testexternals_server_SOURCES = \
  testexternals.cc

testexternals_client_SOURCES = \
  testexternals.cc

testexternals_server_LDADD = \
  libonlhcalmon_server.la

testexternals_client_LDADD = \
  libonlhcalmon_client.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
