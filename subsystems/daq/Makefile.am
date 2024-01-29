AUTOMAKE_OPTIONS = foreign

BUILT_SOURCES = \
  testexternals.cc

SUBDIRS = calib

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(ONLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libonldaqmon_client.la \
  libonldaqmon_server.la

libonldaqmon_server_la_LIBADD = \
  -L$(libdir) \
  -L$(ONLINE_MAIN)/lib \
  -lonlmonserver \
  -lonlmondb

libonldaqmon_client_la_LIBADD = \
  -L$(libdir) \
  -L$(ONLINE_MAIN)/lib \
  -lonlmonclient \
  -lonlmondb


daqincludedir=$(pkgincludedir)/daq

daqinclude_HEADERS = \
  DaqMon.h \
  DaqMonDraw.h

libonldaqmon_server_la_SOURCES = \
  DaqMon.cc

libonldaqmon_client_la_SOURCES = \
  DaqMonDraw.cc

bin_SCRIPTS = \
  DaqMonSetup.csh \
  DaqMonSetup.sh

noinst_PROGRAMS = \
  testexternals_server \
  testexternals_client


testexternals_server_SOURCES = \
  testexternals.cc

testexternals_client_SOURCES = \
  testexternals.cc

testexternals_server_LDADD = \
  libonldaqmon_server.la

testexternals_client_LDADD = \
  libonldaqmon_client.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
