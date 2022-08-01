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
  libonlmbdmon_client.la \
  libonlmbdmon_server.la

libonlmbdmon_server_la_LIBADD = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lonlmonserver \
  -lonlmondb

libonlmbdmon_client_la_LIBADD = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lonlmonclient \
  -lSpectrum \
  -lonlmondb


mbdincludedir=$(pkgincludedir)/mbd

mbdinclude_HEADERS = \
  MbdMonDefs.h \
  MbdMon.h \
  MbdMonDraw.h

libonlmbdmon_server_la_SOURCES = \
  MbdMon.cc

libonlmbdmon_client_la_SOURCES = \
  MbdMonDraw.cc

bin_SCRIPTS = \
  MbdMonSetup.csh \
  MbdMonSetup.sh

noinst_PROGRAMS = \
  testexternals_server \
  testexternals_client


testexternals_server_SOURCES = \
  testexternals.cc

testexternals_client_SOURCES = \
  testexternals.cc

testexternals_server_LDADD = \
  libonlmbdmon_server.la

testexternals_client_LDADD = \
  libonlmbdmon_client.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
