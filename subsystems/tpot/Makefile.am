AUTOMAKE_OPTIONS = foreign

BUILT_SOURCES = \
  testexternals.cc

SUBDIRS = calib

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(ONLINE_MAIN)/include \
  -I$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(ONLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib

lib_LTLIBRARIES = \
  libonltpotmon_client.la \
  libonltpotmon_server.la

libonltpotmon_server_la_LIBADD = \
  -lmicromegas_io \
  -lonlmonserver \
  -lonlmondb

libonltpotmon_client_la_LIBADD = \
  -lmicromegas_io \
  -lonlmonclient \
  -lonlmondb


tpotincludedir=$(pkgincludedir)/tpot

tpotinclude_HEADERS = \
  TpotMon.h \
  TpotMonDraw.h \
  MicromegasGeometry.h

libonltpotmon_server_la_SOURCES = \
  TpotMon.cc \
  MicromegasGeometry.cc

libonltpotmon_client_la_SOURCES = \
  TpotMonDraw.cc \
  MicromegasGeometry.cc

bin_SCRIPTS = \
  TpotMonSetup.csh \
  TpotMonSetup.sh

noinst_PROGRAMS = \
  testexternals_server \
  testexternals_client

testexternals_server_SOURCES = \
  testexternals.cc

testexternals_client_SOURCES = \
  testexternals.cc

testexternals_server_LDADD = \
  libonltpotmon_server.la

testexternals_client_LDADD = \
  libonltpotmon_client.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
