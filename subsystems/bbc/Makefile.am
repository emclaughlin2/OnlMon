AUTOMAKE_OPTIONS = foreign

BUILT_SOURCES = \
  testexternals.cc

SUBDIRS = calib

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(ONLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libonlbbcmon_client.la \
  libonlbbcmon_server.la

libonlbbcmon_server_la_LIBADD = \
  -L$(libdir) \
  -L$(ONLINE_MAIN)/lib \
  -lbbc_io \
  -lonlmonserver

libonlbbcmon_client_la_LIBADD = \
  -L$(libdir) \
  -L$(ONLINE_MAIN)/lib \
  -L$(ROOTSYS)/lib \
  -lonlmonclient \
  -lSpectrum

bbcincludedir=$(pkgincludedir)/bbc

bbcinclude_HEADERS = \
  BbcMonDefs.h \
  BbcMon.h \
  OnlBbcSig.h \
  OnlBbcEvent.h \
  BbcMonDraw.h

libonlbbcmon_server_la_SOURCES = \
  OnlBbcSig.cc \
  OnlBbcEvent.cc \
  BbcMon.cc

libonlbbcmon_client_la_SOURCES = \
  BbcMonDraw.cc

bin_SCRIPTS = \
  BbcMonSetup.csh \
  BbcMonSetup.sh

noinst_PROGRAMS = \
  testexternals_server \
  testexternals_client


testexternals_server_SOURCES = \
  testexternals.cc

testexternals_client_SOURCES = \
  testexternals.cc

testexternals_server_LDADD = \
  libonlbbcmon_server.la

testexternals_client_LDADD = \
  libonlbbcmon_client.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
