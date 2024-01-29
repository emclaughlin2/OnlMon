AUTOMAKE_OPTIONS = foreign

BUILT_SOURCES = \
  testexternals.cc

SUBDIRS = calib

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(ONLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libonlinttmon_client.la \
  libonlinttmon_server.la

libonlinttmon_server_la_LIBADD = \
  -L$(libdir) \
  -L$(ONLINE_MAIN)/lib \
  -lonlmonserver \
  -lonlmondb

libonlinttmon_client_la_LIBADD = \
  -L$(libdir) \
  -L$(ONLINE_MAIN)/lib \
  -lonlmonclient \
  -lonlmondb


inttincludedir=$(pkgincludedir)/intt

inttinclude_HEADERS = \
  InttFelixMap.h \
  InttMon.h \
  InttMonConstants.h \
  InttMonDraw.h

libonlinttmon_server_la_SOURCES = \
  InttFelixMap.cc \
  InttMon.cc \
  InttMonConstants.cc

libonlinttmon_client_la_SOURCES = \
  InttFelixMap.cc \
  InttMonDraw.cc \
  InttMonConstants.cc

bin_SCRIPTS = \
  InttMonSetup.csh \
  InttMonSetup.sh

noinst_PROGRAMS = \
  testexternals_server \
  testexternals_client


testexternals_server_SOURCES = \
  testexternals.cc

testexternals_client_SOURCES = \
  testexternals.cc

testexternals_server_LDADD = \
  libonlinttmon_server.la

testexternals_client_LDADD = \
  libonlinttmon_client.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
