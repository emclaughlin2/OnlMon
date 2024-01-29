AUTOMAKE_OPTIONS = foreign

BUILT_SOURCES = \
  testexternals.cc

SUBDIRS = calib

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(ONLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libonlcemcmon_client.la \
  libonlcemcmon_server.la

libonlcemcmon_server_la_LIBADD = \
  -L$(libdir) \
  -L$(ONLINE_MAIN)/lib \
  -lonlmonserver \
  -lonlmonutils \
  -lcalo_io \
  -lcalo_reco

libonlcemcmon_client_la_LIBADD = \
  -L$(libdir) \
  -L$(ONLINE_MAIN)/lib \
  -lonlmonclient

cemcincludedir=$(pkgincludedir)/cemc

cemcinclude_HEADERS = \
  CemcMon.h \
  CemcMonDraw.h

libonlcemcmon_server_la_SOURCES = \
  CemcMon.cc

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
