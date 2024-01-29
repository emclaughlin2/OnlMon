AUTOMAKE_OPTIONS = foreign

BUILT_SOURCES = \
  testexternals.cc

SUBDIRS = calib

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(ONLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libonlmvtxmon_client.la \
  libonlmvtxmon_server.la

libonlmvtxmon_server_la_LIBADD = \
  -L$(libdir) \
  -L$(ONLINE_MAIN)/lib \
  -lonlmonserver \
  -lonlmondb 

libonlmvtxmon_client_la_LIBADD = \
  -L$(libdir) \
  -L$(ONLINE_MAIN)/lib \
  -lonlmonclient \
  -lonlmondb


mvtxincludedir=$(pkgincludedir)/mvtx

mvtxinclude_HEADERS = \
  MvtxMon.h \
  MvtxMonDraw.h

libonlmvtxmon_server_la_SOURCES = \
  MvtxMon.cc

libonlmvtxmon_client_la_SOURCES = \
  MvtxMonDraw.cc

bin_SCRIPTS = \
  MvtxMonSetup.csh \
  MvtxMonSetup.sh

noinst_PROGRAMS = \
  testexternals_server \
  testexternals_client


testexternals_server_SOURCES = \
  testexternals.cc

testexternals_client_SOURCES = \
  testexternals.cc

testexternals_server_LDADD = \
  libonlmvtxmon_server.la

testexternals_client_LDADD = \
  libonlmvtxmon_client.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
