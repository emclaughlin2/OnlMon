AUTOMAKE_OPTIONS = foreign

BUILT_SOURCES = \
  testexternals.cc

SUBDIRS = calib

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(ONLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libonlll1mon_client.la \
  libonlll1mon_server.la

libonlll1mon_server_la_LIBADD = \
  -L$(libdir) \
  -L$(ONLINE_MAIN)/lib \
  -lonlmonserver \
  -lonlmondb

libonlll1mon_client_la_LIBADD = \
  -L$(libdir) \
  -L$(ONLINE_MAIN)/lib \
  -lonlmonclient \
  -lonlmondb


ll1includedir=$(pkgincludedir)/ll1

ll1include_HEADERS = \
  LL1Mon.h \
  LL1MonDraw.h

libonlll1mon_server_la_SOURCES = \
  LL1Mon.cc

libonlll1mon_client_la_SOURCES = \
  LL1MonDraw.cc

bin_SCRIPTS = \
  ll1MonSetup.csh \
  ll1MonSetup.sh

noinst_PROGRAMS = \
  testexternals_server \
  testexternals_client


testexternals_server_SOURCES = \
  testexternals.cc

testexternals_client_SOURCES = \
  testexternals.cc

testexternals_server_LDADD = \
  libonlll1mon_server.la

testexternals_client_LDADD = \
  libonlll1mon_client.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
