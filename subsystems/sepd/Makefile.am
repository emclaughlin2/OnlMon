AUTOMAKE_OPTIONS = foreign

BUILT_SOURCES = \
  testexternals.cc

SUBDIRS = calib

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(ONLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libonlsepdmon_client.la \
  libonlsepdmon_server.la

libonlsepdmon_server_la_LIBADD = \
  -L$(libdir) \
  -L$(ONLINE_MAIN)/lib \
  -lonlmonserver \
  -lonlmondb \
  -lonlmonutils \
  -lcalo_io \
  -lcalo_reco

libonlsepdmon_client_la_LIBADD = \
  -L$(libdir) \
  -L$(ONLINE_MAIN)/lib \
  -lonlmonclient \
  -lonlmondb


sepdincludedir=$(pkgincludedir)/sepd

sepdinclude_HEADERS = \
  SepdMon.h \
  SepdMonDraw.h

libonlsepdmon_server_la_SOURCES = \
  SepdMon.cc

libonlsepdmon_client_la_SOURCES = \
  SepdMonDraw.cc

bin_SCRIPTS = \
  SepdMonSetup.csh \
  SepdMonSetup.sh

noinst_PROGRAMS = \
  testexternals_server \
  testexternals_client


testexternals_server_SOURCES = \
  testexternals.cc

testexternals_client_SOURCES = \
  testexternals.cc

testexternals_server_LDADD = \
  libonlsepdmon_server.la

testexternals_client_LDADD = \
  libonlsepdmon_client.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
