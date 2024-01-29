AUTOMAKE_OPTIONS = foreign

BUILT_SOURCES = \
  testexternals.cc

SUBDIRS = calib

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(ONLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libonlzdcmon_client.la \
  libonlzdcmon_server.la

libonlzdcmon_server_la_LIBADD = \
  -L$(libdir) \
  -L$(ONLINE_MAIN)/lib \
  -lonlmonserver \
  -lonlmonutils \
  -lcalo_io \
  -lcalo_reco \
  -lonlmondb

libonlzdcmon_client_la_LIBADD = \
  -L$(libdir) \
  -L$(ONLINE_MAIN)/lib \
  -lonlmonclient \
  -lonlmondb


zdcincludedir=$(pkgincludedir)/zdc

zdcinclude_HEADERS = \
  ZdcMon.h \
  ZdcMonDraw.h

libonlzdcmon_server_la_SOURCES = \
  ZdcMon.cc

libonlzdcmon_client_la_SOURCES = \
  ZdcMonDraw.cc

bin_SCRIPTS = \
  ZdcMonSetup.csh \
  ZdcMonSetup.sh

noinst_PROGRAMS = \
  testexternals_server \
  testexternals_client


testexternals_server_SOURCES = \
  testexternals.cc

testexternals_client_SOURCES = \
  testexternals.cc

testexternals_server_LDADD = \
  libonlzdcmon_server.la

testexternals_client_LDADD = \
  libonlzdcmon_client.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
