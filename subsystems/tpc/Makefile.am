AUTOMAKE_OPTIONS = foreign

BUILT_SOURCES = \
  testexternals.cc

SUBDIRS = calib

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(ONLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libonltpcmon_client.la \
  libonltpcmon_server.la

libonltpcmon_server_la_LIBADD = \
  -L$(libdir) \
  -L$(ONLINE_MAIN)/lib \
  -lonlmonserver \
  -lonlmondb \
  -ltpc 

libonltpcmon_client_la_LIBADD = \
  -L$(libdir) \
  -L$(ONLINE_MAIN)/lib \
  -lonlmonclient \
  -lonlmondb


tpcincludedir=$(pkgincludedir)/tpc

tpcinclude_HEADERS = \
  TpcMon.h \
  TpcMonDraw.h

libonltpcmon_server_la_SOURCES = \
  TpcMon.cc

libonltpcmon_client_la_SOURCES = \
  TpcMonDraw.cc

bin_SCRIPTS = \
  TpcMonSetup.csh \
  TpcMonSetup.sh

noinst_PROGRAMS = \
  testexternals_server \
  testexternals_client


testexternals_server_SOURCES = \
  testexternals.cc

testexternals_client_SOURCES = \
  testexternals.cc

testexternals_server_LDADD = \
  libonltpcmon_server.la

testexternals_client_LDADD = \
  libonltpcmon_client.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
