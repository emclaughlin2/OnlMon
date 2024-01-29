AUTOMAKE_OPTIONS = foreign

BUILT_SOURCES = \
  testexternals.cc

SUBDIRS = calib

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(ONLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libonlmymon_client.la \
  libonlmymon_server.la

libonlmymon_server_la_LIBADD = \
  -L$(libdir) \
  -L$(ONLINE_MAIN)/lib \
  -lonlmonserver \
  -lonlmondb

libonlmymon_client_la_LIBADD = \
  -L$(libdir) \
  -L$(ONLINE_MAIN)/lib \
  -lonlmonclient \
  -lonlmondb


exampleincludedir=$(pkgincludedir)/example

exampleinclude_HEADERS = \
  MyMon.h \
  MyMonDraw.h

libonlmymon_server_la_SOURCES = \
  MyMon.cc

libonlmymon_client_la_SOURCES = \
  MyMonDraw.cc

bin_SCRIPTS = \
  MyMonSetup.csh \
  MyMonSetup.sh


noinst_PROGRAMS = \
  testexternals_server \
  testexternals_client


testexternals_server_SOURCES = \
  testexternals.cc

testexternals_client_SOURCES = \
  testexternals.cc

testexternals_server_LDADD = \
  libonlmymon_server.la

testexternals_client_LDADD = \
  libonlmymon_client.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
