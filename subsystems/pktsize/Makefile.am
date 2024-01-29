AUTOMAKE_OPTIONS = foreign

BUILT_SOURCES = \
  testexternals.cc

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(ONLINE_MAIN)/include \
  -I$(OPT_SPHENIX)/include \
  -isystem$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libonlpktsizemon_server.la \
  libonlpktsizemon_client.la

libonlpktsizemon_server_la_LIBADD = \
  -L$(libdir) \
  -L$(ONLINE_MAIN)/lib \
  -lonlmonserver \
  -lonlmondb \
  -lodbc++

libonlpktsizemon_client_la_LIBADD = \
  -L$(libdir) \
  -L$(ONLINE_MAIN)/lib \
  libonlpktsizemon_server.la \
  -lonlmonclient \
  -lonlmondb

pktsizeincludedir=$(pkgincludedir)/pktsize

pktsizeinclude_HEADERS = \
  PktSizeMon.h \
  PktSizeMonDraw.h

libonlpktsizemon_server_la_SOURCES = \
  PktSizeCommon.cc \
  PktSizeDBodbc.cc \
  PktSizeMon.cc

libonlpktsizemon_client_la_SOURCES = \
  PktSizeMonDraw.cc


noinst_PROGRAMS = \
  testexternals_server \
  testexternals_client

testexternals_server_SOURCES = \
  testexternals.cc

testexternals_client_SOURCES = \
  testexternals.cc

testexternals_server_LDADD = \
  libonlpktsizemon_server.la

testexternals_client_LDADD = \
  libonlpktsizemon_client.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@


clean-local:
	rm -f $(BUILT_SOURCES)
