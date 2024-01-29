AUTOMAKE_OPTIONS = foreign

AM_LDFLAGS = -L$(libdir)

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(ONLINE_MAIN)/include \
  -I$(OPT_SPHENIX)/include \
  -isystem$(ROOTSYS)/include 

BUILT_SOURCES = \
  testexternals.cc

pkginclude_HEADERS = \
  OnCalDBodbc.h \
  OnlMonDB.h \
  OnlMonDBodbc.h \
  OnlMonDBReturnCodes.h \
  RunDBodbc.h

lib_LTLIBRARIES = \
  libonlmonodbc.la \
  libonlmondb.la

libonlmondb_la_SOURCES = \
  OnlMonDB.cc \
  OnlMonDBVar.cc

libonlmonodbc_la_SOURCES = \
  OnCalDBodbc.cc \
  OnlMonDBodbc.cc \
  RunDBodbc.cc

libonlmondb_la_LIBADD = \
  libonlmonodbc.la \
  -lonlmonserver

libonlmonodbc_la_LIBADD = \
  -L$(OPT_SPHENIX)/lib \
  -lonlmonserver \
  -lodbc \
  -lodbc++

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = \
  testexternals.cc

testexternals_LDADD = \
  libonlmondb.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
