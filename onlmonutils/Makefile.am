AUTOMAKE_OPTIONS = foreign

BUILT_SOURCES = \
  testexternals.cc

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(ONLMON_MAIN)/include \
  -I$(ONLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libonlmonutils.la

libonlmonutils_la_LIBADD = \
  -L$(libdir) \
  -L$(ONLINE_MAIN)/lib \
  -lonlmonserver \
  -lonlmondb

pkginclude_HEADERS = \
  runningMean.h \
  pseudoRunningMean.h \
  fullRunningMean.h

libonlmonutils_la_SOURCES = \
  runningMean.cc \
  pseudoRunningMean.cc \
  fullRunningMean.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = \
  testexternals.cc

testexternals_LDADD = \
  libonlmonutils.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
