AC_INIT(onlmon,[2.00])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE

AC_PROG_CXX(CC g++)
LT_INIT([disable-static])

INSTALL="/usr/bin/install -D -p"

dnl   no point in suppressing warnings people should 
dnl   at least see them, so here we go for g++: -Wall
if test $ac_cv_cxx_compiler_gnu = yes; then
dnl  CXXFLAGS="$CXXFLAGS -Wall"
  CXXFLAGS="$CXXFLAGS -Wno-deprecated-declarations -Wall -Werror -Wextra -Wshadow -Woverloaded-virtual"
fi

ROOTLIBS="`$ROOTSYS/bin/root-config --glibs`"
AC_SUBST(ROOTLIBS)

dnl PLEASE USE ALPHABETICAL ORDER FOR THE MAKEFILE LIST
dnl otherwise noone will be able to find your system
AC_CONFIG_FILES([
Makefile 
onlmonclient/Makefile 
onlmondatabase/Makefile 
poms/Makefile 
onlmonserver/Makefile 
onlmonutils/Makefile 
scripts/Makefile
macros/Makefile
subsystems/Makefile
subsystems/ll1/Makefile
subsystems/ll1/calib/Makefile
subsystems/cemc/Makefile
subsystems/cemc/calib/Makefile
subsystems/daq/Makefile
subsystems/daq/calib/Makefile
subsystems/sepd/Makefile
subsystems/sepd/calib/Makefile
subsystems/example/Makefile
subsystems/example/calib/Makefile
subsystems/hcal/Makefile
subsystems/hcal/calib/Makefile
subsystems/intt/Makefile
subsystems/intt/calib/Makefile
subsystems/mvtx/Makefile
subsystems/mvtx/calib/Makefile
subsystems/pktsize/Makefile
subsystems/tpc/Makefile
subsystems/tpc/calib/Makefile
subsystems/tpot/Makefile
subsystems/tpot/calib/Makefile
subsystems/zdc/Makefile
subsystems/zdc/calib/Makefile
])
AC_OUTPUT
