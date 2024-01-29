AC_INIT(bbcmon,[2.00])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE

AC_PROG_CXX(CC g++)
LT_INIT([disable-static])

dnl   no point in suppressing warnings people should 
dnl   at least see them, so here we go for gcc: -Wall
dnl if test $ac_cv_prog_gcc = yes; then
  CXXFLAGS="$CXXFLAGS -std=c++17 -Wno-deprecated-declarations -Wall -Werror -Wextra -Wshadow -Woverloaded-virtual"
dnl fi

AC_OUTPUT(Makefile calib/Makefile)
