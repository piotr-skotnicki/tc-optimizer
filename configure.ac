# TC Optimizing Compiler

AC_INIT([tc], [0.4.1], [pskotnicki@zut.edu.pl])

AC_CONFIG_AUX_DIR([build-aux])

AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([foreign])

AC_PROG_CXX
AC_PROG_LIBTOOL

debug="false"
AC_ARG_ENABLE(debug, [  --enable-debug  debug compilation ], debug="true")

AM_CONDITIONAL(TC_DEBUG, test x$debug = xtrue)

if test $debug = "true"; then
    CFLAGS="$CFLAGS -g -O0"
else
    CFLAGS="$CFLAGS -O2"
fi

AC_CONFIG_FILES([Makefile src/Makefile])

AC_CONFIG_SUBDIRS([isl])

AC_CONFIG_SUBDIRS([pet])

AC_CONFIG_SUBDIRS([barvinok])

AC_OUTPUT

echo "Run \`make' to compile tc"
