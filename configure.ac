# TC Optimizing Compiler

AC_INIT([tc], [0.4.1], [pskotnicki@zut.edu.pl])

AC_CONFIG_AUX_DIR([build-aux])

AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([foreign])

AC_PROG_CXX
AC_PROG_LIBTOOL

debug="no"

AC_ARG_ENABLE(
    [debug],
    [AS_HELP_STRING([--enable-debug], [enable debug mode])],
    [debug="yes"],
    [debug="no"]
)

if test "x$debug" = "xyes"; then
    CFLAGS="$CFLAGS -g -O0 -Wall"
    CXXFLAGS="$CXXFLAGS -g -O0 -Wall"
else
    CFLAGS="$CFLAGS -O2 -DNDEBUG"
    CXXFLAGS="$CXXFLAGS -O2 -DNDEBUG"
fi

AC_CONFIG_FILES([Makefile src/Makefile])

AC_CONFIG_SUBDIRS([isl])

AC_CONFIG_SUBDIRS([pet])

AC_CONFIG_SUBDIRS([barvinok])

AC_OUTPUT

echo "Run \`make' to compile tc"
