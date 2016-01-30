# TC Optimizing Compiler

bin_PROGRAMS = tc

tc_SOURCES = \
	main.cpp \
	utility.cpp \
	scop.cpp \
	tiling.cpp \
	tarjan.cpp \
	slicing.cpp \
	stencil.cpp \
	free_schedule_tiling.cpp \
	static_tile_correction.cpp \
	transitive_closure.cpp \
	options.cpp

tc_CPPFLAGS = \
	-I$(top_srcdir)/isl/include \
	-I$(top_srcdir)/pet/include \
	-I$(top_srcdir)/barvinok

tc_LDFLAGS = -static

tc_LDADD = \
	$(top_srcdir)/isl/libisl.la \
	$(top_srcdir)/pet/libpet.la \
	$(top_srcdir)/barvinok/libbarvinok.la \
	-lm
