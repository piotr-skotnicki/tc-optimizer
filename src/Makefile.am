# TC Optimizing Compiler

bin_PROGRAMS = tc

tc_SOURCES = \
	main.cpp \
	utility.cpp \
	scop.cpp \
	options.cpp \
	debug.cpp \
	tiling.cpp \
	codegen.cpp \
	ast.cpp \
	slicing.cpp \
	regular_tiling.cpp \
	stencil_tiling.cpp \
	merge_tiling.cpp \
	correction_tiling.cpp \
	lex_scheduling.cpp \
	sfs_scheduling.cpp \
	free_scheduling.cpp \
	dynamic_free_scheduling.cpp \
	serial_codegen.cpp \
	omp_cpu_codegen.cpp \    
	transitive_closure.cpp \
	tarjan_transitive_closure.cpp \
	floyd_warshall_transitive_closure.cpp \
	iterative_transitive_closure.cpp \
	tile_statistics.cpp \
	for_decorator.cpp \
	timer.cpp

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
