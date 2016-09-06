#ifndef TC_OPTIONS_H
#define	TC_OPTIONS_H

#include <isl/ctx.h>
#include <isl/set.h>

#include <stdio.h>

#include <string>
#include <vector>
#include <map>

enum tc_algorithm_enum
{
    tc_algorithm_enum_regular_tiling,
    tc_algorithm_enum_correction_tiling,
    tc_algorithm_enum_merge_tiling,
    tc_algorithm_enum_stencil_tiling,
    tc_algorithm_enum_unknown
};

enum tc_scheduling_enum
{
    tc_scheduling_enum_lex,
    tc_scheduling_enum_sfs_tile,
    tc_scheduling_enum_free_rk,
    tc_scheduling_enum_free_karl,
    tc_scheduling_enum_free_finite,
    tc_scheduling_enum_free_dynamic,
    tc_scheduling_enum_sfs_single,
    tc_scheduling_enum_sfs_multiple,
    tc_scheduling_enum_unknown
};

enum tc_codegen_enum
{
    tc_codegen_enum_serial,
    tc_codegen_enum_omp_cpu_for,
    tc_codegen_enum_omp_cpu_task,
    tc_codegen_enum_unknown
};

enum tc_transitive_closure_enum
{
    tc_transitive_closure_enum_isl_map,
    tc_transitive_closure_enum_isl_union_map,
    tc_transitive_closure_enum_floyd_warshall,
    tc_transitive_closure_enum_tarjan,
    tc_transitive_closure_enum_unknown
};

struct tc_options
{
    int argc;
    
    char** argv;
    
    FILE* output;
};

struct tc_options* tc_options_alloc(int argc, char* argv[]);

void tc_options_free(struct tc_options* options);

void tc_options_help();

void tc_options_credits();

void tc_options_error(const char* msg, ...);

char* tc_options_get_command_line(struct tc_options* options);

int tc_options_get_int(struct tc_options* options, const char* short_name, const char* long_name);

long tc_options_get_long(struct tc_options* options, const char* short_name, const char* long_name);

const char* tc_options_get_string(struct tc_options* options, const char* short_name, const char* long_name);

int tc_options_is_set(struct tc_options* options, const char* short_name, const char* long_name);

const char* tc_options_source_file(struct tc_options* options);

const char* tc_options_output_file(struct tc_options* options);

int tc_options_is_verbose(struct tc_options* options);

int tc_options_is_report(struct tc_options* options);

int tc_options_cache_line(struct tc_options* options);

enum tc_algorithm_enum tc_options_algorithm(struct tc_options* options);

enum tc_scheduling_enum tc_options_scheduling(struct tc_options* options);

enum tc_codegen_enum tc_options_codegen(struct tc_options* options);

enum tc_transitive_closure_enum tc_options_transitive_closure(struct tc_options* options);

std::map<std::string, std::vector<int> > tc_options_blocks(struct tc_options* options);

std::vector<std::vector<std::string> > tc_options_groups(struct tc_options* options);

__isl_give isl_set* tc_options_collect_values(struct tc_options* options, const char* short_name, const char* long_name, __isl_keep isl_ctx* ctx);

__isl_give isl_set* tc_options_get_defines(struct tc_options* options, __isl_keep isl_ctx* ctx);

__isl_give isl_set* tc_options_get_report_bounds(struct tc_options* options, __isl_keep isl_ctx* ctx);

void tc_options_check_spelling(struct tc_options* options);

#endif // TC_OPTIONS_H
