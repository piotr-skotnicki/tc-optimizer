#ifndef TC_OPTIONS_H
#define	TC_OPTIONS_H

#include <string>
#include <vector>
#include <map>

enum tc_algorithm_enum
{
    tc_algorithm_enum_stencil_tiling,
    tc_algorithm_enum_free_schedule_tiling,
    tc_algorithm_enum_static_correction_tiling,
    tc_algorithm_enum_sfs_tiling,
    tc_algorithm_enum_unknown
};

void tc_options_help();

void tc_options_error(const char* msg, ...);

char* tc_options_source_file(int argc, char* argv[]);

enum tc_algorithm_enum tc_options_algorithm(int argc, char* argv[]);

std::map<std::string, std::vector<int> > tc_options_blocks(int argc, char* argv[]);

#endif // TC_OPTIONS_H
