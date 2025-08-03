#ifndef TC_INPUT_OUTPUT_H
#define TC_INPUT_OUTPUT_H

#include <stdio.h>

struct tc_options;

enum tc_exit_code
{
    tc_exit_code_success   = 0,
    tc_exit_code_failure   = 1,
    tc_exit_code_assertion = 2,
    tc_exit_code_internal  = 3,
    tc_exit_code_usage     = 4,
    tc_exit_code_inexact   = 5,
    tc_exit_code_input     = 6,
    tc_exit_code_nonlex    = 7,
};

bool tc_io_confirm(struct tc_options* options, const char* question);

FILE* tc_io_open_output(const char* path);

void tc_io_close(FILE* file);

#endif // TC_INPUT_OUTPUT_H
