#ifndef TC_INPUT_OUTPUT_H
#define TC_INPUT_OUTPUT_H

struct tc_options;

enum tc_exit_code
{
    tc_exit_code_success = 0,
    tc_exit_code_failure = 1,
    tc_exit_code_inexact = 2,
};

bool tc_io_confirm(struct tc_options* options, const char* question);

#endif
