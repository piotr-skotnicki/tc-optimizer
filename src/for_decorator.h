#ifndef TC_FOR_DECORATOR_H
#define TC_FOR_DECORATOR_H

#include <isl/ast.h>
#include <isl/ast_build.h>
#include <isl/printer.h>

__isl_give isl_printer* tc_for_decorator_omp_task_nested(__isl_take isl_printer* printer, __isl_take isl_ast_print_options* options, __isl_keep isl_ast_node* node, void* user);

__isl_give isl_printer* tc_for_decorator_omp_task_first(__isl_take isl_printer* printer, __isl_take isl_ast_print_options* options, __isl_keep isl_ast_node* node, void* user);

__isl_give isl_printer* tc_for_decorator_omp_parallel_for_nested(__isl_take isl_printer* printer, __isl_take isl_ast_print_options* options, __isl_keep isl_ast_node* node, void* user);

__isl_give isl_printer* tc_for_decorator_omp_parallel_for_first(__isl_take isl_printer* printer, __isl_take isl_ast_print_options* options, __isl_keep isl_ast_node* node, void* user);

__isl_give isl_printer* tc_for_decorator_omp_for_first(__isl_take isl_printer* printer, __isl_take isl_ast_print_options* options, __isl_keep isl_ast_node* node, void* user);

__isl_give isl_printer* tc_for_decorator_omp_teams_distribute_parallel_for(__isl_take isl_printer* printer, __isl_take isl_ast_print_options* options, __isl_keep isl_ast_node* node, void* user);

#endif // TC_FOR_DECORATOR_H
