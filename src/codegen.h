#ifndef TC_CODEGEN_H
#define	TC_CODEGEN_H

#include "scop.h"
#include "options.h"

#include <isl/ctx.h>
#include <isl/id.h>
#include <isl/ast.h>
#include <isl/ast_build.h>
#include <isl/printer.h>

struct tc_codegen_context
{
    isl_bool in_parallel_region;
    
    struct tc_options* options;
};

struct tc_codegen_context* tc_codegen_context_alloc(__isl_keep isl_ctx* ctx);

void tc_codegen_context_free(struct tc_codegen_context* context);

__isl_give isl_printer* tc_codegen_print_statements_macros(struct tc_scop* scop, __isl_take isl_printer* printer);

__isl_give isl_printer* tc_codegen_print_user(__isl_take isl_printer* printer, __isl_take isl_ast_print_options* options, __isl_keep isl_ast_node* node, void* user);

__isl_give isl_printer* tc_codegen_print_prologue(struct tc_scop* scop, struct tc_options* options, __isl_take isl_printer* printer);

__isl_give isl_printer* tc_codegen_print_epilogue(struct tc_scop* scop, struct tc_options* options, __isl_take isl_printer* printer);

#endif // TC_CODEGEN_H
