#include "omp_gpu_codegen.h"
#include "codegen.h"
#include "ast.h"
#include "for_decorator.h"
#include "scop.h"
#include "debug.h"
#include "config.h"

#include <isl/ctx.h>
#include <isl/id.h>
#include <isl/set.h>
#include <isl/union_map.h>
#include <isl/union_set.h>
#include <isl/ast_build.h>
#include <isl/printer.h>

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

void tc_codegen_omp_gpu(struct tc_scop* scop, struct tc_options* options, __isl_take isl_union_map* S, __isl_keep isl_id_list* iterators, __isl_keep isl_id_list* parallel_iterators)
{
    isl_ctx* ctx = isl_union_map_get_ctx(S);
    
    isl_ast_build* ast_build = isl_ast_build_from_context(isl_set_copy(scop->pet->context));

    ast_build = isl_ast_build_set_iterators(ast_build, isl_id_list_copy(iterators));
    
    struct tc_ast_visitor_context* visitor_context = tc_ast_visitor_context_alloc(ctx);
    visitor_context->scop = scop;
    visitor_context->parallel_iterators = parallel_iterators;

    std::vector<struct tc_ast_for_annotation*> annotations_stack;
    visitor_context->annotations_stack = &annotations_stack;
    
    ast_build = isl_ast_build_set_before_each_for(ast_build, &tc_ast_visitor_before_for, visitor_context);
    ast_build = isl_ast_build_set_after_each_for(ast_build, &tc_ast_visitor_after_for, visitor_context);
    ast_build = isl_ast_build_set_at_each_domain(ast_build, &tc_ast_visitor_at_each_domain, visitor_context);
    
    isl_printer* printer = isl_printer_to_str(ctx);
    
    printer = isl_printer_set_output_format(printer, ISL_FORMAT_C);

    isl_ast_print_options* ast_options = isl_ast_print_options_alloc(ctx);

    struct tc_codegen_context* codegen_context = tc_codegen_context_alloc(ctx);
    codegen_context->options = options;
    codegen_context->scop = scop;
    
    ast_options = isl_ast_print_options_set_print_for(ast_options, &tc_for_decorator_omp_teams_distribute_parallel_for, codegen_context);
    ast_options = isl_ast_print_options_set_print_user(ast_options, &tc_codegen_print_user, NULL);
    
    isl_ast_node* ast_tile = isl_ast_build_ast_from_schedule(ast_build, S);
    
    printer = tc_codegen_print_prologue(scop, options, printer);
    
    printer = isl_ast_node_print_macros(ast_tile, printer);

    printer = isl_printer_print_str(printer, "#pragma scop\n");

    if (scop->pet->n_array > 0)
    {
        printer = isl_printer_print_str(printer, "#pragma omp target data");
    }

    for (int i = 0; i < scop->pet->n_array; ++i)
    {
        struct pet_array* array = scop->pet->arrays[i];

        const char* name = isl_set_get_tuple_name(array->extent);
    
        isl_bool is_read = isl_bool_false;
        isl_bool is_write = isl_bool_false;

        isl_map_list* reads = isl_union_map_get_map_list(scop->reads);
        for (int j = 0; j < isl_map_list_n_map(reads) && is_read == isl_bool_false; ++j)
        {
            isl_map* read = isl_map_list_get_map(reads, j);

            if (0 == strcmp(isl_map_get_tuple_name(read, isl_dim_out), name))
            {
                is_read = isl_bool_true;
            }

            isl_map_free(read);
        }
        isl_map_list_free(reads);

        isl_map_list* writes = isl_union_map_get_map_list(scop->writes);
        for (int j = 0; j < isl_map_list_n_map(writes) && is_write == isl_bool_false; ++j)
        {
            isl_map* write = isl_map_list_get_map(writes, j);

            if (0 == strcmp(isl_map_get_tuple_name(write, isl_dim_out), name))
            {
                is_write = isl_bool_true;
            }

            isl_map_free(write);
        }
        isl_map_list_free(writes);

        printer = isl_printer_print_str(printer, " map(");

        if (is_read == isl_bool_true && is_write == isl_bool_true)
        {
            printer = isl_printer_print_str(printer, "alloc");
        }
        else if (is_read == isl_bool_true && is_write == isl_bool_false)
        {
            printer = isl_printer_print_str(printer, "to");
        }
        else if (is_read == isl_bool_false && is_write == isl_bool_true)
        {
            printer = isl_printer_print_str(printer, "alloc");
        }
        else if (is_read == isl_bool_false && is_write == isl_bool_false)
        {
            printer = isl_printer_print_str(printer, "to");
        }

        printer = isl_printer_print_str(printer, ":");
        printer = isl_printer_print_str(printer, name);

        const int n_dim_extent = isl_set_n_dim(array->extent);

        for (int j = 0; j < n_dim_extent; ++j)
        {
            printer = isl_printer_print_str(printer, "[0:");

            isl_pw_aff* ub = isl_set_dim_max(isl_set_copy(array->extent), j);
            ub = isl_pw_aff_add(ub, isl_pw_aff_read_from_str(ctx, "{[(1)]}"));

            char* ub_str = isl_pw_aff_to_str(ub);

            char* begin = strstr(ub_str, "[(") + 2;
            char* end = strstr(ub_str, ")]");

            char buffer[512] = {0};
            strncpy(buffer, begin, end - begin);
            printer = isl_printer_print_str(printer, buffer);

            free(ub_str);
            isl_pw_aff_free(ub);

            printer = isl_printer_print_str(printer, "]");
        }
        printer = isl_printer_print_str(printer, ")");
    }

    if (scop->pet->n_array > 0)
    {
        printer = isl_printer_print_str(printer, "\n");
        printer = isl_printer_print_str(printer, "{\n");
        printer = isl_printer_indent(printer, TC_CONF_INDENT_SIZE);
    }
    
    printer = isl_ast_node_print(ast_tile, printer, ast_options);

    if (scop->pet->n_array > 0)
    {
        printer = isl_printer_indent(printer, -TC_CONF_INDENT_SIZE);
        printer = isl_printer_print_str(printer, "}\n");
    }

    printer = isl_printer_print_str(printer, "#pragma endscop\n");
    printer = tc_codegen_print_epilogue(scop, options, printer);
    
    char* code = isl_printer_get_str(printer);
    
    fprintf(options->output, "%s", code);

    fflush(options->output);
    
    free(code);
    
    isl_printer_free(printer);
    isl_ast_node_free(ast_tile);
    isl_ast_build_free(ast_build);
    
    tc_codegen_context_free(codegen_context);
    tc_ast_visitor_context_free(visitor_context);
}
