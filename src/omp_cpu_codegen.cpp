#include "omp_cpu_codegen.h"
#include "codegen.h"
#include "ast.h"
#include "for_decorator.h"
#include "scop.h"
#include "debug.h"
#include "utility.h"

#include <isl/ctx.h>
#include <isl/id.h>
#include <isl/set.h>
#include <isl/union_map.h>
#include <isl/union_set.h>
#include <isl/ast_build.h>
#include <isl/printer.h>

#include <stdio.h>
#include <stddef.h>

void tc_codegen_omp_parallel_for(struct tc_scop* scop, struct tc_options* options, __isl_take isl_union_map* S, __isl_take isl_set* tile, __isl_keep isl_id_list* iterators, __isl_keep isl_id_list* parallel_iterators, int nested)
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
    
    isl_union_map* S_prim = isl_union_map_intersect_range(S, isl_union_set_from_set(tile));

    tc_debug_umap(S_prim, "S_prim");

    isl_printer* printer = isl_printer_to_str(ctx);    
    
    printer = isl_printer_set_output_format(printer, ISL_FORMAT_C);
            
    isl_ast_print_options* ast_options = isl_ast_print_options_alloc(ctx);
    
    if (!tc_options_is_set(options, NULL, "--use-macros"))
    {
        ast_options = isl_ast_print_options_set_print_user(ast_options, &tc_codegen_print_user, NULL);
    }
        
    struct tc_codegen_context* codegen_context = tc_codegen_context_alloc(ctx);
    codegen_context->options = options;
    
    if (nested)
    {
        ast_options = isl_ast_print_options_set_print_for(ast_options, &tc_for_decorator_omp_parallel_for_nested, codegen_context);
        printer = isl_printer_print_str(printer, "#include <omp.h>\n");
    }
    else
    {
        ast_options = isl_ast_print_options_set_print_for(ast_options, &tc_for_decorator_omp_parallel_for_first, codegen_context);
    }
    
    isl_ast_node* ast_tile = isl_ast_build_ast_from_schedule(ast_build, S_prim);
    
    printer = tc_codegen_print_prologue(scop, options, printer);
    
    printer = isl_ast_node_print_macros(ast_tile, printer);
    
    if (tc_options_is_set(options, NULL, "--use-macros"))
    {
        printer = tc_codegen_print_statements_macros(scop, printer);
    }
    
    printer = isl_printer_print_str(printer, "#pragma scop\n");
    
    if (nested)
    {
        printer = isl_printer_print_str(printer, "omp_set_nested(1);\n");
    }
    
    printer = isl_ast_node_print(ast_tile, printer, ast_options);    
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

void tc_codegen_omp_task_for(struct tc_scop* scop, struct tc_options* options, __isl_take isl_union_map* S, __isl_take isl_set* tile, __isl_keep isl_id_list* iterators, __isl_keep isl_id_list* parallel_iterators, int nested)
{
    isl_ctx* ctx = isl_union_map_get_ctx(S);
    
    isl_ast_build* ast_build = isl_ast_build_from_context(isl_set_copy(scop->pet->context));
                
    ast_build = isl_ast_build_set_iterators(ast_build, isl_id_list_copy(iterators));

    struct tc_ast_visitor_context* visitor_context = tc_ast_visitor_context_alloc(ctx);
    visitor_context->scop = scop;
    
    std::vector<struct tc_ast_for_annotation*> annotations_stack;
    visitor_context->annotations_stack = &annotations_stack;
    
    ast_build = isl_ast_build_set_at_each_domain(ast_build, &tc_ast_visitor_at_each_domain, visitor_context);

    isl_union_map* S_prim = isl_union_map_intersect_range(S, isl_union_set_from_set(tile));
            
    isl_printer* printer = isl_printer_to_str(ctx);    
    
    printer = isl_printer_set_output_format(printer, ISL_FORMAT_C);
        
    isl_ast_print_options* ast_options = isl_ast_print_options_alloc(ctx);
    
    if (!tc_options_is_set(options, NULL, "--use-macros"))
    {
        ast_options = isl_ast_print_options_set_print_user(ast_options, &tc_codegen_print_user, NULL);
    }
    
    if (nested)
    {
        ast_options = isl_ast_print_options_set_print_for(ast_options, &tc_for_decorator_omp_task_nested, parallel_iterators);
    }
    else
    {
        ast_options = isl_ast_print_options_set_print_for(ast_options, &tc_for_decorator_omp_task_first, parallel_iterators);
    }
    
    isl_ast_node* ast_tile = isl_ast_build_ast_from_schedule(ast_build, S_prim);

    printer = tc_codegen_print_prologue(scop, options, printer);
    
    printer = isl_ast_node_print_macros(ast_tile, printer);
    
    if (tc_options_is_set(options, NULL, "--use-macros"))
    {
        printer = tc_codegen_print_statements_macros(scop, printer);
    }
    
    printer = isl_printer_print_str(printer, "#pragma scop\n");
    printer = isl_printer_print_str(printer, "#pragma omp parallel\n");
    printer = isl_printer_print_str(printer, "#pragma omp single\n");
    printer = isl_ast_node_print(ast_tile, printer, ast_options);
    printer = isl_printer_print_str(printer, "#pragma endscop\n");
    
    printer = tc_codegen_print_epilogue(scop, options, printer);
    
    char* code = isl_printer_get_str(printer);
    
    fprintf(options->output, "%s", code);

    fflush(options->output);
    
    free(code);
    
    isl_printer_free(printer);
    isl_ast_node_free(ast_tile);
    isl_ast_build_free(ast_build);
}
