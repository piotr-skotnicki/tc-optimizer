#include "serial_codegen.h"
#include "ast.h"
#include "codegen.h"
#include "scop.h"

#include <isl/ctx.h>
#include <isl/set.h>
#include <isl/id.h>
#include <isl/union_map.h>
#include <isl/union_set.h>
#include <isl/ast_build.h>
#include <isl/printer.h>

#include <stdio.h>

void tc_codegen_serial(struct tc_scop* scop, struct tc_options* options, __isl_take isl_union_map* S, __isl_take isl_set* tile, __isl_keep isl_id_list* iterators)
{
    isl_ctx* ctx = isl_union_map_get_ctx(S);
    
    isl_ast_build* ast_build = isl_ast_build_from_context(isl_set_copy(scop->pet->context));
    
    ast_build = isl_ast_build_set_iterators(ast_build, isl_id_list_copy(iterators));
    
    ast_build = isl_ast_build_set_at_each_domain(ast_build, &tc_ast_visitor_at_each_domain, scop);
                
    isl_union_map* S_prim = isl_union_map_intersect_range(S, isl_union_set_from_set(tile));
            
    isl_printer* printer = isl_printer_to_str(ctx);
    
    printer = isl_printer_set_output_format(printer, ISL_FORMAT_C);
        
    isl_ast_print_options* ast_options = isl_ast_print_options_alloc(ctx);
    
    if (!tc_options_is_set(options, NULL, "--use-macros"))
    {
        ast_options = isl_ast_print_options_set_print_user(ast_options, &tc_codegen_print_user, NULL);
    }

    isl_ast_node* ast_tile = isl_ast_build_ast_from_schedule(ast_build, S_prim);
    
    printer = tc_codegen_print_prologue(scop, options, printer);
    
    printer = isl_ast_node_print_macros(ast_tile, printer);
    
    if (tc_options_is_set(options, NULL, "--use-macros"))
    {
        printer = tc_codegen_print_statements_macros(scop, printer);
    }
    
    printer = isl_printer_print_str(printer, "#pragma scop\n");    
    printer = isl_ast_node_print(ast_tile, printer, ast_options);
    printer = isl_printer_print_str(printer, "#pragma endscop\n");    
    printer = tc_codegen_print_epilogue(scop, options, printer);
    
    char* code = isl_printer_get_str(printer);
    
    fprintf(options->output, "%s", code);
    
    free(code);
    
    isl_printer_free(printer);
    isl_ast_node_free(ast_tile);
    isl_ast_build_free(ast_build);
}
