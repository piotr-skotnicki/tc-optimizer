#include "codegen.h"
#include "ast.h"
#include "utility.h"
#include "options.h"
#include "config.h"

#include <isl/ctx.h>
#include <isl/id.h>
#include <isl/ast.h>
#include <isl/ast_build.h>
#include <isl/printer.h>

#include <pet.h>

#include <stdlib.h>
#include <stddef.h>

struct tc_codegen_context* tc_codegen_context_alloc(__isl_keep isl_ctx* ctx)
{
    struct tc_codegen_context* context = (struct tc_codegen_context*)malloc(sizeof(struct tc_codegen_context));
    
    context->in_parallel_region = isl_bool_false;    
    context->options = NULL;
    
    return context;
}

void tc_codegen_context_free(struct tc_codegen_context* context)
{
    if (NULL == context)
        return;
            
    free(context);
}

__isl_give isl_printer* tc_codegen_print_user(__isl_take isl_printer* printer, __isl_take isl_ast_print_options* options, __isl_keep isl_ast_node* node, void* user)
{
    isl_id* annotation_id = isl_ast_node_get_annotation(node);
    
    struct tc_ast_stmt_annotation* annotation = (struct tc_ast_stmt_annotation*)isl_id_get_user(annotation_id);
    
    printer = pet_stmt_print_body(annotation->stmt, printer, annotation->id2expr);

    isl_ast_print_options_free(options);
    isl_id_free(annotation_id);

    return printer;
}

__isl_give isl_printer* tc_codegen_print_statements_macros(struct tc_scop* scop, __isl_take isl_printer* printer)
{
    struct pet_scop* pet = scop->pet;
    
    char buff[2048];
    
    for (int i = 0; i < pet->n_stmt; ++i) 
    {
        struct pet_stmt* stmt = pet->stmts[i];

        if (!pet_stmt_is_kill(stmt))
        {
            isl_set* statement_domain = isl_set_copy(stmt->domain);
            
            isl_ctx* ctx = isl_set_get_ctx(statement_domain);
            
            if (pet->stmts[i]->n_arg > 0)
            {
                statement_domain = isl_map_domain(isl_set_unwrap(statement_domain));
            }
            
            const char* label = isl_set_get_tuple_name(statement_domain);
            
            int n_set_dim = isl_set_dim(statement_domain, isl_dim_set);
            
            isl_id_list* iterators = isl_id_list_alloc(ctx, n_set_dim);
            
            for (int j = 0; j < n_set_dim; ++j)
            {            
                iterators = isl_id_list_add(iterators, isl_id_alloc(ctx, isl_set_get_dim_name(statement_domain, isl_dim_set, j), NULL));
            }
            
            sprintf(buff, "#define %s_I(%s) ", label, tc_comma(iterators).c_str());            
            
            printer = isl_printer_print_str(printer, buff);
            
            isl_ast_build* ast_stmt_build = isl_ast_build_from_context(isl_set_copy(statement_domain));
            
            isl_id_to_ast_expr* id2ast = pet_stmt_build_ast_exprs(stmt, ast_stmt_build, NULL, NULL, NULL, NULL);
            
            isl_ast_build_free(ast_stmt_build);
                        
            printer = pet_stmt_print_body(stmt, printer, id2ast);
            
            isl_id_to_ast_expr_free(id2ast);
            
            sprintf(buff, "#define %s(%s) %s_I(%s)\n", label, tc_comma(iterators).c_str(), label, tc_parens_comma(iterators).c_str());
            
            printer = isl_printer_print_str(printer, buff);
            
            isl_id_list_free(iterators);            
            isl_set_free(statement_domain);
        }
    }
    
    return printer;
}

__isl_give isl_printer* tc_codegen_print_prologue(struct tc_scop* scop, struct tc_options* options, __isl_take isl_printer* printer)
{
    char* command_line = tc_options_get_command_line(options);

    printer = isl_printer_print_str(printer, "/* TC Optimizing Compiler " TC_CONF_VERSION " */\n");
    printer = isl_printer_print_str(printer, "/* ");
    printer = isl_printer_print_str(printer, command_line);
    printer = isl_printer_print_str(printer, " */\n");
    
    free(command_line);
    
    return printer;
}

__isl_give isl_printer* tc_codegen_print_epilogue(struct tc_scop* scop, struct tc_options* options, __isl_take isl_printer* printer)
{
    //printer = isl_printer_print_str(printer, "/**/\n");
    
    return printer;
}
