#include "for_decorator.h"
#include "codegen.h"
#include "config.h"
#include "debug.h"
#include "utility.h"
#include "ast.h"

#include <isl/ctx.h>
#include <isl/id.h>
#include <isl/printer.h>
#include <isl/ast.h>
#include <isl/ast_build.h>

#include <string.h>
#include <stddef.h>
#include <stdlib.h>

__isl_give isl_printer* tc_for_decorator_omp_task_nested(__isl_take isl_printer* printer, __isl_take isl_ast_print_options* ast_options, __isl_keep isl_ast_node* node, void* user)
{
    isl_ctx* ctx = isl_printer_get_ctx(printer);
            
    isl_ast_expr* expr = isl_ast_node_for_get_iterator(node);
    
    isl_id* iterator = isl_ast_expr_get_id(expr);
    
    isl_id_list* task_id_list = (isl_id_list*)user;
        
    isl_bool is_task_id = isl_bool_false;
    
    for (int i = 0; i < isl_id_list_n_id(task_id_list); ++i)
    {
        isl_id* task_id = isl_id_list_get_id(task_id_list, i);
                        
        if (0 == strcmp(isl_id_get_name(iterator), isl_id_get_name(task_id)))
        {
            is_task_id = isl_bool_true;
            
            const char* name = isl_id_get_name(iterator);
            const char* type = isl_options_get_ast_iterator_type(ctx);
            
            isl_ast_node* body = isl_ast_node_for_get_body(node);
            isl_ast_expr* init = isl_ast_node_for_get_init(node);
            isl_ast_expr* cond = isl_ast_node_for_get_cond(node);
            isl_ast_expr* inc = isl_ast_node_for_get_inc(node);
            
            printer = isl_printer_start_line(printer);
            printer = isl_printer_print_str(printer, "{");
            printer = isl_printer_end_line(printer);
            printer = isl_printer_indent(printer, TC_CONF_INDENT_SIZE);
            
            printer = isl_printer_start_line(printer);
            printer = isl_printer_print_str(printer, "for (");
            printer = isl_printer_print_str(printer, type);
            printer = isl_printer_print_str(printer, " ");
            printer = isl_printer_print_str(printer, name);
            printer = isl_printer_print_str(printer, " = ");
            printer = isl_printer_print_ast_expr(printer, init);
            printer = isl_printer_print_str(printer, "; ");
            printer = isl_printer_print_ast_expr(printer, cond);
            printer = isl_printer_print_str(printer, "; ");
            printer = isl_printer_print_str(printer, name);
            printer = isl_printer_print_str(printer, " += ");
            printer = isl_printer_print_ast_expr(printer, inc);
            printer = isl_printer_print_str(printer, ") {");
            printer = isl_printer_end_line(printer);                                            
            
            printer = isl_printer_indent(printer, TC_CONF_INDENT_SIZE);
                                
            printer = isl_printer_start_line(printer);
            printer = isl_printer_print_str(printer, "#pragma omp task");
            printer = isl_printer_end_line(printer);
            
            printer = isl_ast_node_print(body, printer, ast_options);
                        
            printer = isl_printer_indent(printer, -TC_CONF_INDENT_SIZE);
            
            printer = isl_printer_start_line(printer);
            printer = isl_printer_print_str(printer, "}");
            printer = isl_printer_end_line(printer);
            
            printer = isl_printer_start_line(printer);
            printer = isl_printer_print_str(printer, "#pragma omp taskwait");
            printer = isl_printer_end_line(printer);            
            
            printer = isl_printer_indent(printer, -TC_CONF_INDENT_SIZE);
            printer = isl_printer_start_line(printer);
            printer = isl_printer_print_str(printer, "}");
            printer = isl_printer_end_line(printer);
            
            isl_ast_expr_free(init);
            isl_ast_expr_free(cond);
            isl_ast_expr_free(inc);
            isl_ast_node_free(body);
        }
        
        isl_id_free(task_id);
    }

    if (!is_task_id)
    {
        printer = isl_ast_node_for_print(node, printer, ast_options);
    }
    
    isl_id_free(iterator);
    isl_ast_expr_free(expr);
    
    return printer;
}

__isl_give isl_printer* tc_for_decorator_omp_task_first(__isl_take isl_printer* printer, __isl_take isl_ast_print_options* ast_options, __isl_keep isl_ast_node* node, void* user)
{
    isl_ctx* ctx = isl_printer_get_ctx(printer);
            
    isl_ast_expr* expr = isl_ast_node_for_get_iterator(node);
    
    isl_id* iterator = isl_ast_expr_get_id(expr);
    
    isl_id_list* task_id_list = (isl_id_list*)user;

    isl_bool is_task_id = isl_bool_false;
    
    for (int i = 0; i < isl_id_list_n_id(task_id_list); ++i)
    {
        isl_id* task_id = isl_id_list_get_id(task_id_list, i);
                        
        if (0 == strcmp(isl_id_get_name(iterator), isl_id_get_name(task_id)))
        {
            is_task_id = isl_bool_true;
            
            const char* name = isl_id_get_name(iterator);
            const char* type = isl_options_get_ast_iterator_type(ctx);
            
            isl_ast_node* body = isl_ast_node_for_get_body(node);
            isl_ast_expr* init = isl_ast_node_for_get_init(node);
            isl_ast_expr* cond = isl_ast_node_for_get_cond(node);
            isl_ast_expr* inc = isl_ast_node_for_get_inc(node);
            
            printer = isl_printer_start_line(printer);
            printer = isl_printer_print_str(printer, "{");
            printer = isl_printer_end_line(printer);
            printer = isl_printer_indent(printer, TC_CONF_INDENT_SIZE);
            
            printer = isl_printer_start_line(printer);
            printer = isl_printer_print_str(printer, "for (");
            printer = isl_printer_print_str(printer, type);
            printer = isl_printer_print_str(printer, " ");
            printer = isl_printer_print_str(printer, name);
            printer = isl_printer_print_str(printer, " = ");
            printer = isl_printer_print_ast_expr(printer, init);
            printer = isl_printer_print_str(printer, "; ");
            printer = isl_printer_print_ast_expr(printer, cond);
            printer = isl_printer_print_str(printer, "; ");
            printer = isl_printer_print_str(printer, name);
            printer = isl_printer_print_str(printer, " += ");
            printer = isl_printer_print_ast_expr(printer, inc);
            printer = isl_printer_print_str(printer, ") {");
            printer = isl_printer_end_line(printer);                                            
            
            printer = isl_printer_indent(printer, TC_CONF_INDENT_SIZE);
                                
            printer = isl_printer_start_line(printer);
            printer = isl_printer_print_str(printer, "#pragma omp task");
            printer = isl_printer_end_line(printer);
            
            ast_options = isl_ast_print_options_set_print_for(ast_options, NULL, NULL);
            printer = isl_ast_node_print(body, printer, ast_options);
                        
            printer = isl_printer_indent(printer, -TC_CONF_INDENT_SIZE);
            
            printer = isl_printer_start_line(printer);
            printer = isl_printer_print_str(printer, "}");
            printer = isl_printer_end_line(printer);
            
            printer = isl_printer_start_line(printer);
            printer = isl_printer_print_str(printer, "#pragma omp taskwait");
            printer = isl_printer_end_line(printer);            
            
            printer = isl_printer_indent(printer, -TC_CONF_INDENT_SIZE);
            printer = isl_printer_start_line(printer);
            printer = isl_printer_print_str(printer, "}");
            printer = isl_printer_end_line(printer);
            
            isl_ast_expr_free(init);
            isl_ast_expr_free(cond);
            isl_ast_expr_free(inc);
            isl_ast_node_free(body);
        }
        
        isl_id_free(task_id);
    }

    if (!is_task_id)
    {
        printer = isl_ast_node_for_print(node, printer, ast_options);
    }
    
    isl_id_free(iterator);
    isl_ast_expr_free(expr);
    
    return printer;
}

__isl_give isl_printer* tc_for_decorator_omp_parallel_for_nested(__isl_take isl_printer* printer, __isl_take isl_ast_print_options* ast_options, __isl_keep isl_ast_node* node, void* user)
{            
    isl_ctx* ctx = isl_printer_get_ctx(printer);

    isl_ast_expr* expr = isl_ast_node_for_get_iterator(node);
    
    isl_id* iterator = isl_ast_expr_get_id(expr);
    
    struct tc_codegen_context* context = (struct tc_codegen_context*)user;
        
    isl_id* annotation_id = isl_ast_node_get_annotation(node);
    
    struct tc_ast_for_annotation* annotation = (struct tc_ast_for_annotation*)isl_id_get_user(annotation_id);
    
    const char* name = isl_id_get_name(iterator);
    const char* type = isl_options_get_ast_iterator_type(ctx);
    
    char lb_variable[512], ub_variable[512];
    snprintf(lb_variable, sizeof(lb_variable), "%s_lb", name);
    snprintf(ub_variable, sizeof(ub_variable), "%s_ub", name);
        
    isl_ast_node* body = isl_ast_node_for_get_body(node);
    isl_ast_expr* init = isl_ast_node_for_get_init(node);
    isl_ast_expr* cond = isl_ast_node_for_get_cond(node);
    isl_ast_expr* inc = isl_ast_node_for_get_inc(node);
    
    isl_ast_expr* cond_rhs = isl_ast_expr_get_op_arg(cond, 1);
    
    isl_bool inline_variables = isl_bool_true;
    
    if (!tc_options_is_inline(context->options))
    {
        if (isl_ast_expr_get_type(cond_rhs) == isl_ast_expr_op || isl_ast_expr_get_type(init) == isl_ast_expr_op)
        {
            inline_variables = isl_bool_false;
        }
    }
    
    if (!inline_variables)
    {
        char* init_str = isl_ast_expr_to_C_str(init);
        char* cond_rhs_str = isl_ast_expr_to_C_str(cond_rhs);

        char bounds_declaration[512];
        snprintf(bounds_declaration, sizeof(bounds_declaration), "const %s %s = %s, %s = %s;", type, lb_variable, init_str, ub_variable, cond_rhs_str);

        printer = isl_printer_start_line(printer);
        printer = isl_printer_print_str(printer, bounds_declaration);
        printer = isl_printer_end_line(printer);

        free(init_str);
        free(cond_rhs_str);
    }

    isl_bool is_new_parallel_iterator = isl_bool_false;

    if (annotation->is_parallel)
    {
        is_new_parallel_iterator = isl_bool_true;
    }

    if (is_new_parallel_iterator)
    {                        
        printer = isl_printer_start_line(printer);
        printer = isl_printer_print_str(printer, "#pragma omp parallel for");

        printer = isl_printer_end_line(printer);
    }

    if (!inline_variables)
    {
        isl_ast_expr_free(init);
        init = isl_ast_expr_from_id(isl_id_alloc(ctx, lb_variable, NULL));
        
        cond = isl_ast_expr_set_op_arg(cond, 1, isl_ast_expr_from_id(isl_id_alloc(ctx, ub_variable, NULL)));
    }
    
    char* init_str = isl_ast_expr_to_C_str(init);
    char* cond_str = isl_ast_expr_to_C_str(cond);
    char* inc_str = isl_ast_expr_to_C_str(inc);

    char for_declaration[4096];
    snprintf(for_declaration, sizeof(for_declaration), "for (register %s %s = %s; %s; %s += %s) {", type, name, init_str, cond_str, name, inc_str);
    //snprintf(for_declaration, sizeof(for_declaration), "for (%s %s = %s; %s; %s += %s) {", type, name, init_str, cond_str, name, inc_str);

    free(init_str);
    free(cond_str);
    free(inc_str);

    printer = isl_printer_start_line(printer);
    printer = isl_printer_print_str(printer, for_declaration);
    printer = isl_printer_end_line(printer);
    
    printer = isl_printer_indent(printer, TC_CONF_INDENT_SIZE);
    
    if (isl_ast_node_get_type(body) == isl_ast_node_block)
    {
        isl_ast_node_list* block_statements = isl_ast_node_block_get_children(body);
        
        for (int i = 0; i < isl_ast_node_list_n_ast_node(block_statements); ++i)
        {
            isl_ast_node* block_statement = isl_ast_node_list_get_ast_node(block_statements, i);
         
            printer = isl_ast_node_print(block_statement, printer, isl_ast_print_options_copy(ast_options));
            
            isl_ast_node_free(block_statement);
        }
        
        isl_ast_node_list_free(block_statements);
    }
    else
    {
        printer = isl_ast_node_print(body, printer, isl_ast_print_options_copy(ast_options));
    }
    
    printer = isl_printer_indent(printer, -TC_CONF_INDENT_SIZE);

    printer = isl_printer_start_line(printer);
    printer = isl_printer_print_str(printer, "}");
    printer = isl_printer_end_line(printer);
        
    isl_ast_expr_free(init);
    isl_ast_expr_free(cond);
    isl_ast_expr_free(inc);
    isl_ast_node_free(body);
    isl_ast_expr_free(cond_rhs);
    
    isl_id_free(annotation_id);
    isl_id_free(iterator);
    isl_ast_expr_free(expr);
    
    isl_ast_print_options_free(ast_options);
    
    return printer;
}

__isl_give isl_printer* tc_for_decorator_omp_parallel_for_first(__isl_take isl_printer* printer, __isl_take isl_ast_print_options* ast_options, __isl_keep isl_ast_node* node, void* user)
{
    isl_ctx* ctx = isl_printer_get_ctx(printer);

    isl_ast_expr* expr = isl_ast_node_for_get_iterator(node);
    
    isl_id* iterator = isl_ast_expr_get_id(expr);
    
    struct tc_codegen_context* context = (struct tc_codegen_context*)user;
        
    isl_id* annotation_id = isl_ast_node_get_annotation(node);
    
    struct tc_ast_for_annotation* annotation = (struct tc_ast_for_annotation*)isl_id_get_user(annotation_id);
    
    const char* name = isl_id_get_name(iterator);
    const char* type = isl_options_get_ast_iterator_type(ctx);
    
    char lb_variable[512], ub_variable[512];
    //snprintf(lb_variable, sizeof(lb_variable), "%s_lb", name);
    //snprintf(ub_variable, sizeof(ub_variable), "%s_ub", name);
    snprintf(lb_variable, sizeof(lb_variable), "%s_%d_lb", name, rand() % 100000);
    snprintf(ub_variable, sizeof(ub_variable), "%s_%d_ub", name, rand() % 100000);
        
    isl_ast_node* body = isl_ast_node_for_get_body(node);
    isl_ast_expr* init = isl_ast_node_for_get_init(node);
    isl_ast_expr* cond = isl_ast_node_for_get_cond(node);
    isl_ast_expr* inc = isl_ast_node_for_get_inc(node);
    
    isl_ast_expr* cond_rhs = isl_ast_expr_get_op_arg(cond, 1);
    
    isl_bool inline_variables = isl_bool_true;
    
    if (!tc_options_is_inline(context->options))
    {
        if (isl_ast_expr_get_type(cond_rhs) == isl_ast_expr_op || isl_ast_expr_get_type(init) == isl_ast_expr_op)
        {
            inline_variables = isl_bool_false;
        }
    }
    
    if (!inline_variables)
    {
        char* init_str = isl_ast_expr_to_C_str(init);
        char* cond_rhs_str = isl_ast_expr_to_C_str(cond_rhs);

        char bounds_declaration[512];
        snprintf(bounds_declaration, sizeof(bounds_declaration), "const %s %s = %s, %s = %s;", type, lb_variable, init_str, ub_variable, cond_rhs_str);

        printer = isl_printer_start_line(printer);
        printer = isl_printer_print_str(printer, bounds_declaration);
        printer = isl_printer_end_line(printer);

        free(init_str);
        free(cond_rhs_str);
    }

    isl_bool is_new_parallel_iterator = isl_bool_false;

    if (!context->in_parallel_region && annotation->is_parallel)
    {
        is_new_parallel_iterator = isl_bool_true;
    }

    if (is_new_parallel_iterator)
    {
        context->in_parallel_region = isl_bool_true;
                        
        printer = isl_printer_start_line(printer);
        printer = isl_printer_print_str(printer, "#pragma omp parallel for");
        printer = isl_printer_end_line(printer);
    }

    if (!inline_variables)
    {
        isl_ast_expr_free(init);
        init = isl_ast_expr_from_id(isl_id_alloc(ctx, lb_variable, NULL));
        
        cond = isl_ast_expr_set_op_arg(cond, 1, isl_ast_expr_from_id(isl_id_alloc(ctx, ub_variable, NULL)));
    }
    
    char* init_str = isl_ast_expr_to_C_str(init);
    char* cond_str = isl_ast_expr_to_C_str(cond);
    char* inc_str = isl_ast_expr_to_C_str(inc);

    char for_declaration[4096];
    snprintf(for_declaration, sizeof(for_declaration), "for (register %s %s = %s; %s; %s += %s) {", type, name, init_str, cond_str, name, inc_str);
    //snprintf(for_declaration, sizeof(for_declaration), "for (%s %s = %s; %s; %s += %s) {", type, name, init_str, cond_str, name, inc_str);

    free(init_str);
    free(cond_str);
    free(inc_str);

    printer = isl_printer_start_line(printer);
    printer = isl_printer_print_str(printer, for_declaration);
    printer = isl_printer_end_line(printer);
    
    printer = isl_printer_indent(printer, TC_CONF_INDENT_SIZE);
    
    if (isl_ast_node_get_type(body) == isl_ast_node_block)
    {
        isl_ast_node_list* block_statements = isl_ast_node_block_get_children(body);
        
        for (int i = 0; i < isl_ast_node_list_n_ast_node(block_statements); ++i)
        {
            isl_ast_node* block_statement = isl_ast_node_list_get_ast_node(block_statements, i);
         
            printer = isl_ast_node_print(block_statement, printer, isl_ast_print_options_copy(ast_options));
            
            isl_ast_node_free(block_statement);
        }
        
        isl_ast_node_list_free(block_statements);
    }
    else
    {
        printer = isl_ast_node_print(body, printer, isl_ast_print_options_copy(ast_options));
    }
    
    printer = isl_printer_indent(printer, -TC_CONF_INDENT_SIZE);

    printer = isl_printer_start_line(printer);
    printer = isl_printer_print_str(printer, "}");
    printer = isl_printer_end_line(printer);
    
    if (is_new_parallel_iterator)
    {
        context->in_parallel_region = isl_bool_false;
    }
    
    isl_ast_expr_free(init);
    isl_ast_expr_free(cond);
    isl_ast_expr_free(inc);
    isl_ast_node_free(body);
    isl_ast_expr_free(cond_rhs);
    
    isl_id_free(annotation_id);
    isl_id_free(iterator);
    isl_ast_expr_free(expr);
    
    isl_ast_print_options_free(ast_options);
    
    return printer;
}

__isl_give isl_printer* tc_for_decorator_omp_for_first(__isl_take isl_printer* printer, __isl_take isl_ast_print_options* options, __isl_keep isl_ast_node* node, void* user)
{
    isl_ast_expr* expr = isl_ast_node_for_get_iterator(node);
    
    isl_id* iterator = isl_ast_expr_get_id(expr);
    
    struct tc_codegen_context* context = (struct tc_codegen_context*)user;

    isl_id* annotation_id = isl_ast_node_get_annotation(node);
    
    struct tc_ast_for_annotation* annotation = (struct tc_ast_for_annotation*)isl_id_get_user(annotation_id);

    isl_bool is_new_parallel_iterator = isl_bool_false;

    if (!context->in_parallel_region && annotation->is_parallel)
    {
        is_new_parallel_iterator = isl_bool_true;
    }

    if (is_new_parallel_iterator)
    {
        context->in_parallel_region = isl_bool_true;
        
        printer = isl_printer_start_line(printer);
        printer = isl_printer_print_str(printer, "#pragma omp for");
        printer = isl_printer_end_line(printer);
    }

    printer = isl_ast_node_for_print(node, printer, options);

    if (is_new_parallel_iterator)
    {
        context->in_parallel_region = isl_bool_false;
    }
    
    isl_id_free(annotation_id);
    isl_id_free(iterator);
    isl_ast_expr_free(expr);
    
    return printer;
}

__isl_give isl_printer* tc_for_decorator_omp_teams_distribute_parallel_for(__isl_take isl_printer* printer, __isl_take isl_ast_print_options* options, __isl_keep isl_ast_node* node, void* user)
{
    isl_ctx* ctx = isl_printer_get_ctx(printer);
    
    isl_ast_expr* expr = isl_ast_node_for_get_iterator(node);
    
    isl_id* iterator = isl_ast_expr_get_id(expr);
    
    struct tc_codegen_context* context = (struct tc_codegen_context*)user;

    isl_id* annotation_id = isl_ast_node_get_annotation(node);
    
    struct tc_ast_for_annotation* annotation = (struct tc_ast_for_annotation*)isl_id_get_user(annotation_id);

    isl_bool is_new_parallel_iterator = isl_bool_false;

    if (!context->in_parallel_region && annotation->is_parallel)
    {
        is_new_parallel_iterator = isl_bool_true;
    }

    if (is_new_parallel_iterator)
    {
        context->in_parallel_region = isl_bool_true;

        //////////////

        isl_bool is_first_read_and_write = isl_bool_true;

        for (int i = 0; i < context->scop->pet->n_array; ++i)
        {
            struct pet_array* array = context->scop->pet->arrays[i];

            const char* array_name = isl_set_get_tuple_name(array->extent);

            isl_bool is_read = isl_bool_false;
            isl_bool is_write = isl_bool_false;

            for (int j = 0; j < isl_id_list_n_id(annotation->nested_statements); ++j)            
            {
                isl_id* nested_statement = isl_id_list_get_id(annotation->nested_statements, j);

                isl_union_map* reads = tc_get_union_map_for_input_tuple(context->scop->reads, isl_id_get_name(nested_statement));
                isl_union_map* writes = context->scop->writes;

                if (tc_union_map_has_output_tuple(reads, array_name) == isl_bool_true)
                {
                    is_read = isl_bool_true;
                }

                if (tc_union_map_has_output_tuple(writes, array_name) == isl_bool_true)
                {
                    is_write = isl_bool_true;
                }

                isl_union_map_free(reads);

                isl_id_free(nested_statement);
            }

            if (is_read == isl_bool_true && is_write == isl_bool_true)
            {
                if (is_first_read_and_write == isl_bool_true)
                {
                    is_first_read_and_write = isl_bool_false;
                    printer = isl_printer_start_line(printer);
                    printer = isl_printer_print_str(printer, "#pragma omp target update");
                }

                printer = isl_printer_print_str(printer, " to(");
                printer = isl_printer_print_str(printer, array_name);

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
        }

        if (is_first_read_and_write == isl_bool_false)
        {
            printer = isl_printer_end_line(printer);
        }

        //////////////
        
        printer = isl_printer_start_line(printer);
        printer = isl_printer_print_str(printer, "#pragma omp target teams distribute parallel for schedule(static,1) thread_limit(32)");
        printer = isl_printer_end_line(printer);
    }

    printer = isl_ast_node_for_print(node, printer, options);

    if (is_new_parallel_iterator)
    {
        context->in_parallel_region = isl_bool_false;

        //////////////

        isl_bool is_first_read_and_write = isl_bool_true;

        for (int i = 0; i < context->scop->pet->n_array; ++i)
        {
            struct pet_array* array = context->scop->pet->arrays[i];

            const char* array_name = isl_set_get_tuple_name(array->extent);

            isl_bool is_read = isl_bool_false;
            isl_bool is_write = isl_bool_false;

            for (int j = 0; j < isl_id_list_n_id(annotation->nested_statements); ++j)            
            {
                isl_id* nested_statement = isl_id_list_get_id(annotation->nested_statements, j);

                //isl_union_map* reads = context->scop->reads;
                isl_union_map* writes = tc_get_union_map_for_input_tuple(context->scop->writes, isl_id_get_name(nested_statement));

                /*
                if (tc_union_map_has_output_tuple(reads, array_name) == isl_bool_true)
                {
                    is_read = isl_bool_true;
                }
                 */                

                if (tc_union_map_has_output_tuple(writes, array_name) == isl_bool_true)
                {
                    is_write = isl_bool_true;
                }

                isl_union_map_free(writes);

                isl_id_free(nested_statement);
            }

            if (/*is_read == isl_bool_true && */is_write == isl_bool_true)
            {
                if (is_first_read_and_write == isl_bool_true)
                {
                    is_first_read_and_write = isl_bool_false;
                    printer = isl_printer_start_line(printer);
                    printer = isl_printer_print_str(printer, "#pragma omp target update");
                }

                printer = isl_printer_print_str(printer, " from(");
                printer = isl_printer_print_str(printer, array_name);

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
        }

        if (is_first_read_and_write == isl_bool_false)
        {
            printer = isl_printer_end_line(printer);
        }

        //////////////
    }
    
    isl_id_free(annotation_id);
    isl_id_free(iterator);
    isl_ast_expr_free(expr);
    
    return printer;
}
