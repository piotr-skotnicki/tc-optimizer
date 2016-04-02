#include "for_decorator.h"
#include "config.h"

#include <isl/ctx.h>
#include <isl/id.h>
#include <isl/printer.h>
#include <isl/ast.h>
#include <isl/ast_build.h>

#include <string.h>

__isl_give isl_printer* tc_for_decorator_omp_task_nested(__isl_take isl_printer* printer, __isl_take isl_ast_print_options* options, __isl_keep isl_ast_node* node, void* user)
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
            
            printer = isl_ast_node_print(body, printer, options);
                        
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
        printer = isl_ast_node_for_print(node, printer, options);
    }
    
    isl_id_free(iterator);
    isl_ast_expr_free(expr);
    
    return printer;
}

__isl_give isl_printer* tc_for_decorator_omp_task_first(__isl_take isl_printer* printer, __isl_take isl_ast_print_options* options, __isl_keep isl_ast_node* node, void* user)
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
            
            options = isl_ast_print_options_set_print_for(options, NULL, NULL);
            printer = isl_ast_node_print(body, printer, options);
                        
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
        printer = isl_ast_node_for_print(node, printer, options);
    }
    
    isl_id_free(iterator);
    isl_ast_expr_free(expr);
    
    return printer;
}

__isl_give isl_printer* tc_for_decorator_omp_parallel_for_nested(__isl_take isl_printer* printer, __isl_take isl_ast_print_options* options, __isl_keep isl_ast_node* node, void* user)
{            
    isl_ast_expr* expr = isl_ast_node_for_get_iterator(node);
    
    isl_id* iterator = isl_ast_expr_get_id(expr);
    
    isl_id_list* parallel_id_list = (isl_id_list*)user;
        
    for (int i = 0; i < isl_id_list_n_id(parallel_id_list); ++i)
    {
        isl_id* parallel_id = isl_id_list_get_id(parallel_id_list, i);
                
        if (0 == strcmp(isl_id_get_name(iterator), isl_id_get_name(parallel_id)))
        {
            printer = isl_printer_start_line(printer);
            printer = isl_printer_print_str(printer, "#pragma omp parallel for");
            printer = isl_printer_end_line(printer);
        }
        
        isl_id_free(parallel_id);
    }

    printer = isl_ast_node_for_print(node, printer, options);
    
    isl_id_free(iterator);
    isl_ast_expr_free(expr);
    
    return printer;
}

__isl_give isl_printer* tc_for_decorator_omp_parallel_for_first(__isl_take isl_printer* printer, __isl_take isl_ast_print_options* options, __isl_keep isl_ast_node* node, void* user)
{            
    isl_ast_expr* expr = isl_ast_node_for_get_iterator(node);
    
    isl_id* iterator = isl_ast_expr_get_id(expr);
    
    isl_id_list* parallel_id_list = (isl_id_list*)user;
        
    for (int i = 0; i < isl_id_list_n_id(parallel_id_list); ++i)
    {
        isl_id* parallel_id = isl_id_list_get_id(parallel_id_list, i);
                
        if (0 == strcmp(isl_id_get_name(iterator), isl_id_get_name(parallel_id)))
        {
            printer = isl_printer_start_line(printer);
            printer = isl_printer_print_str(printer, "#pragma omp parallel for");
            printer = isl_printer_end_line(printer);
            
            options = isl_ast_print_options_set_print_for(options, NULL, NULL);
        }
        
        isl_id_free(parallel_id);
    }

    printer = isl_ast_node_for_print(node, printer, options);
    
    isl_id_free(iterator);
    isl_ast_expr_free(expr);
    
    return printer;
}
