#include "debug.h"

#include <isl/printer.h>

#include <stdarg.h>

int tc_debug_flag;

void tc_debug(const char* msg, ...)
{
    if (!tc_debug_flag)
        return;
    
    va_list varargs;
    va_start(varargs, msg);
    vfprintf(stderr, msg, varargs);    
    va_end(varargs);
    
    fprintf(stderr, "\n\n");    
    fflush(stderr);
}

void tc_debug_space(__isl_keep isl_space* set, const char* msg, ...)
{
    if (!tc_debug_flag)
        return;
    
    va_list varargs;
    va_start(varargs, msg);
    vfprintf(stderr, msg, varargs);    
    va_end(varargs);
    
    char* str = isl_space_to_str(set);
    
    fprintf(stderr, " := %s;\n\n", str);
    fflush(stderr);
    
    free(str);
}

void tc_debug_set(__isl_keep isl_set* set, const char* msg, ...)
{
    if (!tc_debug_flag)
        return;
    
    va_list varargs;
    va_start(varargs, msg);
    vfprintf(stderr, msg, varargs);    
    va_end(varargs);
    
    char* str = isl_set_to_str(set);
    
    fprintf(stderr, " := %s;\n\n", str);
    fflush(stderr);
    
    free(str);
}

void tc_debug_set_latex(__isl_keep isl_set* set, const char* msg, ...)
{
    if (!tc_debug_flag)
        return;
    
    va_list varargs;
    va_start(varargs, msg);
    vfprintf(stderr, msg, varargs);
    va_end(varargs);
    
    isl_ctx* ctx = isl_set_get_ctx(set);
            
    isl_printer* printer = isl_printer_to_str(ctx);    

    printer = isl_printer_set_output_format(printer, ISL_FORMAT_LATEX);

    printer = isl_printer_flush(printer);

    printer = isl_printer_print_set(printer, set);

    char* str = isl_printer_get_str(printer);

    fprintf(stderr, " := %s;\n\n", str);    
    fflush(stderr);
    
    isl_printer_free(printer);
    
    free(str);
}

void tc_debug_uset(__isl_keep isl_union_set* uset, const char* msg, ...)
{
    if (!tc_debug_flag)
        return;
    
    va_list varargs;
    va_start(varargs, msg);
    vfprintf(stderr, msg, varargs);
    va_end(varargs);
    
    char* str = isl_union_set_to_str(uset);
    
    fprintf(stderr, " := %s;\n\n", str);    
    fflush(stderr);
    
    free(str);
}

void tc_debug_uset_latex(__isl_keep isl_union_set* uset, const char* msg, ...)
{
    if (!tc_debug_flag)
        return;
    
    va_list varargs;
    va_start(varargs, msg);
    vfprintf(stderr, msg, varargs);
    va_end(varargs);
    
    isl_ctx* ctx = isl_union_set_get_ctx(uset);
            
    isl_printer* printer = isl_printer_to_str(ctx);    

    printer = isl_printer_set_output_format(printer, ISL_FORMAT_LATEX);

    printer = isl_printer_flush(printer);

    printer = isl_printer_print_union_set(printer, uset);

    char* str = isl_printer_get_str(printer);

    fprintf(stderr, " := %s;\n\n", str);    
    fflush(stderr);
    
    isl_printer_free(printer);
    
    free(str);
}

void tc_debug_map(__isl_keep isl_map* map, const char* msg, ...)
{
    if (!tc_debug_flag)
        return;
    
    va_list varargs;
    va_start(varargs, msg);
    vfprintf(stderr, msg, varargs);
    va_end(varargs);
    
    char* str = isl_map_to_str(map);
    
    fprintf(stderr, " := %s;\n\n", str);    
    fflush(stderr);
    
    free(str);
}

void tc_debug_map_latex(__isl_keep isl_map* map, const char* msg, ...)
{
    if (!tc_debug_flag)
        return;
    
    va_list varargs;
    va_start(varargs, msg);
    vfprintf(stderr, msg, varargs);
    va_end(varargs);
    
    isl_ctx* ctx = isl_map_get_ctx(map);
            
    isl_printer* printer = isl_printer_to_str(ctx);    

    printer = isl_printer_set_output_format(printer, ISL_FORMAT_LATEX);

    printer = isl_printer_flush(printer);

    printer = isl_printer_print_map(printer, map);

    char* str = isl_printer_get_str(printer);

    fprintf(stderr, " := %s;\n\n", str);    
    fflush(stderr);
    
    isl_printer_free(printer);
    
    free(str);
}

void tc_debug_umap(__isl_keep isl_union_map* umap, const char* msg, ...)
{
    if (!tc_debug_flag)
        return;
    
    va_list varargs;
    va_start(varargs, msg);
    vfprintf(stderr, msg, varargs);
    va_end(varargs);
    
    char* str = isl_union_map_to_str(umap);
    
    fprintf(stderr, " := %s;\n\n", str);    
    fflush(stderr);
    
    free(str);
}

void tc_debug_umap_latex(__isl_keep isl_union_map* umap, const char* msg, ...)
{
    if (!tc_debug_flag)
        return;
    
    va_list varargs;
    va_start(varargs, msg);
    vfprintf(stderr, msg, varargs);
    va_end(varargs);
    
    isl_ctx* ctx = isl_union_map_get_ctx(umap);
            
    isl_printer* printer = isl_printer_to_str(ctx);    

    printer = isl_printer_set_output_format(printer, ISL_FORMAT_LATEX);

    printer = isl_printer_flush(printer);

    printer = isl_printer_print_union_map(printer, umap);

    char* str = isl_printer_get_str(printer);

    fprintf(stderr, " := %s;\n\n", str);    
    fflush(stderr);
    
    isl_printer_free(printer);
    
    free(str);
}

void tc_debug_val(__isl_keep isl_val* val, const char* msg, ...)
{
    if (!tc_debug_flag)
        return;
    
    va_list varargs;
    va_start(varargs, msg);
    vfprintf(stderr, msg, varargs);
    va_end(varargs);
    
    char* str = isl_val_to_str(val);
    
    fprintf(stderr, " := %s;\n\n", str);    
    fflush(stderr);
    
    free(str);
}

void tc_debug_bool(isl_bool condition, const char* msg, ...)
{
    if (!tc_debug_flag)
        return;
    
    fprintf(stderr, "(");
    
    va_list varargs;
    va_start(varargs, msg);
    vfprintf(stderr, msg, varargs);
    va_end(varargs);
    
    const char* str = NULL;
    
    switch (condition)
    {
        case isl_bool_true: str = "True"; break;
        case isl_bool_false: str = "False"; break;
        case isl_bool_error: str = "Error"; break;
        default: str = "Unknown"; break;
    }
        
    fprintf(stderr, ") = %s;\n\n", str);
    fflush(stderr);
}

void tc_debug_set_card(__isl_keep isl_set* set, const char* msg, ...)
{
    if (!tc_debug_flag)
        return;
        
    fprintf(stderr, "card ");
    
    va_list varargs;
    va_start(varargs, msg);
    vfprintf(stderr, msg, varargs);
    va_end(varargs);
    
    isl_pw_qpolynomial* qpoly = isl_set_card(isl_set_copy(set));
    
    tc_debug_pw_qpolynomial(qpoly, "");
    
    isl_pw_qpolynomial_free(qpoly);
}

void tc_debug_map_card(__isl_keep isl_map* map, const char* msg, ...)
{
    if (!tc_debug_flag)
        return;
        
    fprintf(stderr, "card ");
    
    va_list varargs;
    va_start(varargs, msg);
    vfprintf(stderr, msg, varargs);
    va_end(varargs);
    
    isl_set* wrapped_map = isl_map_wrap(isl_map_copy(map));
    
    isl_pw_qpolynomial* qpoly = isl_set_card(wrapped_map);
    
    tc_debug_pw_qpolynomial(qpoly, "");
    
    isl_pw_qpolynomial_free(qpoly);
}

void tc_debug_pw_qpolynomial(__isl_keep isl_pw_qpolynomial* qpoly, const char* msg, ...)
{
    if (!tc_debug_flag)
        return;
    
    va_list varargs;
    va_start(varargs, msg);
    vfprintf(stderr, msg, varargs);
    va_end(varargs);
    
    isl_ctx* ctx = isl_pw_qpolynomial_get_ctx(qpoly);
            
    isl_printer* printer = isl_printer_to_str(ctx);    

    //printer = isl_printer_set_output_format(printer, ISL_FORMAT_C);

    printer = isl_printer_flush(printer);

    printer = isl_printer_print_pw_qpolynomial(printer, qpoly);

    char* str = isl_printer_get_str(printer);

    fprintf(stderr, " := %s;\n\n", str);    
    fflush(stderr);
    
    isl_printer_free(printer);
    
    free(str);
}

void tc_debug_qpolynomial(__isl_keep isl_qpolynomial* poly, const char* msg, ...)
{
    if (!tc_debug_flag)
        return;
        
    va_list varargs;
    va_start(varargs, msg);
    vfprintf(stderr, msg, varargs);
    va_end(varargs);
    
    isl_ctx* ctx = isl_qpolynomial_get_ctx(poly);
            
    isl_printer* printer = isl_printer_to_str(ctx);    

    //printer = isl_printer_set_output_format(printer, ISL_FORMAT_C);

    printer = isl_printer_flush(printer);

    printer = isl_printer_print_qpolynomial(printer, poly);

    char* str = isl_printer_get_str(printer);
    
    fprintf(stderr, " := %s;\n\n", str);    
    fflush(stderr);
    
    isl_printer_free(printer);
    
    free(str);
}
