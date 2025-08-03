#include "debug.h"
#include "input_output.h"

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

void tc_warn(const char* msg, ...)
{
    fprintf(stderr, "\033[1;33mWarning: \033[0m");

    va_list varargs;
    va_start(varargs, msg);
    vfprintf(stderr, msg, varargs);
    va_end(varargs);

    fprintf(stderr, "\n\n");
    fflush(stderr);
}

void tc_error(const char* msg, ...)
{
    fprintf(stderr, "\033[1;31mError: \033[0m");
    
    va_list varargs;
    va_start(varargs, msg);
    vfprintf(stderr, msg, varargs);
    va_end(varargs);
    
    fprintf(stderr, "\n");
    fflush(stderr);
}

void tc_assert(int condition, const char* msg, ...)
{
    if (condition)
        return;

    fprintf(stderr, "\033[1;31mAssertion failed: \033[0m");

    va_list varargs;
    va_start(varargs, msg);
    vfprintf(stderr, msg, varargs);
    va_end(varargs);

    fprintf(stderr, "\n");
    fflush(stderr);

    tc_die(tc_exit_code_assertion);
}

void tc_die(int status)
{
    exit(status);
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

void tc_debug_bset(__isl_keep isl_basic_set* bset, const char* msg, ...)
{
    if (!tc_debug_flag)
        return;

    va_list varargs;
    va_start(varargs, msg);
    vfprintf(stderr, msg, varargs);
    va_end(varargs);

    char* str = isl_basic_set_to_str(bset);

    fprintf(stderr, " := %s;\n\n", str);
    fflush(stderr);

    free(str);
}

void tc_debug_bset_latex(__isl_keep isl_basic_set* bset, const char* msg, ...)
{
    if (!tc_debug_flag)
        return;

    va_list varargs;
    va_start(varargs, msg);
    vfprintf(stderr, msg, varargs);
    va_end(varargs);

    isl_ctx* ctx = isl_basic_set_get_ctx(bset);

    isl_printer* printer = isl_printer_to_str(ctx);

    printer = isl_printer_set_output_format(printer, ISL_FORMAT_LATEX);

    printer = isl_printer_flush(printer);

    printer = isl_printer_print_basic_set(printer, bset);

    char* str = isl_printer_get_str(printer);

    fprintf(stderr, " := %s;\n\n", str);
    fflush(stderr);

    isl_printer_free(printer);

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

void tc_debug_bmap(__isl_keep isl_basic_map* bmap, const char* msg, ...)
{
    if (!tc_debug_flag)
        return;

    va_list varargs;
    va_start(varargs, msg);
    vfprintf(stderr, msg, varargs);
    va_end(varargs);

    char* str = isl_basic_map_to_str(bmap);

    fprintf(stderr, " := %s;\n\n", str);
    fflush(stderr);

    free(str);
}

void tc_debug_bmap_latex(__isl_keep isl_basic_map* bmap, const char* msg, ...)
{
    if (!tc_debug_flag)
        return;

    va_list varargs;
    va_start(varargs, msg);
    vfprintf(stderr, msg, varargs);
    va_end(varargs);

    isl_ctx* ctx = isl_basic_map_get_ctx(bmap);

    isl_printer* printer = isl_printer_to_str(ctx);

    printer = isl_printer_set_output_format(printer, ISL_FORMAT_LATEX);

    printer = isl_printer_flush(printer);

    printer = isl_printer_print_basic_map(printer, bmap);

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

void tc_debug_id(__isl_keep isl_id* id, const char* msg, ...)
{
    if (!tc_debug_flag)
        return;

    va_list varargs;
    va_start(varargs, msg);
    vfprintf(stderr, msg, varargs);
    va_end(varargs);

    char* str = isl_id_to_str(id);

    fprintf(stderr, " := \"%s\";\n\n", str);
    fflush(stderr);

    free(str);
}

void tc_debug_id_list(__isl_keep isl_id_list* list, const char* msg, ...)
{
    if (!tc_debug_flag)
        return;

    va_list varargs;
    va_start(varargs, msg);
    vfprintf(stderr, msg, varargs);
    va_end(varargs);

    fprintf(stderr, " := [");

    for (int i = 0; i < isl_id_list_n_id(list); ++i)
    {
        isl_id* id = isl_id_list_get_id(list, i);

        char* str = isl_id_to_str(id);

        fprintf(stderr, "%s", str);

        if (i + 1 < isl_id_list_n_id(list))
        {
            fprintf(stderr, ", ");
        }

        free(str);
        isl_id_free(id);
    }

    fprintf(stderr, "];\n\n");

    fflush(stderr);
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

void tc_debug_aff(__isl_keep isl_aff* aff, const char* msg, ...)
{
    if (!tc_debug_flag)
        return;

    va_list varargs;
    va_start(varargs, msg);
    vfprintf(stderr, msg, varargs);
    va_end(varargs);

    char* str = isl_aff_to_str(aff);

    fprintf(stderr, " := %s;\n\n", str);
    fflush(stderr);

    free(str);
}

void tc_debug_multi_aff(__isl_keep isl_multi_aff* aff, const char* msg, ...)
{
    if (!tc_debug_flag)
        return;

    va_list varargs;
    va_start(varargs, msg);
    vfprintf(stderr, msg, varargs);
    va_end(varargs);

    char* str = isl_multi_aff_to_str(aff);

    fprintf(stderr, " := %s;\n\n", str);
    fflush(stderr);

    free(str);
}

void tc_debug_pw_aff(__isl_keep isl_pw_aff* aff, const char* msg, ...)
{
    if (!tc_debug_flag)
        return;

    va_list varargs;
    va_start(varargs, msg);
    vfprintf(stderr, msg, varargs);
    va_end(varargs);

    char* str = isl_pw_aff_to_str(aff);

    fprintf(stderr, " := %s;\n\n", str);
    fflush(stderr);

    free(str);
}

void tc_debug_pw_multi_aff(__isl_keep isl_pw_multi_aff* maff, const char* msg, ...)
{
    if (!tc_debug_flag)
        return;

    va_list varargs;
    va_start(varargs, msg);
    vfprintf(stderr, msg, varargs);
    va_end(varargs);

    char* str = isl_pw_multi_aff_to_str(maff);

    fprintf(stderr, " := %s;\n\n", str);
    fflush(stderr);

    free(str);
}

void tc_debug_union_pw_multi_aff(__isl_keep isl_union_pw_multi_aff* umaff, const char* msg, ...)
{
    if (!tc_debug_flag)
        return;

    va_list varargs;
    va_start(varargs, msg);
    vfprintf(stderr, msg, varargs);
    va_end(varargs);

    char* str = isl_union_pw_multi_aff_to_str(umaff);

    fprintf(stderr, " := %s;\n\n", str);
    fflush(stderr);

    free(str);
}

void tc_debug_multi_union_pw_aff(__isl_keep isl_multi_union_pw_aff* muaff, const char* msg, ...)
{
    if (!tc_debug_flag)
        return;

    va_list varargs;
    va_start(varargs, msg);
    vfprintf(stderr, msg, varargs);
    va_end(varargs);

    char* str = isl_multi_union_pw_aff_to_str(muaff);

    fprintf(stderr, " := %s;\n\n", str);
    fflush(stderr);

    free(str);
}

void tc_debug_schedule(__isl_keep isl_schedule* schedule, const char* msg, ...)
{
    if (!tc_debug_flag)
        return;

    va_list varargs;
    va_start(varargs, msg);
    vfprintf(stderr, msg, varargs);
    va_end(varargs);

    char* str = isl_schedule_to_str(schedule);

    fprintf(stderr, " := %s;\n\n", str);
    fflush(stderr);
    
    free(str);
}

void tc_debug_ast_expr(__isl_keep isl_ast_expr* expr, const char* msg, ...)
{
    if (!tc_debug_flag)
        return;

    va_list varargs;
    va_start(varargs, msg);
    vfprintf(stderr, msg, varargs);
    va_end(varargs);

    char* str = isl_ast_expr_to_str(expr);

    fprintf(stderr, " := %s;\n\n", str);
    fflush(stderr);
    
    free(str);
}
