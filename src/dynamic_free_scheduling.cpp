#include "dynamic_free_scheduling.h"
#include "ast.h"
#include "codegen.h"
#include "config.h"
#include "utility.h"
#include "scop.h"

#include <isl/ctx.h>
#include <isl/id.h>
#include <isl/set.h>
#include <isl/map.h>
#include <isl/union_map.h>
#include <isl/union_set.h>
#include <isl/ast_build.h>
#include <isl/printer.h>

#include <stdio.h>
#include <stddef.h>

void tc_scheduling_dynamic_free_schedule(struct tc_scop* scop, struct tc_options* options, __isl_take isl_union_set* LD, __isl_take isl_union_map* S, __isl_take isl_union_map* R, __isl_take isl_set* ii_set, __isl_take isl_set* tile, __isl_take isl_map* Rtile, __isl_take isl_id_list* II, __isl_take isl_id_list* I)
{
    isl_ctx* ctx = isl_union_set_get_ctx(LD);
    
    isl_ast_build* ast_build = isl_ast_build_from_context(isl_set_copy(scop->pet->context));
    
    ast_build = isl_ast_build_set_at_each_domain(ast_build, &tc_ast_visitor_at_each_domain, scop);

    isl_union_map* S_prim = isl_union_map_intersect_range(isl_union_map_copy(S), isl_union_set_from_set(isl_set_copy(tile)));
        
    char buff[1024];
    
    isl_printer* printer = isl_printer_to_str(ctx);    
    
    printer = isl_printer_set_output_format(printer, ISL_FORMAT_C);
        
    isl_ast_print_options* ast_options = isl_ast_print_options_alloc(ctx);
    
    if (!tc_options_is_set(options, NULL, "--use-macros"))
    {
        ast_options = isl_ast_print_options_set_print_user(ast_options, &tc_codegen_print_user, NULL);
    }
    
    printer = isl_printer_print_str(printer,
        "#include <isl/ctx.h>\n"
        "#include <isl/space.h>\n"
        "#include <isl/point.h>\n"
        "#include <isl/val.h>\n"
        "#include <isl/set.h>\n"
        "#include <isl/map.h>\n"
    );
    
    isl_ast_node* ast_tile = isl_ast_build_ast_from_schedule(ast_build, S_prim);
    
    printer = isl_ast_node_print_macros(ast_tile, printer);
    
    if (tc_options_is_set(options, NULL, "--use-macros"))
    {
        printer = tc_codegen_print_statements_macros(scop, printer);
    }
    
    printer = isl_printer_print_str(printer, "\n");
    printer = isl_printer_print_str(printer, "int create_task(isl_point* point, void* user) {\n");
    
    for (int i = 0; i < isl_id_list_n_id(II); ++i)
    {
        isl_id* II_i = isl_id_list_get_id(II, i);
        
        snprintf(buff, sizeof(buff), "  isl_val* p%d = isl_point_get_coordinate_val(point, isl_dim_set, %d);\n", i, i);
        printer = isl_printer_print_str(printer, buff);
        
        snprintf(buff, sizeof(buff), "  int %s = isl_val_get_num_si(p%d);\n", isl_id_get_name(II_i), i);
        printer = isl_printer_print_str(printer, buff);
        
        snprintf(buff, sizeof(buff), "  isl_val_free(p%d);\n", i);
        printer = isl_printer_print_str(printer, buff);
        
        isl_id_free(II_i);
    }
    
    printer = isl_printer_print_str(printer,
        "  isl_point_free(point);\n"
        "  #pragma omp task\n"
        "  {\n"
    );
    
    printer = isl_printer_set_indent(printer, TC_CONF_INDENT_SIZE*2);
    printer = isl_ast_node_print(ast_tile, printer, ast_options);
    
    printer = isl_printer_print_str(printer,
        "  }\n"
        "  return 0;\n"
        "}\n"
        "\n"
    );
    
    printer = tc_codegen_print_prologue(scop, options, printer);
    printer = isl_printer_print_str(printer,
        "#pragma scop\n"
        "  isl_ctx* ctx = isl_ctx_alloc();\n"
        "  isl_map* rtile = isl_map_read_from_str(ctx, \""
    );    
    printer = isl_printer_set_output_format(printer, ISL_FORMAT_ISL);
    printer = isl_printer_print_map(printer, Rtile);
    printer = isl_printer_set_output_format(printer, ISL_FORMAT_C);
    printer = isl_printer_print_str(printer,
        "\");\n"
        "  isl_set* ii_set = isl_set_read_from_str(ctx, \""
    );
    printer = isl_printer_set_output_format(printer, ISL_FORMAT_ISL);
    printer = isl_printer_print_set(printer, ii_set);
    printer = isl_printer_set_output_format(printer, ISL_FORMAT_C);
    printer = isl_printer_print_str(printer,
        "\");\n"
    );
    
    isl_id_list* Rtile_params = tc_get_map_params_names(Rtile);
    isl_id_list* ii_set_params = tc_get_set_params_names(ii_set);
    for (int i = 0; i < isl_id_list_n_id(Rtile_params); ++i)
    {
        isl_id* Rtile_param = isl_id_list_get_id(Rtile_params, i);
        isl_id* ii_set_param = isl_id_list_get_id(ii_set_params, i);
        
        snprintf(buff, sizeof(buff), "  rtile = isl_map_fix_si(rtile, isl_dim_param, %d, %s);\n", i, isl_id_get_name(Rtile_param));
        printer = isl_printer_print_str(printer, buff);
        
        snprintf(buff, sizeof(buff), "  ii_set = isl_set_fix_si(ii_set, isl_dim_param, %d, %s);\n", i, isl_id_get_name(ii_set_param));
        printer = isl_printer_print_str(printer, buff);
        
        isl_id_free(Rtile_param);
        isl_id_free(ii_set_param);
    }
    isl_id_list_free(Rtile_params);
    isl_id_list_free(ii_set_params);
                
    printer = isl_printer_print_str(printer, 
        "  isl_set* uds = isl_set_subtract(isl_map_domain(isl_map_copy(rtile)), isl_map_range(isl_map_copy(rtile)));\n"
        "  \n"
        "  #pragma omp parallel\n"
        "  #pragma omp single\n"
        "  {\n"
        "    while (!isl_set_is_empty(uds)) {\n"
        "      isl_set_foreach_point(uds, &create_task, NULL);\n"
        "      #pragma omp taskwait\n"
        "\n"
        "      ii_set = isl_set_subtract(ii_set, uds);\n"
        "      rtile = isl_map_intersect_domain(rtile, isl_set_copy(ii_set));\n"
        "      uds = isl_set_subtract(isl_map_domain(isl_map_copy(rtile)), isl_map_range(isl_map_copy(rtile)));\n"
        "    }\n"
        "\n"
        "    isl_set_foreach_point(ii_set, &create_task, NULL);\n"
        "    #pragma omp taskwait\n"
        "  }\n"
        "\n"
        "  isl_set_free(uds);\n"
        "  isl_set_free(ii_set);\n"
        "  isl_map_free(rtile);\n"
        "  isl_ctx_free(ctx);\n"
        "#pragma endscop\n"
    );
    printer = tc_codegen_print_epilogue(scop, options, printer);
    
    char* code = isl_printer_get_str(printer);
    
    printf("%s", code);
    
    free(code);
    
    isl_printer_free(printer);
    isl_ast_node_free(ast_tile);
    isl_ast_build_free(ast_build);
    
    isl_id_list_free(I);
    isl_id_list_free(II);
    isl_set_free(tile);
    isl_set_free(ii_set);
    isl_map_free(Rtile);
    isl_union_set_free(LD);
    isl_union_map_free(S);
    isl_union_map_free(R);
}
