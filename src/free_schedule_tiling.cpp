#include "free_schedule_tiling.h"
#include "tiling.h"
#include "utility.h"
#include "options.h"

#include <isl/ctx.h>
#include <isl/space.h>
#include <isl/id.h>
#include <isl/set.h>
#include <isl/map.h>
#include <isl/union_set.h>
#include <isl/union_map.h>

#include <vector>
#include <map>
#include <string>

void tc_algorithm_free_schedule_tiling(int argc, char* argv[], struct tc_scop* scop)
{
    isl_ctx* ctx = scop->ctx;
    
    isl_union_set* LD = scop->domain;
    
    isl_union_map* S = scop->schedule;
    
    isl_union_map* R = scop->relation;
    
    isl_basic_set* sample = isl_set_sample(tc_normalize_union_set(LD, S));
    
    isl_space* space = isl_basic_set_get_space(sample);
            
    isl_id_list* I = tc_ids_sequence(ctx, "i", isl_space_dim(space, isl_dim_set));
    isl_id_list* II = tc_ids_sequence(ctx, "ii", isl_space_dim(space, isl_dim_set));
    
    isl_id_list* IIprim = tc_ids_prim(II);
    isl_id_list* IIbis = tc_ids_bis(II);
        
    std::map<std::string, std::vector<int> > blocks = tc_options_blocks(argc, argv);
        
    isl_set* tile;
    isl_set* ii_set;
    
    tc_tile_loop_nest(scop->domain, scop->schedule, blocks, II, I, &tile, &ii_set);
        
    isl_map* R_normalized = tc_normalize_union_map(R, S);
        
    isl_map* Rtile = tc_Rtile_map(tile, R_normalized);
    
    isl_map* backward = tc_make_map(ctx, NULL, IIprim, IIbis, tc_tuples_gt(IIprim, IIbis).c_str());
    backward = isl_map_intersect(isl_map_copy(Rtile), backward);
    int is_backward_relation = !isl_map_is_empty(backward);
    isl_map_free(backward);
    
    if (is_backward_relation)
    {
        fprintf(stderr, "Error: Backward relation detected");
        return;
    }
    
    const int INDENT_SIZE = 2;
        
    isl_ast_build* ast_build = isl_ast_build_from_context(isl_set_copy(scop->pet->context));

    isl_union_map* S_prim = isl_union_map_intersect_range(isl_union_map_copy(S), isl_union_set_from_set(isl_set_copy(tile)));
    
    isl_ast_node* ast_tile = isl_ast_build_ast_from_schedule(ast_build, S_prim);
    
    char buff[1024];
    
    isl_printer* printer = isl_printer_to_str(ctx);    
    
    printer = isl_printer_set_output_format(printer, ISL_FORMAT_C);
        
    isl_ast_print_options* ast_options = isl_ast_print_options_alloc(ctx);
    
    printer = isl_printer_print_str(printer,
        "#include <omp.h>\n"
        "#include <isl/ctx.h>\n"
        "#include <isl/space.h>\n"
        "#include <isl/point.h>\n"
        "#include <isl/val.h>\n"
        "#include <isl/set.h>\n"
        "#include <isl/map.h>\n"
    );
    
    printer = isl_ast_node_print_macros(ast_tile, printer);
    
    printer = tc_print_statements_macros(scop, printer, ast_build);
    
    printer = isl_printer_print_str(printer, "\n");
    printer = isl_printer_print_str(printer, "int create_task(isl_point* point, void* user) {\n");
    
    for (int i = 0; i < isl_id_list_n_id(II); ++i)
    {
        isl_id* II_i = isl_id_list_get_id(II, i);
        
        sprintf(buff, "  isl_val* p%d = isl_point_get_coordinate_val(point, isl_dim_set, %d);\n", i, i);
        printer = isl_printer_print_str(printer, buff);
        
        sprintf(buff, "  int %s = isl_val_get_num_si(p%d);\n", isl_id_get_name(II_i), i);
        printer = isl_printer_print_str(printer, buff);
        
        sprintf(buff, "  isl_val_free(p%d);\n", i);
        printer = isl_printer_print_str(printer, buff);
        
        isl_id_free(II_i);
    }
    
    printer = isl_printer_print_str(printer,
        "  isl_point_free(point);\n"
        "  #pragma omp task\n"
        "  {\n"
    );
    
    printer = isl_printer_set_indent(printer, INDENT_SIZE*2);
    printer = isl_ast_node_print(ast_tile, printer, ast_options);
    
    printer = isl_printer_print_str(printer,
        "  }\n"
        "  return 0;\n"
        "}\n"
        "\n"
    );
    
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
        
        sprintf(buff, "  rtile = isl_map_fix_si(rtile, isl_dim_param, %d, %s);\n", i, isl_id_get_name(Rtile_param));
        printer = isl_printer_print_str(printer, buff);
        
        sprintf(buff, "  ii_set = isl_set_fix_si(ii_set, isl_dim_param, %d, %s);\n", i, isl_id_get_name(ii_set_param));
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
    
    char* code = isl_printer_get_str(printer);
    
    printf("%s", code);
    
    free(code);
    
    isl_printer_free(printer);
    isl_ast_node_free(ast_tile);
    isl_ast_build_free(ast_build);
    
    isl_id_list_free(I);
    isl_id_list_free(II);
    isl_id_list_free(IIprim);
    isl_id_list_free(IIbis);
    isl_set_free(tile);
    isl_set_free(ii_set);
    isl_map_free(Rtile);
    isl_map_free(R_normalized);
    isl_basic_set_free(sample);
    isl_space_free(space);
}
