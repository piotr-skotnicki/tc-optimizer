#include "sfs_tiling.h"
#include "utility.h"
#include "tiling.h"
#include "slicing.h"
#include "options.h"

#include <isl/ctx.h>
#include <isl/space.h>
#include <isl/id.h>
#include <isl/set.h>
#include <isl/map.h>
#include <isl/union_set.h>
#include <isl/union_map.h>
#include <isl/val.h>

#include <string.h>

#include <vector>
#include <map>
#include <string>

__isl_give isl_printer* tc_algorithm_sfs_tiling_print_for(__isl_take isl_printer* printer, __isl_take isl_ast_print_options* options, __isl_keep isl_ast_node* node, void* user)
{    
    isl_ast_expr* expr = isl_ast_node_for_get_iterator(node);
    
    isl_id* id = isl_ast_expr_get_id(expr);
    
    isl_id* parallel_id = (isl_id*)user;
    
    if (NULL != parallel_id && 0 == strcmp(isl_id_get_name(id), isl_id_get_name(parallel_id)))
    {
        printer = isl_printer_start_line(printer);
        printer = isl_printer_print_str(printer, "#pragma omp parallel for");
        printer = isl_printer_end_line(printer);
    }    

    printer = isl_ast_node_for_print(node, printer, options);
    
    isl_id_free(id);
    isl_ast_expr_free(expr);
    
    return printer;
}

void tc_algorithm_sfs_tiling(int argc, char* argv[], struct tc_scop* scop)
{
    isl_ctx* ctx = scop->ctx;
    
    isl_union_set* LD = scop->domain;
    
    isl_union_map* S = scop->schedule;
    
    isl_union_map* R = scop->relation;
    
    isl_basic_set* sample = isl_set_sample(tc_normalize_union_set(LD, S));
    
    isl_space* space = isl_basic_set_get_space(sample);
            
    isl_id_list* I = tc_ids_sequence(ctx, "i", isl_space_dim(space, isl_dim_set));
    isl_id_list* II = tc_ids_sequence(ctx, "ii", isl_space_dim(space, isl_dim_set));
                
    std::map<std::string, std::vector<int> > blocks = tc_options_blocks(argc, argv);
        
    isl_set* tile;
    isl_set* ii_set;
    
    tc_tile_loop_nest(scop->domain, scop->schedule, blocks, II, I, &tile, &ii_set);
        
    isl_map* R_normalized = tc_normalize_union_map(R, S);
        
    isl_map* Rtile = tc_Rtile_map(tile, R_normalized);
    
    isl_map* Rusc = tc_Rusc_map(Rtile);
    
    isl_set* Suds = tc_uds_set(Rtile);
    
    isl_set* Srepr = isl_set_subtract(isl_set_copy(Suds), isl_map_range(isl_map_copy(Rusc)));
    
    Srepr = tc_parameterize_all(Srepr, II);
        
    int exact;
    isl_map* Rtile_plus = isl_map_transitive_closure(isl_map_copy(Rtile), &exact);
    isl_map* Rtile_star = isl_map_union(isl_map_copy(Rtile_plus), tc_make_identity(isl_map_copy(Rtile)));
                
    isl_set* Sslice = NULL;
    
    if (isl_map_is_empty(Rusc))
    {
        // S_slice = R*(e)
        Sslice = isl_set_apply(isl_set_copy(Srepr), isl_map_copy(Rtile_star));
    }
    else
    {
        isl_map* Rusc_plus = isl_map_transitive_closure(isl_map_copy(Rusc), &exact);
        isl_map* Rusc_star = isl_map_union(isl_map_copy(Rusc_plus), tc_make_identity(isl_map_copy(Rusc)));
    
        // S_slice = R*((R_USC*)(e))
        Sslice = isl_set_apply(isl_set_apply(isl_set_copy(Srepr), isl_map_copy(Rusc_star)), isl_map_copy(Rtile_star));
        
        isl_map_free(Rusc_plus);
        isl_map_free(Rusc_star);
    }
        
    Sslice = tc_lift_up_set_params(Sslice, II);
    
    tile = tc_lift_up_set_params(tile, II);
    
    Sslice = isl_set_add_dims(Sslice, isl_dim_set, isl_id_list_n_id(II));
    tile = isl_set_insert_dims(tile, isl_dim_set, 0, isl_id_list_n_id(II));
    
    isl_set* tile_ext = isl_set_intersect(isl_set_copy(Sslice), isl_set_copy(tile));
    
    tile_ext = isl_set_coalesce(tile_ext);
            
    isl_union_map* S_ext = NULL;
    isl_map_list* S_maps = tc_collect_maps(S);
    
    for (int i = 0; i < isl_map_list_n_map(S_maps); ++i)
    {
        isl_map* map = isl_map_list_get_map(S_maps, i);
        
        map = isl_map_insert_dims(map, isl_dim_out, 0, 2*isl_id_list_n_id(II));
    
        if (NULL == S_ext)
        {
            S_ext = isl_union_map_from_map(map);
        }
        else
        {
            S_ext = isl_union_map_add_map(S_ext, map);
        }
    }
    
    isl_map_list_free(S_maps);
    
    const int INDENT_SIZE = 2;
        
    isl_ast_build* ast_build = isl_ast_build_from_context(isl_set_copy(scop->pet->context));
    
    isl_id_list* iterators = tc_ids_sequence(ctx, "ii", isl_set_n_dim(tile_ext));
        
    isl_id* parallel_iterator = NULL;
    for (int i = 0; i < isl_space_dim(space, isl_dim_set); ++i)
    {
        isl_val* val = isl_set_plain_get_val_if_fixed(tile_ext, isl_dim_set, i);
        isl_bool is_nan = isl_val_is_nan(val);
        isl_val_free(val);
        if (is_nan)
        {
            parallel_iterator = isl_id_list_get_id(iterators, i);
            break;
        }
    }
    
    ast_build = isl_ast_build_set_iterators(ast_build, iterators);

    isl_union_map* S_ext_prim = isl_union_map_intersect_range(S_ext, isl_union_set_from_set(isl_set_copy(tile_ext)));
    
    isl_ast_node* ast_tile = isl_ast_build_ast_from_schedule(ast_build, S_ext_prim);
        
    isl_printer* printer = isl_printer_to_str(ctx);    
    
    printer = isl_printer_set_output_format(printer, ISL_FORMAT_C);
        
    isl_ast_print_options* ast_options = isl_ast_print_options_alloc(ctx);
    
    ast_options = isl_ast_print_options_set_print_for(ast_options, &tc_algorithm_sfs_tiling_print_for, parallel_iterator);
        
    printer = isl_ast_node_print_macros(ast_tile, printer);
    
    printer = tc_print_statements_macros(scop, printer, ast_build);
    
    printer = isl_printer_print_str(printer, "#pragma scop\n");
    printer = isl_printer_set_indent(printer, INDENT_SIZE);
    printer = isl_ast_node_print(ast_tile, printer, ast_options);
    printer = isl_printer_print_str(printer, "#pragma endscop\n");
    
    char* code = isl_printer_get_str(printer);
    
    printf("%s", code);
    
    free(code);
    
    isl_id_free(parallel_iterator);
    isl_id_list_free(I);
    isl_id_list_free(II);
    isl_printer_free(printer);
    isl_ast_node_free(ast_tile);
    isl_ast_build_free(ast_build);
    isl_map_free(R_normalized);    
    isl_map_free(Rtile);
    isl_map_free(Rtile_plus);
    isl_map_free(Rtile_star);
    isl_map_free(Rusc);
    isl_set_free(Suds);
    isl_set_free(Srepr);
    isl_set_free(tile);
    isl_set_free(tile_ext);
    isl_set_free(ii_set);
    isl_basic_set_free(sample);
    isl_set_free(Sslice);
    isl_space_free(space);
}
