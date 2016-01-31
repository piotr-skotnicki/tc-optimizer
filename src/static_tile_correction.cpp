#include "static_tile_correction.h"
#include "tiling.h"
#include "utility.h"
#include "scop.h"
#include "options.h"

#include <string>
#include <vector>
#include <map>

void tc_algorithm_static_tile_correction(int argc, char* argv[], struct tc_scop* scop)
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
    
    int is_exact;
    isl_union_map* R_plus = isl_union_map_transitive_closure(isl_union_map_copy(R), &is_exact);
    isl_map* R_plus_normalized = tc_normalize_union_map(R_plus, S);
    
    isl_set* tile_lt = tc_tile_lt_set(tile, ii_set, II);
    isl_set* tile_gt = tc_tile_gt_set(tile, ii_set, II);
    
    // TILE_ITR = TILE - R+(TILE_GT)
    isl_set* tile_itr = isl_set_subtract(isl_set_copy(tile), isl_set_apply(isl_set_copy(tile_gt), isl_map_copy(R_plus_normalized)));
    
    // TVLD_LT = (R+(TILE_ITR) * TILE_LT) - R+(TILE_GT)
    isl_set* tvld_lt = isl_set_apply(isl_set_copy(tile_itr), isl_map_copy(R_plus_normalized));
    tvld_lt = isl_set_intersect(tvld_lt, tile_lt);
    tvld_lt = isl_set_subtract(tvld_lt, isl_set_apply(tile_gt, isl_map_copy(R_plus_normalized)));
    
    // TILE_VLD = TILE_ITR + TVLD_LT
    isl_set* tile_vld = isl_set_union(tile_itr, tvld_lt);
    
    tile_vld = isl_set_coalesce(tile_vld);
        
    isl_set* tile_vld_ext = tc_lift_up_set_params(tile_vld, II);
    
    tile_vld_ext = isl_set_coalesce(tile_vld_ext);
    
    isl_union_map* S_ext = NULL;
    isl_map_list* S_maps = tc_collect_maps(S);
    
    for (int i = 0; i < isl_map_list_n_map(S_maps); ++i)
    {
        isl_map* map = isl_map_list_get_map(S_maps, i);
        
        map = isl_map_insert_dims(map, isl_dim_out, 0, isl_id_list_n_id(II));
    
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

    isl_union_map* S_ext_prim = isl_union_map_intersect_range(S_ext, isl_union_set_from_set(isl_set_copy(tile_vld_ext)));
    
    isl_ast_node* ast_tile = isl_ast_build_ast_from_schedule(ast_build, S_ext_prim);
        
    isl_printer* printer = isl_printer_to_str(ctx);    
    
    printer = isl_printer_set_output_format(printer, ISL_FORMAT_C);
        
    isl_ast_print_options* ast_options = isl_ast_print_options_alloc(ctx);
        
    printer = isl_ast_node_print_macros(ast_tile, printer);
    
    printer = tc_print_statements_macros(scop, printer, ast_build);
    
    printer = isl_printer_print_str(printer, "#pragma scop\n");
    printer = isl_printer_set_indent(printer, INDENT_SIZE);
    printer = isl_ast_node_print(ast_tile, printer, ast_options);
    printer = isl_printer_print_str(printer, "#pragma endscop\n");
    
    char* code = isl_printer_get_str(printer);
    
    printf("%s", code);
    
    free(code);
    
    isl_id_list_free(I);
    isl_id_list_free(II);
    isl_printer_free(printer);
    isl_ast_node_free(ast_tile);
    isl_ast_build_free(ast_build);
    isl_set_free(tile);
    isl_set_free(tile_vld_ext);
    isl_set_free(ii_set);
    isl_map_free(R_plus_normalized);
    isl_union_map_free(R_plus);
    isl_basic_set_free(sample);
    isl_space_free(space);
}
