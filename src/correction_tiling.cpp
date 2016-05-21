#include "correction_tiling.h"
#include "lex_scheduling.h"
#include "sfs_scheduling.h"
#include "free_scheduling.h"
#include "dynamic_free_scheduling.h"
#include "tiling.h"
#include "utility.h"
#include "scop.h"
#include "options.h"
#include "tile_statistics.h"
#include "debug.h"
#include "transitive_closure.h"

#include <stdio.h>

#include <string>
#include <vector>
#include <map>

void tc_algorithm_correction_tiling(struct tc_scop* scop, struct tc_options* options)
{
    isl_ctx* ctx = scop->ctx;
    
    isl_union_set* LD = isl_union_set_copy(scop->domain);
    tc_debug_uset(LD, "LD");
    
    isl_union_map* S = isl_union_map_copy(scop->schedule);
    tc_debug_umap(S, "S");
    
    isl_union_map* R = isl_union_map_copy(scop->relation);
    tc_debug_umap(R, "R");
    
    isl_set* LD_normalized = tc_normalize_union_set(LD, S);
    
    tc_debug_set_card(LD_normalized, "LD");
    
    isl_basic_set* sample = isl_set_sample(isl_set_copy(LD_normalized));
    
    isl_space* space = isl_basic_set_get_space(sample);
            
    isl_id_list* I = tc_ids_sequence(ctx, "i", isl_space_dim(space, isl_dim_set));
    isl_id_list* II = tc_ids_sequence(ctx, "ii", isl_space_dim(space, isl_dim_set));
            
    std::map<std::string, std::vector<int> > blocks = tc_options_blocks(options);
    
    std::vector<std::vector<std::string> > groups = tc_options_groups(options);
    
    //R = tc_remove_loop_independent_dependences(R, S, groups);
    
    isl_set* tile;
    isl_set* ii_set;
    
    tc_tile_loop_nest(LD, S, II, I, &tile, &ii_set, blocks, groups);
    
    tc_debug_set(tile, "TILE");
    tc_debug_set(ii_set, "II_SET");
    
    isl_map* R_normalized = tc_normalize_union_map(R, S);
    
    int exact;
    isl_map* R_plus_normalized = tc_transitive_closure(isl_map_copy(R_normalized), S, &exact);
        
    tc_debug_map(R_plus_normalized, "R^+ (exact=%d)", exact);    
    
    isl_set* tile_lt = tc_tile_lt_set(tile, ii_set, II);
    isl_set* tile_gt = tc_tile_gt_set(tile, ii_set, II);
    
    // TILE_ITR = TILE - R+(TILE_GT)
    isl_set* tile_itr = isl_set_subtract(isl_set_copy(tile), isl_set_apply(isl_set_copy(tile_gt), isl_map_copy(R_plus_normalized)));
    
    tc_debug_set(tile_itr, "TILE_ITR");
    
    // TVLD_LT = (R+(TILE_ITR) * TILE_LT) - R+(TILE_GT)
    isl_set* tvld_lt = isl_set_apply(isl_set_copy(tile_itr), isl_map_copy(R_plus_normalized));
    tvld_lt = isl_set_intersect(tvld_lt, tile_lt);
    tvld_lt = isl_set_subtract(tvld_lt, isl_set_apply(tile_gt, isl_map_copy(R_plus_normalized)));
    
    tc_debug_set(tvld_lt, "TVLD_LT");
    
    // TILE_VLD = TILE_ITR + TVLD_LT
    isl_set* tile_vld = isl_set_union(tile_itr, tvld_lt);    
    tile_vld = isl_set_coalesce(tile_vld);
    
    tc_debug_set(tile_vld, "TILE_VLD");

    tc_debug_bool(isl_set_is_equal(tile, tile_vld), "TILE = TILE_VLD");
    
    if (tc_options_is_report(options))
    {
        isl_set* bounds = tc_options_get_report_bounds(options, ctx);
        
        struct tc_tile_statistics* stats = tc_compute_tile_statistics(tile_vld, ii_set, II, bounds, LD, S, scop, blocks);
        
        tc_tile_statistics_print(options->output, stats);
        
        tc_tile_statistics_free(stats);
        isl_set_free(bounds);
    }
    
    isl_map* Rtile = tc_Rtile_map(II, tile_vld, R_normalized);
    
    tc_debug_map(Rtile, "R_TILE");
    
    enum tc_scheduling_enum scheduling = tc_options_scheduling(options);
    
    if (tc_scheduling_enum_lex == scheduling)
    {
        tc_scheduling_lex(scop, options, LD, S, R, ii_set, tile_vld, Rtile, II);
    }
    else if (tc_scheduling_enum_sfs_tile == scheduling)
    {
        tc_scheduling_sfs_tiles(scop, options, LD, S, R, ii_set, tile_vld, Rtile, II);
    }
    else if (tc_scheduling_enum_free_rk == scheduling)
    {
        tc_scheduling_free_schedule_rk(scop, options, LD, S, R, ii_set, tile_vld, Rtile, II);
    }
    else if (tc_scheduling_enum_free_karl == scheduling)
    {
        tc_scheduling_free_schedule_karl(scop, options, LD, S, R, ii_set, tile_vld, Rtile, II);
    }
    else if (tc_scheduling_enum_free_dynamic == scheduling)
    {
        tc_scheduling_dynamic_free_schedule(scop, options, LD, S, R, ii_set, tile_vld, Rtile, II);
    }
    
    isl_id_list_free(I);
    isl_set_free(tile);
    isl_map_free(R_normalized);
    isl_map_free(R_plus_normalized);
    isl_basic_set_free(sample);
    isl_space_free(space);
    isl_set_free(LD_normalized);
}
