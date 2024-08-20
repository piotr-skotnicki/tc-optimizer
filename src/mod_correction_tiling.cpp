#include "mod_correction_tiling.h"
#include "scheduling.h"
#include "tiling.h"
#include "utility.h"
#include "scop.h"
#include "options.h"
#include "tile_statistics.h"
#include "debug.h"
#include "transitive_closure.h"

#include <isl/ctx.h>
#include <isl/space.h>
#include <isl/id.h>
#include <isl/set.h>
#include <isl/map.h>
#include <isl/union_set.h>
#include <isl/union_map.h>

#include <barvinok/isl.h>

#include <stdio.h>
#include <stddef.h>

#include <string>
#include <vector>
#include <map>

void tc_algorithm_mod_correction_tiling(struct tc_scop* scop, struct tc_options* options)
{
    isl_ctx* ctx = scop->ctx;
    
    isl_union_set* LD = isl_union_set_copy(scop->domain);
    tc_debug_uset(LD, "LD");
    
    isl_union_map* S = isl_union_map_copy(scop->schedule);
    tc_debug_umap(S, "S");
    
    isl_union_map* R = isl_union_map_copy(scop->relation);
    tc_debug_umap(R, "R");
    
    isl_union_map* RA = scop->reads;
    tc_debug_umap(RA, "RA");
    
    isl_union_map* WA = scop->writes;
    tc_debug_umap(WA, "WA");
    
    isl_set* LD_normalized = tc_normalize_union_set(LD, S);
    
    tc_debug_set_card(LD_normalized, "LD");
    
    isl_basic_set* sample = isl_set_sample(isl_set_copy(LD_normalized));
    
    isl_space* space = isl_basic_set_get_space(sample);
            
    isl_id_list* I = tc_ids_sequence(ctx, "i", isl_space_dim(space, isl_dim_set));
    isl_id_list* II = tc_ids_sequence(ctx, "ii", isl_space_dim(space, isl_dim_set));
            
    std::map<std::string, std::vector<int> > blocks = tc_options_blocks(options);
    
    std::vector<std::vector<std::string> > groups = tc_options_groups(options);
        
    isl_set* tile;
    isl_set* ii_set;
    
    tc_tile_loop_nest(LD, S, II, I, &tile, &ii_set, blocks, groups);
    
    tc_debug_set(tile, "TILE");
    tc_debug_set(ii_set, "II_SET");
        
    isl_map* R_normalized = tc_normalize_union_map(R, S);
    tc_debug_map(R_normalized, "R_norm");
    
    isl_bool exact = isl_bool_false;
    isl_map* R_plus_normalized = tc_transitive_closure(isl_map_copy(R_normalized), S, &exact);
        
    tc_debug_map(R_plus_normalized, "R^+ (exact=%d)", exact);

    if (exact != isl_bool_true)
    {
        tc_error("Inexact R+");
        //tc_warn("Inexact R+");
    }

    isl_map* Rtile = tc_Rtile_map(II, tile, R_normalized);
    tc_debug_map(Rtile, "R_TILE");

    isl_map* Rtile_plus = tc_transitive_closure(isl_map_copy(Rtile), S, &exact);
    //isl_map* Rtile_plus = tc_transitive_closure_adapter_isl_map(isl_map_copy(Rtile), S, &exact);
    //isl_map* Rtile_plus = tc_transitive_closure_adapter_isl_union_map(isl_map_copy(Rtile), S, &exact);
    //isl_map* Rtile_plus = tc_transitive_closure_adapter_iterative(isl_map_copy(Rtile), S, &exact);
    //isl_map* Rtile_plus = tc_transitive_closure_adapter_floyd_warshall(isl_map_copy(Rtile), S, &exact);    
    
    Rtile_plus = isl_map_coalesce(Rtile_plus);    
    tc_debug_map(Rtile_plus, "R_TILE^+ (exact=%d)", exact);

    if (exact == 0)
    {
        //tc_error("Inexact R_TILE^+");
        tc_warn("Inexact R_TILE^+");
    }
        
    isl_map* Tcycle = tc_Tcycle_map(II, Rtile_plus);
    Tcycle = isl_map_coalesce(Tcycle);
    tc_debug_map(Tcycle, "T_CYCLE");

    isl_set* Incycles = isl_set_union(isl_map_domain(isl_map_copy(Tcycle)), isl_map_range(isl_map_copy(Tcycle)));
    Incycles = isl_set_coalesce(Incycles);
    
    isl_map_free(Tcycle);

    // II_SET_A = II_SET - IN_CYCLES
    isl_set* ii_set_a = isl_set_subtract(isl_set_copy(ii_set), isl_set_copy(Incycles));
    tc_debug_set(ii_set_a, "II_SET_A");

    isl_set* tile_a = isl_set_intersect_params(isl_set_copy(tile), tc_make_set_constraints(isl_set_copy(ii_set_a), II));
    tile_a = isl_set_coalesce(tile_a);
    tc_debug_set(tile_a, "TILE_A");

    // II_SET_C = II_SET * IN_CYCLES
    isl_set* ii_set_c = isl_set_copy(Incycles);
    //isl_set* ii_set_c = isl_set_union(isl_map_domain(isl_map_copy(Incycles)), isl_map_range(isl_map_copy(Incycles)));
    tc_debug_set(ii_set_c, "II_SET_C");

    isl_set* tile_c = isl_set_intersect_params(isl_set_copy(tile), tc_make_set_constraints(isl_set_copy(ii_set_c), II));
    tile_c = isl_set_coalesce(tile_c);
    tc_debug_set(tile_c, "TILE_C");
    
    // CORRECTION

    isl_set* tile_lt = tc_tile_lt_set(tile_c, ii_set_c, II);
    isl_set* tile_gt = tc_tile_gt_set(tile_c, ii_set_c, II);
    
    tc_debug_set(tile_lt, "TILE_LT");
    
    tc_debug_set(tile_gt, "TILE_GT");
    
    // TILE_ITR = TILE - R+(TILE_GT)
    isl_set* tile_itr = isl_set_subtract(isl_set_copy(tile_c), isl_set_apply(isl_set_copy(tile_gt), isl_map_copy(R_plus_normalized)));
    
    tc_debug_set(tile_itr, "TILE_ITR");
    
    // TVLD_LT = (R+(TILE_ITR) * TILE_LT) - R+(TILE_GT)
    isl_set* tvld_lt = isl_set_apply(isl_set_copy(tile_itr), isl_map_copy(R_plus_normalized));
    tvld_lt = isl_set_intersect(tvld_lt, isl_set_copy(tile_lt));
    tvld_lt = isl_set_subtract(tvld_lt, isl_set_apply(isl_set_copy(tile_gt), isl_map_copy(R_plus_normalized)));
    
    tc_debug_set(tvld_lt, "TVLD_LT");
    
    // TILE_VLD = TILE_ITR + TVLD_LT
    isl_set* tile_vld = isl_set_union(tile_itr, tvld_lt);
    tile_vld = isl_set_coalesce(tile_vld);
    
    tc_debug_set(tile_vld, "TILE_VLD");

    // CORRECTION END

    isl_set* tile_mc = isl_set_union(tile_a, tile_vld);
    isl_set* ii_set_mc = isl_set_union(ii_set_a, ii_set_c);
    //isl_set* ii_set_mc = tc_lift_up_set_params(isl_set_params(isl_set_copy(tile_mc)), II);

    tile_mc = isl_set_coalesce(tile_mc);
    ii_set_mc = isl_set_coalesce(ii_set_mc);

    tc_debug_set(tile_mc, "TILE_MC");
    tc_debug_set(ii_set_mc, "II_SET_MC");

    tc_debug_bool(isl_set_is_equal(tile, tile_mc), "TILE = TILE_MC");

    if (tc_options_is_report(options))
    {
        isl_set* bounds = tc_options_get_report_bounds(options, ctx);
        
        struct tc_tile_statistics* stats = tc_compute_tile_statistics(tile_mc, ii_set_mc, II, bounds, LD, S, scop->reads, scop->writes, scop, options, blocks);
        
        tc_tile_statistics_print(options->output, stats);
        
        tc_tile_statistics_free(stats);

        isl_set_free(bounds);
    }

    isl_map* Rtile_mc = tc_Rtile_map(II, tile_mc, R_normalized);
    tc_debug_map(Rtile_mc, "R_TILE_MC");

    tc_scheduling(scop, options, LD, S, R, ii_set_mc, tile_mc, Rtile_mc, II, I);

    isl_set_free(tile_c);

    isl_set_free(tile_lt);
    isl_set_free(tile_gt);

    isl_set_free(tile);
    isl_set_free(ii_set);

    isl_map_free(Rtile);
    isl_map_free(Rtile_plus);

    isl_set_free(Incycles);

    isl_map_free(R_normalized);
    isl_map_free(R_plus_normalized);

    isl_basic_set_free(sample);
    isl_space_free(space);
    isl_set_free(LD_normalized);
}
