#include "correction_inv_tiling.h"
#include "scheduling.h"
#include "tiling.h"
#include "utility.h"
#include "scop.h"
#include "options.h"
#include "tile_statistics.h"
#include "debug.h"
#include "transitive_closure.h"
#include "slicing.h"

#include <isl/ctx.h>
#include <isl/space.h>
#include <isl/id.h>
#include <isl/set.h>
#include <isl/map.h>
#include <isl/union_set.h>
#include <isl/union_map.h>

#include <barvinok/isl.h>

#include <stdio.h>

#include <string>
#include <vector>
#include <map>

void tc_algorithm_correction_inv_tiling(struct tc_scop* scop, struct tc_options* options)
{
    isl_ctx* ctx = scop->ctx;
    
    isl_union_set* LD = isl_union_set_copy(scop->domain);
    tc_debug_uset(LD, "LD");
    
    isl_union_map* S = isl_union_map_copy(scop->schedule);
    //S = tc_loop_interchange(S, "S3", "j", "k");
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
    
    int exact;
    isl_map* R_plus_normalized = tc_transitive_closure(isl_map_copy(R_normalized), S, &exact);
    tc_debug_map(R_plus_normalized, "R^+ (exact=%d)", exact);

    if (exact == 0)
    {
        tc_error("Inexact R+");
    }

    isl_set* tile_lt = tc_tile_lt_set(tile, ii_set, II);
    isl_set* tile_gt = tc_tile_gt_set(tile, ii_set, II);
    
    tc_debug_set(tile_lt, "TILE_LT");
    tc_debug_set(tile_gt, "TILE_GT");

    isl_map* R_plus_normalized_inv = isl_map_reverse(isl_map_copy(R_plus_normalized));

    // TILE_ITR = TILE - (R+)^-1(TILE_LT)
    isl_set* tile_itr = isl_set_subtract(isl_set_copy(tile), isl_set_apply(isl_set_copy(tile_lt), isl_map_copy(R_plus_normalized_inv)));
    
    tc_debug_set(tile_itr, "TILE_ITR");
    
    // TVLD_GT = ((R+)^-1(TILE_ITR)) * TILE_GT - (R+)^-1(TILE_LT)
    isl_set* tvld_gt = isl_set_apply(isl_set_copy(tile_itr), isl_map_copy(R_plus_normalized_inv));
    tvld_gt = isl_set_intersect(tvld_gt, isl_set_copy(tile_gt));
    tvld_gt = isl_set_subtract(tvld_gt, isl_set_apply(isl_set_copy(tile_lt), isl_map_copy(R_plus_normalized_inv)));
    
    tc_debug_set(tvld_gt, "TVLD_GT");
    
    // TILE_TRG = TILE_ITR + TVLD_GT
    isl_set* tile_trg = isl_set_union(tile_itr, tvld_gt);
    tile_trg = isl_set_coalesce(tile_trg);

    tc_debug_set(tile_trg, "TILE_TRG");

    tc_debug_bool(isl_set_is_equal(tile, tile_trg), "TILE = TILE_TRG");

    if (tc_options_is_report(options))
    {
        isl_set* bounds = tc_options_get_report_bounds(options, ctx);

        struct tc_tile_statistics* stats = tc_compute_tile_statistics(tile_trg, ii_set, II, bounds, LD, S, scop->reads, scop->writes, scop, options, blocks);

        tc_tile_statistics_print(options->output, stats);

        tc_tile_statistics_free(stats);

        isl_set_free(bounds);
    }

    isl_map* Rtile = tc_Rtile_map(II, tile_trg, R_normalized);
    tc_debug_map(Rtile, "R_TILE");
    
    tc_scheduling(scop, options, LD, S, R, ii_set, tile_trg, Rtile, II, I);

    isl_set_free(tile_lt);
    isl_set_free(tile_gt);
    isl_set_free(tile);
    isl_map_free(R_normalized);
    isl_map_free(R_plus_normalized);
    isl_map_free(R_plus_normalized_inv);
    isl_basic_set_free(sample);
    isl_space_free(space);
    isl_set_free(LD_normalized);
}
