#include "inv_correction_tiling.h"
#include "scheduling.h"
#include "tiling.h"
#include "utility.h"
#include "scop.h"
#include "options.h"
#include "tile_statistics.h"
#include "debug.h"
#include "transitive_closure.h"
#include "slicing.h"
#include "input_output.h"

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

void tc_algorithm_inv_correction_tiling(struct tc_scop* scop, struct tc_options* options)
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

    isl_map* R_inv_normalized = isl_map_reverse(isl_map_copy(R_normalized));
    tc_debug_map(R_inv_normalized, "R^(-1)_norm");

    isl_bool exact = isl_bool_false;
    isl_map* R_inv_plus_normalized = tc_transitive_closure(R_inv_normalized, S, &exact);
    tc_debug_map(R_inv_plus_normalized, "(R^(-1))+ (exact=%d)", exact);

    if (exact != isl_bool_true)
    {
        tc_warn("Inexact R+. The results can be non-optimal. Restart TC with a different transitive closure method.");
        if (!tc_io_confirm(options, "Continue?"))
        {
            tc_die(tc_exit_code_inexact);
        }
    }

    isl_set* tile_lt = tc_tile_lt_set(tile, ii_set, II);
    isl_set* tile_gt = tc_tile_gt_set(tile, ii_set, II);
    
    tc_debug_set(tile_lt, "TILE_LT");
    tc_debug_set(tile_gt, "TILE_GT");

    // TILE_ISG = ((R^-1)+)(TILE) * TILE_GT
    isl_set* tile_isg = isl_set_apply(isl_set_copy(tile), isl_map_copy(R_inv_plus_normalized));
    tile_isg = isl_set_intersect(tile_isg, isl_set_copy(tile_gt));
    tile_isg = isl_set_coalesce(tile_isg);

    tc_debug_set(tile_isg, "TILE_ISG");

    // TILE_CORR = TILE + TILE_ISG - ((R^-1)+)(TILE_LT)
    isl_set* tile_corr = isl_set_union(isl_set_copy(tile), tile_isg);
    tile_corr = isl_set_subtract(tile_corr, isl_set_apply(isl_set_copy(tile_lt), isl_map_copy(R_inv_plus_normalized)));
    tile_corr = isl_set_coalesce(tile_corr);

    tc_debug_set(tile_corr, "TILE_CORR");

    tc_debug_bool(isl_set_is_equal(tile, tile_corr), "TILE = TILE_CORR");

    if (tc_options_is_report(options))
    {
        isl_set* bounds = tc_options_get_report_bounds(options, ctx);

        struct tc_tile_statistics* stats = tc_compute_tile_statistics(tile_corr, ii_set, II, bounds, LD, S, scop->reads, scop->writes, scop, options, blocks);

        tc_tile_statistics_print(options->output, stats);

        tc_tile_statistics_free(stats);

        isl_set_free(bounds);
    }

    isl_map* Rtile = tc_Rtile_map(II, tile_corr, R_normalized);
    tc_debug_map(Rtile, "R_TILE");

    tc_scheduling(scop, options, LD, S, R, ii_set, tile_corr, Rtile, II, I);

    isl_set_free(tile_lt);
    isl_set_free(tile_gt);
    isl_set_free(tile);
    isl_map_free(R_normalized);
    isl_map_free(R_inv_plus_normalized);
    isl_basic_set_free(sample);
    isl_space_free(space);
    isl_set_free(LD_normalized);
}
