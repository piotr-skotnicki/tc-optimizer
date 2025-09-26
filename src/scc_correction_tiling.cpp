#include "scc_correction_tiling.h"
#include "scheduling.h"
#include "tiling.h"
#include "utility.h"
#include "scop.h"
#include "options.h"
#include "tile_statistics.h"
#include "debug.h"
#include "transitive_closure.h"
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
#include <stddef.h>

#include <string>
#include <vector>
#include <map>

void tc_algorithm_scc_correction_tiling(struct tc_scop* scop, struct tc_options* options)
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
    
    tc_tile_loop_nest(options, LD, S, II, I, &tile, &ii_set, blocks, groups);
    
    tc_debug_set(tile, "TILE");
    tc_debug_set(ii_set, "II_SET");
        
    isl_map* R_normalized = tc_normalize_union_map(R, S);
    tc_debug_map(R_normalized, "R_norm");
    
    isl_bool exact = isl_bool_false;
    isl_map* R_plus_normalized = tc_transitive_closure(isl_map_copy(R_normalized), S, &exact);

    tc_debug_map(R_plus_normalized, "R^+ (exact=%d)", exact);

    if (exact != isl_bool_true)
    {
        tc_warn("Inexact R^+. The results can be non-optimal. Restart TC with a different transitive closure method.");
        if (!tc_io_confirm(options, "Continue?"))
        {
            tc_die(tc_exit_code_inexact);
        }
    }

    isl_map* Rtile = tc_Rtile_map(II, tile, R_normalized);
    tc_debug_map(Rtile, "R_TILE");

    isl_map* Rtile_plus = tc_transitive_closure(isl_map_copy(Rtile), S, &exact);

    Rtile_plus = isl_map_coalesce(Rtile_plus);    
    tc_debug_map(Rtile_plus, "R_TILE^+ (exact=%d)", exact);

    if (exact != isl_bool_true)
    {
        tc_warn("Inexact R_TILE^+. The results can be non-optimal. Restart TC with a different transitive closure method.");
        if (!tc_io_confirm(options, "Continue?"))
        {
            tc_die(tc_exit_code_inexact);
        }
    }
        
    isl_map* Tcycle = tc_Tcycle_map(II, Rtile_plus);
    Tcycle = isl_map_coalesce(Tcycle);
    tc_debug_map(Tcycle, "T_CYCLE");

    // CORRECTION

    isl_id_list* JJ = tc_ids_sequence(ctx, "jj", isl_space_dim(space, isl_dim_set));
    isl_id_list* II_JJ = isl_id_list_concat(isl_id_list_copy(II), isl_id_list_copy(JJ));

    // TILE_C_LT = [II] -> { [I] | exists JJ: II,JJ in II_SET and I in TILE(JJ) and II in T_CYCLE(JJ) }
    // TILE_C_GT = [II] -> { [I] | exists JJ: II,JJ in II_SET and I in TILE(JJ) and JJ in T_CYCLE(II) }

    isl_set* Tcycle_II_JJ = tc_make_map_constraints(isl_map_copy(Tcycle), II, JJ);
    isl_set* Tcycle_JJ_II = tc_make_map_constraints(isl_map_copy(Tcycle), JJ, II);

    isl_set* ii_set_constraints = tc_make_set_constraints(isl_set_copy(ii_set), II);
    isl_set* ii_set_constraints_JJ = tc_make_set_constraints(isl_set_copy(ii_set), JJ);

    isl_set* tile_JJ = tc_rename_params(isl_set_copy(tile), II, JJ);

    isl_set* tile_lt = tc_make_set(ctx, II_JJ, I, NULL);
    tile_lt = isl_set_intersect_params(tile_lt, isl_set_copy(ii_set_constraints));
    tile_lt = isl_set_intersect_params(tile_lt, isl_set_copy(ii_set_constraints_JJ));
    tile_lt = isl_set_intersect_params(tile_lt, Tcycle_JJ_II);
    tile_lt = isl_set_intersect(tile_lt, isl_set_copy(tile_JJ));
    tile_lt = tc_project_out_params(tile_lt, JJ);
    tile_lt = isl_set_coalesce(tile_lt);

    isl_set* tile_gt = tc_make_set(ctx, II_JJ, I, NULL);
    tile_gt = isl_set_intersect_params(tile_gt, ii_set_constraints);
    tile_gt = isl_set_intersect_params(tile_gt, ii_set_constraints_JJ);
    tile_gt = isl_set_intersect_params(tile_gt, Tcycle_II_JJ);
    tile_gt = isl_set_intersect(tile_gt, tile_JJ);
    tile_gt = tc_project_out_params(tile_gt, JJ);
    tile_gt = isl_set_coalesce(tile_gt);

    tc_debug_set(tile_lt, "TILE_LT_SCC");
    
    tc_debug_set(tile_gt, "TILE_GT_SCC");

    // TILE_ITL = R+(TILE) * TILE_LT_SCC
    isl_set* tile_itl = isl_set_apply(isl_set_copy(tile), isl_map_copy(R_plus_normalized));
    tile_itl = isl_set_intersect(tile_itl, isl_set_copy(tile_lt));
    tile_itl = isl_set_coalesce(tile_itl);

    tc_debug_set(tile_itl, "TILE_ITL");

    // TILE_CORR = TILE + TILE_ITL - R+(TILE_GT_SCC)
    isl_set* tile_corr = isl_set_union(isl_set_copy(tile), tile_itl);
    tile_corr = isl_set_subtract(tile_corr, isl_set_apply(isl_set_copy(tile_gt), isl_map_copy(R_plus_normalized)));
    tile_corr = isl_set_coalesce(tile_corr);

    tc_debug_set(tile_corr, "TILE_CORR_SCC");

    // CORRECTION END

    tc_debug_set(ii_set, "II_SET_SCC");

    tc_debug_bool(isl_set_is_equal(tile, tile_corr), "TILE = TILE_CORR_SCC");

    if (tc_options_is_report(options))
    {
        isl_set* bounds = tc_options_get_report_bounds(options, ctx);
        
        struct tc_tile_statistics* stats = tc_compute_tile_statistics(tile_corr, ii_set, II, bounds, LD, S, scop->reads, scop->writes, scop, options, blocks);
        
        tc_tile_statistics_print(options->output, stats);
        
        tc_tile_statistics_free(stats);

        isl_set_free(bounds);
    }

    isl_map* Rtile_corr_scc = tc_Rtile_map(II, tile_corr, R_normalized);
    tc_debug_map(Rtile_corr_scc, "R_TILE_CORR_SCC");

    tc_scheduling(scop, options, LD, S, R, ii_set, tile_corr, Rtile_corr_scc, II, I);

    isl_set_free(tile_lt);
    isl_set_free(tile_gt);

    isl_set_free(tile);

    isl_map_free(Rtile);
    isl_map_free(Rtile_plus);

    isl_map_free(Tcycle);

    isl_id_list_free(JJ);
    isl_id_list_free(II_JJ);

    isl_map_free(R_normalized);
    isl_map_free(R_plus_normalized);

    isl_basic_set_free(sample);
    isl_space_free(space);
    isl_set_free(LD_normalized);
}
