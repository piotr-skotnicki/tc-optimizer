#include "split_tiling.h"
#include "tiling.h"
#include "scop.h"
#include "utility.h"
#include "options.h"
#include "debug.h"
#include "serial_codegen.h"
#include "omp_cpu_codegen.h"
#include "transitive_closure.h"
#include "tile_statistics.h"
#include "scheduling.h"
#include "tuples.h"
#include "input_output.h"

#include <isl/ctx.h>
#include <isl/set.h>
#include <isl/map.h>
#include <isl/union_set.h>
#include <isl/union_map.h>
#include <isl/val.h>
#include <isl/polynomial.h>

#include <barvinok/isl.h>

#include <stddef.h>
#include <stdio.h>

#include <vector>
#include <map>
#include <string>

void tc_algorithm_split_tiling(struct tc_scop* scop, struct tc_options* options)
{
    isl_ctx* ctx = scop->ctx;
        
    isl_union_set* LD = isl_union_set_copy(scop->domain);
    tc_debug_uset(LD, "LD");
    
    isl_union_map* S = isl_union_map_copy(scop->schedule);
    //S = tc_loop_interchange(S, "S3", "j", "k");
    tc_debug_umap(S, "S");
    
    isl_union_map* R = isl_union_map_copy(scop->relation);
    tc_debug_umap(R, "R");
    
    isl_set* LD_normalized = tc_normalize_union_set(LD, S);
    
    isl_basic_set* sample = isl_set_sample(isl_set_copy(LD_normalized));
    
    isl_space* space = isl_basic_set_get_space(sample);
    
    isl_id_list* params = tc_get_set_params_names(LD_normalized);
    
    isl_id_list* I = tc_ids_sequence(ctx, "i", isl_space_dim(space, isl_dim_set));
    isl_id_list* II = tc_ids_sequence(ctx, "ii", isl_space_dim(space, isl_dim_set));
    isl_id_list* II_k = isl_id_list_insert(isl_id_list_copy(II), 1, isl_id_alloc(ctx, "z", NULL));
    isl_id_list* II_1_to_n = tc_ids_sub(II, 1, isl_space_dim(space, isl_dim_set));
        
    std::map<std::string, std::vector<int> > blocks = tc_options_blocks(options);
    
    int max = tc_options_get_int(options, "-m", "--max");
    
    isl_set* tile;
    isl_set* ii_set;
    
    tc_tile_loop_nest(LD, S, II, I, &tile, &ii_set, blocks);
    
    tc_debug_set(tile, "TILE");
    tc_debug_set(ii_set, "II_SET");
    
    isl_map* R_normalized = tc_normalize_union_map(R, S);
        
    tc_debug_map(R_normalized, "R_norm");

    isl_bool exact = isl_bool_false;
    isl_map* R_plus_normalized = tc_transitive_closure(isl_map_copy(R_normalized), S, &exact);        
    tc_debug_map(R_plus_normalized, "R^+ (exact=%d)", exact);

    isl_union_map* R_plus = tc_denormalize_map(R_plus_normalized, S);
    tc_debug_umap(R_plus, "R^+ denorm");
    isl_union_map_free(R_plus);

    if (exact != isl_bool_true)
    {
        tc_warn("Inexact R^+. The results can be non-optimal. Restart TC with a different transitive closure method.");
        if (!tc_io_confirm(options, "Continue?"))
        {
            tc_die(tc_exit_code_inexact);
        }
    }
    
    int k = 1;
    
    isl_set_list* tiles_k = isl_set_list_alloc(ctx, 1);
    
    while (!isl_set_is_empty(tile))
    {
        tc_debug("# k = %d", k);
        
        //isl_set* tile_gt = tc_tile_gt_set(tile, ii_set, II);
        isl_set* tile_gt = tc_tile_set_of(tile, ii_set, II, &tc_tuples_tgt);
        tc_debug_set(tile_gt, "TILE_GT_%d", k);

        // If one or more elements of the card(TVLD_LT) are parametric, or varied, or greater than max
        if (k <= max)
        {
            // INVALID =  R+(TILE_GT) * TILE
            isl_set* invalid = isl_set_intersect(isl_set_apply(isl_set_copy(tile_gt), isl_map_copy(R_plus_normalized)), isl_set_copy(tile));
            invalid = isl_set_coalesce(invalid);
            tc_debug_set(invalid, "INVALID_%d", k);

            // PROBLEMATIC = R*(INVALID) = INVALID u (R+)(INVALID)
            isl_set* problematic = isl_set_union(invalid, isl_set_apply(isl_set_copy(invalid), isl_map_copy(R_plus_normalized)));
            tc_debug_set(problematic, "PROBLEMATIC_%d", k);

            isl_set* problematic_T = tc_project_out_params(problematic, II_1_to_n);
            //isl_set* problematic_T = tc_project_out_params(problematic, II);
            tc_debug_set(problematic_T, "PROBLEMATIC_T_%d", k);

            // TILEk = TILE - PROBLEMATIC_T
            isl_set* tile_k = isl_set_subtract(isl_set_copy(tile), problematic_T);
            tile_k = isl_set_coalesce(tile_k);

            tc_debug_set(tile_k, "TILE_k_%d", k);

            // TILE = TILE - TILEk
            tile = isl_set_subtract(tile, isl_set_copy(tile_k));

            tc_debug_set(tile, "TILE_%d", k);
            
            tiles_k = isl_set_list_add(tiles_k, tile_k);

            k = k + 1;

            isl_set_free(tile_gt);
        }
        else
        {
            //isl_set* tile_lt = tc_tile_lt_set(tile, ii_set, II);
            isl_set* tile_lt = tc_tile_set_of(tile, ii_set, II, &tc_tuples_tlt);
            tc_debug_set(tile_lt, "TILE_LT_%d", k);

            // TILE_ITR = TILE - R+(TILE_GT)
            isl_set* tile_itr = isl_set_subtract(isl_set_copy(tile), isl_set_apply(isl_set_copy(tile_gt), isl_map_copy(R_plus_normalized)));
            tile_itr = isl_set_coalesce(tile_itr);
            
            tc_debug_set(tile_itr, "TILE_ITR_%d", k);

            // TVLD_LT = (R+(TILE_ITR) * TILE_LT) - R+(TILE_GT)
            isl_set* tvld_lt = isl_set_apply(isl_set_copy(tile_itr), isl_map_copy(R_plus_normalized));
            tvld_lt = isl_set_intersect(tvld_lt, isl_set_copy(tile_lt));
            tvld_lt = isl_set_subtract(tvld_lt, isl_set_apply(isl_set_copy(tile_gt), isl_map_copy(R_plus_normalized)));
            tvld_lt = isl_set_coalesce(tvld_lt);
            
            tc_debug_set(tvld_lt, "TVLD_LT_%d", k);

            // TILEk = TILE_ITR + TVLD_LT
            isl_set* tile_k = isl_set_union(tile_itr, tvld_lt);
            tile_k = isl_set_coalesce(tile_k);

            tc_debug_set(tile_k, "TILE_k_%d", k);
            
            tiles_k = isl_set_list_add(tiles_k, tile_k);

            k = k + 1;

            isl_set_free(tile_lt);
            isl_set_free(tile_gt);

            break;
        }
    }
        
    isl_set* tile_ext = NULL;
    
    for (int i = 0; i < isl_set_list_n_set(tiles_k); ++i)
    {
        isl_set* tile_k = isl_set_list_get_set(tiles_k, i);

        // 7. Enlarge the tuple of set II:  [ii1, ii2, ..., iid] to the tuple [ii1, k, ii2, ...,iid]
        isl_set* tile_k_ext = tc_lift_up_set_params(tile_k, II);        
        tile_k_ext = isl_set_insert_dims(tile_k_ext, isl_dim_set, 1, 1);
        tile_k_ext = isl_set_fix_si(tile_k_ext, isl_dim_set, 1, i);
        tile_k_ext = tc_lift_down_set_vars(tile_k_ext, II_k);

        tc_debug_set(tile_k_ext, "TILE_K_EXT_%d", i);
        
        if (NULL == tile_ext)
        {
            tile_ext = tile_k_ext;
        }
        else
        {
            tile_ext = isl_set_union(tile_ext, tile_k_ext);
            tile_ext = isl_set_coalesce(tile_ext);
        }
    }

    // tile_ext = isl_set_coalesce(tile_ext);

    isl_set_list_free(tiles_k);

    isl_set* tile_ext_params = isl_set_params(isl_set_copy(tile_ext));
    isl_set* ii_set_ext = tc_lift_up_set_params(tile_ext_params, II_k);

    tc_debug_set(tile_ext, "TILE_EXT");
    tc_debug_set(ii_set_ext, "II_SET_EXT");

    if (tc_options_is_report(options))
    {
        isl_set* bounds = tc_options_get_report_bounds(options, ctx);

        struct tc_tile_statistics* stats = tc_compute_tile_statistics(tile_ext, ii_set_ext, II_k, bounds, LD, S, scop->reads, scop->writes, scop, options, blocks);

        tc_tile_statistics_print(options->output, stats);

        tc_tile_statistics_free(stats);

        isl_set_free(bounds);
    }

    isl_map* Rtile = tc_Rtile_map(II_k, tile_ext, R_normalized);    
    tc_debug_map(Rtile, "R_TILE");

    tc_scheduling(scop, options, LD, S, R, ii_set_ext, tile_ext, Rtile, II_k, I);
    
    isl_set_free(LD_normalized);
    isl_id_list_free(params);
    isl_map_free(R_normalized);
    isl_map_free(R_plus_normalized);
    isl_set_free(tile);
    isl_set_free(ii_set);
    //isl_set_free(tile_lt);
    //isl_set_free(tile_gt);
    isl_id_list_free(II_1_to_n);
    isl_id_list_free(II);
    isl_space_free(space);
    isl_basic_set_free(sample);
}
