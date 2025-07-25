#include "correction_tiling.h"
#include "scheduling.h"
#include "tiling.h"
#include "utility.h"
#include "scop.h"
#include "options.h"
#include "tile_statistics.h"
#include "debug.h"
#include "transitive_closure.h"
#include "slicing.h"
#include "serial_codegen.h"
#include "tuples.h"

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

void tc_algorithm_correction_tiling(struct tc_scop* scop, struct tc_options* options)
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
    
    //R = tc_remove_loop_independent_dependences(R, S, groups);
    
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
    
    isl_set* tile_lt = tc_tile_lt_set(tile, ii_set, II);
    isl_set* tile_gt = tc_tile_gt_set(tile, ii_set, II);
    
    tc_debug_set(tile_lt, "TILE_LT");
    
    tc_debug_set(tile_gt, "TILE_GT");
    
    // TILE_ITR = TILE - R+(TILE_GT)
    isl_set* tile_itr = isl_set_subtract(isl_set_copy(tile), isl_set_apply(isl_set_copy(tile_gt), isl_map_copy(R_plus_normalized)));
    
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

    tc_debug_bool(isl_set_is_equal(tile, tile_vld), "TILE = TILE_VLD");
    
    if (tc_options_is_report(options))
    {
        isl_set* bounds = tc_options_get_report_bounds(options, ctx);
        
        struct tc_tile_statistics* stats = tc_compute_tile_statistics(tile_vld, ii_set, II, bounds, LD, S, scop->reads, scop->writes, scop, options, blocks);
        
        tc_tile_statistics_print(options->output, stats);
        
        // Tileability
        
        long n_inst = stats->n_statement_instances;
                
        // INVALID = R+(TILE_GT) * TILE
        isl_set* invalid = isl_set_intersect(isl_set_apply(isl_set_copy(tile_gt), isl_map_copy(R_plus_normalized)), isl_set_copy(tile));
        tc_debug_set(invalid, "INVALID");
                
        double per_tileable = 0.0;
        
        if (isl_set_is_empty(invalid))
        {
            per_tileable = 100.0;
            fprintf(options->output, "Tileability: %.8f %%\n", per_tileable);
        }
        else
        {
            isl_set* invalid_bounded = tc_set_fix_params_bounds(isl_set_copy(invalid), isl_set_copy(bounds));
            long n_invalid = tc_set_card_value(invalid_bounded);
            double per_invalid = (100.0 * n_invalid) / n_inst;
            fprintf(options->output, "Invalid: %.8f %% (%ld)\n", per_invalid, n_invalid);

            /*
            // SPACE_LT, SPACE_GT
            {
                isl_id_list* e = tc_ids_sequence(ctx, "e", isl_set_dim(invalid, isl_dim_set));    
                isl_id_list* e_prim = tc_ids_prim(e);
                isl_id_list* e_e_prim = isl_id_list_concat(isl_id_list_copy(e), isl_id_list_copy(e_prim));

                //isl_set* invalid_T = tc_project_out_params(isl_set_copy(invalid), II);
                isl_set* invalid_bounded = tc_set_fix_params_bounds(isl_set_copy(invalid), isl_set_copy(bounds));
                isl_set* invalid_T = tc_project_out_params(invalid_bounded, II);
                tc_debug_set(invalid_T, "INVALID_T");

                isl_set* invalid_T_lexmin = isl_set_lexmin(isl_set_copy(invalid_T));
                isl_set* invalid_T_lexmin_params = tc_make_set_constraints(invalid_T_lexmin, e);

                isl_set* invalid_T_lexmax = isl_set_lexmax(isl_set_copy(invalid_T));
                isl_set* invalid_T_lexmax_params = tc_make_set_constraints(invalid_T_lexmax, e);

                isl_set* tile_bounded = tc_set_fix_params_bounds(isl_set_copy(tile), isl_set_copy(bounds));

                isl_set* tile_T = tc_project_out_params(tile_bounded, II);

                isl_set* tile_T_params = tc_make_set_constraints(tile_T, e_prim);

                isl_set* space_lt = tc_make_set(ctx, e_e_prim, e_prim, tc_tuples_lt(e_prim, e).c_str());

                space_lt = isl_set_intersect_params(space_lt, isl_set_copy(tile_T_params));
                space_lt = isl_set_intersect_params(space_lt, invalid_T_lexmin_params);

                space_lt = tc_project_out_params(space_lt, e_e_prim);
                space_lt = isl_set_coalesce(space_lt);
                tc_debug_set(space_lt, "SPACE_LT");

                isl_set* space_gt = tc_make_set(ctx, e_e_prim, e_prim, tc_tuples_gt(e_prim, e).c_str());

                space_gt = isl_set_intersect_params(space_gt, isl_set_copy(tile_T_params));
                space_gt = isl_set_intersect_params(space_gt, invalid_T_lexmax_params);

                space_gt = tc_project_out_params(space_gt, e_e_prim);
                space_gt = isl_set_coalesce(space_gt);
                tc_debug_set(space_gt, "SPACE_GT");

                long n_space_lt_inst = tc_set_card_value(space_lt);
                double per_space_lt_tileable = (100.0 * n_space_lt_inst) / n_inst;
                fprintf(options->output, "Space LT: %.8f %% (%ld)\n", per_space_lt_tileable, n_space_lt_inst);

                long n_space_gt_inst = tc_set_card_value(space_gt);
                double per_space_gt_tileable = (100.0 * n_space_gt_inst) / n_inst;
                fprintf(options->output, "Space GT: %.8f %% (%ld)\n", per_space_gt_tileable, n_space_gt_inst);

                double per_space_tileable = (100.0 * (n_space_lt_inst + n_space_gt_inst)) / n_inst;
                fprintf(options->output, "Space LT+GT: %.8f %% (%ld)\n", per_space_tileable, n_space_lt_inst + n_space_gt_inst);

                isl_set_free(invalid_T);
                isl_set_free(tile_T_params);
                isl_id_list_free(e);
                isl_id_list_free(e_prim);
                isl_id_list_free(e_e_prim);
            }
            */

            {
                isl_map* R_plus_normalized_inv = isl_map_reverse(isl_map_copy(R_plus_normalized));

                // PROBLEMATIC = R*^-1(INVALID) = INVALID u (R+^-1)(INVALID)            
                //isl_set* problematic = isl_set_union(isl_set_copy(invalid), isl_set_apply(isl_set_copy(invalid), isl_map_copy(R_plus_normalized_inv)));

                // PROBLEMATIC = R+^-1(INVALID)
                isl_set* problematic = isl_set_apply(isl_set_copy(invalid), isl_map_copy(R_plus_normalized_inv));

                isl_set* problematic_bounded = tc_set_fix_params_bounds(isl_set_copy(problematic), isl_set_copy(bounds));
                
                long n_prob = tc_set_card_value(problematic_bounded);

                per_tileable = (100.0 * (n_inst - n_prob)) / n_inst;

                isl_map_free(R_plus_normalized_inv);

                /*
                isl_set* problematic_T = tc_project_out_params(isl_set_copy(problematic), II);

                isl_set* tile_np = isl_set_subtract(isl_set_copy(tile), isl_set_copy(problematic_T));

                tc_codegen_serial(scop, options, isl_union_map_copy(S), isl_set_copy(problematic_T), I);

                tc_scheduling_lex(scop, options, isl_union_set_copy(LD), isl_union_map_copy(S), isl_union_map_copy(R), isl_set_copy(ii_set), isl_set_copy(tile_np), isl_map_copy(Rtile), isl_id_list_copy(II), isl_id_list_copy(I));

                isl_set_free(problematic_T);
                isl_set_free(tile_np);
                */
                isl_set_free(problematic);

                fprintf(options->output, "Tileability after: %.8f %%\n", per_tileable);
            }

            // ----------------------------------------

            {
                // PROBLEMATIC = R*(INVALID) = INVALID u (R+)(INVALID)            
                isl_set* problematic = isl_set_union(isl_set_copy(invalid), isl_set_apply(isl_set_copy(invalid), isl_map_copy(R_plus_normalized)));

                isl_set* problematic_bounded = tc_set_fix_params_bounds(isl_set_copy(problematic), isl_set_copy(bounds));
                
                long n_prob = tc_set_card_value(problematic_bounded);

                per_tileable = (100.0 * (n_inst - n_prob)) / n_inst;

                isl_set_free(problematic);

                fprintf(options->output, "Tileability before: %.8f %%\n", per_tileable);
            }

            // ----------------------------------------

            /*
            isl_map* R_plus_normalized_inv = isl_map_reverse(isl_map_copy(R_plus_normalized));
            
            isl_map* R_plus_U_R_plus_inv_normalized = isl_map_union(isl_map_copy(R_plus_normalized), isl_map_copy(R_plus_normalized_inv));
            
            // PROBLEMATIC = INVALID + (R+ U R+^-1)(INVALID)            
            isl_set* problematic = isl_set_union(isl_set_copy(invalid), isl_set_apply(isl_set_copy(invalid), isl_map_copy(R_plus_U_R_plus_inv_normalized)));
            
            isl_set* problematic_bounded = tc_set_fix_params_bounds(problematic, isl_set_copy(bounds));
                            
            isl_pw_qpolynomial* problematic_card = isl_set_card(problematic_bounded);

            isl_val* problematic_card_val = isl_pw_qpolynomial_max(problematic_card);

            long n_prob = isl_val_get_num_si(problematic_card_val);
            
            isl_val_free(problematic_card_val);
            
            per_tileable = (100.0 * (n_inst - n_prob)) / n_inst;
            
            isl_map_free(R_plus_normalized_inv);
            isl_map_free(R_plus_U_R_plus_inv_normalized);
            */

            // ----------------------------------------
        }
        
        //fprintf(options->output, "Tileability: %.8f %%\n", per_tileable);
        
        isl_set_free(invalid);
        
        tc_tile_statistics_free(stats);
        isl_set_free(bounds);
    }
    
    isl_map* Rtile = tc_Rtile_map(II, tile_vld, R_normalized);
    
    tc_debug_map(Rtile, "R_TILE");
    
    tc_scheduling(scop, options, LD, S, R, ii_set, tile_vld, Rtile, II, I);
    
    isl_set_free(tile_lt);
    isl_set_free(tile_gt);
    isl_set_free(tile);
    isl_map_free(R_normalized);
    isl_map_free(R_plus_normalized);
    isl_basic_set_free(sample);
    isl_space_free(space);
    isl_set_free(LD_normalized);
}
