#include "diamond_tiling.h"
#include "scop.h"
#include "tiling.h"
#include "utility.h"
#include "options.h"
#include "debug.h"
#include "transitive_closure.h"
#include "omp_cpu_codegen.h"
#include "omp_gpu_codegen.h"
#include "serial_codegen.h"
#include "tile_statistics.h"
#include "input_output.h"
#include "tuples.h"

#include <isl/id.h>

#include <string.h>
#include <stddef.h>

#include <vector>
#include <map>
#include <string>
#include <list>

void tc_algorithm_diamond_tiling(struct tc_scop* scop, struct tc_options* options, bool expansion)
{
    isl_ctx* ctx = scop->ctx;

    isl_union_set* LD = isl_union_set_copy(scop->domain);
    tc_debug_uset(LD, "LD");

    isl_union_map* R = isl_union_map_copy(scop->relation);
    tc_debug_umap(R, "R");

    isl_union_map* S = isl_union_map_copy(scop->schedule);
    tc_debug_umap(S, "S");

    isl_set* LD_normalized = tc_normalize_union_set(LD, S);
    tc_debug_set(LD_normalized, "LDnorm");

    isl_id_list* PARAMS = tc_get_params_ids(LD_normalized);

    isl_map* R_normalized = tc_normalize_union_map(R, S);
    tc_debug_map(R_normalized, "Rnorm");

    isl_union_set* LD_unbounded = tc_union_set_remove_params(isl_union_set_copy(LD), PARAMS);
    //LD_unbounded = isl_union_set_universe(LD_unbounded);
    tc_debug_uset(LD_unbounded, "LD_unbounded");

    isl_set* LD_unbounded_norm = tc_remove_params(isl_set_copy(LD_normalized), PARAMS);
    //LD_unbounded_norm = isl_set_universe(isl_set_get_space(isl_set_copy(LD_unbounded_norm)));
    tc_debug_set(LD_unbounded_norm, "LD_unbounded_norm");

    const bool drop_bounds = (tc_options_is_set(options, NULL, "--drop-bounds") == 1);

    if (drop_bounds)
    {
        R_normalized = tc_map_project_out_params(R_normalized, PARAMS);
        tc_debug_map(R_normalized, "Rnorm_T");
    }

    isl_basic_set* sample = isl_set_sample(isl_set_copy(LD_normalized));
    tc_debug_bset(sample, "SAMPLE");

    isl_space* space = isl_basic_set_get_space(sample);
    tc_debug_space(space, "SPACE");

    isl_id_list* I = tc_ids_sequence(ctx, "i", isl_space_dim(space, isl_dim_set));
    isl_id_list* II = tc_ids_sequence(ctx, "ii", isl_space_dim(space, isl_dim_set));
    isl_id_list* IItime = tc_ids_sub(II, 0, 1);
    isl_id_list* IIspace = tc_ids_sub(II, 1, isl_id_list_n_id(II));
    isl_id_list* IIspace0 = tc_ids_sub(IIspace, 0, 1);

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

    std::map<std::string, std::vector<int> > blocks = tc_options_blocks(options);

    isl_set* tile;
    isl_set* ii_set;
    //tc_tile_loop_nest(LD, S, II, I, &tile, &ii_set, blocks);
    //tc_tile_loop_nest2(LD_unbounded, S, II, I, &tile, &ii_set, blocks);
    //tc_tile_loop_nest3(LD, S, II, I, &tile, &ii_set, blocks);
    tc_tile_normalized_loop_nest(LD_normalized, II, I, &tile, &ii_set, blocks);
    //isl_set* tile_orig = isl_set_copy(tile);
    //isl_set* ii_set_orig = isl_set_copy(ii_set);

    if (drop_bounds)
    {
        tile = tc_project_out_params(tile, PARAMS);
        ii_set = tc_project_out_params(ii_set, PARAMS);
    }

    tc_debug_set(tile, "TILE");
    tc_debug_set(ii_set, "II_SET");

    isl_map* R_plus_normalized_inv = isl_map_reverse(isl_map_copy(R_plus_normalized));

    isl_map* R_star_normalized = isl_map_union(isl_map_copy(R_plus_normalized), tc_make_identity(isl_map_copy(R_plus_normalized)));

    isl_map* R_star_normalized_inv = isl_map_reverse(isl_map_copy(R_star_normalized));

    // TSLICE_CURR = [t,ii1,...,iin] -> { [e] : exists [t,ii1',...,iin] in II_SET and e TILE([t,ii1',...,iin]) }
    isl_set* tslice_curr = tc_tile_set_of(tile, ii_set, II, &tc_tuples_eq_params_except_1);
    //isl_set* tslice_curr = tc_project_out_params(isl_set_copy(tile), IIspace);
    tslice_curr = isl_set_remove_redundancies(tslice_curr);
    tslice_curr = isl_set_coalesce(tslice_curr);
    tc_debug_set(tslice_curr, "TSLICE_CURR");

    // TILE_FIX = (R*)^-1(TILE) * TSLICE_CURR
    isl_set* tile_fix = isl_set_intersect(isl_set_apply(isl_set_copy(tile), isl_map_copy(R_star_normalized_inv)),
                                          isl_set_copy(tslice_curr));
    tile_fix = isl_set_remove_redundancies(tile_fix);
    tile_fix = isl_set_coalesce(tile_fix);
    tc_debug_set(tile_fix, "TILE_FIX");

    tc_debug("# Overlapping regions extraction");

    isl_set_list* tiles_phases = isl_set_list_alloc(ctx, 1);
    tiles_phases = isl_set_list_add(tiles_phases, isl_set_copy(tile_fix));
    int k = 2;
    isl_id_list* IIspace0_k_prev = IIspace0;
    while (true)
    {
        // TILE_OVER_k = [II,JJ,KK,...] -> { [e] | II,JJ,KK,... in II_SET and e in TILE_FIX(II) and e in TILE_FIX(JJ) e in TILE_FIX(KK) and e in ... and IItime = JJtime and JJtime = KKtime ... and IIspace << JJspace and JJspace << KKspace and ... }

        // TILE_OVER_k = [II,JJ,KK,...] -> { [e] | II,JJ,KK... in II_SET and e in (TILE_FIX(II) * TILE_FIX(JJ) * TILE_FIX(KK) * ...) and ii0 = jj0 and jj0 = kk0 and ii2 = jj2 and jj2 = kk2 and ... and ii1 < jj1 and jj1 < kk1 and ...  }

        isl_set* tile_over_prev = isl_set_list_get_set(tiles_phases, isl_set_list_n_set(tiles_phases) - 1);

        char k_str[20];
        snprintf(k_str, sizeof(k_str), "_%d", k);
        isl_id_list* IIspace0_k = tc_ids_add_suffix(IIspace0, k_str);

        isl_set* ii_set_constraints = tc_make_set_constraints(isl_set_copy(ii_set), II);
        isl_set* ii_set_constraints_k = tc_rename_params(isl_set_copy(ii_set_constraints), IIspace0, IIspace0_k);

        isl_set* tile_fix_k = tc_rename_params(isl_set_copy(tile_fix), IIspace0, IIspace0_k);

        isl_id_list* IIall = isl_id_list_concat(isl_id_list_copy(IIspace0_k_prev), isl_id_list_copy(IIspace0_k));
        isl_set* constraints_params = tc_make_params(ctx, IIall, tc_tuples_lt(IIspace0_k_prev, IIspace0_k).c_str());

        isl_set* tile_over_k = isl_set_intersect_params(isl_set_copy(tile_over_prev), constraints_params);

        tile_over_k = isl_set_intersect(tile_over_k, tile_fix_k);
        tile_over_k = isl_set_intersect_params(tile_over_k, ii_set_constraints);
        tile_over_k = isl_set_intersect_params(tile_over_k, ii_set_constraints_k);

        tile_over_k = isl_set_remove_redundancies(tile_over_k);
        tile_over_k = isl_set_coalesce(tile_over_k); // ?

        tc_debug_set(tile_over_k, "TILE_OVER_%d", k);

        isl_bool tile_over_is_empty = isl_set_is_empty(tile_over_k);

        if (tile_over_is_empty == isl_bool_false)
        {
            isl_id_list* IIspace_all = isl_id_list_copy(IIspace0);
            for (int i = 2; i <= k; ++i)
            {
                char k_str[20];
                snprintf(k_str, sizeof(k_str), "_%d", i);
                IIspace_all = isl_id_list_concat(IIspace_all, tc_ids_add_suffix(IIspace0, k_str));
            }
            isl_id_list* IIspace_all_prev = tc_ids_sub(IIspace_all, 0, isl_id_list_n_id(IIspace_all) - 1);

            tc_debug_id_list(IIspace_all, "IIspace_all");

            isl_set* tile_over_k_T = tc_project_out_params(isl_set_copy(tile_over_k), IIspace_all);
            tile_over_k_T = tc_project_out_params(tile_over_k_T, PARAMS);
            tc_debug_set(tile_over_k_T, "TILE_OVER_K_T");

            isl_set* tile_over_prev_T = tc_project_out_params(isl_set_copy(tile_over_prev), IIspace_all_prev);
            tile_over_prev_T = tc_project_out_params(tile_over_prev_T, PARAMS);
            tc_debug_set(tile_over_prev_T, "TILE_OVER_PREV_T");

            isl_bool tile_over_k_equals_prev = isl_set_is_equal(tile_over_k_T, tile_over_prev_T);
            //isl_set_free(tile_over_k_T);
            //isl_set_free(tile_over_prev_T);

            if (tile_over_k_equals_prev == isl_bool_true)
            {
                tc_warn("TILE_OVER_k equals TILE_OVER_k-1, breaking.");

                isl_set* tile_over_k_fixed = isl_set_copy(tile_over_prev);
                tile_over_k_fixed = tc_project_out_params(tile_over_k_fixed, IIspace_all_prev);
                tile_over_k_fixed = isl_set_intersect_params(tile_over_k_fixed, tc_make_params(ctx, IIspace_all_prev, ""));
                for (int i = 0; i < isl_id_list_n_id(IIspace_all_prev); ++i)
                {
                    tile_over_k_fixed = tc_set_fix_param_value(tile_over_k_fixed, isl_id_list_get_id(IIspace_all_prev, i), 0);
                }

                tc_debug_set(tile_over_k_fixed, "TILE_OVER_%d_fixed", k);

                int n_set = isl_set_list_n_set(tiles_phases);
                tc_debug("n_set =%d", n_set);
                tiles_phases = isl_set_list_set_set(tiles_phases, n_set - 1, tile_over_k_fixed);
                break;
            }
            else
            {
                tiles_phases = isl_set_list_add(tiles_phases, tile_over_k);
                IIspace0_k_prev = IIspace0_k;
                ++k;
            }
        }
        else
        {
            isl_set_free(tile_over_k);
            break;
        }
    }

    tc_debug("# Overlapping regions removal phase");

    isl_set* space_over = NULL;
    for (int k = isl_set_list_n_set(tiles_phases); k >= 1; --k)
    {
        isl_set* tile_over = isl_set_list_get_set(tiles_phases, k - 1);

        if (space_over != NULL)
        {
            tc_debug_set(space_over, "SPACE_OVER");
        }

        isl_id_list* IIspace0_all = isl_id_list_copy(IIspace0);
        for (int i = 2; i <= k; ++i)
        {
            char k_str[20];
            snprintf(k_str, sizeof(k_str), "_%d", i);
            isl_id_list* IIspace0_k = tc_ids_add_suffix(IIspace0, k_str);

            IIspace0_all = isl_id_list_concat(IIspace0_all, isl_id_list_copy(IIspace0_k));
        }

        if (space_over != NULL)
        {
            tile_over = isl_set_subtract(tile_over, isl_set_copy(space_over));
            tile_over = isl_set_remove_redundancies(tile_over);
            tile_over = isl_set_coalesce(tile_over);
        }
        tc_debug_set(tile_over, "TILE_PHASE_%d", k);

        tiles_phases = isl_set_list_set_set(tiles_phases, k - 1, isl_set_copy(tile_over));

        if (space_over == NULL)
        {
            //!space_over = isl_set_copy(tile_over);
            //space_over = tc_project_out_params(isl_set_copy(tile_over), IIspace0_all);
            space_over = tc_project_out_params(isl_set_copy(tile_over), IIspace0_all);
            //space_over = isl_set_compute_divs(space_over);
            space_over = isl_set_remove_redundancies(space_over);
            space_over = isl_set_coalesce(space_over);
            //space_over = isl_set_detect_equalities(space_over);
        }
        else
        {
            //!space_over = isl_set_union(space_over, isl_set_copy(tile_over));
            space_over = isl_set_union(space_over, tc_project_out_params(isl_set_copy(tile_over), IIspace0_all));
            //space_over = isl_set_compute_divs(space_over);
            space_over = isl_set_remove_redundancies(space_over);
            space_over = isl_set_coalesce(space_over);
            //space_over = isl_set_detect_equalities(space_over);
        }
    }

    if (expansion)
    {
        tc_debug("# Expansion phase");

        isl_set* tile_k_last = isl_set_list_get_set(tiles_phases, 0);
        tc_debug_set(tile_k_last, "TILE_K_LAST");

        isl_set* ii_set_k_last = isl_set_params(isl_set_copy(tile_k_last));
        ii_set_k_last = tc_lift_up_set_params(ii_set_k_last, II);
        ii_set_k_last = isl_set_remove_redundancies(ii_set_k_last);
        ii_set_k_last = isl_set_coalesce(ii_set_k_last);
        tc_debug_set(ii_set_k_last, "II_SET_K_LAST");

        isl_set* tslice_curr_for_last = tc_tile_set_of(tile_k_last, ii_set_k_last, II, &tc_tuples_eq_params_except_1);
        //isl_set* tslice_curr_for_last = tc_tile_set_of(tile_k_last, ii_set_k_last, II, &tc_tuples_eq_param_0);
        tslice_curr_for_last = isl_set_remove_redundancies(tslice_curr_for_last);
        tslice_curr_for_last = isl_set_coalesce(tslice_curr_for_last);
        tc_debug_set(tslice_curr_for_last, "TSLICE_CURR_FOR_LAST");

        // TSLICE_NEXT = [t,ii1,...,iin] -> { [e] : exists [t',ii1',...,iin'] in II_SET : t' = t + 1 and e in TILE([t',ii1',...,iin']) }
        isl_set* tslice_next_for_last = tc_tile_set_of(tile, ii_set, II, &tc_tuples_next_param_0_eq_other_params_except_1);
        //tslice_next_for_last = tc_project_out_params(tslice_next_for_last, IIspace);
        tslice_next_for_last = isl_set_remove_redundancies(tslice_next_for_last);
        tslice_next_for_last = isl_set_coalesce(tslice_next_for_last);
        tc_debug_set(tslice_next_for_last, "TSLICE_NEXT");

        // TILE_K_LAST_COMPL = TSLICE_CURR_FOR_LAST - TILE_K_LAST
        isl_set* tile_k_last_compl = isl_set_subtract(isl_set_copy(tslice_curr_for_last), isl_set_copy(tile_k_last));
        tile_k_last_compl = isl_set_remove_redundancies(tile_k_last_compl);
        tile_k_last_compl = isl_set_coalesce(tile_k_last_compl);
        tc_debug_set(tile_k_last_compl, "TILE_K_LAST_COMPL");

        // TILE_RANGE = R+(TILE_K_LAST) * TSLICE_NEXT
        isl_set* tile_range = isl_set_intersect(isl_set_apply(isl_set_copy(tile_k_last), isl_map_copy(R_plus_normalized)),
                                                isl_set_copy(tslice_next_for_last));
        tile_range = isl_set_remove_redundancies(tile_range);
        tile_range = isl_set_coalesce(tile_range);
        tc_debug_set(tile_range, "TILE_RANGE");

        // TILE_RANGE_COMPL = TSLICE_NEXT - TILE_RANGE
        isl_set* tile_range_compl = isl_set_subtract(isl_set_copy(tslice_next_for_last), isl_set_copy(tile_range));
        tile_range_compl = isl_set_remove_redundancies(tile_range_compl);
        tile_range_compl = isl_set_coalesce(tile_range_compl);
        tc_debug_set(tile_range_compl, "TILE_RANGE_COMPL");

        // TILE_EXP = TILE_K_LAST + TILE_RANGE - R+(TILE_K_LAST_COMPL) - R+(TILE_RANGE_COMPL)
        isl_set* tile_exp = isl_set_subtract(
                                isl_set_subtract(
                                    isl_set_union(isl_set_copy(tile_k_last),
                                                  isl_set_copy(tile_range)),
                                    isl_set_apply(isl_set_copy(tile_k_last_compl), isl_map_copy(R_plus_normalized))),
                                isl_set_apply(isl_set_copy(tile_range_compl), isl_map_copy(R_plus_normalized)));
        tile_exp = isl_set_remove_redundancies(tile_exp);
        tile_exp = isl_set_coalesce(tile_exp);
        tc_debug_set(tile_exp, "TILE_EXP");

        // TILE_EXP == TILE_k_last !
        tiles_phases = isl_set_list_set_set(tiles_phases, 0, isl_set_copy(tile_exp));

        isl_set* space_exp = tc_project_out_params(isl_set_copy(tile_exp), IIspace);
        space_exp = isl_set_remove_redundancies(space_exp);
        space_exp = isl_set_coalesce(space_exp);
        tc_debug_set(space_exp, "SPACE_EXP");

        isl_set* ii_set_space_exp = isl_set_params(isl_set_copy(space_exp));
        ii_set_space_exp = tc_lift_up_set_params(ii_set_space_exp, IItime);

        // TSLICE_PREV = [t] -> { [e] : exists t : t' + 1 = t and e in SPACE_EXP([t']) };
        //isl_set* tslice_prev = tc_tile_set_of(space_exp, ii_set_space_exp, IItime, &tc_tuples_prev_param_0);
        isl_set* tslice_prev = tc_tile_set_of(space_exp, ii_set_space_exp, IItime, &tc_tuples_prev_param_0);//_eq_other_params_except_1);
        tslice_prev = isl_set_remove_redundancies(tslice_prev);
        tslice_prev = isl_set_coalesce(tslice_prev);
        tc_debug_set(tslice_prev, "TSLICE_PREV");

        // TILE_K := TILE_K - TSLICE_PREV
        for (int k = isl_set_list_n_set(tiles_phases); k >= 1; --k)
        {
            isl_set* tile_phase = isl_set_list_get_set(tiles_phases, k - 1);

            tile_phase = isl_set_subtract(tile_phase, isl_set_copy(tslice_prev));
            tile_phase = isl_set_remove_redundancies(tile_phase);
            tile_phase = isl_set_coalesce(tile_phase);

            tiles_phases = isl_set_list_set_set(tiles_phases, k - 1, tile_phase);
        }
    }

    tc_debug("# Correction phase");

    for (int k = isl_set_list_n_set(tiles_phases); k >= 1; --k)
    {
        isl_set* tile_phase = isl_set_list_get_set(tiles_phases, k - 1);
        tc_debug_set(tile_phase, "TILE_PHASE");

        isl_id_list* IIspace0_all = isl_id_list_copy(IIspace0);
        isl_id_list* II_all = isl_id_list_copy(II);
        for (int i = 2; i <= k; ++i)
        {
            char k_str[20];
            snprintf(k_str, sizeof(k_str), "_%d", i);
            isl_id_list* IIspace0_k = tc_ids_add_suffix(IIspace0, k_str);

            IIspace0_all = isl_id_list_concat(IIspace0_all, isl_id_list_copy(IIspace0_k));
            II_all = isl_id_list_concat(II_all, IIspace0_k);
        }
        tc_debug_id_list(IIspace0_all, "IIspace0_all(%d)", k);

        isl_set* ii_set_phase = tc_ii_set(isl_set_copy(tile_phase), II_all);
        tc_debug_set(ii_set_phase, "II_SET_ALL");

        isl_set* tile_lt = tc_tile_set_of(tile_phase, ii_set_phase, II_all, &tc_tuples_t_ii_lt);
        //tile_lt = isl_set_remove_redundancies(tile_lt);
        tc_debug_set(tile_lt, "TILE_LT");

        isl_set* tile_gt = tc_tile_set_of(tile_phase, ii_set_phase, II_all, &tc_tuples_t_ii_gt);
        //tile_gt = isl_set_remove_redundancies(tile_gt);
        tc_debug_set(tile_gt, "TILE_GT");

        isl_set* tslice_curr_subspace = tc_tile_set_of(tile_phase, ii_set_phase, II_all, &tc_tuples_eq_t_ii);
        //isl_set* tslice_curr = tc_project_out_params(isl_set_copy(tile), IIspace);
        tslice_curr_subspace = isl_set_remove_redundancies(tslice_curr_subspace);
        tslice_curr_subspace = isl_set_coalesce(tslice_curr_subspace);
        tc_debug_set(tslice_curr_subspace, "TSLICE_CURR_SUBSPACE");

        isl_set* tile_gt_subspace = isl_set_intersect(isl_set_apply(isl_set_copy(tile_gt), isl_map_copy(R_plus_normalized)),
                                                                    isl_set_copy(tslice_curr_subspace));
        tile_gt_subspace = isl_set_remove_redundancies(tile_gt_subspace);
        tile_gt_subspace = isl_set_coalesce(tile_gt_subspace);

        // TILE_ITR = TILE_ALL - R+(TILE_GT) * TSLICE_CURR_SUBSPACE
        isl_set* tile_itr = isl_set_subtract(isl_set_copy(tile_phase),
                                            isl_set_copy(tile_gt_subspace));
        tile_itr = isl_set_remove_redundancies(tile_itr);
        tile_itr = isl_set_coalesce(tile_itr);
        tc_debug_set(tile_itr, "TILE_ITR");

        // TVLD_LT = (R+(TILE_ITR) * TILE_LT) - (R+(TILE_GT) * TSLICE_CURR_SUBSPACE)
        isl_set* tvld_lt = isl_set_apply(isl_set_copy(tile_itr), isl_map_copy(R_plus_normalized));
        tvld_lt = isl_set_intersect(tvld_lt, isl_set_copy(tile_lt));
        tvld_lt = isl_set_subtract(tvld_lt, isl_set_copy(tile_gt_subspace));
        tvld_lt = isl_set_remove_redundancies(tvld_lt);
        tvld_lt = isl_set_coalesce(tvld_lt);

        tc_debug_set(tvld_lt, "TVLD_LT");

        // TILE_VLD = TILE_ITR + TVLD_LT
        isl_set* tile_vld = isl_set_union(tile_itr, tvld_lt);
        tile_vld = isl_set_remove_redundancies(tile_vld);
        tile_vld = isl_set_coalesce(tile_vld);

        tc_debug_set(tile_vld, "TILE_VLD");

        tc_debug_bool(isl_set_is_equal(tile_phase, tile_vld), "TILE_PHASE = TILE_VLD");

        tiles_phases = isl_set_list_set_set(tiles_phases, k - 1, tile_vld);
    }

    tc_debug("# Id extension phase");

    isl_id* k_id = isl_id_alloc(ctx, "k", NULL);
    isl_id_list* K = isl_id_list_from_id(isl_id_copy(k_id));
    isl_id_list* II_k = isl_id_list_copy(II);
    II_k = isl_id_list_insert(II_k, 1, isl_id_copy(k_id));
    tc_debug_id_list(II_k, "IIk");

    isl_set* tile_all = NULL;
    for (int k = isl_set_list_n_set(tiles_phases); k >= 1; --k)
    {
        // TILE_PHASE [t,ii1,ii2] -> TILE_PHASE_EXT [t,k,ii1,ii2]
        isl_set* tile_phase = isl_set_list_get_set(tiles_phases, k - 1);

        for (int i = k; i < isl_set_list_n_set(tiles_phases); ++i)
        {
            char k_str[20];
            snprintf(k_str, sizeof(k_str), "_%d", i + 1);
            isl_id_list* IIspace0_k = tc_ids_add_suffix(IIspace0, k_str);

            isl_set* IIspace0_k_param = tc_make_params(ctx, IIspace0_k, "");

            IIspace0_k_param = isl_set_fix_si(IIspace0_k_param, isl_dim_param, 0, 0);
            //tc_debug_set(tile_phase, "tile_phase");
            //tc_debug_set(IIspace0_k_param, "IIspace0_k_param");

            tile_phase = isl_set_intersect_params(tile_phase, IIspace0_k_param);
            //tc_debug_set(tile_phase, "tile_phase_intersected");
            //tile_phase = tc_lift_up_set_params(tile_phase, IIspace0_k);
            //tile_phase = isl_set_coalesce(tile_phase);
            //tc_debug_set(tile_phase, "tile_phase_parameterized");
        }

        isl_set* k_param = tc_make_params(ctx, K, "");
        k_param = isl_set_fix_si(k_param, isl_dim_param, 0, isl_set_list_n_set(tiles_phases) - k);

        isl_set* tile_phase_ext = isl_set_intersect_params(tile_phase, k_param);
        tile_phase_ext = isl_set_coalesce(tile_phase_ext);

        if (tile_all == NULL)
        {
            tile_all = tile_phase_ext;
        }
        else
        {
            tile_all = isl_set_union(tile_all, tile_phase_ext);
        }
    }

    tc_debug_set(tile_all, "TILE_ALL");

    for (int i = 2; i <= isl_set_list_n_set(tiles_phases); ++i)
    {
        char k_str[20];
        snprintf(k_str, sizeof(k_str), "_%d", i);
        isl_id_list* IIspace0_k = tc_ids_add_suffix(IIspace0, k_str);
        II_k = isl_id_list_concat(II_k, IIspace0_k);
    }

    tc_debug("# Tuple extension phase");

    isl_set* tile_vld = isl_set_copy(tile_all);
    isl_set* ii_set_all = tc_ii_set(isl_set_copy(tile_vld), II_k);

    tc_debug_set(tile_vld, "TILE_PLOT");
    tc_debug_set(ii_set_all, "II_SET_PLOT");

    isl_set* tile_ext = tc_lift_up_set_params(isl_set_copy(tile_vld), II_k);

    //tile_ext = isl_set_compute_divs(tile_ext);
    tile_ext = isl_set_remove_redundancies(tile_ext);
    tile_ext = isl_set_coalesce(tile_ext);
    tile_ext = isl_set_detect_equalities(tile_ext);
    tc_debug_set(tile_ext, "TILE_EXT");

    if (drop_bounds)
    {
        isl_set* LD_normalized_ext = isl_set_insert_dims(isl_set_copy(LD_normalized), isl_dim_set, 0, isl_id_list_n_id(II_k));
        tc_debug_set(LD_normalized_ext, "LD_normalized_ext");
        tile_ext = isl_set_intersect(tile_ext, LD_normalized_ext);
        tile_ext = isl_set_coalesce(tile_ext);
        tc_debug_set(tile_ext, "TILE_EXT_domain");

        isl_set* tile_vld_LD = isl_set_intersect(isl_set_copy(tile_vld), isl_set_copy(LD_normalized));
        tile_vld_LD = isl_set_coalesce(tile_vld_LD);
        tc_debug_set(tile_vld_LD, "TILE_VLD_LD");

        isl_set* ii_set_ext = tc_ii_set(isl_set_copy(tile_vld_LD), II_k);
        tc_debug_set(ii_set_ext, "II_SET_LD");
    }


    if (tc_options_is_report(options))
    {
        isl_set* bounds = tc_options_get_report_bounds(options, ctx);

        struct tc_tile_statistics* stats = tc_compute_tile_statistics(tile_vld, ii_set_all, II_k, bounds, LD, S, scop->reads, scop->writes, scop, options, blocks);

        tc_tile_statistics_print(options->output, stats);

        tc_tile_statistics_free(stats);
        isl_set_free(bounds);
    }

    isl_union_map* S_ext = tc_extend_schedule(isl_union_map_copy(S), isl_id_list_n_id(II_k));
    tc_debug_umap(S_ext, "S_ext");

    isl_id_list* iterators = isl_id_list_copy(II_k);
    iterators = isl_id_list_concat(iterators, isl_id_list_copy(I));

    //isl_id_list* parallel_iterators = isl_id_list_drop(isl_id_list_copy(II), 0, 1);
    isl_id_list* parallel_iterators = tc_ids_sub(II, 1, 2);

    tc_debug_id_list(iterators, "iterators");
    tc_debug_id_list(parallel_iterators, "parallel_iterators");

    enum tc_codegen_enum codegen = tc_options_codegen(options);

    if (tc_codegen_enum_serial == codegen)
    {
        tc_codegen_serial(scop, options, S_ext, tile_ext, iterators);
    }
    else if (tc_codegen_enum_omp_cpu_for == codegen)
    {
        tc_codegen_omp_parallel_for(scop, options, S_ext, tile_ext, iterators, parallel_iterators, 0);
    }
    else if (tc_codegen_enum_omp_cpu_task == codegen)
    {
        tc_codegen_omp_task_for(scop, options, S_ext, tile_ext, iterators, parallel_iterators, 0);
    }
    else if (tc_codegen_enum_omp_gpu == codegen)
    {
        tc_codegen_omp_gpu(scop, options, S_ext, tile_ext, iterators, II);
    }

    isl_union_set_free(LD);
    isl_union_map_free(R);
    isl_union_map_free(S);
    isl_map_free(R_normalized);
    isl_space_free(space);
    isl_basic_set_free(sample);
    isl_id_list_free(I);
    isl_id_list_free(II);
    isl_set_free(tile);
    isl_set_free(ii_set);
    isl_map_free(R_star_normalized);
    isl_map_free(R_plus_normalized);
    isl_map_free(R_star_normalized_inv);
    isl_map_free(R_plus_normalized_inv);
    isl_set_free(tslice_curr);
    //isl_set_free(tslice_next);
    isl_set_free(tile_fix);
    //isl_set_free(tile_exp);
}
