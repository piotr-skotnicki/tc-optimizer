#include "isl_scheduling.h"
#include "serial_codegen.h"
#include "omp_cpu_codegen.h"
#include "omp_gpu_codegen.h"
#include "debug.h"
#include "utility.h"
#include "scop.h"
#include "options.h"
#include "tuples.h"
#include "input_output.h"

#include <isl/ctx.h>
#include <isl/map.h>
#include <isl/set.h>
#include <isl/schedule_node.h>
#include <isl/options.h>

#include <string.h>
#include <stddef.h>

void tc_scheduling_isl(struct tc_scop* scop, struct tc_options* options, __isl_take isl_union_set* LD, __isl_take isl_union_map* S, __isl_take isl_union_map* R, __isl_take isl_set* ii_set, __isl_take isl_set* tile, __isl_take isl_map* Rtile, __isl_take isl_id_list* II, __isl_take isl_id_list* I, enum tc_scheduling_enum scheduling)
{
    isl_ctx* ctx = isl_union_set_get_ctx(LD);

    if (tc_scheduling_enum_isl == scheduling || tc_scheduling_enum_isl_wavefronting == scheduling)
    {
        isl_options_set_schedule_algorithm(ctx, ISL_SCHEDULE_ALGORITHM_ISL);
    }
    else if (tc_scheduling_enum_feautrier == scheduling)
    {
        isl_options_set_schedule_algorithm(ctx, ISL_SCHEDULE_ALGORITHM_FEAUTRIER);
    }

    isl_union_set* ii_set_denorm = tc_denormalize_set(ii_set, S);
    tc_debug_uset(ii_set_denorm, "II_SET_MC_denorm");

    isl_union_map* Rtile_denorm = tc_denormalize_map(Rtile, S);
    tc_debug_umap(Rtile_denorm, "R_TILE_MC_denorm");

    isl_schedule* schedule = tc_schedule_compute(ii_set_denorm, Rtile_denorm);
    tc_debug_schedule(schedule, "SCHED");

    if (schedule == NULL)
    {
        tc_error("Cannot compute ISL schedule");
        tc_die(tc_exit_code_failure);
    }

    isl_union_map* F = tc_schedule_umap(schedule);
    tc_debug_umap(F, "F");

    isl_union_map* F_best = NULL;
    int parallel_dim = -1;

    if (tc_scheduling_enum_isl == scheduling || tc_scheduling_enum_feautrier == scheduling)
    {
        F_best = isl_union_map_copy(F);
    }
    else if (tc_scheduling_enum_isl_wavefronting == scheduling)
    {
        isl_union_pw_multi_aff* F_umaff = isl_union_pw_multi_aff_from_union_map(isl_union_map_copy(F));

        isl_pw_multi_aff_list* F_maffs = isl_union_pw_multi_aff_get_pw_multi_aff_list(F_umaff);

        isl_union_pw_multi_aff* F_wave = NULL;

        isl_schedule_node* node = isl_schedule_get_root(schedule);

        while (node)
        {
            if (isl_schedule_node_has_parent(node) == isl_bool_true && isl_schedule_node_get_parent_type(node) == isl_schedule_node_band)
            {
                parallel_dim = 1;

                isl_union_pw_multi_aff* node_umaff = isl_schedule_node_get_prefix_schedule_union_pw_multi_aff(node);

                isl_pw_multi_aff_list* node_maffs = isl_union_pw_multi_aff_get_pw_multi_aff_list(node_umaff);

                isl_pw_aff_list* affs = isl_pw_aff_list_alloc(ctx, 10);

                for (int i = 0; i < isl_pw_multi_aff_list_n_pw_multi_aff(node_maffs); ++i)
                {
                    isl_pw_multi_aff* node_maff = isl_pw_multi_aff_list_get_pw_multi_aff(node_maffs, i);

                    isl_pw_aff* aff = NULL;

                    for (int j = 0; j < isl_pw_multi_aff_dim(node_maff, isl_dim_out); ++j)
                    {
                        if (aff == NULL)
                        {
                            aff = isl_pw_multi_aff_get_pw_aff(node_maff, j);
                        }
                        else
                        {
                            aff = isl_pw_aff_add(aff, isl_pw_multi_aff_get_pw_aff(node_maff, j));
                        }
                    }

                    affs = isl_pw_aff_list_add(affs, aff);

                    isl_pw_multi_aff_free(node_maff);
                }

                for (int i = 0; i < isl_pw_multi_aff_list_n_pw_multi_aff(F_maffs); ++i)
                {
                    isl_pw_multi_aff* F_aff = isl_pw_multi_aff_list_get_pw_multi_aff(F_maffs, i);

                    for (int j = 0; j < isl_pw_aff_list_n_pw_aff(affs); ++j)
                    {
                        isl_pw_aff* aff = isl_pw_aff_list_get_pw_aff(affs, j);

                        isl_id* aff_in_id = isl_pw_aff_get_tuple_id(aff, isl_dim_in);

                        if (0 == strcmp(isl_pw_multi_aff_get_tuple_name(F_aff, isl_dim_in), isl_id_get_name(aff_in_id)))
                        {
                            isl_pw_multi_aff* F_new_aff = isl_pw_multi_aff_zero(isl_space_alloc(ctx, 0, isl_pw_multi_aff_dim(F_aff, isl_dim_in)
                                                                                                      , isl_pw_multi_aff_dim(F_aff, isl_dim_out) + 1));
                            F_new_aff = isl_pw_multi_aff_set_tuple_id(F_new_aff, isl_dim_in, aff_in_id);

                            F_new_aff = isl_pw_multi_aff_set_pw_aff(F_new_aff, 0, aff);

                            for (int k = 0; k < isl_pw_multi_aff_dim(F_aff, isl_dim_out); ++k)
                            {
                                F_new_aff = isl_pw_multi_aff_set_pw_aff(F_new_aff, k + 1, isl_pw_multi_aff_get_pw_aff(F_aff, k));
                            }

                            if (F_wave == NULL)
                            {
                                F_wave = isl_union_pw_multi_aff_from_pw_multi_aff(F_new_aff);
                            }
                            else
                            {
                                F_wave = isl_union_pw_multi_aff_add_pw_multi_aff(F_wave, F_new_aff);
                            }

                            break;
                        }
                        else
                        {
                            isl_id_free(aff_in_id);
                            isl_pw_aff_free(aff);
                        }
                    }

                    isl_pw_multi_aff_free(F_aff);
                }

                isl_pw_aff_list_free(affs);
                isl_pw_multi_aff_list_free(node_maffs);
                isl_union_pw_multi_aff_free(node_umaff);

                isl_schedule_node_free(node);
                node = NULL;
            }
            else
            {
                if (isl_schedule_node_has_children(node) == isl_bool_true)
                {
                    if (isl_schedule_node_n_children(node) > 1)
                    {
                        tc_error("Unable to process the scheduling tree");
                        tc_die(tc_exit_code_failure);
                    }
                    else
                    {
                        isl_schedule_node* child = isl_schedule_node_get_child(node, 0);

                        isl_schedule_node_free(node);
                        node = child;
                    }
                }
                else
                {
                    isl_schedule_node_free(node);
                    node = NULL;
                }
            }
        }

        isl_union_pw_multi_aff_free(F_umaff);
        isl_pw_multi_aff_list_free(F_maffs);

        F_best = isl_union_map_from_union_pw_multi_aff(F_wave);
    }

    tc_debug_umap(F_best, "F_best");

    //F_best = isl_union_map_read_from_str(ctx, "[_PB_N] -> { S3[i, j, k] -> [i + j + 2k, j + k, k, i, 1]; S1[i, j, k] -> [i + j + 2k, i + k, j, k, 2]; S2[i, j] -> [3j + i, i + j, j, j, 0] }");

    isl_map* F_best_norm = isl_map_from_union_map(isl_union_map_apply_domain(F_best, isl_union_map_copy(S)));
    tc_debug_map(F_best_norm, "F_best_norm");

    isl_set* deltas_V = tc_schedule_validity(F_best_norm, Rtile);

    if (tc_schedule_validity_is_valid(deltas_V) == isl_bool_true)
    {
        tc_debug("# Schedule F_best is valid");
    }
    else
    {
        tc_error("# Schedule F_best is not valid");
        tc_die(tc_exit_code_failure);
    }    
    
    // TILE_EXT := { [k, II, I] : k = F(II) and II in II_SET and I in TILE(II) }

    isl_id_list* k = tc_ids_sequence(ctx, "k", isl_map_dim(F_best_norm, isl_dim_out));

    isl_id_list* k_II_I = isl_id_list_concat(isl_id_list_concat(isl_id_list_copy(k), isl_id_list_copy(II)), isl_id_list_copy(I));
    isl_id_list* k_I = isl_id_list_concat(isl_id_list_copy(k), isl_id_list_copy(I));

    isl_set* tile_ext = tc_make_set(ctx, k_II_I, k_I, NULL);

    isl_set* F_constraints = tc_make_map_constraints(isl_map_copy(F_best_norm), II, k);
    
    isl_set* tile_constraints = tc_make_set_constraints(isl_set_copy(tile), I);
    
    isl_set* ii_set_constraints = tc_make_set_constraints(isl_set_copy(ii_set), II);

    tile_ext = isl_set_intersect_params(tile_ext, F_constraints);
    tile_ext = isl_set_intersect_params(tile_ext, tile_constraints);
    tile_ext = isl_set_intersect_params(tile_ext, ii_set_constraints);

    tile_ext = tc_project_out_params(tile_ext, k_II_I);
    tc_debug_set(tile_ext, "tile_ext");

    tile_ext = isl_set_coalesce(tile_ext);

    isl_union_map* S_prim = tc_extend_schedule(isl_union_map_copy(S), isl_id_list_n_id(k));
    tc_debug_umap(S_prim, "S_prim");
    S_prim = isl_union_map_intersect_range(S_prim, isl_union_set_from_set(tile_ext));

    isl_id_list* parallel_iterators = NULL;

    if (parallel_dim != -1)
    {
        parallel_iterators = tc_ids_get(k, parallel_dim);
    }

    tc_debug_id_list(parallel_iterators, "parallel_iterators");

    enum tc_codegen_enum codegen = tc_options_codegen(options);

    if (tc_codegen_enum_serial == codegen)
    {
        tc_codegen_serial(scop, options, S_prim, k_I);
    }
    else if (tc_codegen_enum_omp_cpu_for == codegen)
    {
        tc_codegen_omp_parallel_for(scop, options, S_prim, k_I, parallel_iterators, 0);
    }
    else if (tc_codegen_enum_omp_cpu_task == codegen)
    {
        tc_codegen_omp_task_for(scop, options, S_prim, k_I, parallel_iterators, 0);
    }
    else if (tc_codegen_enum_omp_gpu == codegen)
    {
        tc_codegen_omp_gpu(scop, options, S_prim, k_I, parallel_iterators);
    }

    isl_id_list_free(k);
    isl_id_list_free(k_I);
    isl_id_list_free(k_II_I);
    isl_id_list_free(parallel_iterators);

    isl_map_free(Rtile);
    isl_set_free(tile);
    isl_set_free(ii_set);
    isl_id_list_free(I);
    isl_id_list_free(II);
    isl_union_set_free(ii_set_denorm);
    isl_union_map_free(Rtile_denorm);

    isl_schedule_free(schedule);
    isl_map_free(F_best_norm);
    isl_union_map_free(F);

    isl_union_set_free(LD);
    isl_union_map_free(S);
    isl_union_map_free(R);
}

void tc_scheduling_adapter_isl(struct tc_scop* scop, struct tc_options* options, __isl_take isl_union_set* LD, __isl_take isl_union_map* S, __isl_take isl_union_map* R, __isl_take isl_set* ii_set, __isl_take isl_set* tile, __isl_take isl_map* Rtile, __isl_take isl_id_list* II, __isl_take isl_id_list* I)
{
    tc_scheduling_isl(scop, options, LD, S, R, ii_set, tile, Rtile, II, I, tc_scheduling_enum_isl);
}

void tc_scheduling_adapter_feautrier(struct tc_scop* scop, struct tc_options* options, __isl_take isl_union_set* LD, __isl_take isl_union_map* S, __isl_take isl_union_map* R, __isl_take isl_set* ii_set, __isl_take isl_set* tile, __isl_take isl_map* Rtile, __isl_take isl_id_list* II, __isl_take isl_id_list* I)
{
    tc_scheduling_isl(scop, options, LD, S, R, ii_set, tile, Rtile, II, I, tc_scheduling_enum_feautrier);
}

void tc_scheduling_adapter_isl_wavefronting(struct tc_scop* scop, struct tc_options* options, __isl_take isl_union_set* LD, __isl_take isl_union_map* S, __isl_take isl_union_map* R, __isl_take isl_set* ii_set, __isl_take isl_set* tile, __isl_take isl_map* Rtile, __isl_take isl_id_list* II, __isl_take isl_id_list* I)
{
    tc_scheduling_isl(scop, options, LD, S, R, ii_set, tile, Rtile, II, I, tc_scheduling_enum_isl_wavefronting);
}

__isl_give isl_schedule* tc_schedule_compute(__isl_keep isl_union_set* IS, __isl_keep isl_union_map* R)
{
    isl_schedule* schedule = isl_union_set_compute_schedule(
                                 isl_union_set_copy(IS)
                               , isl_union_map_copy(R)
                               , isl_union_map_copy(R));

    return schedule;
}

__isl_give isl_schedule* tc_schedule_compute_normalized(__isl_keep isl_set* IS, __isl_keep isl_map* R)
{
    isl_schedule* schedule = isl_union_set_compute_schedule(
                                 isl_union_set_from_set(isl_set_copy(IS))
                               , isl_union_map_from_map(isl_map_copy(R))
                               , isl_union_map_from_map(isl_map_copy(R)));

    return schedule;
}

__isl_give isl_schedule* tc_schedule_compute_from_map(__isl_keep isl_map* R)
{
    isl_set* IS = isl_set_union(isl_map_domain(isl_map_copy(R)), isl_map_range(isl_map_copy(R)));

    isl_schedule* schedule = tc_schedule_compute_normalized(IS, R);

    isl_set_free(IS);

    return schedule;
}

__isl_give isl_schedule* tc_schedule_compute_from_umap(__isl_keep isl_union_map* R)
{
    isl_union_set* IS = isl_union_set_union(isl_union_map_domain(isl_union_map_copy(R)), isl_union_map_range(isl_union_map_copy(R)));

    isl_schedule* schedule = tc_schedule_compute(IS, R);

    isl_union_set_free(IS);

    return schedule;
}

__isl_give isl_map* tc_schedule_map(__isl_keep isl_schedule* schedule)
{
    return isl_map_from_union_map(isl_schedule_get_map(schedule));
}

__isl_give isl_union_map* tc_schedule_umap(__isl_keep isl_schedule* schedule)
{
    return isl_schedule_get_map(schedule);
}

__isl_give isl_map* tc_schedule_validity_map(__isl_keep isl_id_list* II, __isl_keep isl_map* R, __isl_keep isl_map* S)
{
    // VALIDITY := { [II] -> [JJ] : II in domain(R) and JJ in R(II) and S(II) >= S(JJ) }

    // VALIDITY := { [II] -> [JJ] : II in domain(R) and JJ in R(II) and exists II', JJ' : II' = S(II) and JJ' = S(JJ) and II' >= JJ' }
   
    isl_ctx* ctx = isl_map_get_ctx(R);

    isl_id_list* JJ = tc_ids_sequence(ctx, "jj", isl_id_list_n_id(II));

    isl_id_list* II_S = tc_ids_sequence(ctx, "ii", isl_map_dim(S, isl_dim_out));
    isl_id_list* JJ_S = tc_ids_sequence(ctx, "jj", isl_map_dim(S, isl_dim_out));

    isl_id_list* IIprim = tc_ids_prim(II_S);
    isl_id_list* JJprim = tc_ids_prim(JJ_S);

    isl_id_list* II_JJ_IIprim_JJprim = isl_id_list_concat(isl_id_list_concat(isl_id_list_copy(II), isl_id_list_copy(JJ))
                                                        , isl_id_list_concat(isl_id_list_copy(IIprim), isl_id_list_copy(JJprim)));

    isl_map* validity = tc_make_map(ctx, II_JJ_IIprim_JJprim, II, JJ, tc_disjunction(tc_tuples_gt(IIprim, JJprim), tc_tuples_eq(IIprim, JJprim)).c_str());

    isl_set* S_constraints_II = tc_make_map_constraints(isl_map_copy(S), II, IIprim);
    isl_set* S_constraints_JJ = tc_make_map_constraints(isl_map_copy(S), JJ, JJprim);
    isl_set* R_constraints = tc_make_map_constraints(isl_map_copy(R), II, JJ);
    
    validity = isl_map_intersect_params(validity, R_constraints);
    validity = isl_map_intersect_params(validity, S_constraints_II);
    validity = isl_map_intersect_params(validity, S_constraints_JJ);

    validity = tc_map_project_out_params(validity, II_JJ_IIprim_JJprim);

    validity = isl_map_coalesce(validity);

    isl_id_list_free(II_S);
    isl_id_list_free(JJ_S);
    isl_id_list_free(JJ);
    isl_id_list_free(IIprim);
    isl_id_list_free(JJprim);
    isl_id_list_free(II_JJ_IIprim_JJprim);

    return validity;
}

__isl_give isl_set* tc_schedule_validity(__isl_keep isl_map* m, __isl_keep isl_map* R)
{
    // C = (m^-1) . R . m
    isl_map* C = isl_map_apply_range(isl_map_apply_range(isl_map_reverse(isl_map_copy(m)), isl_map_copy(R)), isl_map_copy(m));
    //tc_debug_map(C, "C");

    isl_set* deltas = isl_map_deltas(C);
    //tc_debug_set(deltas, "deltas");

    isl_ctx* ctx = isl_map_get_ctx(R);

    int n_out = isl_map_dim(m, isl_dim_out);

    isl_id_list* I = tc_ids_sequence(ctx, "i", n_out);
    isl_id_list* J = tc_ids_sequence(ctx, "j", n_out);

    isl_set* m_range = isl_map_range(isl_map_copy(m));

    isl_set* zero = isl_set_from_point(isl_point_zero(isl_set_get_space(m_range)));
    zero = isl_set_remove_dims(zero, isl_dim_param, 0, isl_set_n_param(zero));

    isl_set* zero_constraints_J = tc_make_set_constraints(zero, J);

    isl_set* V = tc_make_set(ctx, J, I, tc_disjunction(tc_tuples_lt(I, J), tc_tuples_eq(I, J)).c_str());

    V = isl_set_intersect_params(V, zero_constraints_J);

    V = tc_project_out_params(V, J);

    V = isl_set_coalesce(V);

    isl_set_free(m_range);
    isl_id_list_free(I);
    isl_id_list_free(J);

    return isl_set_intersect(deltas, V);
}

__isl_give isl_set* tc_schedule_validity_bounded(__isl_keep isl_map* m, __isl_keep isl_map* R, __isl_keep isl_set* bounds)
{
    return tc_set_fix_params_bounds(tc_schedule_validity(m, R), bounds);
}

isl_bool tc_schedule_validity_is_valid(__isl_take isl_set* validity)
{
    const isl_bool is_valid = isl_set_is_empty(validity);

    isl_set_free(validity);

    return is_valid;
}

isl_bool tc_schedule_validity_map_is_valid(__isl_take isl_map* validity)
{
    const isl_bool is_valid = isl_map_is_empty(validity);

    isl_map_free(validity);

    return is_valid;
}
