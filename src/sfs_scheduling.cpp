#include "sfs_scheduling.h"
#include "omp_cpu_codegen.h"
#include "omp_gpu_codegen.h"
#include "utility.h"
#include "scop.h"
#include "slicing.h"
#include "serial_codegen.h"
#include "options.h"
#include "debug.h"
#include "transitive_closure.h"

#include <isl/ctx.h>
#include <isl/id.h>
#include <isl/union_map.h>
#include <isl/union_set.h>

#include <stdio.h>
#include <stddef.h>

void tc_scheduling_sfs_tiles(struct tc_scop* scop, struct tc_options* options, __isl_take isl_union_set* LD, __isl_take isl_union_map* S, __isl_take isl_union_map* R, __isl_take isl_set* ii_set, __isl_take isl_set* tile, __isl_take isl_map* Rtile, __isl_take isl_id_list* II, __isl_take isl_id_list* I)
{
    isl_ctx* ctx = isl_union_set_get_ctx(LD);

    if (!tc_is_lex_forward(Rtile))
    {
        tc_error("Backward relation detected.");
        return;
    }
    
    isl_set* Suds = tc_uds_set(Rtile);
    
    tc_debug_set(Suds, "UDS");
    
    isl_set* ind = tc_ind_set(ii_set, Rtile);
    
    tc_debug_set(ind, "IND");
    
    isl_map* Rusc = tc_Rusc_map(Rtile, S);
    
    tc_debug_map(Rusc, "R_USC");
        
    isl_set* repr = isl_set_subtract(isl_set_copy(Suds), isl_map_range(isl_map_copy(Rusc)));
    
    tc_debug_set(repr, "REPR");
    
    tc_debug_set_card(repr, "REPR");
    
    isl_set* repr_ind = isl_set_union(repr, ind);
    
    tc_debug_set(repr_ind, "REPR_IND");
    
    tc_debug_set_card(repr_ind, "REPR_IND");
    
    repr_ind = tc_parameterize_all(repr_ind, II);
        
    int exact;
    isl_map* Rtile_plus = tc_transitive_closure(isl_map_copy(Rtile), S, &exact);
    isl_map* Rtile_star = isl_map_union(isl_map_copy(Rtile_plus), tc_make_identity(isl_map_copy(Rtile)));
    
    tc_debug_map(Rtile_star, "R_TILE* (exact=%d)", exact);
                    
    isl_map* Rusc_plus = tc_transitive_closure(isl_map_copy(Rusc), S, &exact);
    isl_map* Rusc_star = isl_map_union(Rusc_plus, tc_make_identity(isl_map_copy(Rusc)));
    
    tc_debug_map(Rusc_star, "R_USC* (exact=%d)", exact);

    // S_slice = R*((R_USC*)(e))
    isl_set* Sslice = isl_set_apply(isl_set_apply(isl_set_copy(repr_ind), Rusc_star), isl_map_copy(Rtile_star));
    
    tc_debug_set(Sslice, "SLICE");
    
    Sslice = tc_lift_up_set_params(Sslice, II);
    
    tile = tc_lift_up_set_params(tile, II);
    
    Sslice = isl_set_add_dims(Sslice, isl_dim_set, isl_id_list_n_id(II));
    tile = isl_set_insert_dims(tile, isl_dim_set, 0, isl_id_list_n_id(II));
    
    isl_set* tile_ext = isl_set_intersect(isl_set_copy(Sslice), isl_set_copy(tile));
    
    tile_ext = isl_set_coalesce(tile_ext);
    
    isl_union_map* S_ext = tc_extend_schedule(isl_union_map_copy(S), 2 * isl_id_list_n_id(II));
    
    isl_id_list* IR = tc_ids_sequence(ctx, "ir", isl_id_list_n_id(II));
    isl_id_list* iterators = isl_id_list_concat(isl_id_list_concat(isl_id_list_copy(IR), isl_id_list_copy(II)), isl_id_list_copy(I));
                
    enum tc_codegen_enum codegen = tc_options_codegen(options);
    
    if (tc_codegen_enum_serial == codegen)
    {
        tc_codegen_serial(scop, options, S_ext, tile_ext, iterators);
    }
    else if (tc_codegen_enum_omp_cpu_for == codegen)
    {
        tc_codegen_omp_parallel_for(scop, options, S_ext, tile_ext, iterators, IR, 0);
    }
    else if (tc_codegen_enum_omp_cpu_task == codegen)
    {
        tc_codegen_omp_task_for(scop, options, S_ext, tile_ext, iterators, IR, 0);
    }
    else if (tc_codegen_enum_omp_gpu == codegen)
    {
        tc_codegen_omp_gpu(scop, options, S_ext, tile_ext, iterators, IR);
    }

    isl_map_free(Rtile);
    isl_map_free(Rtile_plus);
    isl_map_free(Rtile_star);
    isl_map_free(Rusc);
    isl_set_free(Suds);
    isl_set_free(repr_ind);
    isl_set_free(tile);
    isl_set_free(ii_set);
    isl_set_free(Sslice);
    isl_union_set_free(LD);
    isl_union_map_free(S);
    isl_union_map_free(R);
    isl_id_list_free(IR);
    isl_id_list_free(iterators);
    isl_id_list_free(II);
    isl_id_list_free(I);
}

void tc_scheduling_sfs_single(struct tc_scop* scop, struct tc_options* options, __isl_take isl_union_set* LD, __isl_take isl_union_map* S, __isl_take isl_union_map* R, __isl_take isl_set* ii_set, __isl_take isl_set* tile_vld, __isl_take isl_map* Rtile, __isl_take isl_id_list* II, __isl_take isl_id_list* I)
{
    isl_ctx* ctx = isl_union_set_get_ctx(LD);
    
    if (!tc_is_lex_forward(Rtile))
    {
        tc_error("Backward relation detected.");
        return;
    }
        
    isl_set* LD_normalized = tc_normalize_union_set(LD, S);
    
    isl_map* R_normalized = tc_normalize_union_map(R, S);
    
    int exact;
    
    isl_map* R_plus_normalized = tc_transitive_closure(isl_map_copy(R_normalized), S, &exact);
    
    isl_map* R_star_normalized = isl_map_union(isl_map_copy(R_plus_normalized), tc_make_identity(isl_map_copy(R_normalized)));
        
    tc_debug_map(R_star_normalized, "R* (exact=%d)", exact);
    
    isl_id_list* IR = tc_ids_sequence(ctx, "ir", isl_set_n_dim(tile_vld));
    isl_id_list* IIprim = tc_ids_prim(II);
    
    isl_set* uds = tc_uds_set(R_normalized);
    
    tc_debug_set(uds, "UDS");
            
    isl_map* Rusc = tc_Rusc_map(R_normalized, S);
    
    tc_debug_map(Rusc, "R_USC");
                    
    isl_set* repr = isl_set_subtract(uds, isl_map_range(isl_map_copy(Rusc)));
    
    tc_debug_set(repr, "REPR");
    
    tc_debug_set_card(repr, "REPR");
    
    isl_set* ind = tc_ind_set(LD_normalized, R_normalized);
    
    tc_debug_set(ind, "IND");
    
    isl_set* repr_ind = isl_set_union(repr, ind);
    
    tc_debug_set(repr_ind, "REPR_IND");
    
    tc_debug_set_card(repr_ind, "REPR_IND");
    
    isl_set* Srepr_ind = tc_parameterize_all(repr_ind, IR);
    
    isl_set* tile_repr_ind = isl_set_intersect(isl_set_copy(tile_vld), Srepr_ind);
    tile_repr_ind = isl_set_coalesce(tile_repr_ind);    
    
    tc_debug_set(tile_repr_ind, "TILE_REPR_IND");
    
    isl_map* Rusc_plus = tc_transitive_closure(isl_map_copy(Rusc), S, &exact);
    isl_map* Rusc_star = isl_map_union(Rusc_plus, tc_make_identity(isl_map_copy(Rusc)));
    Rusc_star = isl_map_coalesce(Rusc_star);
    
    tc_debug_map(Rusc_star, "R_USC* (exact=%d)", exact);

    // SFS := R^∗(R_USC^∗(TILE_REPR_IND))
    isl_set* sfs = isl_set_apply(isl_set_apply(tile_repr_ind, Rusc_star), R_star_normalized);
    sfs = isl_set_coalesce(sfs);
    
    tc_debug_set(sfs, "SFS");
        
    // II_SFS := [IR,II] -> { [II'] : II' in II_SET AND exists I s.t. I in TILE_VLD(II') AND I in SFS(IR,II) }
        
    isl_id_list* IR_II = isl_id_list_concat(isl_id_list_copy(IR), isl_id_list_copy(II));
    isl_id_list* IIprim_I = isl_id_list_concat(isl_id_list_copy(IIprim), isl_id_list_copy(I));
    isl_id_list* IR_II_IIprim_I = isl_id_list_concat(isl_id_list_copy(IR_II), isl_id_list_copy(IIprim_I));
    
    isl_set* ii_sfs = tc_make_set(ctx, IR_II_IIprim_I, IIprim, "");
    
    isl_set* tile_vld_constraints_prim = tc_rename_params(isl_set_copy(tile_vld), II, IIprim);
    tile_vld_constraints_prim = tc_make_set_constraints(tile_vld_constraints_prim, I);
    
    isl_set* sfs_constraints = tc_make_set_constraints(isl_set_copy(sfs), I);
    
    isl_set* ii_set_constraints_prim = tc_make_set_constraints(isl_set_copy(ii_set), IIprim);
    
    ii_sfs = isl_set_intersect_params(ii_sfs, isl_set_copy(tile_vld_constraints_prim));
    ii_sfs = isl_set_intersect_params(ii_sfs, isl_set_copy(sfs_constraints));
    ii_sfs = isl_set_intersect_params(ii_sfs, ii_set_constraints_prim);
    
    ii_sfs = tc_project_out_params(ii_sfs, IIprim_I);
    ii_sfs = isl_set_coalesce(ii_sfs);

    tc_debug_set(ii_sfs, "II_SFS");
    
    //tc_debug_set_card(ii_sfs, "II_SFS");
        
    // S := { [IR, II, II', I] : [IR,II] in params SFS AND II' in II_SFS([IR,II]) AND I in TILE_VLD(II') AND I in SFS([IR,II]) }
    
    isl_set* slice = tc_make_params(ctx, IR_II_IIprim_I, "");
    
    isl_set* sfs_params = isl_set_params(isl_set_copy(sfs));
    
    tc_debug_set(sfs_params, "params SFS");
        
    isl_set* ii_sfs_constraints = tc_make_set_constraints(ii_sfs, IIprim);
        
    slice = isl_set_intersect_params(slice, sfs_params);
    slice = isl_set_intersect_params(slice, sfs_constraints);
    slice = isl_set_intersect_params(slice, ii_sfs_constraints);
    slice = isl_set_intersect_params(slice, tile_vld_constraints_prim);
        
    slice = tc_lift_up_set_params(slice, I);
    slice = tc_lift_up_set_params(slice, IIprim);
    slice = tc_lift_up_set_params(slice, II);
    slice = tc_lift_up_set_params(slice, IR);
    
    slice = isl_set_coalesce(slice);
    
    tc_debug_set(slice, "SLICE");
        
    isl_union_map* S_ext = tc_extend_schedule(isl_union_map_copy(S), isl_id_list_n_id(IR) + isl_id_list_n_id(II) + isl_id_list_n_id(IIprim));
    
    isl_id_list* iterators = isl_id_list_concat(isl_id_list_concat(isl_id_list_copy(IR_II), isl_id_list_copy(IIprim)), isl_id_list_copy(I));
    
    enum tc_codegen_enum codegen = tc_options_codegen(options);
    
    if (tc_codegen_enum_serial == codegen)
    {
        tc_codegen_serial(scop, options, S_ext, slice, iterators);
    }
    else if (tc_codegen_enum_omp_cpu_for == codegen)
    {
        tc_codegen_omp_parallel_for(scop, options, S_ext, slice, iterators, IR_II, 0);
    }
    else if (tc_codegen_enum_omp_cpu_task == codegen)
    {
        tc_codegen_omp_task_for(scop, options, S_ext, slice, iterators, IR_II, 1);
    }
    else if (tc_codegen_enum_omp_gpu == codegen)
    {
        tc_codegen_omp_gpu(scop, options, S_ext, slice, iterators, IR_II);
    }
            
    isl_id_list_free(I);
    isl_id_list_free(II);
    isl_id_list_free(IIprim);
    isl_id_list_free(IR_II);
    isl_id_list_free(iterators);
    isl_id_list_free(IIprim_I);
    isl_id_list_free(IR_II_IIprim_I);
    isl_map_free(Rusc);
    isl_set_free(LD_normalized);
    isl_map_free(R_normalized);
    isl_map_free(R_plus_normalized);
    isl_set_free(tile_vld);
    isl_set_free(ii_set);
    isl_set_free(sfs);
    isl_union_set_free(LD);
    isl_union_map_free(S);
    isl_union_map_free(R);
    isl_id_list_free(IR);
    isl_map_free(Rtile);
}

void tc_scheduling_sfs_multiple(struct tc_scop* scop, struct tc_options* options, __isl_take isl_union_set* LD, __isl_take isl_union_map* S, __isl_take isl_union_map* R, __isl_take isl_set* ii_set, __isl_take isl_set* tile_vld, __isl_take isl_map* Rtile, __isl_take isl_id_list* II, __isl_take isl_id_list* I)
{
    isl_ctx* ctx = isl_union_set_get_ctx(LD);
    
    if (!tc_is_lex_forward(Rtile))
    {
        tc_error("Backward relation detected.");
        return;
    }
        
    isl_set* LD_normalized = tc_normalize_union_set(LD, S);
    
    isl_map* R_normalized = tc_normalize_union_map(R, S);    
    
    int exact;
    
    isl_map* R_plus_normalized = tc_transitive_closure(isl_map_copy(R_normalized), S, &exact);
    
    isl_map* R_star_normalized = isl_map_union(isl_map_copy(R_plus_normalized), tc_make_identity(isl_map_copy(R_normalized)));
        
    tc_debug_map(R_star_normalized, "R* (exact=%d)", exact);
    
    isl_id_list* IIprim = tc_ids_prim(II);
    
    isl_set* uds = tc_uds_set(R_normalized);
    
    tc_debug_set(uds, "UDS");
    
    isl_map* Rusc = tc_Rusc_map(R_normalized, S);
    
    tc_debug_map(Rusc, "R_USC");
            
    isl_set* repr = isl_set_subtract(uds, isl_map_range(isl_map_copy(Rusc)));
    
    tc_debug_set(repr, "REPR");
    
    tc_debug_set_card(repr, "REPR");
    
    isl_set* ind = tc_ind_set(LD_normalized, R_normalized);
    
    tc_debug_set(ind, "IND");
    
    isl_set* repr_ind = isl_set_union(repr, ind);
    
    tc_debug_set(repr_ind, "REPR_IND");
    
    tc_debug_set_card(repr_ind, "REPR_IND");
        
    isl_set* tile_repr_ind = isl_set_intersect(isl_set_copy(tile_vld), repr_ind);
    tile_repr_ind = isl_set_coalesce(tile_repr_ind);
    
    tc_debug_set(tile_repr_ind, "TILE_REPR_IND");
    
    if (tc_debug_flag)
    {
        isl_set* tile_repr_ind_params = tc_get_params_set(isl_set_copy(tile_repr_ind), II);
        tc_debug_set_card(tile_repr_ind_params, "params TILE_REPR_IND");
        isl_set_free(tile_repr_ind_params);
    }
    
    isl_map* Rusc_plus = tc_transitive_closure(isl_map_copy(Rusc), S, &exact);
    isl_map* Rusc_star = isl_map_union(Rusc_plus, tc_make_identity(isl_map_copy(Rusc)));
    
    tc_debug_map(Rusc_star, "R_USC* (exact=%d)", exact);

    // SFS := R^∗(R_USC^∗(TILE_REPR_IND))
    isl_set* sfs = isl_set_apply(isl_set_apply(tile_repr_ind, Rusc_star), R_star_normalized);
    sfs = isl_set_coalesce(sfs);
    
    tc_debug_set(sfs, "SFS");
        
    // II_SFS := [II] -> { [II'] : II' in II_SET AND exists I s.t. I in TILE_VLD(II') AND I in SFS(II) }
        
    isl_id_list* IIprim_I = isl_id_list_concat(isl_id_list_copy(IIprim), isl_id_list_copy(I));
    isl_id_list* II_IIprim_I = isl_id_list_concat(isl_id_list_copy(II), isl_id_list_copy(IIprim_I));
    
    isl_set* ii_sfs = tc_make_set(ctx, II_IIprim_I, IIprim, "");
    
    isl_set* tile_vld_constraints_prim = tc_rename_params(isl_set_copy(tile_vld), II, IIprim);
    tile_vld_constraints_prim = tc_make_set_constraints(tile_vld_constraints_prim, I);
    
    isl_set* sfs_constraints = tc_make_set_constraints(isl_set_copy(sfs), I);
    
    isl_set* ii_set_constraints_prim = tc_make_set_constraints(isl_set_copy(ii_set), IIprim);
    
    ii_sfs = isl_set_intersect_params(ii_sfs, isl_set_copy(tile_vld_constraints_prim));
    ii_sfs = isl_set_intersect_params(ii_sfs, isl_set_copy(sfs_constraints));
    ii_sfs = isl_set_intersect_params(ii_sfs, ii_set_constraints_prim);
    
    ii_sfs = tc_project_out_params(ii_sfs, IIprim_I);
    ii_sfs = isl_set_coalesce(ii_sfs);
    
    tc_debug_set(ii_sfs, "II_SFS");
    
    //tc_debug_set_card(ii_sfs, "II_SFS");
        
    // S := { [II, II', I] : [II] in params SFS AND II' in II_SFS([II]) AND I in TILE_VLD(II') AND I in SFS([II]) }
    
    isl_set* slice = tc_make_params(ctx, II_IIprim_I, "");
    
    isl_set* sfs_params = isl_set_params(isl_set_copy(sfs));
    
    tc_debug_set(sfs_params, "params SFS");
        
    isl_set* ii_sfs_constraints = tc_make_set_constraints(ii_sfs, IIprim);
        
    slice = isl_set_intersect_params(slice, sfs_params);
    slice = isl_set_intersect_params(slice, sfs_constraints);
    slice = isl_set_intersect_params(slice, ii_sfs_constraints);
    slice = isl_set_intersect_params(slice, tile_vld_constraints_prim);
        
    slice = tc_lift_up_set_params(slice, I);
    slice = tc_lift_up_set_params(slice, IIprim);
    slice = tc_lift_up_set_params(slice, II);
    
    slice = isl_set_coalesce(slice);
    
    tc_debug_set(slice, "SLICE");
        
    isl_union_map* S_ext = tc_extend_schedule(isl_union_map_copy(S), isl_id_list_n_id(II) + isl_id_list_n_id(IIprim));
    
    isl_id_list* iterators = isl_id_list_concat(isl_id_list_concat(isl_id_list_copy(II), isl_id_list_copy(IIprim)), isl_id_list_copy(I));
            
    enum tc_codegen_enum codegen = tc_options_codegen(options);
    
    if (tc_codegen_enum_serial == codegen)
    {
        tc_codegen_serial(scop, options, S_ext, slice, iterators);
    }
    else if (tc_codegen_enum_omp_cpu_for == codegen)
    {
        tc_codegen_omp_parallel_for(scop, options, S_ext, slice, iterators, II, 0);
    }
    else if (tc_codegen_enum_omp_cpu_task == codegen)
    {
        tc_codegen_omp_task_for(scop, options, S_ext, slice, iterators, II, 1);
    }
    else if (tc_codegen_enum_omp_gpu == codegen)
    {
        tc_codegen_omp_gpu(scop, options, S_ext, slice, iterators, II);
    }

    isl_id_list_free(I);
    isl_id_list_free(II);
    isl_id_list_free(IIprim);
    isl_id_list_free(IIprim_I);
    isl_id_list_free(II_IIprim_I);
    isl_id_list_free(iterators);
    isl_map_free(Rusc);
    isl_set_free(LD_normalized);
    isl_map_free(R_normalized);
    isl_map_free(R_plus_normalized);
    isl_set_free(tile_vld);
    isl_set_free(ii_set);
    isl_set_free(sfs);
    isl_union_set_free(LD);
    isl_union_map_free(S);
    isl_union_map_free(R);
    isl_map_free(Rtile);
}
