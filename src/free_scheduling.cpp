#include "free_scheduling.h"
#include "omp_cpu_codegen.h"
#include "slicing.h"
#include "utility.h"
#include "scop.h"
#include "serial_codegen.h"
#include "debug.h"
#include "transitive_closure.h"

#include <isl/ctx.h>
#include <isl/id.h>
#include <isl/set.h>
#include <isl/map.h>
#include <isl/union_map.h>
#include <isl/union_set.h>

#include <stdio.h>
#include <stddef.h>

void tc_scheduling_free_schedule_rk(struct tc_scop* scop, struct tc_options* options, __isl_take isl_union_set* LD, __isl_take isl_union_map* S, __isl_take isl_union_map* R, __isl_take isl_set* ii_set, __isl_take isl_set* tile, __isl_take isl_map* Rtile, __isl_take isl_id_list* II, __isl_take isl_id_list* I)
{
    int exact;    
    
    isl_map* Rtile_plus = tc_transitive_closure(isl_map_copy(Rtile), S, &exact);
    Rtile_plus = isl_map_coalesce(Rtile_plus);
    tc_debug_map(Rtile_plus, "R_TILE_M^+ (exact=%d)", exact);
        
    isl_id* k_param;
    isl_set* ii_set_k = tc_Sk_set(ii_set, Rtile, Rtile_plus, S, &k_param);
    
    isl_set* ii_set_k_ext = isl_set_add_dims(ii_set_k, isl_dim_set, isl_set_n_dim(tile));
    
    isl_set* tile_ext = tc_lift_up_set_params(tile, II);        
    tile_ext = isl_set_intersect(tile_ext, ii_set_k_ext);
    tile_ext = isl_set_coalesce(tile_ext);
    
    tc_debug_set(tile_ext, "TILE_M_K_EXT");

    isl_id_list* k_param_list = isl_id_list_from_id(isl_id_copy(k_param));
    tile_ext = tc_lift_up_set_params(tile_ext, k_param_list);
    tile_ext = isl_set_coalesce(tile_ext);
    
    isl_union_map* S_ext = tc_extend_schedule(isl_union_map_copy(S), isl_id_list_n_id(II) + 1);
    
    isl_id_list* iterators = isl_id_list_concat(isl_id_list_concat(isl_id_list_copy(k_param_list), isl_id_list_copy(II)), isl_id_list_copy(I));
    
    enum tc_codegen_enum codegen = tc_options_codegen(options);
    
    if (tc_codegen_enum_serial == codegen)
    {
        tc_codegen_serial(scop, options, S_ext, tile_ext, iterators);
    }
    else if (tc_codegen_enum_omp_cpu_for == codegen)
    {
        tc_codegen_omp_parallel_for(scop, options, S_ext, tile_ext, iterators, II, 0);
    }
    else if (tc_codegen_enum_omp_cpu_task == codegen)
    {
        tc_codegen_omp_task_for(scop, options, S_ext, tile_ext, iterators, II, 0);
    }
    
    isl_id_list_free(k_param_list);
    isl_union_set_free(LD);
    isl_union_map_free(S);
    isl_union_map_free(R);
    isl_id_free(k_param);
    isl_map_free(Rtile);
    isl_map_free(Rtile_plus);
    isl_set_free(ii_set);
    isl_id_list_free(I);
    isl_id_list_free(II);
    isl_id_list_free(iterators);
}

void tc_scheduling_free_schedule_karl(struct tc_scop* scop, struct tc_options* options, __isl_take isl_union_set* LD, __isl_take isl_union_map* S, __isl_take isl_union_map* R, __isl_take isl_set* ii_set, __isl_take isl_set* tile, __isl_take isl_map* Rtile, __isl_take isl_id_list* II, __isl_take isl_id_list* I)
{
    isl_ctx* ctx = isl_union_set_get_ctx(LD);
    
    isl_set* LD_normalized = tc_normalize_union_set(LD, S);
    
    isl_set* FS = tc_FS_set(LD_normalized, Rtile, S);
    
    isl_set_free(LD_normalized);
    
    FS = isl_set_add_dims(FS, isl_dim_set, isl_set_n_dim(tile));
    
    isl_id* k_param = isl_id_alloc(ctx, "k", NULL);
    isl_id_list* k_param_list = isl_id_list_from_id(isl_id_copy(k_param));
    
    isl_set* tile_ext = tc_lift_up_set_params(tile, II);    
    tile_ext = isl_set_insert_dims(tile_ext, isl_dim_set, 0, 1);
    tile_ext = isl_set_coalesce(tile_ext);
        
    tile_ext = isl_set_intersect(tile_ext, FS);
                    
    isl_union_map* S_ext = tc_extend_schedule(isl_union_map_copy(S), isl_id_list_n_id(II) + 1);
    
    isl_id_list* iterators = isl_id_list_concat(isl_id_list_concat(isl_id_list_copy(k_param_list), isl_id_list_copy(II)), isl_id_list_copy(I));
    
    enum tc_codegen_enum codegen = tc_options_codegen(options);
    
    if (tc_codegen_enum_serial == codegen)
    {
        tc_codegen_serial(scop, options, S_ext, tile_ext, iterators);
    }
    else if (tc_codegen_enum_omp_cpu_for == codegen)
    {
        tc_codegen_omp_parallel_for(scop, options, S_ext, tile_ext, iterators, II, 0);
    }
    else if (tc_codegen_enum_omp_cpu_task == codegen)
    {
        tc_codegen_omp_task_for(scop, options, S_ext, tile_ext, iterators, II, 0);
    }
    
    isl_id_list_free(k_param_list);
    isl_union_set_free(LD);
    isl_union_map_free(S);
    isl_union_map_free(R);
    isl_id_free(k_param);
    isl_map_free(Rtile);
    isl_set_free(ii_set);
    isl_id_list_free(I);
    isl_id_list_free(II);
    isl_id_list_free(iterators);
}

void tc_scheduling_free_schedule_finite(struct tc_scop* scop, struct tc_options* options, __isl_take isl_union_set* LD, __isl_take isl_union_map* S, __isl_take isl_union_map* R, __isl_take isl_set* ii_set, __isl_take isl_set* tile, __isl_take isl_map* Rtile, __isl_take isl_id_list* II, __isl_take isl_id_list* I)
{
    isl_ctx* ctx = isl_union_set_get_ctx(LD);
            
    isl_id* k_param = isl_id_alloc(ctx, "k", NULL);
    isl_id_list* k_param_list = isl_id_list_from_id(isl_id_copy(k_param));
    
    isl_set* ind = tc_ind_set(ii_set, Rtile);    
    isl_set* ind_parameterized = isl_set_insert_dims(ind, isl_dim_param, 0, 1);
    ind_parameterized = isl_set_set_dim_id(ind_parameterized, isl_dim_param, 0, isl_id_copy(k_param));
    ind_parameterized = isl_set_fix_si(ind_parameterized, isl_dim_param, 0, 0);
    
    isl_set* FS = ind_parameterized;
    
    int k = 0;
    isl_set* uds = tc_uds_set(Rtile);
    while (!isl_set_is_empty(uds))
    {   
        isl_set* uds_parameterized = isl_set_insert_dims(isl_set_copy(uds), isl_dim_param, 0, 1);
        uds_parameterized = isl_set_set_dim_id(uds_parameterized, isl_dim_param, 0, isl_id_copy(k_param));
        uds_parameterized = isl_set_fix_si(uds_parameterized, isl_dim_param, 0, k);
        
        FS = isl_set_union(FS, uds_parameterized);
        FS = isl_set_coalesce(FS);
        
        ii_set = isl_set_subtract(ii_set, isl_set_copy(uds));
        Rtile = isl_map_subtract_domain(Rtile, uds);
        uds = tc_uds_set(Rtile);
        
        k = k + 1;
    }
    
    isl_set_free(uds);
    
    isl_set* ii_set_parameterized = isl_set_insert_dims(isl_set_copy(ii_set), isl_dim_param, 0, 1);
    ii_set_parameterized = isl_set_set_dim_id(ii_set_parameterized, isl_dim_param, 0, isl_id_copy(k_param));
    ii_set_parameterized = isl_set_fix_si(ii_set_parameterized, isl_dim_param, 0, k);
    
    FS = isl_set_union(FS, ii_set_parameterized);
    FS = isl_set_coalesce(FS);
    
    FS = isl_set_add_dims(FS, isl_dim_set, isl_set_n_dim(tile));
    
    FS = tc_lift_up_set_params(FS, k_param_list);
    
    tc_debug_set(FS, "FS");
        
    isl_set* tile_ext = tc_lift_up_set_params(tile, II);    
    tile_ext = isl_set_insert_dims(tile_ext, isl_dim_set, 0, 1);
    tile_ext = isl_set_coalesce(tile_ext);
            
    tile_ext = isl_set_intersect(tile_ext, FS);
    tile_ext = isl_set_coalesce(tile_ext);
    
    tc_debug_set(tile_ext, "TILE_EXT");
    
    isl_union_map* S_ext = tc_extend_schedule(isl_union_map_copy(S), isl_id_list_n_id(II) + 1);
    
    isl_id_list* iterators = isl_id_list_concat(isl_id_list_concat(isl_id_list_copy(k_param_list), isl_id_list_copy(II)), isl_id_list_copy(I));
    
    enum tc_codegen_enum codegen = tc_options_codegen(options);
    
    if (tc_codegen_enum_serial == codegen)
    {
        tc_codegen_serial(scop, options, S_ext, tile_ext, iterators);
    }
    else if (tc_codegen_enum_omp_cpu_for == codegen)
    {
        tc_codegen_omp_parallel_for(scop, options, S_ext, tile_ext, iterators, II, 0);
    }
    else if (tc_codegen_enum_omp_cpu_task == codegen)
    {
        tc_codegen_omp_task_for(scop, options, S_ext, tile_ext, iterators, II, 0);
    }
    
    isl_id_list_free(k_param_list);
    isl_union_set_free(LD);
    isl_union_map_free(S);
    isl_union_map_free(R);
    isl_id_free(k_param);
    isl_map_free(Rtile);
    isl_set_free(ii_set);
    isl_id_list_free(I);
    isl_id_list_free(II);
    isl_id_list_free(iterators);
}
