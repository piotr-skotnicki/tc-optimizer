#include "stencil_tiling.h"
#include "scop.h"
#include "tiling.h"
#include "for_decorator.h"
#include "utility.h"
#include "options.h"
#include "omp_cpu_codegen.h"
#include "serial_codegen.h"
#include "debug.h"
#include "transitive_closure.h"

#include <isl/id.h>

#include <string.h>
#include <stddef.h>

#include <vector>
#include <map>
#include <string>

void tc_algorithm_stencil_tiling(struct tc_scop* scop, struct tc_options* options)
{
    isl_ctx* ctx = scop->ctx;
    
    isl_union_set* LD = isl_union_set_copy(scop->domain);
    tc_debug_uset(LD, "LD");
    
    isl_union_map* S = isl_union_map_copy(scop->schedule);
    tc_debug_umap(S, "S");
    
    isl_union_map* R = isl_union_map_copy(scop->relation);
    tc_debug_umap(R, "R");
    
    isl_basic_set* sample = isl_set_sample(tc_normalize_union_set(LD, S));
    
    isl_space* space = isl_basic_set_get_space(sample);
            
    isl_id_list* I = tc_ids_sequence(ctx, "i", isl_space_dim(space, isl_dim_set));
    isl_id_list* II = tc_ids_sequence(ctx, "ii", isl_space_dim(space, isl_dim_set));
    
    isl_id_list* Iprim = tc_ids_prim(I);
    isl_id_list* Ibis = tc_ids_bis(I);
    
    isl_id_list* IIprim = tc_ids_prim(II);
    isl_id_list* IIbis = tc_ids_bis(II);
    
    isl_map* R_normalized = tc_normalize_union_map(R, S);

    isl_bool exact = isl_bool_false;
    isl_map* R_plus_normalized = tc_transitive_closure(R_normalized, S, &exact);
    
    tc_debug_map(R_plus_normalized, "R^+ (exact=%d)", exact);
        
    std::map<std::string, std::vector<int> > blocks = tc_options_blocks(options);
    
    isl_set* tile;
    isl_set* ii_set;
    
    tc_tile_loop_nest(LD, S, II, I, &tile, &ii_set, blocks);
    
    tc_debug_set(tile, "TILE");
    tc_debug_set(ii_set, "II_SET");
        
    int k = 0;
    
    isl_set* subtiles_ext = NULL;
        
    while (!isl_set_is_empty(tile))
    {
        k = k + 1;
        
        isl_set* tile_lt = tc_tile_set_of(tile, ii_set, II, &tc_tuples_tlt);
        
        isl_set* tile_gt = tc_tile_set_of(tile, ii_set, II, &tc_tuples_tgt);
        
        isl_set* subtile = isl_set_copy(tile);
        
        subtile = isl_set_subtract(subtile, isl_set_apply(tile_lt, isl_map_copy(R_plus_normalized)));
        
        subtile = isl_set_subtract(subtile, isl_set_apply(tile_gt, isl_map_copy(R_plus_normalized)));
        
        tile = isl_set_subtract(tile, isl_set_copy(subtile));
        
        tc_debug_set(tile, "TILE_%d", k);
        
        // Extend the tuple of set II: [ii1, ii2, ..., iid] to the tuple [ii1, k, ii2, ...,iid]; 
        isl_set* subtile_ext = tc_lift_up_set_params(subtile, II);
        
        subtile_ext = isl_set_insert_dims(subtile_ext, isl_dim_set, 1, 1);
        
        subtile_ext = isl_set_fix_si(subtile_ext, isl_dim_set, 1, k - 1);
        
        if (NULL == subtiles_ext)
        {
            subtiles_ext = subtile_ext;
        }
        else
        {
            subtiles_ext = isl_set_union(subtiles_ext, subtile_ext);
        }
    }
    
    subtiles_ext = isl_set_coalesce(subtiles_ext);
    
    isl_union_map* S_ext = tc_extend_schedule(isl_union_map_copy(S), isl_id_list_n_id(II) + 1);
                            
    isl_id_list* iterators = isl_id_list_copy(II);
    iterators = isl_id_list_concat(iterators, isl_id_list_copy(I));
    iterators = isl_id_list_insert(iterators, 1, isl_id_alloc(ctx, "k", NULL));
    
    isl_id_list* parallel_iterators = isl_id_list_drop(isl_id_list_copy(II), 0, 1);
        
    enum tc_codegen_enum codegen = tc_options_codegen(options);
    
    if (tc_codegen_enum_serial == codegen)
    {
        tc_codegen_serial(scop, options, S_ext, subtiles_ext, iterators);
    }
    else if (tc_codegen_enum_omp_cpu_for == codegen)
    {
        tc_codegen_omp_parallel_for(scop, options, S_ext, subtiles_ext, iterators, parallel_iterators, 0);
    }
    else if (tc_codegen_enum_omp_cpu_task == codegen)
    {
        tc_codegen_omp_task_for(scop, options, S_ext, subtiles_ext, iterators, parallel_iterators, 0);
    }
    
    isl_set_free(ii_set);
    isl_set_free(tile);
    isl_basic_set_free(sample);
    isl_map_free(R_plus_normalized);
    isl_id_list_free(I);
    isl_id_list_free(Iprim);
    isl_id_list_free(Ibis);
    isl_id_list_free(II);
    isl_id_list_free(IIprim);
    isl_id_list_free(IIbis);
    isl_space_free(space);
    isl_union_set_free(LD);
    isl_union_map_free(S);
    isl_union_map_free(R);
    isl_id_list_free(iterators);
    isl_id_list_free(parallel_iterators);
}
