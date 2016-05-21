#include "sfs_scheduling.h"
#include "omp_cpu_codegen.h"
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

void tc_scheduling_sfs_tiles(struct tc_scop* scop, struct tc_options* options, __isl_take isl_union_set* LD, __isl_take isl_union_map* S, __isl_take isl_union_map* R, __isl_take isl_set* ii_set, __isl_take isl_set* tile, __isl_take isl_map* Rtile, __isl_take isl_id_list* II)
{
    isl_ctx* ctx = isl_union_set_get_ctx(LD);

    if (!tc_is_lex_forward(Rtile))
    {
        tc_options_error("Backward relation detected");
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

    isl_id_list* IR = tc_ids_sequence(ctx, "ir", isl_id_list_n_id(II));
    isl_id_list* iterators = isl_id_list_concat(isl_id_list_copy(IR), isl_id_list_copy(II));
    
    isl_union_map* S_ext = tc_extend_schedule(isl_union_map_copy(S), 2 * isl_id_list_n_id(II));
                
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
}

