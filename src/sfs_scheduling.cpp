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
#include "input_output.h"

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
        tc_die(tc_exit_code_nonlex);
    }

    isl_set* uds = tc_uds_set(Rtile);
    
    tc_debug_set(uds, "UDS");

    isl_map* Rusc = tc_Rusc_map(Rtile, S);
    
    tc_debug_map(Rusc, "R_USC");

    isl_set* repr = isl_set_subtract(uds, isl_map_range(isl_map_copy(Rusc)));

    tc_debug_set(repr, "REPR");

    tc_debug_set_card(repr, "REPR");
    
    isl_set* ind = tc_ind_set(ii_set, Rtile);

    tc_debug_set(ind, "IND");

    isl_set* repr_ind = isl_set_union(repr, ind);
    
    tc_debug_set(repr_ind, "REPR_IND");
    
    tc_debug_set_card(repr_ind, "REPR_IND");

    isl_map* repr_ind_map = tc_make_identity_from_set(isl_set_copy(repr_ind));

    repr_ind = tc_parameterize_all(repr_ind, II);

    isl_bool exact = isl_bool_false;
    isl_map* Rtile_plus = tc_transitive_closure(isl_map_copy(Rtile), S, &exact);
    isl_map* Rtile_star = isl_map_union(isl_map_copy(Rtile_plus), tc_make_identity(isl_map_copy(Rtile)));
    
    tc_debug_map(Rtile_star, "R_TILE* (exact=%d)", exact);

    if (exact != isl_bool_true)
    {
        tc_warn("Inexact R_TILE*. The results can be non-optimal. Restart TC with a different transitive closure method.");
    }

    isl_map* Rusc_star = NULL;

    // SFS = R*(REPR_IND)
    if (isl_map_is_empty(Rusc))
    {
        Rusc_star = repr_ind_map;
    }
    else
    {
        isl_map* Rusc_plus = tc_transitive_closure(isl_map_copy(Rusc), S, &exact);

        Rusc_star = isl_map_union(Rusc_plus, tc_make_identity(isl_map_copy(Rusc)));
        Rusc_star = isl_map_coalesce(Rusc_star);

        tc_debug_map(Rusc_star, "R_USC* (exact=%d)", exact);

        if (exact != isl_bool_true)
        {
            tc_warn("Inexact R_USC*. The results can be non-optimal. Restart TC with a different transitive closure method.");
        }

        Rusc_star = isl_map_union(Rusc_star, repr_ind_map);
    }

    // SFS = R*(R_USC*(REPR_IND))
    isl_set* sfs = isl_set_apply(isl_set_apply(repr_ind, Rusc_star), isl_map_copy(Rtile_star));
    sfs = isl_set_coalesce(sfs);

    tc_debug_set(sfs, "SFS");
    
    sfs = tc_lift_up_set_params(sfs, II);
    
    tile = tc_lift_up_set_params(tile, II);
    
    sfs = isl_set_add_dims(sfs, isl_dim_set, isl_id_list_n_id(II));
    tile = isl_set_insert_dims(tile, isl_dim_set, 0, isl_id_list_n_id(II));
    
    isl_set* tile_ext = isl_set_intersect(isl_set_copy(sfs), isl_set_copy(tile));
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
    isl_set_free(tile);
    isl_set_free(ii_set);
    isl_set_free(sfs);
    isl_union_set_free(LD);
    isl_union_map_free(S);
    isl_union_map_free(R);
    isl_id_list_free(IR);
    isl_id_list_free(iterators);
    isl_id_list_free(II);
    isl_id_list_free(I);
}

void tc_scheduling_sfs_single(struct tc_scop* scop, struct tc_options* options, __isl_take isl_union_set* LD, __isl_take isl_union_map* S, __isl_take isl_union_map* R, __isl_take isl_set* ii_set, __isl_take isl_set* tile, __isl_take isl_map* Rtile, __isl_take isl_id_list* II, __isl_take isl_id_list* I)
{
    isl_ctx* ctx = isl_union_set_get_ctx(LD);
    
    if (!tc_is_lex_forward(Rtile))
    {
        tc_error("Backward relation detected.");
        tc_die(tc_exit_code_nonlex);
    }

    isl_set* LD_normalized = tc_normalize_union_set(LD, S);
    
    isl_map* R_normalized = tc_normalize_union_map(R, S);
    
    isl_bool exact = isl_bool_false;
    
    isl_map* R_plus_normalized = tc_transitive_closure(isl_map_copy(R_normalized), S, &exact);
    
    isl_map* R_star_normalized = isl_map_union(isl_map_copy(R_plus_normalized), tc_make_identity(isl_map_copy(R_normalized)));

    tc_debug_map(R_star_normalized, "R* (exact=%d)", exact);

    if (exact != isl_bool_true)
    {
        tc_warn("Inexact R*. The results can be non-optimal. Restart TC with a different transitive closure method.");
    }

    isl_id_list* IR = tc_ids_sequence(ctx, "ir", isl_set_n_dim(tile));
    
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

    isl_map* repr_ind_map = tc_make_identity_from_set(isl_set_copy(repr_ind));
    
    repr_ind = tc_parameterize_all(repr_ind, IR);

    isl_map* Rusc_star = NULL;

    // SFS = R∗(REPR_IND)
    if (isl_map_is_empty(Rusc) == isl_bool_true)
    {
        Rusc_star = repr_ind_map;
    }
    else
    {
        isl_map* Rusc_plus = tc_transitive_closure(isl_map_copy(Rusc), S, &exact);

        Rusc_star = isl_map_union(Rusc_plus, tc_make_identity(isl_map_copy(Rusc)));
        Rusc_star = isl_map_coalesce(Rusc_star);

        tc_debug_map(Rusc_star, "R_USC* (exact=%d)", exact);

        if (exact != isl_bool_true)
        {
            tc_warn("Inexact R_USC*. The results can be non-optimal. Restart TC with a different transitive closure method.");
        }

        Rusc_star = isl_map_union(Rusc_star, repr_ind_map);
    }

    // SFS = R∗(R_USC∗(REPR_IND))
    isl_set* sfs = isl_set_apply(isl_set_apply(repr_ind, Rusc_star), R_star_normalized);
    sfs = isl_set_coalesce(sfs);

    tc_debug_set(sfs, "SFS");

    isl_id_list* IR_II = isl_id_list_concat(isl_id_list_copy(IR), isl_id_list_copy(II));
    isl_id_list* IR_II_I = isl_id_list_concat(isl_id_list_copy(IR_II), isl_id_list_copy(I));

    // S_SLICE = { [I] -> [IR] | I in SFS(IR) and IR in REPR_IND }

    isl_set* Sslice = tc_make_set_constraints(isl_set_copy(sfs), I);
    tc_debug_set(Sslice, "Sslice");

    // S_TILE = { [I] -> [II] | I in TILE(II) and II in II_SET }

    isl_set* Stile = tc_make_set_constraints(isl_set_copy(tile), I);
    tc_debug_set(Stile, "Stile");
    
    // TILE_EXT = { [IR, II, I] | IR = S_SLICE(I) and II in S_TILE(I) }

    isl_set* tile_ext = tc_make_set(ctx, IR_II_I, I, "");
    tile_ext = isl_set_intersect(tile_ext, isl_set_copy(tile));
    tile_ext = isl_set_intersect_params(tile_ext, Sslice);
    tile_ext = isl_set_intersect_params(tile_ext, Stile);
    
    tile_ext = tc_lift_up_set_params(tile_ext, II);
    tile_ext = tc_lift_up_set_params(tile_ext, IR);

    tile_ext = tc_project_out_params(tile_ext, IR_II_I);
    tile_ext = isl_set_coalesce(tile_ext);
    
    tc_debug_set(tile_ext, "TILE_EXT");
        
    isl_union_map* S_ext = tc_extend_schedule(isl_union_map_copy(S), isl_id_list_n_id(IR) + isl_id_list_n_id(II));
    
    isl_id_list* iterators = isl_id_list_copy(IR_II_I);
    
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
        tc_codegen_omp_task_for(scop, options, S_ext, tile_ext, iterators, IR, 1);
    }
    else if (tc_codegen_enum_omp_gpu == codegen)
    {
        tc_codegen_omp_gpu(scop, options, S_ext, tile_ext, iterators, IR);
    }

    isl_id_list_free(I);
    isl_id_list_free(IR);
    isl_id_list_free(II);
    isl_id_list_free(IR_II);
    isl_id_list_free(IR_II_I);
    isl_id_list_free(iterators);
    isl_map_free(Rusc);
    isl_set_free(LD_normalized);
    isl_map_free(R_normalized);
    isl_map_free(R_plus_normalized);
    isl_set_free(tile);
    isl_set_free(ii_set);
    isl_set_free(sfs);
    isl_union_set_free(LD);
    isl_union_map_free(S);
    isl_union_map_free(R);
    isl_map_free(Rtile);
}

void tc_scheduling_sfs_multiple(struct tc_scop* scop, struct tc_options* options, __isl_take isl_union_set* LD, __isl_take isl_union_map* S, __isl_take isl_union_map* R, __isl_take isl_set* ii_set, __isl_take isl_set* tile, __isl_take isl_map* Rtile, __isl_take isl_id_list* II, __isl_take isl_id_list* I)
{
    isl_ctx* ctx = isl_union_set_get_ctx(LD);
    
    if (!tc_is_lex_forward(Rtile))
    {
        tc_error("Backward relation detected.");
        tc_die(tc_exit_code_nonlex);
    }
        
    isl_set* LD_normalized = tc_normalize_union_set(LD, S);
    
    isl_map* R_normalized = tc_normalize_union_map(R, S);
    
    isl_bool exact = isl_bool_false;
    
    isl_map* R_plus_normalized = tc_transitive_closure(isl_map_copy(R_normalized), S, &exact);
    
    isl_map* R_star_normalized = isl_map_union(isl_map_copy(R_plus_normalized), tc_make_identity(isl_map_copy(R_normalized)));

    tc_debug_map(R_star_normalized, "R* (exact=%d)", exact);

    if (exact != isl_bool_true)
    {
        tc_warn("Inexact R*. The results can be non-optimal. Restart TC with a different transitive closure method.");
    }
    
    isl_id_list* IR = tc_ids_sequence(ctx, "ir", isl_set_n_dim(tile));
    isl_id_list* IIR = tc_ids_sequence(ctx, "iir", isl_set_n_dim(tile));
    
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

    isl_map* repr_ind_map = tc_make_identity_from_set(isl_set_copy(repr_ind));

    repr_ind = tc_parameterize_all(repr_ind, IR);

    isl_map* Rusc_star = NULL;

    // SFS = R∗(REPR_IND)
    if (isl_map_is_empty(Rusc) == isl_bool_true)
    {
        Rusc_star = repr_ind_map;
    }
    else
    {
        isl_map* Rusc_plus = tc_transitive_closure(isl_map_copy(Rusc), S, &exact);

        Rusc_star = isl_map_union(Rusc_plus, tc_make_identity(isl_map_copy(Rusc)));
        Rusc_star = isl_map_coalesce(Rusc_star);

        tc_debug_map(Rusc_star, "R_USC* (exact=%d)", exact);

        if (exact != isl_bool_true)
        {
            tc_warn("Inexact R_USC*. The results can be non-optimal. Restart TC with a different transitive closure method.");
        }

        Rusc_star = isl_map_union(Rusc_star, repr_ind_map);
    }

    // SFS = R∗(R_USC∗(REPR_IND))
    isl_set* sfs = isl_set_apply(isl_set_apply(repr_ind, Rusc_star), R_star_normalized);
    sfs = isl_set_coalesce(sfs);

    tc_debug_set(sfs, "SFS");

    isl_id_list* IIR_II = isl_id_list_concat(isl_id_list_copy(IIR), isl_id_list_copy(II));
    isl_id_list* IIR_II_I = isl_id_list_concat(isl_id_list_copy(IIR_II), isl_id_list_copy(I));
    isl_id_list* IR_IIR_II_I = isl_id_list_concat(isl_id_list_copy(IR), isl_id_list_copy(IIR_II_I));

    // S_SLICE = { [I] -> [IR] | I in SFS(IR) and IR in REPR_IND }

    isl_set* Sslice = tc_make_set_constraints(isl_set_copy(sfs), I);
    tc_debug_set(Sslice, "Sslice");

    // S_TILE = { [I] -> [II] | I in TILE(II) and II in II_SET }

    isl_set* Stile_ii_i = tc_make_set_constraints(isl_set_copy(tile), I);
    tc_debug_set(Stile_ii_i, "Stile");

    isl_set* Stile_iir_ir = tc_make_set_constraints(isl_set_copy(tile), I);
    Stile_iir_ir = tc_rename_params(Stile_iir_ir, II, IIR);
    Stile_iir_ir = tc_rename_params(Stile_iir_ir, I, IR);

    // TILE_EXT = [IR] -> { [IIR, II, I] | IIR = S_TILE(IR) and IR = S_SLICE(I) and II in S_TILE(I) }

    isl_set* tile_ext = tc_make_set(ctx, IR_IIR_II_I, I, "");
    tile_ext = isl_set_intersect(tile_ext, isl_set_copy(tile));
    tile_ext = isl_set_intersect_params(tile_ext, Sslice);
    tile_ext = isl_set_intersect_params(tile_ext, Stile_ii_i);
    tile_ext = isl_set_intersect_params(tile_ext, Stile_iir_ir);

    tile_ext = tc_lift_up_set_params(tile_ext, II);
    tile_ext = tc_lift_up_set_params(tile_ext, IIR);

    tile_ext = tc_project_out_params(tile_ext, IR_IIR_II_I);
    tile_ext = isl_set_coalesce(tile_ext);

    tc_debug_set(tile_ext, "TILE_EXT");

    isl_union_map* S_ext = tc_extend_schedule(isl_union_map_copy(S), isl_id_list_n_id(IIR) + isl_id_list_n_id(II));
    
    isl_id_list* iterators = isl_id_list_copy(IIR_II_I);

    enum tc_codegen_enum codegen = tc_options_codegen(options);
    
    if (tc_codegen_enum_serial == codegen)
    {
        tc_codegen_serial(scop, options, S_ext, tile_ext, iterators);
    }
    else if (tc_codegen_enum_omp_cpu_for == codegen)
    {
        tc_codegen_omp_parallel_for(scop, options, S_ext, tile_ext, iterators, IIR, 0);
    }
    else if (tc_codegen_enum_omp_cpu_task == codegen)
    {
        tc_codegen_omp_task_for(scop, options, S_ext, tile_ext, iterators, IIR, 1);
    }
    else if (tc_codegen_enum_omp_gpu == codegen)
    {
        tc_codegen_omp_gpu(scop, options, S_ext, tile_ext, iterators, IIR);
    }

    isl_id_list_free(I);
    isl_id_list_free(IR);
    isl_id_list_free(II);
    isl_id_list_free(IIR);
    isl_id_list_free(IIR_II);
    isl_id_list_free(IIR_II_I);
    isl_id_list_free(IR_IIR_II_I);
    isl_id_list_free(iterators);
    isl_map_free(Rusc);
    isl_set_free(LD_normalized);
    isl_map_free(R_normalized);
    isl_map_free(R_plus_normalized);
    isl_set_free(tile);
    isl_set_free(ii_set);
    isl_set_free(sfs);
    isl_union_set_free(LD);
    isl_union_map_free(S);
    isl_union_map_free(R);
    isl_map_free(Rtile);
}
