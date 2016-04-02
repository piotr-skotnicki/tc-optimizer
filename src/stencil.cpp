#include "stencil.h"
#include "scop.h"
#include "tiling.h"
#include "for_decorator.h"
#include "utility.h"
#include "options.h"

#include <isl/id.h>
#include <isl/ast.h>
#include <isl/printer.h>

#include <string.h>
#include <stddef.h>

#include <vector>
#include <map>
#include <string>

std::string tc_algorithm_stencil_tuples_tlt(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs)
{
    isl_id_list* ids_0_lhs = isl_id_list_from_id(isl_id_list_get_id(lhs, 0));
    isl_id_list* ids_0_rhs = isl_id_list_from_id(isl_id_list_get_id(rhs, 0));
    
    std::string str = tc_conjunction(tc_tuples_lt(lhs, rhs), tc_tuples_eq(ids_0_lhs, ids_0_rhs));
    
    isl_id_list_free(ids_0_lhs);
    isl_id_list_free(ids_0_rhs);
    
    return str;
}

std::string tc_algorithm_stencil_tuples_tgt(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs)
{
    isl_id_list* ids_0_lhs = isl_id_list_from_id(isl_id_list_get_id(lhs, 0));
    isl_id_list* ids_0_rhs = isl_id_list_from_id(isl_id_list_get_id(rhs, 0));
    
    std::string str = tc_conjunction(tc_tuples_gt(lhs, rhs), tc_tuples_eq(ids_0_lhs, ids_0_rhs));
    
    isl_id_list_free(ids_0_lhs);
    isl_id_list_free(ids_0_rhs);
    
    return str;
}

void tc_algorithm_stencil(int argc, char* argv[], struct tc_scop* scop)
{
    isl_ctx* ctx = scop->ctx;
    
    isl_union_set* LD = scop->domain;
    
    isl_union_map* S = scop->schedule;
    
    isl_union_map* R = scop->relation;
    
    isl_basic_set* sample = isl_set_sample(tc_normalize_union_set(LD, S));
    
    isl_space* space = isl_basic_set_get_space(sample);
            
    isl_id_list* I = tc_ids_sequence(ctx, "i", isl_space_dim(space, isl_dim_set));
    isl_id_list* II = tc_ids_sequence(ctx, "ii", isl_space_dim(space, isl_dim_set));
    
    isl_id_list* Iprim = tc_ids_prim(I);
    isl_id_list* Ibis = tc_ids_bis(I);
    
    isl_id_list* IIprim = tc_ids_prim(II);
    isl_id_list* IIbis = tc_ids_bis(II);
    
    int is_exact;
    isl_union_map* R_plus = isl_union_map_transitive_closure(isl_union_map_copy(R), &is_exact);
            
    isl_map* R_plus_normalized = tc_normalize_union_map(R_plus, S);
        
    std::map<std::string, std::vector<int> > blocks = tc_options_blocks(argc, argv);
    
    isl_set* tile;
    isl_set* ii_set;
    
    tc_tile_loop_nest(scop->domain, scop->schedule, blocks, II, I, &tile, &ii_set);
        
    int k = 0;
    
    isl_set* subtiles_ext = NULL;
        
    while (!isl_set_is_empty(tile))
    {
        k = k + 1;
        
        isl_set* tile_lt = tc_tile_set_of(tile, ii_set, II, &tc_algorithm_stencil_tuples_tlt);
        
        isl_set* tile_gt = tc_tile_set_of(tile, ii_set, II, &tc_algorithm_stencil_tuples_tgt);
        
        isl_set* subtile = isl_set_copy(tile);
        
        subtile = isl_set_subtract(subtile, isl_set_apply(tile_lt, isl_map_copy(R_plus_normalized)));
        
        subtile = isl_set_subtract(subtile, isl_set_apply(tile_gt, isl_map_copy(R_plus_normalized)));
        
        tile = isl_set_subtract(tile, isl_set_copy(subtile));
        
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
    
    isl_union_map* S_ext = NULL;
    isl_map_list* S_maps = tc_collect_maps(S);
    
    for (int i = 0; i < isl_map_list_n_map(S_maps); ++i)
    {
        isl_map* map = isl_map_list_get_map(S_maps, i);
        
        map = isl_map_insert_dims(map, isl_dim_out, 0, isl_id_list_n_id(II) + 1);
    
        if (NULL == S_ext)
        {
            S_ext = isl_union_map_from_map(map);
        }
        else
        {
            S_ext = isl_union_map_add_map(S_ext, map);
        }
    }
    
    isl_map_list_free(S_maps);
                    
    isl_ast_build* ast_build = isl_ast_build_from_context(isl_set_copy(scop->pet->context));
        
    isl_id_list* iterators = isl_id_list_copy(II);
    iterators = isl_id_list_concat(iterators, isl_id_list_copy(I));
    iterators = isl_id_list_insert(iterators, 1, isl_id_alloc(ctx, "k", NULL));
    
    isl_id_list* parallel_iterators = isl_id_list_drop(isl_id_list_copy(II), 0, 1);
    
    ast_build = isl_ast_build_set_iterators(ast_build, iterators);
    //ast_build = isl_ast_build_set_at_each_domain(ast_build, &at_each_domain, scop);
        
    isl_union_map* S_ext_prim = isl_union_map_intersect_range(S_ext, isl_union_set_from_set(subtiles_ext));
    
    isl_ast_node* ast_tile = isl_ast_build_ast_from_schedule(ast_build, S_ext_prim);
    
    isl_printer* printer = isl_printer_to_str(ctx);
    
    printer = isl_printer_set_output_format(printer, ISL_FORMAT_C);
    
    isl_ast_print_options* ast_options = isl_ast_print_options_alloc(ctx);
    
    ast_options = isl_ast_print_options_set_print_for(ast_options, &tc_for_decorator_omp_parallel_for_first, parallel_iterators);
    
    printer = isl_printer_print_str(printer, "#include <omp.h>\n");
    
    printer = isl_ast_node_print_macros(ast_tile, printer);
    
    printer = tc_print_statements_macros(scop, printer, ast_build);
    
    printer = isl_printer_print_str(printer, "#pragma scop\n");
    printer = isl_ast_node_print(ast_tile, printer, ast_options);
    printer = isl_printer_print_str(printer, "#pragma endscop\n");
    
    char* code = isl_printer_get_str(printer);
    
    printf("%s", code);
    
    free(code);
    
    isl_id_list_free(parallel_iterators);
    isl_printer_free(printer);
    isl_ast_node_free(ast_tile);
    isl_ast_build_free(ast_build);
    isl_set_free(ii_set);
    isl_set_free(tile);
    isl_basic_set_free(sample);
    isl_union_map_free(R_plus);
    isl_map_free(R_plus_normalized);
    isl_id_list_free(I);
    isl_id_list_free(Iprim);
    isl_id_list_free(Ibis);
    isl_id_list_free(II);
    isl_id_list_free(IIprim);
    isl_id_list_free(IIbis);
    isl_space_free(space);
}
