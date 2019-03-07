#include "merge_tiling.h"
#include "scheduling.h"
#include "config.h"
#include "tiling.h"
#include "for_decorator.h"
#include "utility.h"
#include "options.h"
#include "tile_statistics.h"
#include "slicing.h"
#include "transitive_closure.h"
#include "debug.h"
#include "slicing.h"

#include <isl/ctx.h>
#include <isl/space.h>
#include <isl/id.h>
#include <isl/set.h>
#include <isl/map.h>
#include <isl/union_set.h>
#include <isl/union_map.h>

#include <string.h>

#include <vector>
#include <map>
#include <string>

void tc_algorithm_merge_tiling(struct tc_scop* scop, struct tc_options* options)
{
    isl_ctx* ctx = scop->ctx;
    
    isl_union_set* LD = isl_union_set_copy(scop->domain);
    tc_debug_uset(LD, "LD");
    
    isl_union_map* S = isl_union_map_copy(scop->schedule);
    tc_debug_umap(S, "S");

    isl_union_map* R = isl_union_map_copy(scop->relation);
    tc_debug_umap(R, "R");
    
    tc_debug_umap(scop->reads, "RA");
    tc_debug_umap(scop->writes, "WA");
        
    isl_basic_set* sample = isl_set_sample(tc_normalize_union_set(LD, S));
    
    isl_space* space = isl_basic_set_get_space(sample);
            
    isl_id_list* I = tc_ids_sequence(ctx, "i", isl_space_dim(space, isl_dim_set));
    isl_id_list* II = tc_ids_sequence(ctx, "ii", isl_space_dim(space, isl_dim_set));
                
    std::map<std::string, std::vector<int> > blocks = tc_options_blocks(options);
            
    isl_set* tile;
    isl_set* ii_set;
    
    tc_tile_loop_nest(LD, S, II, I, &tile, &ii_set, blocks);
    
    tc_debug_set(tile, "TILE");
    tc_debug_set(ii_set, "II_SET");
    //tc_debug_set_card(ii_set, "II_SET");
    
    tile = isl_set_coalesce(tile);
    ii_set = isl_set_coalesce(ii_set);
                
    isl_map* R_normalized = tc_normalize_union_map(R, S);
    R_normalized = isl_map_coalesce(R_normalized);
    
    tc_debug_map(R_normalized, "R_norm");
        
    isl_map* Rtile = tc_Rtile_map(II, tile, R_normalized);
    //Rtile = isl_map_compute_divs(Rtile);
    //Rtile = isl_map_remove_redundancies(Rtile);
    Rtile = isl_map_coalesce(Rtile);
    
    tc_debug_map(Rtile, "R_TILE");
    //tc_debug_map_card(Rtile, "R_TILE");
    
    isl_union_map* Rtile_denormalized = tc_denormalize_map(Rtile, S);
    tc_debug_umap(Rtile_denormalized, "R_TILE_denorm");
    isl_union_map_free(Rtile_denormalized);
                    
    int exact;    
    isl_map* Rtile_plus = tc_transitive_closure(isl_map_copy(Rtile), S, &exact);
    //Rtile_plus = isl_map_compute_divs(Rtile_plus);
    Rtile_plus = isl_map_coalesce(Rtile_plus);
    
    tc_debug_map(Rtile_plus, "R_TILE^+ (exact=%d)", exact);
        
    isl_map* Tcycle = tc_Tcycle_map(II, Rtile_plus);
    Tcycle = isl_map_coalesce(Tcycle);
    tc_debug_map(Tcycle, "T_CYCLE");
    
    isl_union_map* Tcycle_denormalized = tc_denormalize_map(Tcycle, S);
    tc_debug_umap(Tcycle_denormalized, "T_CYCLE_denorm");
    isl_union_map_free(Tcycle_denormalized);
    
    // II_SET_M = II_SET - range(T_CYCLE)
    isl_set* ii_set_m = isl_set_subtract(isl_set_copy(ii_set), isl_map_range(isl_map_copy(Tcycle)));
    ii_set_m = isl_set_coalesce(ii_set_m);
    
    tc_debug_set(ii_set_m, "II_SET_M");
    //tc_debug_set_card(ii_set_m, "II_SET_M");
                    
    isl_set* tile_m = tc_tile_m_set(II, tile, ii_set_m, Tcycle);
    tile_m = isl_set_coalesce(tile_m);
    
    tc_debug_set(tile_m, "TILE_M");
        
    isl_map* Rtile_m = tc_Rtile_map(II, tile_m, R_normalized);
    Rtile_m = isl_map_coalesce(Rtile_m);
    
    tc_debug_map(Rtile_m, "R_TILE_M");
    //tc_debug_map_card(Rtile_m, "R_TILE_M");
    
    tc_debug_bool(isl_set_is_equal(tile, tile_m), "TILE = TILE_M");
    
    if (tc_options_is_report(options))
    {
        isl_set* bounds = tc_options_get_report_bounds(options, ctx);
        
        struct tc_tile_statistics* stats = tc_compute_tile_statistics(tile_m, ii_set_m, II, bounds, LD, S, scop->reads, scop->writes, scop, options, blocks);
        
        tc_tile_statistics_print(options->output, stats);
        
        tc_tile_statistics_free(stats);
        isl_set_free(bounds);
    }
    
    tc_scheduling(scop, options, LD, S, R, ii_set_m, tile_m, Rtile_m, II, I);
    
    isl_set_free(tile);
    isl_set_free(ii_set);
    isl_map_free(Rtile);
    isl_map_free(Rtile_plus);
    isl_map_free(Tcycle);
    isl_map_free(R_normalized);
    isl_basic_set_free(sample);
    isl_space_free(space); 
}
