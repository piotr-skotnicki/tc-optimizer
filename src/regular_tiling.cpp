#include "regular_tiling.h"
#include "lex_scheduling.h"
#include "sfs_scheduling.h"
#include "free_scheduling.h"
#include "dynamic_free_scheduling.h"
#include "utility.h"
#include "tiling.h"
#include "slicing.h"
#include "options.h"
#include "debug.h"
#include "tile_statistics.h"
#include "transitive_closure.h"

#include <isl/ctx.h>
#include <isl/space.h>
#include <isl/id.h>
#include <isl/set.h>
#include <isl/map.h>
#include <isl/union_set.h>
#include <isl/union_map.h>

#include <vector>
#include <map>
#include <string>

void tc_algorithm_regular_tiling(struct tc_scop* scop, struct tc_options* options)
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
                
    std::map<std::string, std::vector<int> > blocks = tc_options_blocks(options);
        
    isl_set* tile;
    isl_set* ii_set;
    
    tc_tile_loop_nest(LD, S, II, I, &tile, &ii_set, blocks);
    
    tc_debug_set(tile, "TILE");
    tc_debug_set(ii_set, "II_SET");
    
    isl_map* R_normalized = tc_normalize_union_map(R, S);
        
    if (tc_options_is_report(options))
    {
        isl_set* bounds = tc_options_get_report_bounds(options, ctx);
        
        struct tc_tile_statistics* stats = tc_compute_tile_statistics(tile, ii_set, II, bounds, LD, S, scop, blocks);
        
        tc_tile_statistics_print(options->output, stats);
        
        tc_tile_statistics_free(stats);
        isl_set_free(bounds);
    }
    
    isl_map* Rtile = tc_Rtile_map(II, tile, R_normalized);

    enum tc_scheduling_enum scheduling = tc_options_scheduling(options);
    
    if (tc_scheduling_enum_lex == scheduling)
    {
        tc_scheduling_lex(scop, options, LD, S, R, ii_set, tile, Rtile, II);
    }
    else if (tc_scheduling_enum_sfs_tile == scheduling)
    {
        tc_scheduling_sfs_tiles(scop, options, LD, S, R, ii_set, tile, Rtile, II);
    }
    else if (tc_scheduling_enum_free_rk == scheduling)
    {
        tc_scheduling_free_schedule_rk(scop, options, LD, S, R, ii_set, tile, Rtile, II);
    }
    else if (tc_scheduling_enum_free_karl == scheduling)
    {
        tc_scheduling_free_schedule_karl(scop, options, LD, S, R, ii_set, tile, Rtile, II);
    }
    else if (tc_scheduling_enum_free_dynamic == scheduling)
    {
        tc_scheduling_dynamic_free_schedule(scop, options, LD, S, R, ii_set, tile, Rtile, II);
    }
    
    isl_map_free(R_normalized);
    isl_id_list_free(I);
    isl_basic_set_free(sample);
    isl_space_free(space);
}