#include "regular_tiling.h"
#include "scheduling.h"
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
            
    tc_debug_umap(scop->reads, "RA");
    tc_debug_umap(scop->writes, "WA");
    
    isl_basic_set* sample = isl_set_sample(tc_normalize_union_set(LD, S));
    
    isl_space* space = isl_basic_set_get_space(sample);
            
    isl_id_list* I = tc_ids_sequence(ctx, "i", isl_space_dim(space, isl_dim_set));
    isl_id_list* II = tc_ids_sequence(ctx, "ii", isl_space_dim(space, isl_dim_set));
                
    std::map<std::string, std::vector<int> > blocks = tc_options_blocks(options);
    
    std::vector<std::vector<std::string> > groups = tc_options_groups(options);
        
    isl_set* tile;
    isl_set* ii_set;
    
    tc_tile_loop_nest(LD, S, II, I, &tile, &ii_set, blocks, groups);
    
    tc_debug_set(tile, "TILE");
    tc_debug_set(ii_set, "II_SET");
    
    isl_map* R_normalized = tc_normalize_union_map(R, S);
    tc_debug_map(R_normalized, "R_norm");
        
    if (tc_options_is_report(options))
    {
        isl_set* bounds = tc_options_get_report_bounds(options, ctx);
        
        struct tc_tile_statistics* stats = tc_compute_tile_statistics(tile, ii_set, II, bounds, LD, S, scop->reads, scop->writes, scop, options, blocks);
        
        tc_tile_statistics_print(options->output, stats);
        
        tc_tile_statistics_free(stats);
        isl_set_free(bounds);
    }
    
    isl_map* Rtile = tc_Rtile_map(II, tile, R_normalized);
    
    tc_debug_map(Rtile, "R_TILE");

    enum tc_scheduling_enum scheduling = tc_options_scheduling(options);
    
    tc_scheduling(scop, options, LD, S, R, ii_set, tile, Rtile, II, I);
    
    isl_map_free(R_normalized);
    isl_basic_set_free(sample);
    isl_space_free(space);
}
