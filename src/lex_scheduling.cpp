#include "lex_scheduling.h"
#include "serial_codegen.h"
#include "utility.h"
#include "scop.h"
#include "options.h"

#include <isl/ctx.h>
#include <isl/id.h>
#include <isl/set.h>
#include <isl/union_map.h>
#include <isl/union_set.h>

void tc_scheduling_lex(struct tc_scop* scop, struct tc_options* options, __isl_take isl_union_set* LD, __isl_take isl_union_map* S, __isl_take isl_union_map* R, __isl_take isl_set* ii_set, __isl_take isl_set* tile, __isl_take isl_map* Rtile, __isl_take isl_id_list* II, __isl_take isl_id_list* I)
{   
    if (!tc_is_lex_forward(Rtile))
    {
        tc_options_error("Backward relation detected");
        return;
    }
    
    isl_set* tile_ext = tc_lift_up_set_params(tile, II);
    
    tile_ext = isl_set_coalesce(tile_ext);
    
    isl_union_map* S_ext = tc_extend_schedule(isl_union_map_copy(S), isl_id_list_n_id(II));
    
    isl_id_list* iterators = isl_id_list_concat(isl_id_list_copy(II), isl_id_list_copy(I));
    
    tc_codegen_serial(scop, options, S_ext, tile_ext, iterators);
    
    isl_id_list_free(iterators);
    isl_id_list_free(I);
    isl_id_list_free(II);
    isl_union_set_free(LD);
    isl_union_map_free(S);
    isl_union_map_free(R);
    isl_set_free(ii_set);
    isl_map_free(Rtile);
}
