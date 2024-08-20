#include "iterative_transitive_closure.h"
#include "utility.h"
#include "debug.h"

#include <isl/ctx.h>
#include <isl/map.h>
#include <isl/union_map.h>

static isl_stat tc_iterative_union_closure(__isl_take isl_basic_map* bmap, void* user)
{
    isl_map** union_map = (isl_map**)user;
    isl_bool exact = isl_bool_false;

    isl_map* R = isl_map_from_basic_map(isl_basic_map_copy(bmap));

    tc_debug_map(R, "R");

    isl_map* R_plus = isl_map_transitive_closure(R, &exact);
    tc_debug_map(R_plus, "R+");

    R_plus = isl_map_compute_divs(R_plus);
    R_plus = isl_map_remove_redundancies(R_plus);
    R_plus = isl_map_coalesce(R_plus);

    if (*union_map)
    {
        *union_map = isl_map_union(*union_map, (exact == isl_bool_true) ? isl_map_copy(R_plus) : isl_map_from_basic_map(isl_basic_map_copy(bmap)));
    }
    else
    {
        *union_map = (exact == isl_bool_true) ? isl_map_copy(R_plus) : isl_map_from_basic_map(isl_basic_map_copy(bmap));
    }

    *union_map = isl_map_compute_divs(*union_map);
    *union_map = isl_map_remove_redundancies(*union_map);
    *union_map = isl_map_coalesce(*union_map);

    isl_basic_map_free(bmap);
    isl_map_free(R_plus);

    return isl_stat_ok;
}

__isl_give isl_map* tc_iterative_transitive_closure(__isl_take isl_map* map, int max_iterations, int max_disjunctions, isl_bool* exact)
{
    isl_map* out = NULL;
    int n = 0;
    
    isl_map* identity = isl_map_identity(isl_map_get_space(map));

    isl_map_foreach_basic_map(map, &tc_iterative_union_closure, &out);
    out = isl_map_union(out, isl_map_copy(identity));
        
    isl_map* compose = isl_map_apply_range(isl_map_copy(out), isl_map_copy(out));
    compose = isl_map_compute_divs(compose);
    compose = isl_map_remove_redundancies(compose);
    compose = isl_map_coalesce(compose);

    *exact = isl_map_is_equal(compose, out);
            
    while (!*exact && /*compose->n < max_disjunctions &&*/ n < max_iterations) 
    {
        tc_debug("n=%d", n);

        compose = isl_map_subtract(compose, isl_map_copy(identity));
        compose = isl_map_compute_divs(compose);
        compose = isl_map_remove_redundancies(compose);
        compose = isl_map_coalesce(compose);

        out = isl_map_free(out);
        isl_map_foreach_basic_map(compose, &tc_iterative_union_closure, &out);
        out = isl_map_union(out, isl_map_copy(identity));

        compose = isl_map_union(compose, isl_map_copy(identity));
        compose = isl_map_apply_range(compose, isl_map_copy(out));
        compose = isl_map_remove_redundancies(compose);
        compose = isl_map_coalesce(compose);

        *exact = isl_map_is_equal(compose, out);

        ++n;
    }
    
    out = isl_map_subtract(out, identity);
    out = isl_map_compute_divs(out);
    out = isl_map_remove_redundancies(out);
    out = isl_map_coalesce(out);
    
    isl_map_free(map);
    isl_map_free(compose);
    
    return out;
}
