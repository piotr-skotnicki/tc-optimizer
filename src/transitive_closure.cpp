#include "transitive_closure.h"
#include "utility.h"

#include <isl/set.h>
#include <isl/map.h>

#include <stddef.h>

__isl_give isl_map* tc_transitive_closure(__isl_take isl_map* R, int max)
{
    isl_ctx* ctx = isl_map_get_ctx(R);
    
    isl_id_list* I = tc_ids_sequence(ctx, "i", isl_map_n_in(R));
    isl_id_list* J = tc_ids_sequence(ctx, "j", isl_map_n_out(R));
        
    isl_map* RP = tc_make_map(ctx, NULL, I, J, tc_tuples_op(I, J, "<=").c_str());
    
    isl_map* R1 = isl_map_intersect(isl_map_copy(R), RP);
    isl_map* R2 = isl_map_subtract(isl_map_copy(R), isl_map_copy(R1));

    R1 = isl_map_project_out(R1, isl_dim_param, 0, isl_map_n_param(R1));
    R2 = isl_map_project_out(R2, isl_dim_param, 0, isl_map_n_param(R2));
    
    int exact;
    isl_map* R1_plus = isl_map_transitive_closure(isl_map_copy(R1), &exact);
    isl_map* R2_plus = isl_map_transitive_closure(isl_map_copy(R2), &exact);
    
    isl_map* R_plus_old = isl_map_union(isl_map_union(R1, R2), isl_map_union(R1_plus, R2_plus));
        
    int k = 1;
    while (1)
    {
        isl_map* R_plus_new = isl_map_union(isl_map_copy(R_plus_old), isl_map_apply_range(isl_map_copy(R_plus_old), isl_map_copy(R_plus_old)));
        
        isl_map* delta = isl_map_subtract(isl_map_copy(R_plus_new), R_plus_old);
        
        R_plus_old = R_plus_new;
        
        isl_bool empty = isl_map_is_empty(delta);
        isl_map_free(delta);
        
        if (empty || k == max)
        {
            break;
        }
        k = k + 1;
    }
            
    isl_set* IS = isl_set_union(isl_map_domain(isl_map_copy(R)), isl_map_range(isl_map_copy(R)));
    
    isl_map* R_IS = tc_make_map(ctx, NULL, I, J, tc_tuples_gt(J, I).c_str());
        
    R_IS = isl_map_intersect_domain(R_IS, isl_set_copy(IS));
    R_IS = isl_map_intersect_range(R_IS, IS);
    
    isl_map* R_plus = isl_map_intersect(R_plus_old, R_IS);
    R_plus = isl_map_coalesce(R_plus);
    
    isl_map_free(R);
    
    isl_id_list_free(I);
    isl_id_list_free(J);
    
    return R_plus;
}
