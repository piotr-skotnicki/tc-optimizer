#include "slicing.h"
#include "utility.h"

#include <isl/set.h>
#include <isl/map.h>
#include <isl/id.h>

__isl_give isl_set* tc_uds_set(__isl_keep isl_map* R)
{
    isl_set* uds = isl_set_subtract(isl_map_domain(isl_map_copy(R)), isl_map_range(isl_map_copy(R)));
    return uds;
}

__isl_give isl_set* tc_udd_set(__isl_keep isl_map* R)
{
    isl_set* udd = isl_set_subtract(isl_map_range(isl_map_copy(R)), isl_map_domain(isl_map_copy(R)));
    return udd;
}

__isl_give isl_set* tc_ind_set(__isl_keep isl_set* LD, __isl_keep isl_map* R)
{
    isl_set* ind = isl_set_subtract(isl_set_subtract(isl_set_copy(LD), isl_map_domain(isl_map_copy(R))), isl_map_range(isl_map_copy(R)));
    return ind;
}

__isl_give isl_map* tc_Rusc_map(__isl_keep isl_map* R)
{
    isl_ctx* ctx = isl_map_get_ctx(R);
    
    int n_in_dim = isl_map_n_in(R);
    
    isl_id_list* e = tc_ids_sequence(ctx, "e", n_in_dim);
    isl_id_list* e_prim = tc_ids_prim(e);
    isl_id_list* e_bis = tc_ids_bis(e);
    
    int exact;
    isl_map* R_plus = isl_map_transitive_closure(isl_map_copy(R), &exact);
    isl_map* R_star = isl_map_union(isl_map_copy(R_plus), tc_make_identity(isl_map_copy(R_plus)));
    
    isl_set* R_star_constraints_prim = tc_make_map_constraints(isl_map_copy(R_star), e_prim, e);
    isl_set* R_star_constraints_bis = tc_make_map_constraints(isl_map_copy(R_star), e_bis, e);
    
    isl_set* uds = tc_uds_set(R);
    
    isl_id_list* e_eprim_ebis = isl_id_list_concat(isl_id_list_copy(e), isl_id_list_concat(isl_id_list_copy(e_prim), isl_id_list_copy(e_bis)));
    isl_map* Rusc = tc_make_map(ctx, e_eprim_ebis, e_prim, e_bis, tc_tuples_lt(e_prim, e_bis).c_str());
        
    Rusc = isl_map_intersect_domain(Rusc, isl_set_copy(uds));
    Rusc = isl_map_intersect_range(Rusc, isl_set_copy(uds));
    Rusc = isl_map_intersect_params(Rusc, R_star_constraints_prim);
    Rusc = isl_map_intersect_params(Rusc, R_star_constraints_bis);
    
    Rusc = tc_map_project_out_params(Rusc, e_eprim_ebis);
    
    isl_id_list_free(e);
    isl_id_list_free(e_prim);
    isl_id_list_free(e_bis);
    isl_id_list_free(e_eprim_ebis);
    isl_set_free(uds);
    isl_map_free(R_plus);
    isl_map_free(R_star);    
    
    return Rusc;
}
