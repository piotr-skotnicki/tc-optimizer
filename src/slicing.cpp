#include "slicing.h"
#include "utility.h"

#include <isl/ctx.h>
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

__isl_give isl_set* tc_cds_set(__isl_keep isl_map* R)
{
    // CDS := { [e] : exists e', e'' : e = R^-1(e') = R^-1(e'') and e',e'' in range(R) and e' <> e'' }
    
    isl_ctx* ctx = isl_map_get_ctx(R);
    
    isl_map* R_inv = isl_map_reverse(isl_map_copy(R));
    
    int n_in_dim = isl_map_n_in(R);
    int n_out_dim = isl_map_n_out(R);
    
    isl_id_list* e = tc_ids_sequence(ctx, "e", n_in_dim);
    isl_id_list* e_out = tc_ids_sequence(ctx, "e", n_out_dim);
    isl_id_list* e_prim = tc_ids_prim(e_out);
    isl_id_list* e_bis = tc_ids_bis(e_out);
        
    // CDS := [e,e',e''] -> { [e] : [e',e]->R^-1 and [e'',e]->R^-1 and e',e'' in range(R) and e' <> e'' }
    
    isl_set* R_inv_constraints_prim = tc_make_map_constraints(isl_map_copy(R_inv), e_prim, e);
    isl_set* R_inv_constraints_bis = tc_make_map_constraints(R_inv, e_bis, e);
    
    isl_set* R_range_constraints_prim = tc_make_set_constraints(isl_map_range(isl_map_copy(R)), e_prim);
    isl_set* R_range_constraints_bis = tc_make_set_constraints(isl_map_range(isl_map_copy(R)), e_bis);
            
    isl_id_list* e_e_prim_e_bis = isl_id_list_concat(isl_id_list_concat(isl_id_list_copy(e), isl_id_list_copy(e_prim)), isl_id_list_copy(e_bis));
    
    isl_set* cds = tc_make_set(ctx, e_e_prim_e_bis, e, tc_tuples_neq(e_prim, e_bis).c_str());
    
    cds = isl_set_intersect_params(cds, R_inv_constraints_prim);
    cds = isl_set_intersect_params(cds, R_inv_constraints_bis);
    cds = isl_set_intersect_params(cds, R_range_constraints_prim);
    cds = isl_set_intersect_params(cds, R_range_constraints_bis);
        
    cds = tc_project_out_params(cds, e_e_prim_e_bis);    
    
    isl_id_list_free(e);
    isl_id_list_free(e_out);
    isl_id_list_free(e_prim);
    isl_id_list_free(e_bis);
    isl_id_list_free(e_e_prim_e_bis);
    
    return cds;    
}

__isl_give isl_set* tc_cdd_set(__isl_keep isl_map* R)
{
    // CDD := { [e] : exists e', e'' : e = R(e') = R(e'') and e',e'' in domain(R) and e' <> e'' }
    
    isl_ctx* ctx = isl_map_get_ctx(R);
        
    int n_in_dim = isl_map_n_in(R);
    int n_out_dim = isl_map_n_out(R);
    
    isl_id_list* e = tc_ids_sequence(ctx, "e", n_out_dim);
    isl_id_list* e_in = tc_ids_sequence(ctx, "e", n_in_dim);
    isl_id_list* e_prim = tc_ids_prim(e_in);
    isl_id_list* e_bis = tc_ids_bis(e_in);
    
    // CDD := [e,e',e''] -> { [e] : [e',e]->R and [e'',e]->R and e',e'' in domain(R) and e' <> e'' }
    
    isl_set* R_constraints_prim = tc_make_map_constraints(isl_map_copy(R), e_prim, e);
    isl_set* R_constraints_bis = tc_make_map_constraints(isl_map_copy(R), e_bis, e);
    
    isl_set* R_domain_constraints_prim = tc_make_set_constraints(isl_map_domain(isl_map_copy(R)), e_prim);
    isl_set* R_domain_constraints_bis = tc_make_set_constraints(isl_map_domain(isl_map_copy(R)), e_bis);
        
    isl_id_list* e_e_prim_e_bis = isl_id_list_concat(isl_id_list_concat(isl_id_list_copy(e), isl_id_list_copy(e_prim)), isl_id_list_copy(e_bis));
    
    isl_set* cdd = tc_make_set(ctx, e_e_prim_e_bis, e, tc_tuples_neq(e_prim, e_bis).c_str());
    
    cdd = isl_set_intersect_params(cdd, R_constraints_prim);
    cdd = isl_set_intersect_params(cdd, R_constraints_bis);
    cdd = isl_set_intersect_params(cdd, R_domain_constraints_prim);
    cdd = isl_set_intersect_params(cdd, R_domain_constraints_bis);
    
    cdd = tc_project_out_params(cdd, e_e_prim_e_bis);
    
    isl_id_list_free(e);
    isl_id_list_free(e_in);
    isl_id_list_free(e_prim);
    isl_id_list_free(e_bis);
    isl_id_list_free(e_e_prim_e_bis);
    
    return cdd;    
}
    
isl_bool tc_topology_is_chain(__isl_keep isl_set* cds, __isl_keep isl_set* cdd)
{
    return isl_set_is_empty(cds) && isl_set_is_empty(cdd) ? isl_bool_true : isl_bool_false;
}

isl_bool tc_topology_is_tree(__isl_keep isl_set* cds, __isl_keep isl_set* cdd)
{
    return !isl_set_is_empty(cds) && isl_set_is_empty(cdd) ? isl_bool_true : isl_bool_false;
}

isl_bool tc_topology_is_graph(__isl_keep isl_set* cds, __isl_keep isl_set* cdd)
{
    return !isl_set_is_empty(cdd) ? isl_bool_true : isl_bool_false;
}
