#include "slicing.h"
#include "utility.h"
#include "debug.h"
#include "transitive_closure.h"
#include "options.h"

#include <isl/ctx.h>
#include <isl/set.h>
#include <isl/map.h>
#include <isl/id.h>

#include <cstddef>

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

__isl_give isl_map* tc_Rusc_map(__isl_keep isl_map* R, __isl_keep isl_union_map* S)
{
    // R_UCS := { [e] -> [e'] : e,e' in UDS and e << e' and e' in (R + R^-1)^*(e) }
    
    isl_ctx* ctx = isl_map_get_ctx(R);
    
    int n_in_dim = isl_map_n_in(R);
    
    isl_set* uds = tc_uds_set(R); 
    
    isl_id_list* e = tc_ids_sequence(ctx, "e", n_in_dim);
    isl_id_list* e_prim = tc_ids_prim(e);
    
    isl_map* R_inv = isl_map_reverse(isl_map_copy(R));
    
    isl_map* R_U_R_inv = isl_map_union(isl_map_copy(R), R_inv);
    
    int exact;
    isl_map* R_U_R_inv_plus = tc_transitive_closure(isl_map_copy(R_U_R_inv), S, &exact);
    isl_map* R_U_R_inv_star = isl_map_union(R_U_R_inv_plus, tc_make_identity(R_U_R_inv));    
    R_U_R_inv_star = isl_map_coalesce(R_U_R_inv_star);
    
    tc_debug_map(R_U_R_inv_star, "(R + R^-1)* (exact=%d)", exact);
       
    isl_set* R_U_R_inv_star_constraints_e_e_prim = tc_make_map_constraints(R_U_R_inv_star, e, e_prim);
        
    isl_id_list* e_eprim = isl_id_list_concat(isl_id_list_copy(e), isl_id_list_copy(e_prim));
    isl_map* Rusc = tc_make_map(ctx, e_eprim, e, e_prim, tc_tuples_lt(e, e_prim).c_str());
        
    Rusc = isl_map_intersect_domain(Rusc, isl_set_copy(uds));
    Rusc = isl_map_intersect_range(Rusc, uds);
    Rusc = isl_map_intersect_params(Rusc, R_U_R_inv_star_constraints_e_e_prim);
    
    Rusc = tc_map_project_out_params(Rusc, e_eprim);
    Rusc = isl_map_coalesce(Rusc);
    
    isl_id_list_free(e);
    isl_id_list_free(e_prim);
    isl_id_list_free(e_eprim);
    
    return Rusc;
}

__isl_give isl_map* tc_Rusc3_map(__isl_keep isl_set* uds, __isl_keep isl_map* R, __isl_keep isl_map* R_plus, __isl_keep isl_union_map* S)
{
    // R_UCS := { [e] -> [e'] : e,e' in UDS and e << e' and (R^*(e) * R^*(e') <> empty) }
    
    isl_ctx* ctx = isl_map_get_ctx(R);
    
    int n_in_dim = isl_map_n_in(R);
    
    isl_id_list* e = tc_ids_sequence(ctx, "e", n_in_dim);
    isl_id_list* e_prim = tc_ids_prim(e);
    isl_id_list* e_bis = tc_ids_bis(e);
    
    if (NULL == R_plus)
    {
        int exact;
        R_plus = tc_transitive_closure(isl_map_copy(R), S, &exact);
    }
    else
    {
        R_plus = isl_map_copy(R_plus);
    }
    
    isl_map* R_star = isl_map_union(R_plus, tc_make_identity(isl_map_copy(R)));
   
    isl_set* R_star_constraints_e_e_bis = tc_make_map_constraints(isl_map_copy(R_star), e, e_bis);
    isl_set* R_star_constraints_e_prim_e_bis = tc_make_map_constraints(isl_map_copy(R_star), e_prim, e_bis);
        
    isl_id_list* e_eprim_ebis = isl_id_list_concat(isl_id_list_copy(e), isl_id_list_concat(isl_id_list_copy(e_prim), isl_id_list_copy(e_bis)));
    isl_map* Rusc = tc_make_map(ctx, e_eprim_ebis, e, e_prim, tc_tuples_lt(e, e_prim).c_str());
        
    Rusc = isl_map_intersect_domain(Rusc, isl_set_copy(uds));
    Rusc = isl_map_intersect_range(Rusc, isl_set_copy(uds));
    Rusc = isl_map_intersect_params(Rusc, R_star_constraints_e_e_bis);
    Rusc = isl_map_intersect_params(Rusc, R_star_constraints_e_prim_e_bis);
    
    Rusc = tc_map_project_out_params(Rusc, e_eprim_ebis);
    
    isl_id_list_free(e);
    isl_id_list_free(e_prim);
    isl_id_list_free(e_bis);
    isl_id_list_free(e_eprim_ebis);
    isl_map_free(R_star);
    
    return Rusc;
}

__isl_give isl_map* tc_Rusc2_map(__isl_keep isl_set* uds, __isl_keep isl_map* R, __isl_keep isl_map* R_plus, __isl_keep isl_union_map* S)
{
    // R_UCS := { [e] -> [e'] : e,e' in UDS and e << e' and (R^+(e) * R^+(e') <> empty) }
    
    isl_ctx* ctx = isl_map_get_ctx(R);
    
    int n_in_dim = isl_map_n_in(R);
    
    isl_id_list* e = tc_ids_sequence(ctx, "e", n_in_dim);
    isl_id_list* e_prim = tc_ids_prim(e);
    isl_id_list* e_bis = tc_ids_bis(e);
    
    if (NULL == R_plus)
    {
        int exact;
        R_plus = tc_transitive_closure(isl_map_copy(R), S, &exact);
    }
    else
    {
        R_plus = isl_map_copy(R_plus);
    }
        
    isl_set* R_plus_constraints_e_e_bis = tc_make_map_constraints(isl_map_copy(R_plus), e, e_bis);
    isl_set* R_plus_constraints_e_prim_e_bis = tc_make_map_constraints(isl_map_copy(R_plus), e_prim, e_bis);
        
    isl_id_list* e_eprim_ebis = isl_id_list_concat(isl_id_list_copy(e), isl_id_list_concat(isl_id_list_copy(e_prim), isl_id_list_copy(e_bis)));
    isl_map* Rusc = tc_make_map(ctx, e_eprim_ebis, e, e_prim, tc_tuples_lt(e, e_prim).c_str());
        
    Rusc = isl_map_intersect_domain(Rusc, isl_set_copy(uds));
    Rusc = isl_map_intersect_range(Rusc, isl_set_copy(uds));
    Rusc = isl_map_intersect_params(Rusc, R_plus_constraints_e_e_bis);
    Rusc = isl_map_intersect_params(Rusc, R_plus_constraints_e_prim_e_bis);
    
    Rusc = tc_map_project_out_params(Rusc, e_eprim_ebis);
    
    isl_id_list_free(e);
    isl_id_list_free(e_prim);
    isl_id_list_free(e_bis);
    isl_id_list_free(e_eprim_ebis);
    isl_map_free(R_plus);
    
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

__isl_give isl_set* tc_Sk_set(__isl_keep isl_set* LD, __isl_keep isl_map* R, __isl_keep isl_map* R_plus, __isl_keep isl_union_map* S, __isl_give isl_id** k_param)
{
    int exact;
    
    isl_map* R_k = tc_map_power(isl_map_copy(R), S, &exact);
    //isl_map* R_k = isl_map_power(isl_map_copy(R), &exact);

    R_k = isl_map_move_dims(R_k, isl_dim_param, 0, isl_dim_in, 0, 1);
    
    *k_param = isl_map_get_dim_id(R_k, isl_dim_param, 0);
        
    R_k = isl_set_unwrap(isl_map_range(R_k));
    
    tc_debug_map(R_k, "R^k (exact=%d)", exact);
        
    isl_set* uds = tc_uds_set(R);
    isl_set* ind = tc_ind_set(LD, R);
    
    isl_set* uds_ind = isl_set_union(uds, ind);
    uds_ind = isl_set_coalesce(uds_ind);
        
    // Sk = R^k(UDS) - (R^+ . R^k)(UDS)
    isl_set* Sk = isl_set_apply(isl_set_copy(uds_ind), isl_map_copy(R_k));
    Sk = isl_set_coalesce(Sk);
    
    isl_set* uds_ind_apply = isl_set_apply(isl_set_copy(uds_ind), isl_map_copy(R_k));
    uds_ind_apply = isl_set_coalesce(uds_ind_apply);
    Sk = isl_set_subtract(Sk, isl_set_apply(uds_ind_apply, isl_map_copy(R_plus)));
    Sk = isl_set_coalesce(Sk);
            
    uds_ind = isl_set_insert_dims(uds_ind, isl_dim_param, 0, 1);
    uds_ind = isl_set_set_dim_id(uds_ind, isl_dim_param, 0, isl_id_copy(*k_param));
    uds_ind = isl_set_fix_si(uds_ind, isl_dim_param, 0, 0);
    
    Sk = isl_set_union(Sk, uds_ind);
                    
    isl_map_free(R_k);
    
    return Sk;
}

__isl_give isl_map* tc_Rprim_map(__isl_keep isl_map* R)
{
    isl_map* Rprim = isl_map_copy(R);
    
    Rprim = isl_map_insert_dims(Rprim, isl_dim_in, 0, 1);
    Rprim = isl_map_insert_dims(Rprim, isl_dim_out, 0, 1);
        
    isl_constraint* equality = isl_constraint_alloc_equality(isl_local_space_from_space(isl_map_get_space(Rprim)));
    
    equality = isl_constraint_set_coefficient_si(equality, isl_dim_in, 0, 1);
    equality = isl_constraint_set_coefficient_si(equality, isl_dim_out, 0, -1);
    equality = isl_constraint_set_constant_si(equality, 1);
    
    Rprim = isl_map_add_constraint(Rprim, equality);
    
    isl_constraint* inequality = isl_constraint_alloc_inequality(isl_local_space_from_space(isl_map_get_space(Rprim)));
    
    inequality = isl_constraint_set_coefficient_si(inequality, isl_dim_in, 0, 1);
    inequality = isl_constraint_set_constant_si(inequality, 0);
    
    Rprim = isl_map_add_constraint(Rprim, inequality);
    
    return Rprim;
}

__isl_give isl_set* tc_FS_set(__isl_keep isl_set* LD, __isl_keep isl_map* R, __isl_keep isl_union_map* S)
{
    isl_ctx* ctx = isl_map_get_ctx(R);
    
    isl_map* Rprim = tc_Rprim_map(R);
    
    isl_union_map* Sprim = tc_extend_union_map(isl_union_map_copy(S), 1);
    
    tc_debug_umap(Sprim, "S'");
    
    int exact;
    isl_map* Rprim_plus = tc_transitive_closure(isl_map_copy(Rprim), Sprim, &exact);
    isl_map* Rprim_star = isl_map_union(isl_map_copy(Rprim_plus), tc_make_identity(Rprim));
    
    tc_debug_map(Rprim_plus, "R'+ (exact=%d)", exact);
    
    isl_set* uds = tc_uds_set(R);
    isl_set* ind = tc_ind_set(LD, R);
    
    isl_set* uds_ind = isl_set_union(uds, ind);
        
    isl_set* uds_ind_ext = isl_set_insert_dims(uds_ind, isl_dim_set, 0, 1);
    uds_ind_ext = isl_set_fix_dim_si(uds_ind_ext, 0, 0);
    
    int n_dim = isl_set_n_dim(uds_ind_ext);
    
    isl_set* Rprim_plus_uds_ind_ext = isl_set_apply(isl_set_copy(uds_ind_ext), Rprim_plus);
    isl_set* Rprim_star_uds_ind_ext = isl_set_apply(uds_ind_ext, Rprim_star);
            
    // FS = { [X] -> [k, Y] : X in UDS_IND and [k, Y] in Range((R')* \ {[0, X]}) and not exists k' > k s.t. [k', Y] in Range((R')+ \ {[0, X]}) }
    
    // FS = { [k, Y] : [k, Y] in (R')*(UDS_IND) and not exists k' > k s.t. [k', Y] in (R')+(UDS_IND) }
    
    isl_id_list* Y = tc_ids_sequence(ctx, "y", n_dim);
    isl_id_list* Yprim = tc_ids_prim(Y);
    isl_id_list* Y_0_1 = tc_ids_sub(Y, 0, 1);
    isl_id_list* Yprim_0_1 = tc_ids_sub(Yprim, 0, 1);
    isl_id_list* Y_1_n = tc_ids_sub(Y, 1, n_dim);
    isl_id_list* Yprim_1_n = tc_ids_sub(Yprim, 1, n_dim);
        
    isl_set* FS = isl_set_copy(Rprim_star_uds_ind_ext);
    
    // FS' = { [k, Y] : [k, Y] in (R')+(UDS_IND) and exists k' : k' > k and [k',Y] in (R')*(UDS_IND) }
    
    // FS' = [k,k',Y] -> { [k, Y] : [k, Y] in (R')+(UDS_IND) and k' > k and [k',Y] in (R')*(UDS_IND) }
    
    isl_id_list* Y_Yprim = isl_id_list_concat(isl_id_list_copy(Y), isl_id_list_copy(Yprim));
    
    isl_set* FS_uds_range = tc_make_set(ctx, Y_Yprim, Y, tc_conjunction(tc_tuples_gt(Yprim_0_1, Y_0_1), tc_tuples_eq(Y_1_n, Yprim_1_n)).c_str());
        
    isl_set* Rprim_plus_uds_ind_ext_constraints = tc_make_set_constraints(Rprim_plus_uds_ind_ext, Y);
    isl_set* Rprim_star_uds_ind_ext_constraints = tc_make_set_constraints(Rprim_star_uds_ind_ext, Yprim);
    
    FS_uds_range = isl_set_intersect_params(FS_uds_range, Rprim_plus_uds_ind_ext_constraints);
    FS_uds_range = isl_set_intersect_params(FS_uds_range, Rprim_star_uds_ind_ext_constraints);
    
    FS_uds_range = tc_project_out_params(FS_uds_range, Y_Yprim);
    FS_uds_range = isl_set_coalesce(FS_uds_range);
    
    FS = isl_set_subtract(FS, FS_uds_range);
    FS = isl_set_coalesce(FS);
            
    isl_id_list_free(Y);
    isl_id_list_free(Yprim);
    isl_id_list_free(Y_0_1);
    isl_id_list_free(Yprim_0_1);
    isl_id_list_free(Y_1_n);
    isl_id_list_free(Yprim_1_n);
    isl_id_list_free(Y_Yprim);
    isl_union_map_free(Sprim);
    
    return FS;
}
