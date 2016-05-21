#include "tiling.h"
#include "utility.h"
#include "debug.h"

#include <isl/id.h>
#include <isl/set.h>
#include <isl/map.h>

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include <vector>
#include <string>
#include <map>

__isl_give isl_set* tc_tile_set(__isl_keep isl_id_list* II, __isl_keep isl_id_list* I, const std::vector<int>& BLOCK, __isl_keep isl_set* set, __isl_keep isl_map* statement_schedule, __isl_keep isl_union_set* LD, __isl_keep isl_union_map* S)
{
    isl_ctx* ctx = isl_set_get_ctx(set);
    
    isl_set* LD_normalized = tc_normalize_union_set(LD, S);
    isl_set* set_normalized = isl_set_apply(isl_set_copy(set), isl_map_copy(statement_schedule));
        
    int n_LD_normalized_dim = isl_set_dim(set_normalized, isl_dim_set);
            
    isl_id_list* LB = tc_ids_sequence(ctx, "LB", n_LD_normalized_dim);
    isl_id_list* UB = tc_ids_sequence(ctx, "UB", n_LD_normalized_dim);
    
    isl_set* statement_schedule_range = isl_map_range(isl_map_copy(statement_schedule));
            
    char buff[2048];   
    std::string tile_set_str;
    for (int i = 0, j = 0; i < n_LD_normalized_dim; ++i)
    {
        isl_val* val = isl_set_plain_get_val_if_fixed(statement_schedule_range, isl_dim_set, i);
        
        if (isl_bool_true == isl_val_is_nan(val))
        {
            isl_id* LB_i = isl_id_list_get_id(LB, i);
            isl_id* UB_i = isl_id_list_get_id(UB, i);
            isl_id* II_i = isl_id_list_get_id(II, i);
            isl_id* I_i = isl_id_list_get_id(I, i);

            sprintf(buff, "%d * %s + %s <= %s <= min( %d * (%s + 1) + %s - 1, %s ) and %s >= 0 and ",
                       BLOCK[j], isl_id_get_name(II_i), isl_id_get_name(LB_i), isl_id_get_name(I_i), BLOCK[j], isl_id_get_name(II_i), isl_id_get_name(LB_i), isl_id_get_name(UB_i), isl_id_get_name(II_i));
            //sprintf(buff, "%d * %s <= %s <= min( %d * (%s + 1) - 1, %s ) and %s >= 0 and ",
            //           BLOCK[j], isl_id_get_name(II_i), isl_id_get_name(I_i), BLOCK[j], isl_id_get_name(II_i), isl_id_get_name(UB_i), isl_id_get_name(II_i));

            tile_set_str += buff;
            
            j = j + 1;

            isl_id_free(LB_i);
            isl_id_free(UB_i);
            isl_id_free(II_i);
            isl_id_free(I_i);
        }
        
        isl_val_free(val);
    }
    tile_set_str += "true";
    
    statement_schedule_range = tc_parameterize_all(statement_schedule_range, II);
    statement_schedule_range = isl_set_params(statement_schedule_range);
        
    isl_set* bounds = tc_get_set_bounds(LD_normalized, LB, UB);
    //isl_set* bounds = tc_get_set_bounds(set_normalized, LB, UB);    
        
    isl_id_list* LB_UB = isl_id_list_concat(LB, UB);
    
    isl_id_list* LB_UB_II = isl_id_list_concat(isl_id_list_copy(LB_UB), isl_id_list_copy(II));
        
    isl_set* tile = tc_make_set(ctx, LB_UB_II, I, tile_set_str.c_str());
    
    tile = isl_set_intersect(tile, set_normalized);
    tile = isl_set_intersect_params(tile, statement_schedule_range);
    tile = isl_set_intersect_params(tile, bounds);
        
    tile = tc_project_out_params(tile, LB_UB);
        
    tile = isl_set_coalesce(tile);

    isl_set_free(LD_normalized);
    isl_id_list_free(LB_UB);
    isl_id_list_free(LB_UB_II);
    
    return tile;
}
    
__isl_give isl_set* tc_ii_set_set(__isl_keep isl_id_list* II, const std::vector<int>& BLOCK, __isl_keep isl_set* set, __isl_keep isl_map* tile_schedule, __isl_keep isl_union_set* LD, __isl_keep isl_union_map* S)
{    
    isl_ctx* ctx = isl_set_get_ctx(set);
    
    isl_set* LD_normalized = tc_normalize_union_set(LD, S);
    isl_set* set_normalized = isl_set_apply(isl_set_copy(set), isl_map_copy(tile_schedule));
        
    int n_LD_normalized_dim = isl_set_dim(LD_normalized, isl_dim_set);
        
    isl_id_list* LB = tc_ids_sequence(ctx, "LB", n_LD_normalized_dim);
    isl_id_list* UB = tc_ids_sequence(ctx, "UB", n_LD_normalized_dim);
    
    isl_set* statement_schedule_range = isl_map_range(isl_map_copy(tile_schedule));
    
    char buff[2048];   
    std::string ii_set_str;
    for (int i = 0, j = 0; i < n_LD_normalized_dim; ++i)
    {
        isl_val* val = isl_set_plain_get_val_if_fixed(statement_schedule_range, isl_dim_set, i);
        
        if (isl_bool_true == isl_val_is_nan(val))
        {
            isl_id* LB_i = isl_id_list_get_id(LB, i);
            isl_id* UB_i = isl_id_list_get_id(UB, i);
            isl_id* II_i = isl_id_list_get_id(II, i);
        
            sprintf(buff, "%s >= 0 and %d * %s + %s <= %s and ",
                       isl_id_get_name(II_i), BLOCK[j], isl_id_get_name(II_i), isl_id_get_name(LB_i), isl_id_get_name(UB_i));
            //sprintf(buff, "%s >= 0 and %d * %s <= %s and ",
            //           isl_id_get_name(II_i), BLOCK[j], isl_id_get_name(II_i), isl_id_get_name(UB_i));

            ii_set_str += buff;
            
            j = j + 1;

            isl_id_free(LB_i);
            isl_id_free(UB_i);
            isl_id_free(II_i);
        }
        
        isl_val_free(val);
    }
    ii_set_str += "true";
        
    isl_set* bounds = tc_get_set_bounds(LD_normalized, LB, UB);
    //isl_set* bounds = tc_get_set_bounds(set_normalized, LB, UB);

    isl_id_list* LB_UB = isl_id_list_concat(LB, UB);
    
    isl_set* ii_set = tc_make_set(ctx, LB_UB, II, ii_set_str.c_str());
            
    ii_set = isl_set_intersect(ii_set, statement_schedule_range);
    ii_set = isl_set_intersect_params(ii_set, bounds);
    
    ii_set = tc_project_out_params(ii_set, LB_UB);
    
    ii_set = isl_set_coalesce(ii_set);
    
    isl_id_list_free(LB_UB);
    isl_set_free(set_normalized);
    isl_set_free(LD_normalized);
    
    return ii_set;    
}
    
__isl_give isl_set* tc_tile_set_of(__isl_keep isl_set* set, __isl_keep isl_set* ii_set, __isl_keep isl_id_list* II, std::string(*f)(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs))
{
    isl_ctx* ctx = isl_set_get_ctx(set);
    
    isl_id_list* IIprim = tc_ids_prim(II);
    
    isl_set* ii_set_constraints = tc_make_set_constraints(isl_set_copy(ii_set), II);
    
    isl_set* set_prim = tc_rename_params(isl_set_copy(set), II, IIprim);
    isl_set* ii_set_constraints_prim = tc_rename_params(isl_set_copy(ii_set_constraints), II, IIprim);
    
    isl_set* new_set = isl_set_universe(isl_set_get_space(set));
    
    isl_id_list* IIprim_II = isl_id_list_concat(isl_id_list_copy(IIprim), isl_id_list_copy(II));
    
    isl_set* transform_set = tc_make_params(ctx, IIprim_II, f(IIprim, II).c_str());
    
    new_set = isl_set_intersect_params(new_set, transform_set);
    new_set = isl_set_intersect_params(new_set, ii_set_constraints);
    new_set = isl_set_intersect_params(new_set, ii_set_constraints_prim);
    new_set = isl_set_intersect(new_set, set_prim);
    
    new_set = tc_project_out_params(new_set, IIprim);
    
    new_set = isl_set_coalesce(new_set);
    
    isl_id_list_free(IIprim);
    isl_id_list_free(IIprim_II);
    
    return new_set;
}
    
__isl_give isl_set* tc_tile_lt_set(__isl_keep isl_set* tile, __isl_keep isl_set* ii_set, __isl_keep isl_id_list* II)
{
    return tc_tile_set_of(tile, ii_set, II, &tc_tuples_lt);
}

__isl_give isl_set* tc_tile_gt_set(__isl_keep isl_set* tile, __isl_keep isl_set* ii_set, __isl_keep isl_id_list* II)
{
    return tc_tile_set_of(tile, ii_set, II, &tc_tuples_gt);
}

__isl_give isl_map* tc_get_tile_schedule(const char* statement_label, __isl_keep isl_union_map* S, const std::vector<std::vector<std::string> >& groups)
{
    isl_map* tile_schedule = tc_get_map_for_input_tuple(S, statement_label);
    
    isl_id* id = isl_map_get_tuple_id(tile_schedule, isl_dim_in);
    
    for (int i = 0; i < groups.size(); ++i)
    {
        for (int j = 0; j < groups[i].size(); ++j)
        {
            if (0 == strcmp(groups[i][j].c_str(), statement_label))
            {
                isl_map_free(tile_schedule);
                tile_schedule = tc_get_map_for_input_tuple(S, groups[i][0].c_str());
                tile_schedule = isl_map_set_tuple_id(tile_schedule, isl_dim_in, isl_id_copy(id));
                break;
            }
        }
    }
    
    isl_id_free(id);
    
    return tile_schedule;
}

__isl_give isl_union_map* tc_remove_loop_independent_dependences(__isl_take isl_union_map* R, __isl_keep isl_union_map* S, const std::vector<std::vector<std::string> >& groups)
{
    isl_space* space = isl_union_map_get_space(R);
    
    isl_union_map* loop_independent_dependences = isl_union_map_empty(space);
    
    for (int i = 0; i < groups.size(); ++i)
    {
        const std::vector<std::string>& group = groups[i];
        for (int j = 0; j < group.size(); ++j)
        {
            for (int k = j + 1; k < group.size(); ++k)
            {                
                isl_map* tuple_map1 = tc_get_map_for_input_tuple(S, group[j].c_str());                
                isl_map* tuple_map2 = tc_get_map_for_input_tuple(S, group[k].c_str());
                
                isl_map* loop_independent_dependence1 = isl_map_from_domain_and_range(isl_map_domain(isl_map_copy(tuple_map1)), isl_map_domain(isl_map_copy(tuple_map2)));
                isl_map* loop_independent_dependence2 = isl_map_from_domain_and_range(isl_map_domain(isl_map_copy(tuple_map2)), isl_map_domain(isl_map_copy(tuple_map1)));

                for (int l = 0; l < isl_map_n_in(tuple_map1); ++l)
                {
                    loop_independent_dependence1 = isl_map_equate(loop_independent_dependence1, isl_dim_in, l, isl_dim_out, l);
                    loop_independent_dependence2 = isl_map_equate(loop_independent_dependence2, isl_dim_in, l, isl_dim_out, l);
                }
                                
                loop_independent_dependences = isl_union_map_add_map(loop_independent_dependences, loop_independent_dependence1);
                loop_independent_dependences = isl_union_map_add_map(loop_independent_dependences, loop_independent_dependence2);
                
                isl_map_free(tuple_map1);
                isl_map_free(tuple_map2);
            }
        }
    }
    
    R = isl_union_map_subtract(R, loop_independent_dependences);
    
    return R;
}

__isl_give isl_set* tc_normalize_params(__isl_take isl_set* tile, __isl_keep isl_map* tile_schedule, __isl_keep isl_id_list* IIprim, __isl_keep isl_id_list* II)
{            
    isl_ctx* ctx = isl_set_get_ctx(tile);
        
    isl_set* params_input = tc_make_set(ctx, IIprim, IIprim, "");
    
    if (isl_map_has_tuple_id(tile_schedule, isl_dim_in))
    {
        params_input = isl_set_set_tuple_id(params_input, isl_map_get_tuple_id(tile_schedule, isl_dim_in));
    }
    
    params_input = isl_set_apply(params_input, isl_map_copy(tile_schedule));
    
    isl_set* params_rename = tc_make_set(ctx, II, II, "");
    
    params_rename = isl_set_intersect(params_rename, params_input);
    params_rename = isl_set_params(params_rename);
        
    tile = isl_set_intersect_params(tile, params_rename);
    
    tile = tc_project_out_params(tile, IIprim);
    
    tile = isl_set_coalesce(tile);
    
    return tile;
}
    
__isl_give isl_map* tc_Rtile_map(__isl_keep isl_id_list* II, __isl_keep isl_set* tile, __isl_keep isl_map* R)
{
    isl_ctx* ctx = isl_set_get_ctx(tile);
    
    isl_id_list* I = tc_ids_sequence(ctx, "i", isl_set_dim(tile, isl_dim_set));
    isl_id_list* Iprim = tc_ids_prim(I);
    isl_id_list* Ibis = tc_ids_bis(I);
    
    isl_id_list* IIprim = tc_ids_prim(II);
    isl_id_list* IIbis = tc_ids_bis(II);
    
    isl_id_list* Iprim_Ibis_IIprim_IIbis = isl_id_list_concat(isl_id_list_concat(isl_id_list_copy(Iprim), isl_id_list_copy(Ibis)), isl_id_list_concat(isl_id_list_copy(IIprim), isl_id_list_copy(IIbis)));

    isl_map* Rtile = tc_make_map(ctx, Iprim_Ibis_IIprim_IIbis, IIprim, IIbis, "");
    
    isl_set* R_constraints = tc_make_map_constraints(isl_map_copy(R), Iprim, Ibis);
    
    isl_set* tile_prim = tc_rename_params(isl_set_copy(tile), II, IIprim);
    tile_prim = tc_make_set_constraints(tile_prim, Iprim);
    
    isl_set* tile_bis = tc_rename_params(isl_set_copy(tile), II, IIbis);
    tile_bis = tc_make_set_constraints(tile_bis, Ibis);
    
    Rtile = isl_map_intersect_params(Rtile, R_constraints);
    Rtile = isl_map_intersect_params(Rtile, tile_prim);
    Rtile = isl_map_intersect_params(Rtile, tile_bis);
    
    Rtile = tc_map_project_out_params(Rtile, Iprim_Ibis_IIprim_IIbis);
    
    Rtile = isl_map_subtract(Rtile, tc_make_identity(isl_map_copy(Rtile)));
        
    Rtile = isl_map_coalesce(Rtile);
    
    isl_id_list_free(I);
    isl_id_list_free(Iprim);
    isl_id_list_free(Ibis);
    isl_id_list_free(IIprim);
    isl_id_list_free(IIbis);
    isl_id_list_free(Iprim_Ibis_IIprim_IIbis);
    
    return Rtile;
}

__isl_give isl_map* tc_Tcycle_map(__isl_keep isl_id_list* II, __isl_keep isl_map* R_plus)
{
    // T_CYCLE := { [II] -> [JJ] : II in R_TILE^+(JJ) and JJ in R_TILE^+(II) and II << JJ  }
    
    isl_ctx* ctx = isl_map_get_ctx(R_plus);
        
    isl_id_list* IIprim = tc_ids_prim(II);
    isl_id_list* IIbis = tc_ids_bis(II);
    
    isl_id_list* IIprim_IIbis = isl_id_list_concat(isl_id_list_copy(IIprim), isl_id_list_copy(IIbis));

    isl_map* Tcycle = tc_make_map(ctx, IIprim_IIbis, IIprim, IIbis, tc_tuples_lt(IIprim, IIbis).c_str());
    
    isl_set* R_plus_constraints_prim_bis = tc_make_map_constraints(isl_map_copy(R_plus), IIprim, IIbis);
    isl_set* R_plus_constraints_bis_prim = tc_make_map_constraints(isl_map_copy(R_plus), IIbis, IIprim);
        
    Tcycle = isl_map_intersect_params(Tcycle, R_plus_constraints_prim_bis);
    Tcycle = isl_map_intersect_params(Tcycle, R_plus_constraints_bis_prim);
    
    Tcycle = tc_map_project_out_params(Tcycle, IIprim_IIbis);
            
    //Tcycle = isl_map_coalesce(Tcycle);
    
    isl_id_list_free(IIprim);
    isl_id_list_free(IIbis);
    isl_id_list_free(IIprim_IIbis);
    
    return Tcycle;
}

__isl_give isl_set* tc_tile_m_set(__isl_keep isl_id_list* II, __isl_keep isl_set* tile, __isl_keep isl_set* ii_set_m, __isl_keep isl_map* Tcycle)
{
    // TILE_M := [II] -> { [I] : (exists JJ : II in II_SET_M and JJ in T_CYCLE^*(II) and I in TILE(JJ)) }

    isl_ctx* ctx = isl_set_get_ctx(tile);
        
    isl_id_list* I = tc_ids_sequence(ctx, "i", isl_set_n_dim(tile));
    
    isl_id_list* IIprim = tc_ids_prim(II);
    
    isl_id_list* II_IIprim = isl_id_list_concat(isl_id_list_copy(II), isl_id_list_copy(IIprim));
    
    isl_set* tile_m = tc_make_set(ctx, II_IIprim, I, NULL);
    
    isl_set* tile_prim = tc_rename_params(isl_set_copy(tile), II, IIprim);
    
    isl_set* ii_set_m_constraints = tc_make_set_constraints(isl_set_copy(ii_set_m), II);
        
    int exact;
    isl_map* Tcycle_plus = isl_map_transitive_closure(isl_map_copy(Tcycle), &exact);
    isl_map* Tcycle_star = isl_map_union(Tcycle_plus, tc_make_identity(isl_map_copy(Tcycle)));
    
    isl_set* Tcycle_star_constraints = tc_make_map_constraints(Tcycle_star, II, IIprim);
    
    tile_m = isl_set_intersect_params(tile_m, ii_set_m_constraints);
    tile_m = isl_set_intersect_params(tile_m, Tcycle_star_constraints);
    tile_m = isl_set_intersect(tile_m, tile_prim);
    
    tile_m = tc_project_out_params(tile_m, IIprim);
    
    //tile_m = isl_set_coalesce(tile_m);
    
    isl_id_list_free(I);
    isl_id_list_free(IIprim);
    isl_id_list_free(II_IIprim);
    
    return tile_m;
}

isl_bool tc_tile_check_vld(__isl_keep isl_set* tile, __isl_keep isl_set* ii_set, __isl_keep isl_id_list* II, __isl_keep isl_map* R_plus)
{
    isl_set* tile_gt = tc_tile_gt_set(tile, ii_set, II);
    
    // CHECK_VLD = TILE * R+(TILE_GT)
    isl_set* check_vld = isl_set_intersect(isl_set_copy(tile), isl_set_apply(tile_gt, isl_map_copy(R_plus)));
    
    isl_bool is_valid = isl_set_is_empty(check_vld);
    
    isl_set_free(check_vld);
    
    return is_valid;
}

void tc_tile_loop_nest(__isl_keep isl_union_set* LD, __isl_keep isl_union_map* S, __isl_keep isl_id_list* II, __isl_keep isl_id_list* I, __isl_give isl_set** tiles, __isl_give isl_set** ii_sets, const std::map<std::string, std::vector<int> >& blocks, const std::vector<std::vector<std::string> >& groups)
{
    *tiles = NULL;
    *ii_sets = NULL;
    
    isl_set_list* statements_domains = tc_collect_sets(LD);
    
    for (int i = 0; i < isl_set_list_n_set(statements_domains); ++i)
    {
        isl_set* statement_domain = isl_set_list_get_set(statements_domains, i);
                  
        isl_id* statement_id = isl_set_get_tuple_id(statement_domain);
        const char* statement_label = isl_id_get_name(statement_id);

        isl_map* statement_schedule = tc_get_map_for_input_tuple(S, statement_label);

        std::vector<int> block;
        if (blocks.count(statement_label) > 0)
        {
            block = blocks.at(statement_label);
        }
        else
        {
            block = blocks.at("__DEFAULT__");
        }

        isl_map* tile_schedule = tc_get_tile_schedule(statement_label, S, groups);

        isl_set* tile = tc_tile_set(II, I, block, statement_domain, statement_schedule, LD, S);

        if (NULL == *tiles)
        {
            *tiles = tile;
        }
        else
        {
            *tiles = isl_set_union(*tiles, tile);
        }

        isl_set* ii_set = tc_ii_set_set(II, block, statement_domain, tile_schedule, LD, S);

        if (NULL == *ii_sets)
        {
            *ii_sets = ii_set;
        }
        else
        {
            *ii_sets = isl_set_union(*ii_sets, ii_set);
        }

        isl_map_free(statement_schedule);
        isl_map_free(tile_schedule);
        isl_id_free(statement_id);
        isl_set_free(statement_domain);
    }
    
    isl_set_list_free(statements_domains);
}
