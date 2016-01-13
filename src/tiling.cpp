#include "tiling.h"
#include "utility.h"

#include <isl/id.h>
#include <isl/set.h>
#include <isl/map.h>
#include <isl/union_map.h>

#include <stdio.h>
#include <stddef.h>

#include <vector>
#include <string>
#include <map>

__isl_give isl_set* tc_tile_set(__isl_keep isl_id_list* II, __isl_keep isl_id_list* I, const std::vector<int>& BLOCK, __isl_keep isl_set* set)
{
    int n_set_dim = isl_set_dim(set, isl_dim_set);
    
    isl_ctx* ctx = isl_set_get_ctx(set);
    
    isl_id_list* LB = tc_ids_sequence(ctx, "LB", n_set_dim);
    isl_id_list* UB = tc_ids_sequence(ctx, "UB", n_set_dim);
    
    char buff[2048];   
    std::string tile_set_str;
    for (int i = 0; i < isl_id_list_n_id(I); ++i)
    {
        isl_id* LB_i = isl_id_list_get_id(LB, i);
        isl_id* UB_i = isl_id_list_get_id(UB, i);
        isl_id* II_i = isl_id_list_get_id(II, i);
        isl_id* I_i = isl_id_list_get_id(I, i);
        
        sprintf(buff, "%d * %s + %s <= %s <= min( %d * (%s + 1) + %s - 1, %s ) and %s >= 0 and ",
                   BLOCK[i], isl_id_get_name(II_i), isl_id_get_name(LB_i), isl_id_get_name(I_i), BLOCK[i], isl_id_get_name(II_i), isl_id_get_name(LB_i), isl_id_get_name(UB_i), isl_id_get_name(II_i));

        tile_set_str += buff;
        
        isl_id_free(LB_i);
        isl_id_free(UB_i);
        isl_id_free(II_i);
        isl_id_free(I_i);
    }
    tile_set_str += "true";
    
    isl_set* bounds = tc_get_set_bounds(set, LB, UB);
    
    isl_id_list* LB_UB = isl_id_list_concat(LB, UB);
    isl_id_list* LB_UB_II = isl_id_list_concat(isl_id_list_copy(LB_UB), isl_id_list_copy(II));
    
    isl_set* tile = tc_make_set(ctx, LB_UB_II, I, tile_set_str.c_str());
    
    if (isl_set_has_tuple_id(set))
    {
        tile = isl_set_set_tuple_id(tile, isl_set_get_tuple_id(set));
    }

    tile = isl_set_intersect(tile, isl_set_copy(set));
    tile = isl_set_intersect_params(tile, bounds);
    
    tile = tc_project_out_params(tile, LB_UB);
    
    tile = isl_set_coalesce(tile);
        
    isl_id_list_free(LB_UB);
    isl_id_list_free(LB_UB_II);
    
    return tile;
}
    
__isl_give isl_set* tc_ii_set_set(__isl_keep isl_id_list* II, const std::vector<int>& BLOCK, __isl_keep isl_set* set)
{
    int n_set_dim = isl_set_dim(set, isl_dim_set);
    
    isl_ctx* ctx = isl_set_get_ctx(set);
    
    isl_id_list* LB = tc_ids_sequence(ctx, "LB", n_set_dim);
    isl_id_list* UB = tc_ids_sequence(ctx, "UB", n_set_dim);
    
    char buff[2048];   
    std::string ii_set_str;
    for (int i = 0; i < isl_id_list_n_id(II); ++i)
    {
        isl_id* LB_i = isl_id_list_get_id(LB, i);
        isl_id* UB_i = isl_id_list_get_id(UB, i);
        isl_id* II_i = isl_id_list_get_id(II, i);
        
        sprintf(buff, "%s >= 0 and %d * %s + %s <= %s and ",
                   isl_id_get_name(II_i), BLOCK[i], isl_id_get_name(II_i), isl_id_get_name(LB_i), isl_id_get_name(UB_i));

        ii_set_str += buff;
        
        isl_id_free(LB_i);
        isl_id_free(UB_i);
        isl_id_free(II_i);
    }
    ii_set_str += "true";
    
    isl_set* bounds = tc_get_set_bounds(set, LB, UB);
    
    isl_id_list* LB_UB = isl_id_list_concat(LB, UB);
    
    isl_set* ii_set = tc_make_set(ctx, LB_UB, II, ii_set_str.c_str());
    
    if (isl_set_has_tuple_id(set))
    {
        ii_set = isl_set_set_tuple_id(ii_set, isl_set_get_tuple_id(set));
    }
        
    ii_set = isl_set_intersect_params(ii_set, bounds);
    
    ii_set = tc_project_out_params(ii_set, LB_UB);
    
    ii_set = isl_set_coalesce(ii_set);
    
    isl_id_list_free(LB_UB);
    
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

__isl_give isl_map* tc_get_tile_schedule(const char* statement_label, /*const std::map<std::string, std::vector<std::string> >& groups,*/ __isl_keep isl_union_map* S)
{
    isl_map* tile_schedule = tc_get_map_for_input_tuple(S, statement_label);
    return tile_schedule;
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
    
__isl_give isl_map* tc_Rtile_map(__isl_keep isl_set* tile, __isl_keep isl_map* R)
{
    isl_ctx* ctx = isl_set_get_ctx(tile);
    
    isl_id_list* I = tc_ids_sequence(ctx, "i", isl_set_dim(tile, isl_dim_set));
    isl_id_list* Iprim = tc_ids_prim(I);
    isl_id_list* Ibis = tc_ids_bis(I);
    
    isl_id_list* II = tc_ids_sequence(ctx, "ii", isl_set_dim(tile, isl_dim_set));
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
    isl_id_list_free(II);
    isl_id_list_free(IIprim);
    isl_id_list_free(IIbis);
    isl_id_list_free(Iprim_Ibis_IIprim_IIbis);
    
    return Rtile;
}

struct tc_tile_loop_nest_user
{    
    isl_union_map* S;
    
    isl_id_list* I;
    
    isl_id_list* II;
    
    isl_set* tile;
    
    isl_set* ii_set;
    
    std::map<std::string, std::vector<int> > blocks;
};

static isl_stat tc_tile_loop_nest_callback(__isl_take isl_set* statement_domain, void* user)
{
    tc_tile_loop_nest_user* data = (tc_tile_loop_nest_user*)user;
         
    isl_union_map* S = data->S;    
    isl_id_list* I = data->I;
    isl_id_list* II = data->II;
    isl_id_list* Iprim = tc_ids_prim(I);
    isl_id_list* IIprim = tc_ids_prim(II);    
    std::map<std::string, std::vector<int> >& blocks = data->blocks;
    
    isl_id* statement_id = isl_set_get_tuple_id(statement_domain);
    const char* statement_label = isl_id_get_name(statement_id);
    
    isl_map* statement_schedule = tc_get_map_for_input_tuple(S, statement_label);
    
    int n_set_dim = isl_set_dim(statement_domain, isl_dim_set);
    
    std::vector<int> block;
    if (blocks.count(statement_label) > 0)
    {
        block = blocks[statement_label];
    }
    else
    {
        block = blocks["__DEFAULT__"];
    }
    
    isl_map* tile_schedule = tc_get_tile_schedule(statement_label, S);
    
    isl_id_list* Iprim_0_n_set_dim = tc_ids_sub(Iprim, 0, n_set_dim);
    isl_id_list* II_0_n_set_dim = tc_ids_sub(II, 0, n_set_dim);
    isl_id_list* IIprim_0_n_set_dim = tc_ids_sub(IIprim, 0, n_set_dim);
    
    isl_set* tile = tc_tile_set(IIprim_0_n_set_dim, Iprim_0_n_set_dim, block, statement_domain);
    
    tile = isl_set_apply(tile, statement_schedule);

    tile = tc_normalize_params(tile, tile_schedule, IIprim_0_n_set_dim, II);
    
    tile = isl_set_coalesce(tile);
    
    if (NULL == data->tile)
    {
        data->tile = tile;
    }
    else
    {
        data->tile = isl_set_union(data->tile, tile);
    }    
            
    isl_set* ii_set = tc_ii_set_set(II_0_n_set_dim, block, statement_domain);
    
    ii_set = isl_set_apply(ii_set, tile_schedule);
    
    ii_set = isl_set_coalesce(ii_set);
    
    if (NULL == data->ii_set)
    {
        data->ii_set = ii_set;
    }
    else
    {
        data->ii_set = isl_set_union(data->ii_set, ii_set);
    }
    
    isl_id_list_free(Iprim);
    isl_id_list_free(IIprim);
    isl_id_list_free(Iprim_0_n_set_dim);
    isl_id_list_free(II_0_n_set_dim);
    isl_id_list_free(IIprim_0_n_set_dim);
    isl_id_free(statement_id);
    isl_set_free(statement_domain);
}
    
void tc_tile_loop_nest(__isl_keep isl_union_set* LD, __isl_keep isl_union_map* S, const std::map<std::string, std::vector<int> >& blocks, __isl_keep isl_id_list* II, __isl_keep isl_id_list* I, __isl_give isl_set** tile, __isl_give isl_set** ii_set)
{
    tc_tile_loop_nest_user data;
    data.S = S;
    data.I = I;
    data.II = II;
    data.tile = NULL;
    data.ii_set = NULL;    
    data.blocks = blocks;
    
    isl_union_set_foreach_set(LD, &tc_tile_loop_nest_callback, &data);
    
    *tile = data.tile;
    *ii_set = data.ii_set;
}
