#include "utility.h"
#include "debug.h"

#include <isl/ctx.h>
#include <isl/id.h>
#include <isl/set.h>
#include <isl/map.h>
#include <isl/union_set.h>
#include <isl/union_map.h>
#include <isl/polynomial.h>
#include <isl/point.h>
#include <isl/val.h>

#include <barvinok/isl.h>

#include <stdio.h>
#include <stddef.h>
#include <string.h>
#include <stdlib.h>

#include <string>

__isl_give isl_set* tc_rename_dim(__isl_take isl_set* set, isl_dim_type dim, __isl_keep isl_id_list* from, __isl_keep isl_id_list* to)
{    
    for (int i = 0; i < isl_id_list_n_id(from); ++i)
    {
        isl_id* id = isl_id_list_get_id(from, i);
        
        int pos = isl_set_find_dim_by_id(set, dim, id);
        
        set = isl_set_set_dim_id(set, dim, pos, isl_id_list_get_id(to, i));
        
        isl_id_free(id);
    }
    
    return set;
}

__isl_give isl_set* tc_rename_params(__isl_take isl_set* set, __isl_keep isl_id_list* from, __isl_keep isl_id_list* to)
{
    return tc_rename_dim(set, isl_dim_param, from, to);
}

__isl_give isl_map* tc_map_rename_dim(__isl_take isl_map* map, isl_dim_type dim, __isl_keep isl_id_list* from, __isl_keep isl_id_list* to)
{    
    for (int i = 0; i < isl_id_list_n_id(from); ++i)
    {
        isl_id* id = isl_id_list_get_id(from, i);
        
        int pos = isl_map_find_dim_by_id(map, dim, id);
        
        map = isl_map_set_dim_id(map, dim, pos, isl_id_list_get_id(to, i));
        
        isl_id_free(id);
    }
    
    return map;
}

__isl_give isl_map* tc_map_rename_params(__isl_take isl_map* map, __isl_keep isl_id_list* from, __isl_keep isl_id_list* to)
{
    return tc_map_rename_dim(map, isl_dim_param, from, to);
}

__isl_give isl_set* tc_remove_dim_names(__isl_take isl_set* set, isl_dim_type dim, __isl_keep isl_id_list* names)
{    
    for (int i = 0; i < isl_id_list_n_id(names); ++i)
    {
        isl_id* id = isl_id_list_get_id(names, i);
        
        int pos = isl_set_find_dim_by_id(set, dim, id);
        
        set = isl_set_remove_dims(set, dim, pos, 1);
        
        isl_id_free(id);
    }
    
    return set;
}

__isl_give isl_map* tc_map_remove_dim_names(__isl_take isl_map* map, isl_dim_type dim, __isl_keep isl_id_list* names)
{
    for (int i = 0; i < isl_id_list_n_id(names); ++i)
    {
        isl_id* id = isl_id_list_get_id(names, i);

        int pos = isl_map_find_dim_by_id(map, dim, id);

        map = isl_map_remove_dims(map, dim, pos, 1);

        isl_id_free(id);
    }

    return map;
}

__isl_give isl_set* tc_remove_params(__isl_take isl_set* set, __isl_keep isl_id_list* names)
{
    return tc_remove_dim_names(set, isl_dim_param, names);
}

__isl_give isl_map* tc_map_remove_params(__isl_take isl_map* map, __isl_keep isl_id_list* names)
{
    return tc_map_remove_dim_names(map, isl_dim_param, names);
}

__isl_give isl_union_set* tc_union_set_remove_params(__isl_take isl_union_set* uset, __isl_keep isl_id_list* names)
{
    isl_set_list* sets = tc_collect_sets(uset);

    isl_union_set* result = NULL;

    for (int i = 0; i < isl_set_list_n_set(sets); ++i)
    {
        isl_set* set = isl_set_list_get_set(sets, i);
        set = tc_remove_params(set, names);

        if (NULL == result)
        {
            result = isl_union_set_from_set(set);
        }
        else
        {
            result = isl_union_set_add_set(result, set);
        }
    }

    isl_union_set_free(uset);
    isl_set_list_free(sets);

    return result;
}

__isl_give isl_union_map* tc_union_map_remove_params(__isl_take isl_union_map* umap, __isl_keep isl_id_list* names)
{
    isl_map_list* maps = tc_collect_maps(umap);

    isl_union_map* result = NULL;

    for (int i = 0; i < isl_map_list_n_map(maps); ++i)
    {
        isl_map* map = isl_map_list_get_map(maps, i);
        map = tc_map_remove_params(map, names);

        if (NULL == result)
        {
            result = isl_union_map_from_map(map);
        }
        else
        {
            result = isl_union_map_add_map(result, map);
        }
    }

    isl_map_list_free(maps);

    return result;
}

__isl_give isl_set* tc_project_out_dim_names(__isl_take isl_set* set, isl_dim_type dim, __isl_keep isl_id_list* names)
{
    for (int i = 0; i < isl_id_list_n_id(names); ++i)
    {
        isl_id* id = isl_id_list_get_id(names, i);

        int pos = isl_set_find_dim_by_id(set, dim, id);

        if (pos >= 0)
        {
            set = isl_set_project_out(set, dim, pos, 1);
        }

        isl_id_free(id);
    }

    return set;
}

__isl_give isl_set* tc_project_out_params(__isl_take isl_set* set, __isl_keep isl_id_list* names)
{
    return tc_project_out_dim_names(set, isl_dim_param, names);
}

__isl_give isl_set* tc_project_out_dims_except_pos(__isl_take isl_set* set, isl_dim_type dim, int pos)
{    
    set = isl_set_project_out(set, dim, 0, pos);
    set = isl_set_project_out(set, dim, 1, isl_set_dim(set, dim) - 1);
    return set;
}

__isl_give isl_set* tc_project_out_params_except_pos(__isl_take isl_set* set, int pos)
{    
    return tc_project_out_dims_except_pos(set, isl_dim_param, pos);
}

__isl_give isl_set* tc_project_out_dims_except(__isl_take isl_set* set, isl_dim_type dim, __isl_take isl_id* id)
{    
    int pos = isl_set_find_dim_by_id(set, dim, id);
    isl_id_free(id);
    return tc_project_out_dims_except_pos(set, dim, pos);
}

__isl_give isl_set* tc_project_out_params_except(__isl_take isl_set* set, __isl_take isl_id* id)
{    
    return tc_project_out_dims_except(set, isl_dim_param, id);
}

__isl_give isl_union_set* tc_union_set_project_out_dim_names(__isl_take isl_union_set* uset, isl_dim_type dim, __isl_keep isl_id_list* names)
{
    isl_set_list* sets = tc_collect_sets(uset);
    
    isl_union_set* new_uset = NULL;
    
    for (int i = 0; i < isl_set_list_n_set(sets); ++i)
    {
        isl_set* set = isl_set_list_get_set(sets, i);
        
        set = tc_project_out_dim_names(set, dim, names);
        
        if (NULL == new_uset)
        {
            new_uset = isl_union_set_from_set(set);
        }
        else
        {
            new_uset = isl_union_set_add_set(new_uset, set);
        }
    }
    
    isl_union_set_free(uset);
    isl_set_list_free(sets);
    
    return new_uset;
}

__isl_give isl_union_set* tc_union_set_project_out_params(__isl_take isl_union_set* uset, __isl_keep isl_id_list* names)
{
    return tc_union_set_project_out_dim_names(uset, isl_dim_param, names);
}

__isl_give isl_map* tc_map_project_out_dim_names(__isl_take isl_map* map, isl_dim_type dim, __isl_keep isl_id_list* names)
{    
    for (int i = 0; i < isl_id_list_n_id(names); ++i)
    {
        isl_id* id = isl_id_list_get_id(names, i);
        
        int pos = isl_map_find_dim_by_id(map, dim, id);
        
        map = isl_map_project_out(map, dim, pos, 1);
        
        isl_id_free(id);
    }
    
    return map;
}

__isl_give isl_map* tc_map_project_out_params(__isl_take isl_map* map, __isl_keep isl_id_list* names)
{
    return tc_map_project_out_dim_names(map, isl_dim_param, names);
}

__isl_give isl_union_map* tc_union_map_project_out_dim_names(__isl_take isl_union_map* umap, isl_dim_type dim, __isl_keep isl_id_list* names)
{
    isl_map_list* maps = tc_collect_maps(umap);
    
    isl_union_map* new_umap = NULL;
    
    for (int i = 0; i < isl_map_list_n_map(maps); ++i)
    {
        isl_map* map = isl_map_list_get_map(maps, i);
        
        map = tc_map_project_out_dim_names(map, dim, names);
        
        if (NULL == new_umap)
        {
            new_umap = isl_union_map_from_map(map);
        }
        else
        {
            new_umap = isl_union_map_add_map(new_umap, map);
        }
    }
    
    isl_union_map_free(umap);
    isl_map_list_free(maps);
    
    return new_umap;
}

__isl_give isl_union_map* tc_union_map_project_out_params(__isl_take isl_union_map* umap, __isl_keep isl_id_list* names)
{
    return tc_union_map_project_out_dim_names(umap, isl_dim_param, names);
}

__isl_give isl_map* tc_map_set_dim_names(__isl_take isl_map* map, isl_dim_type dim, __isl_keep isl_id_list* names)
{
    for (int i = 0; i < isl_map_dim(map, dim); ++i)
    {
        map = isl_map_set_dim_id(map, dim, i, isl_id_list_get_id(names, i));
    }
    
    return map;
}

isl_stat tc_flatten_union_set_callback(__isl_take isl_set* set, void* user) 
{
    isl_set** result = (isl_set**)user;
    
    if (NULL == *result)
    {
        *result = set;
    }
    else
    {
        *result = isl_set_union(*result, set);
    }
    
    return isl_stat_ok;
}

__isl_give isl_set* tc_flatten_union_set(__isl_take isl_union_set* uset)
{
    isl_set* result = NULL;
    
    if (isl_union_set_is_empty(uset))
    {
        result = isl_set_empty(isl_union_set_get_space(uset));
    }
    else
    {
        isl_union_set_foreach_set(uset, &tc_flatten_union_set_callback, &result);
    }
    
    isl_union_set_free(uset);
    
    return result;
}

isl_stat tc_flatten_union_map_callback(__isl_take isl_map* map, void* user) 
{
    isl_map** result = (isl_map**)user;
    
    if (NULL == *result)
    {
        *result = map;
    }
    else
    {
        *result = isl_map_union(*result, map);
    }
    
    return isl_stat_ok;
}

__isl_give isl_map* tc_flatten_union_map(__isl_take isl_union_map* umap)
{
    isl_map* result = NULL;
    
    if (isl_union_map_is_empty(umap))
    {
        result = isl_map_empty(isl_union_map_get_space(umap));
    }
    else
    {
        isl_union_map_foreach_map(umap, &tc_flatten_union_map_callback, &result);
    }
    
    isl_union_map_free(umap);
    
    return result;
}

__isl_give isl_set* tc_lift_up_set_params(__isl_take isl_set* set, __isl_keep isl_id_list* params)
{    
    for (int i = 0; i < isl_id_list_n_id(params); ++i)
    {
        isl_id* id = isl_id_list_get_id(params, i);
        
        int pos = isl_set_find_dim_by_id(set, isl_dim_param, id);
        
        set = isl_set_move_dims(set, isl_dim_set, i, isl_dim_param, pos, 1);
        
        isl_id_free(id);
    }
    
    return set;
}

__isl_give isl_map* tc_lift_up_map_params(__isl_take isl_map* map, __isl_keep isl_id_list* params, isl_dim_type dim)
{    
    for (int i = 0; i < isl_id_list_n_id(params); ++i)
    {
        isl_id* id = isl_id_list_get_id(params, i);
        
        int pos = isl_map_find_dim_by_id(map, isl_dim_param, id);
        
        map = isl_map_move_dims(map, dim, i, isl_dim_param, pos, 1);
        
        isl_id_free(id);
    }
    
    return map;
}

__isl_give isl_set* tc_lift_down_all_set_vars(__isl_take isl_set* set, __isl_keep isl_id_list* vars)
{    
    for (int i = isl_id_list_n_id(vars) - 1; i >= 0; --i)
    {
        set = isl_set_set_dim_id(set, isl_dim_set, i, isl_id_list_get_id(vars, i));
        
        set = isl_set_move_dims(set, isl_dim_param, 0, isl_dim_set, i, 1);
    }
    
    return set;
}

__isl_give isl_set* tc_lift_down_set_vars(__isl_take isl_set* set, __isl_keep isl_id_list* vars)
{    
    for (int i = 0; i < isl_id_list_n_id(vars); ++i)
    {
        set = isl_set_set_dim_id(set, isl_dim_set, 0, isl_id_list_get_id(vars, i));
        
        set = isl_set_move_dims(set, isl_dim_param, i, isl_dim_set, 0, 1);
    }
    
    return set;
}

__isl_give isl_map* tc_lift_down_all_map_vars(__isl_take isl_map* map
                                            , __isl_keep isl_id_list* in
                                            , __isl_keep isl_id_list* out)
{    
    for (int i = isl_id_list_n_id(in) - 1; i >= 0; --i)
    {
        map = isl_map_set_dim_id(map, isl_dim_in, i, isl_id_list_get_id(in, i));
        
        map = isl_map_move_dims(map, isl_dim_param, 0, isl_dim_in, i, 1);
    }
        
    for (int i = isl_id_list_n_id(out) - 1; i >= 0; --i)
    {
        map = isl_map_set_dim_id(map, isl_dim_out, i, isl_id_list_get_id(out, i));
        
        map = isl_map_move_dims(map, isl_dim_param, 0, isl_dim_out, i, 1);
    }
    
    return map;
}

std::string tc_comma(__isl_keep isl_id_list* list)
{
    std::string str;
    
    for (int i = 0; i < isl_id_list_n_id(list); ++i)
    {
        isl_id* id = isl_id_list_get_id(list, i);
        str += std::string(isl_id_get_name(id)) + ',';
        isl_id_free(id);
    }
    
    if (str.length() > 0)
    {
        str.erase(str.size() - 1);
    }
    
    return str;
}

std::string tc_parens_comma(__isl_keep isl_id_list* list)
{
    std::string str;
    
    for (int i = 0; i < isl_id_list_n_id(list); ++i)
    {
        isl_id* id = isl_id_list_get_id(list, i);
        str += "(" + std::string(isl_id_get_name(id)) + "),";
        isl_id_free(id);
    }
    
    if (str.length() > 0)
    {
        str.erase(str.size() - 1);
    }
    
    return str;
}

__isl_give isl_set* tc_make_set(__isl_keep isl_ctx* ctx
                              , __isl_keep isl_id_list* params
                              , __isl_keep isl_id_list* vars
                              , const char* constraints)
{
    char set_str[2048];
    
    if (NULL == constraints || 0 == strlen(constraints))
    {
        snprintf(set_str, sizeof(set_str), "[%s] -> { [%s] }", NULL != params ? tc_comma(params).c_str() : "", tc_comma(vars).c_str());
    }
    else
    {
        snprintf(set_str, sizeof(set_str), "[%s] -> { [%s] : %s }", NULL != params ? tc_comma(params).c_str() : "", tc_comma(vars).c_str(), constraints);
    }
    
    isl_set* set = isl_set_read_from_str(ctx, set_str);
    
    if (NULL != params)
    {
        for (int i = 0; i < isl_id_list_n_id(params); ++i)
        {
            set = isl_set_set_dim_id(set, isl_dim_param, i, isl_id_list_get_id(params, i));
        }
    }
    
    for (int i = 0; i < isl_id_list_n_id(vars); ++i)
    {
        set = isl_set_set_dim_id(set, isl_dim_set, i, isl_id_list_get_id(vars, i));
    }
    
    return set;
}

__isl_give isl_union_set* tc_make_union_set(__isl_keep isl_ctx* ctx
                                          , __isl_keep isl_id_list* params
                                          , __isl_keep isl_id_list* vars
                                          , const char* constraints)
{
    return isl_union_set_from_set(tc_make_set(ctx, params, vars, constraints));
}

__isl_give isl_map* tc_make_map(__isl_keep isl_ctx* ctx
                              , __isl_keep isl_id_list* params
                              , __isl_keep isl_id_list* in
                              , __isl_keep isl_id_list* out
                              , const char* constraints)
{
    char map_str[2048];
    
    if (NULL == constraints || 0 == strlen(constraints))
    {
        snprintf(map_str, sizeof(map_str), "[%s] -> { [%s] -> [%s] }", NULL != params ? tc_comma(params).c_str() : "", tc_comma(in).c_str(), tc_comma(out).c_str());
    }
    else
    {
        snprintf(map_str, sizeof(map_str), "[%s] -> { [%s] -> [%s] : %s }", NULL != params ? tc_comma(params).c_str() : "", tc_comma(in).c_str(), tc_comma(out).c_str(), constraints);
    }
    
    isl_map* map = isl_map_read_from_str(ctx, map_str);
    
    if (NULL != params)
    {
        for (int i = 0; i < isl_id_list_n_id(params); ++i)
        {
            map = isl_map_set_dim_id(map, isl_dim_param, i, isl_id_list_get_id(params, i));
        }
    }
    
    for (int i = 0; i < isl_id_list_n_id(in); ++i)
    {
        map = isl_map_set_dim_id(map, isl_dim_in, i, isl_id_list_get_id(in, i));
    }
    
    for (int i = 0; i < isl_id_list_n_id(out); ++i)
    {
        map = isl_map_set_dim_id(map, isl_dim_out, i, isl_id_list_get_id(out, i));
    }
    
    return map;
}

__isl_give isl_union_map* tc_make_union_map(__isl_keep isl_ctx* ctx
                                          , __isl_keep isl_id_list* params
                                          , __isl_keep isl_id_list* in
                                          , __isl_keep isl_id_list* out
                                          , const char* constraints)
{
    return isl_union_map_from_map(tc_make_map(ctx, params, in, out, constraints));
}

__isl_give isl_set* tc_make_params(__isl_keep isl_ctx* ctx
                                 , __isl_keep isl_id_list* params
                                 , const char* constraints)
{
    char set_str[2048];
    
    snprintf(set_str, sizeof(set_str), "[%s] -> { : %s }", tc_comma(params).c_str(), constraints);
    
    isl_set* set = isl_set_read_from_str(ctx, set_str);    
    
    for (int i = 0; i < isl_id_list_n_id(params); ++i)
    {
        set = isl_set_set_dim_id(set, isl_dim_param, i, isl_id_list_get_id(params, i));
    }
    
    return set;
}

__isl_give isl_set* tc_parameterize(__isl_take isl_set* set, int pos, __isl_take isl_id* name)
{
    set = isl_set_insert_dims(set, isl_dim_param, 0, 1);
    
    set = isl_set_set_dim_id(set, isl_dim_param, 0, name);
    
    set = isl_set_equate(set, isl_dim_param, 0, isl_dim_set, pos);
    
    return set;
}

__isl_give isl_set* tc_parameterize_all(__isl_take isl_set* set, __isl_keep isl_id_list* names)
{
    int n_len = isl_id_list_n_id(names);
    
    set = isl_set_insert_dims(set, isl_dim_param, 0, n_len);
    
    for (int i = 0; i < n_len; ++i)
    {
        isl_id* id = isl_id_list_get_id(names, i);
        
        set = isl_set_set_dim_id(set, isl_dim_param, i, id);
        
        set = isl_set_equate(set, isl_dim_param, i, isl_dim_set, i);
    }
    
    return set;
}

__isl_give isl_map* tc_parameterize_map_all_in(__isl_take isl_map* map, __isl_keep isl_id_list* names)
{
    int n_len = isl_id_list_n_id(names);
    
    map = isl_map_insert_dims(map, isl_dim_param, 0, n_len);
    
    for (int i = 0; i < n_len; ++i)
    {
        isl_id* id = isl_id_list_get_id(names, i);
        
        map = isl_map_set_dim_id(map, isl_dim_param, i, id);
        
        map = isl_map_equate(map, isl_dim_param, i, isl_dim_in, i);
    }
    
    return map;
}

__isl_give isl_map* tc_parameterize_map_all_out(__isl_take isl_map* map, __isl_keep isl_id_list* names)
{
    int n_len = isl_id_list_n_id(names);
    
    map = isl_map_insert_dims(map, isl_dim_param, 0, n_len);
    
    for (int i = 0; i < n_len; ++i)
    {
        isl_id* id = isl_id_list_get_id(names, i);
        
        map = isl_map_set_dim_id(map, isl_dim_param, i, id);
        
        map = isl_map_equate(map, isl_dim_param, i, isl_dim_out, i);
    }
    
    return map;
}

__isl_give isl_id_list* tc_get_set_params_names(__isl_keep isl_set* set)
{
    int n_params_dim = isl_set_dim(set, isl_dim_param);
    
    isl_id_list* list = isl_id_list_alloc(isl_set_get_ctx(set), n_params_dim);
    
    for (int i = 0; i < n_params_dim; ++i)
    {
        list = isl_id_list_add(list, isl_set_get_dim_id(set, isl_dim_param, i));
    }
    
    return list;
}

__isl_give isl_id_list* tc_get_union_set_params_names(__isl_keep isl_union_set* uset)
{
    isl_set* params = isl_union_set_params(isl_union_set_copy(uset));
    
    isl_id_list* list = tc_get_set_params_names(params);
    
    isl_set_free(params);
    
    return list;
}

__isl_give isl_id_list* tc_get_map_params_names(__isl_keep isl_map* map)
{
    int n_params_dim = isl_map_dim(map, isl_dim_param);
    
    isl_id_list* list = isl_id_list_alloc(isl_map_get_ctx(map), n_params_dim);
    
    for (int i = 0; i < n_params_dim; ++i)
    {
        list = isl_id_list_add(list, isl_map_get_dim_id(map, isl_dim_param, i));
    }
    
    return list;
}

__isl_give isl_id_list* tc_get_union_map_params_names(__isl_keep isl_union_map* umap)
{
    isl_set* params = isl_union_map_params(isl_union_map_copy(umap));
    
    isl_id_list* list = tc_get_set_params_names(params);
    
    isl_set_free(params);
    
    return list;
}

__isl_give isl_set* tc_make_set_constraints(__isl_take isl_set* set, __isl_keep isl_id_list* names)
{
    set = tc_parameterize_all(set, names);
    
    isl_set* params = isl_set_params(set);
    
    return params;
}

__isl_give isl_set* tc_make_set_constraints_for(__isl_take isl_set* set, int pos, __isl_take isl_id* name)
{
    set = tc_parameterize(set, pos, name);
    
    isl_set* params = isl_set_params(set);
    
    return params;
}

__isl_give isl_set* tc_make_map_constraints(__isl_take isl_map* map, __isl_keep isl_id_list* in, __isl_keep isl_id_list* out)
{
    map = tc_lift_down_all_map_vars(map, in, out);
    
    isl_set* params = isl_map_params(map);
    
    return params;
}

__isl_give isl_set* tc_set_fix_param_value(__isl_take isl_set* set, __isl_take isl_id* name, int value)
{
    int pos = isl_set_find_dim_by_id(set, isl_dim_param, name);
    
    set = isl_set_fix_si(set, isl_dim_param, pos, value);
    
    isl_id_free(name);
    
    return set;
}

__isl_give isl_map* tc_map_fix_param_value(__isl_take isl_map* map, __isl_take isl_id* name, int value)
{    
    int pos = isl_map_find_dim_by_id(map, isl_dim_param, name);
    
    map = isl_map_fix_si(map, isl_dim_param, pos, value);
    
    isl_id_free(name);
    
    return map;
}

__isl_give isl_set* tc_set_fix_params_bounds(__isl_take isl_set* set, __isl_take isl_set* bounds)
{
    if (NULL == bounds)
    {
        return set;
    }
    
    isl_id_list* params = tc_get_set_params_names(bounds);
    
    set = isl_set_intersect_params(set, bounds);
    
    set = tc_project_out_params(set, params);
    
    isl_id_list_free(params);
    
    return set;
}

__isl_give isl_map* tc_map_fix_params_bounds(__isl_take isl_map* map, __isl_take isl_set* bounds)
{
    if (NULL == bounds)
    {
        return map;
    }
    
    isl_id_list* params = tc_get_set_params_names(bounds);
    
    map = isl_map_intersect_params(map, bounds);
    
    map = tc_map_project_out_params(map, params);
    
    isl_id_list_free(params);
    
    return map;
}

__isl_give isl_union_set* tc_union_set_fix_params_bounds(__isl_take isl_union_set* uset, __isl_take isl_set* bounds)
{
    if (NULL == bounds)
    {
        return uset;
    }
    
    isl_id_list* params = tc_get_set_params_names(bounds);
    
    uset = isl_union_set_intersect_params(uset, bounds);
    
    uset = tc_union_set_project_out_params(uset, params);
    
    isl_id_list_free(params);
    
    return uset;
}

__isl_give isl_union_map* tc_union_map_fix_params_bounds(__isl_take isl_union_map* umap, __isl_take isl_set* bounds)
{
    if (NULL == bounds)
    {
        return umap;
    }
    
    isl_id_list* params = tc_get_set_params_names(bounds);
    
    umap = isl_union_map_intersect_params(umap, bounds);
    
    umap = tc_union_map_project_out_params(umap, params);
    
    isl_id_list_free(params);
    
    return umap;
}

__isl_give isl_map* tc_make_identity(__isl_take isl_map* map)
{
    isl_map* result = isl_map_identity(isl_map_get_space(map));
    
    isl_map_free(map);
    
    return result;
}

__isl_give isl_union_map* tc_union_map_make_identity(__isl_take isl_union_map* umap)
{
    isl_map_list* maps = tc_collect_maps(umap);
    
    isl_union_map* result = NULL;
    
    for (int i = 0; i < isl_map_list_n_map(maps); ++i)
    {
        isl_map* map = isl_map_list_get_map(maps, i);
        map = tc_make_identity(map);
        
        if (NULL == result)
        {
            result = isl_union_map_from_map(map);
        }
        else
        {
            result = isl_union_map_add_map(result, map);
        }
    }
    
    isl_map_list_free(maps);
    
    return result;    
}

isl_stat tc_scan_set_callback(__isl_take isl_point* point, void* user)
{
    isl_point_dump(point);
    
    isl_point_free(point);
    
    return isl_stat_ok;
}

void tc_scan_set(__isl_keep isl_set* set)
{
    isl_set_foreach_point(set, &tc_scan_set_callback, NULL);
}

void tc_scan_map(__isl_keep isl_map* map)
{
    isl_set* set = isl_map_wrap(isl_map_copy(map));
    
    tc_scan_set(set);
    
    isl_set_free(set);
}

void tc_scan_union_set(__isl_keep isl_union_set* uset)
{
    isl_union_set_foreach_point(uset, &tc_scan_set_callback, NULL);
}

void tc_scan_union_map(__isl_keep isl_union_map* umap)
{
    isl_union_set* uset = isl_union_map_wrap(isl_union_map_copy(umap));
    
    tc_scan_union_set(uset);
    
    isl_union_set_free(uset);
}

__isl_give isl_map* tc_map_closure(__isl_take isl_map* map, int k)
{
    isl_map* k_closure = isl_map_copy(map); 
    
    isl_map* k_power = isl_map_copy(map);
    
    for (int i = 2; i <= k; ++i)
    {
        k_power = isl_map_apply_range(k_power, isl_map_copy(map));
        
        k_closure = isl_map_union(k_closure, isl_map_copy(k_power));
    }
    
    isl_map_free(k_power);
    isl_map_free(map);
    
    return k_closure;
}

__isl_give isl_union_map* tc_union_map_closure(__isl_take isl_union_map* umap, int k)
{
    isl_union_map* k_closure = isl_union_map_copy(umap); 
    
    isl_union_map* k_power = isl_union_map_copy(umap);
    
    for (int i = 2; i <= k; ++i)
    {
        k_power = isl_union_map_apply_range(k_power, isl_union_map_copy(umap));
        
        k_power = isl_union_map_coalesce(k_power);
        
        k_closure = isl_union_map_union(k_closure, isl_union_map_copy(k_power));
        
        k_closure = isl_union_map_coalesce(k_closure);
    }
    
    isl_union_map_free(k_power);
    isl_union_map_free(umap);
    
    return k_closure;
}

__isl_give isl_id_list* tc_ids_sequence(__isl_keep isl_ctx* ctx, const char* id, int size)
{
    char buff[2048];
    
    isl_id_list* result = isl_id_list_alloc(ctx, size);
    
    for (int i = 0; i < size; ++i)
    {
        snprintf(buff, sizeof(buff), "%s%d", id, i);
        
        result = isl_id_list_add(result, isl_id_alloc(ctx, buff, NULL));
    }
    
    return result;
}

__isl_give isl_id_list* tc_ids_sequence_offset(__isl_keep isl_ctx* ctx, const char* id, int size, int offset)
{
    char buff[2048];
    
    isl_id_list* result = isl_id_list_alloc(ctx, size);
    
    for (int i = 0; i < size; ++i)
    {
        snprintf(buff, sizeof(buff), "%s%d", id, i + offset);
        
        result = isl_id_list_add(result, isl_id_alloc(ctx, buff, NULL));
    }
    
    return result;
}

__isl_give isl_id_list* tc_ids_add_suffix(__isl_keep isl_id_list* list, const char* suffix)
{
    char buff[2048];
    
    isl_ctx* ctx = isl_id_list_get_ctx(list);
    
    int n_len = isl_id_list_n_id(list);
    
    isl_id_list* result = isl_id_list_alloc(ctx, n_len);
    
    for (int i = 0; i < n_len; ++i)
    {
        isl_id* id = isl_id_list_get_id(list, i);
        
        snprintf(buff, sizeof(buff), "%s%s", isl_id_get_name(id), suffix);
        
        result = isl_id_list_add(result, isl_id_alloc(ctx, buff, NULL));
        
        isl_id_free(id);
    }    
    
    return result;
}

__isl_give isl_id_list* tc_ids_double(__isl_keep isl_id_list* list)
{
    char buff[2048];
    
    isl_ctx* ctx = isl_id_list_get_ctx(list);
    
    int n_len = isl_id_list_n_id(list);
    
    isl_id_list* result = isl_id_list_alloc(ctx, n_len);
    
    for (int i = 0; i < n_len; ++i)
    {
        isl_id* id = isl_id_list_get_id(list, i);
        
        snprintf(buff, sizeof(buff), "%s%s", isl_id_get_name(id), isl_id_get_name(id));
        
        result = isl_id_list_add(result, isl_id_alloc(ctx, buff, NULL));
        
        isl_id_free(id);
    }
    
    return result;
}

__isl_give isl_id_list* tc_ids_prim(__isl_keep isl_id_list* list)
{
    return tc_ids_add_suffix(list, "_prim");
}

__isl_give isl_id_list* tc_ids_bis(__isl_keep isl_id_list* list)
{
    return tc_ids_add_suffix(list, "_bis");
}

__isl_give isl_id_list* tc_ids_ter(__isl_keep isl_id_list* list)
{
    return tc_ids_add_suffix(list, "_ter");
}

__isl_give isl_id_list* tc_ids_sub(__isl_keep isl_id_list* list, int begin, int end)
{
    isl_id_list* sublist = isl_id_list_alloc(isl_id_list_get_ctx(list), end - begin);

    for (int i = begin; i < end; ++i)
    {
        sublist = isl_id_list_add(sublist, isl_id_list_get_id(list, i));
    }

    return sublist;
}

__isl_give isl_id_list* tc_ids_get(__isl_keep isl_id_list* list, int pos)
{
    return tc_ids_sub(list, pos, pos + 1);
}

__isl_give isl_id_list* tc_ids_single(__isl_keep isl_ctx* ctx, const char* id)
{
    isl_id_list* result = isl_id_list_alloc(ctx, 1);

    result = isl_id_list_add(result, isl_id_alloc(ctx, id, NULL));

    return result;
}

std::string tc_conjunction(const std::string& lhs, const std::string& rhs)
{
    return "(" + lhs + ") and (" + rhs + ")";
}

std::string tc_disjunction(const std::string& lhs, const std::string& rhs)
{
    return "(" + lhs + ") or (" + rhs + ")";
}

std::string tc_negation(const std::string& arg)
{
    return "not (" + arg + ")";
}

std::string tc_tuples_eq(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs)
{
    std::string str;
    
    for (int i = 0; i < isl_id_list_n_id(lhs); ++i)
    {
        isl_id* lhs_id = isl_id_list_get_id(lhs, i);
        isl_id* rhs_id = isl_id_list_get_id(rhs, i);
        
        std::string eq = std::string(isl_id_get_name(lhs_id)) + " = " + std::string(isl_id_get_name(rhs_id));
        
        if (str.size() > 0)
        {
            str = tc_conjunction(str, eq);
        }
        else
        {
            str = eq;
        }
        
        isl_id_free(lhs_id);
        isl_id_free(rhs_id);
    }
        
    return str;
}

std::string tc_tuples_neq(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs)
{
    return tc_negation(tc_tuples_eq(lhs, rhs));
}

std::string tc_tuples_lt(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs)
{
    std::string str;    
    
    for (int i = 0; i < isl_id_list_n_id(lhs); ++i)
    {
        isl_id* lhs_id_i = isl_id_list_get_id(lhs, i);
        isl_id* rhs_id_i = isl_id_list_get_id(rhs, i);
        
        str += "(";
        
        for (int j = 0; j < i; ++j)
        {
            isl_id* lhs_id_j = isl_id_list_get_id(lhs, j);
            isl_id* rhs_id_j = isl_id_list_get_id(rhs, j);
        
            str += std::string(isl_id_get_name(lhs_id_j)) + " = " + std::string(isl_id_get_name(rhs_id_j)) + " and ";
            
            isl_id_free(lhs_id_j);
            isl_id_free(rhs_id_j);
        }
        
        str += std::string(isl_id_get_name(lhs_id_i)) + " < " + std::string(isl_id_get_name(rhs_id_i)) + ") or ";
        
        isl_id_free(lhs_id_i);
        isl_id_free(rhs_id_i);
    }    
    
    str += "false";
    
    return str;
}

std::string tc_tuples_gt(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs)
{
    std::string str; 
    
    for (int i = 0; i < isl_id_list_n_id(lhs); ++i)
    {
        isl_id* lhs_id_i = isl_id_list_get_id(lhs, i);
        isl_id* rhs_id_i = isl_id_list_get_id(rhs, i);
        
        str += "(";
        
        for (int j = 0; j < i; ++j)
        {
            isl_id* lhs_id_j = isl_id_list_get_id(lhs, j);
            isl_id* rhs_id_j = isl_id_list_get_id(rhs, j);
        
            str += std::string(isl_id_get_name(lhs_id_j)) + " = " + std::string(isl_id_get_name(rhs_id_j)) + " and ";
            
            isl_id_free(lhs_id_j);
            isl_id_free(rhs_id_j);
        }
        
        str += std::string(isl_id_get_name(lhs_id_i)) + " > " + std::string(isl_id_get_name(rhs_id_i)) + ") or ";
        
        isl_id_free(lhs_id_i);
        isl_id_free(rhs_id_i);
    }    
    
    str += "false";
    
    return str;
}

std::string tc_tuples_tlt(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs)
{
    isl_id_list* ids_0_lhs = isl_id_list_from_id(isl_id_list_get_id(lhs, 0));
    isl_id_list* ids_0_rhs = isl_id_list_from_id(isl_id_list_get_id(rhs, 0));
    
    std::string str = tc_conjunction(tc_tuples_lt(lhs, rhs), tc_tuples_eq(ids_0_lhs, ids_0_rhs));
    
    isl_id_list_free(ids_0_lhs);
    isl_id_list_free(ids_0_rhs);
    
    return str;
}

std::string tc_tuples_tgt(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs)
{
    isl_id_list* ids_0_lhs = isl_id_list_from_id(isl_id_list_get_id(lhs, 0));
    isl_id_list* ids_0_rhs = isl_id_list_from_id(isl_id_list_get_id(rhs, 0));
    
    std::string str = tc_conjunction(tc_tuples_gt(lhs, rhs), tc_tuples_eq(ids_0_lhs, ids_0_rhs));
    
    isl_id_list_free(ids_0_lhs);
    isl_id_list_free(ids_0_rhs);
    
    return str;
}

std::string tc_tuples_op(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs, const char* op)
{
    std::string str;    
    
    for (int i = 0; i < isl_id_list_n_id(lhs); ++i)
    {
        isl_id* lhs_id_i = isl_id_list_get_id(lhs, i);
        isl_id* rhs_id_i = isl_id_list_get_id(rhs, i);
                
        str += std::string(isl_id_get_name(lhs_id_i)) + " " + std::string(op) + " " + std::string(isl_id_get_name(rhs_id_i)) + " and ";
        
        isl_id_free(lhs_id_i);
        isl_id_free(rhs_id_i);
    }    
    
    str += "true";
    
    return str;
}

__isl_give isl_set* tc_lexmin_set_pos(__isl_take isl_set* set, int pos)
{
    set = tc_project_out_dims_except_pos(set, isl_dim_set, pos);
    
    set = isl_set_lexmin(set);
    
    return set;
}

__isl_give isl_set* tc_lexmax_set_pos(__isl_take isl_set* set, int pos)
{
    set = tc_project_out_dims_except_pos(set, isl_dim_set, pos);
        
    set = isl_set_lexmax(set);
    
    return set;
}

__isl_give isl_set* tc_lexmin_set_var(__isl_take isl_set* set, __isl_take isl_id* id)
{
    int pos = isl_set_find_dim_by_id(set, isl_dim_set, id);
    
    isl_id_free(id);
    
    return tc_lexmin_set_pos(set, pos);
}

__isl_give isl_set* tc_lexmax_set_var(__isl_take isl_set* set, __isl_take isl_id* id)
{
    int pos = isl_set_find_dim_by_id(set, isl_dim_set, id);
    
    isl_id_free(id);
    
    return tc_lexmax_set_pos(set, pos);
}

long tc_lexmin_set_pos_value(__isl_keep isl_set* set, int pos)
{
    isl_set* lexmin = tc_lexmin_set_pos(isl_set_copy(set), pos);
    
    isl_val* value = isl_set_plain_get_val_if_fixed(lexmin, isl_dim_set, 0);
    
    long result = isl_val_get_num_si(value);
    
    isl_val_free(value);    
    isl_set_free(lexmin);
    
    return result;
}

long tc_lexmax_set_pos_value(__isl_keep isl_set* set, int pos)
{
    isl_set* lexmax = tc_lexmax_set_pos(isl_set_copy(set), pos);
    
    isl_val* value = isl_set_plain_get_val_if_fixed(lexmax, isl_dim_set, 0);
    
    long result = isl_val_get_num_si(value);
    
    isl_val_free(value);
    isl_set_free(lexmax);
    
    return result;
}

long tc_lexmin_param_value(__isl_keep isl_set* set, __isl_take isl_id* param)
{
    isl_id_list* params = isl_id_list_from_id(param);
    
    isl_set* lexmin = tc_lift_up_set_params(isl_set_copy(set), params);
    
    long result = tc_lexmin_set_pos_value(lexmin, 0);
    
    isl_id_list_free(params);    
    isl_set_free(lexmin);
    
    return result;
}

long tc_lexmax_param_value(__isl_keep isl_set* set, __isl_take isl_id* param)
{
    isl_id_list* params = isl_id_list_from_id(param);
    
    isl_set* lexmax = tc_lift_up_set_params(isl_set_copy(set), params);
    
    long result = tc_lexmax_set_pos_value(lexmax, 0);
    
    isl_id_list_free(params);
    isl_set_free(lexmax);
    
    return result;
}
 
__isl_give isl_set* tc_get_set_bounds(__isl_keep isl_set* set, __isl_keep isl_id_list* LB, __isl_keep isl_id_list* UB)
{    
    isl_ctx* ctx = isl_set_get_ctx(set);
    
    isl_id_list* LB_UB = isl_id_list_concat(isl_id_list_copy(LB), isl_id_list_copy(UB));
    
    isl_set* bounds = tc_make_params(ctx, LB_UB, "");
    
    isl_id_list_free(LB_UB);
    
    for (int i = 0; i < isl_set_dim(set, isl_dim_set); ++i)
    {
        isl_set* lexmin_set = tc_lexmin_set_pos(isl_set_copy(set), i);
        
        isl_id_list* LB_i = isl_id_list_from_id(isl_id_list_get_id(LB, i));
        
        lexmin_set = tc_make_set_constraints(lexmin_set, LB_i);
        
        isl_id_list_free(LB_i);
        
        bounds = isl_set_intersect(bounds, lexmin_set);
        
        isl_set* lexmax_set = tc_lexmax_set_pos(isl_set_copy(set), i);
        
        isl_id_list* UB_i = isl_id_list_from_id(isl_id_list_get_id(UB, i));
        
        lexmax_set = tc_make_set_constraints(lexmax_set, UB_i);
        
        isl_id_list_free(UB_i);
        
        bounds = isl_set_intersect(bounds, lexmax_set);
    }
    
    return bounds;
}

isl_bool tc_union_map_has_input_tuple(__isl_keep isl_union_map* umap, const char* name)
{
    isl_bool result = isl_bool_false;
    
    isl_map_list* maps = tc_collect_maps(umap);
    
    for (int i = 0; i < isl_map_list_n_map(maps) && result == isl_bool_false; ++i)
    {
        isl_map* map = isl_map_list_get_map(maps, i);
        
        if (0 == strcmp(name, isl_map_get_tuple_name(map, isl_dim_in)))
        {
            result = isl_bool_true;
        }
        
        isl_map_free(map);
    }
    
    isl_map_list_free(maps);
    
    return result;
}

isl_bool tc_union_map_has_output_tuple(__isl_keep isl_union_map* umap, const char* name)
{
    isl_bool result = isl_bool_false;
    
    isl_map_list* maps = tc_collect_maps(umap);
    
    for (int i = 0; i < isl_map_list_n_map(maps) && result == isl_bool_false; ++i)
    {
        isl_map* map = isl_map_list_get_map(maps, i);
        
        if (0 == strcmp(name, isl_map_get_tuple_name(map, isl_dim_out)))
        {
            result = isl_bool_true;
        }
        
        isl_map_free(map);
    }
    
    isl_map_list_free(maps);
    
    return result;
}

__isl_give isl_map* tc_get_map_for_input_tuple(__isl_keep isl_union_map* umap, const char* name)
{
    isl_map* result = NULL;
    
    isl_map_list* maps = tc_collect_maps(umap);
    
    for (int i = 0; i < isl_map_list_n_map(maps); ++i)
    {
        isl_map* map = isl_map_list_get_map(maps, i);
        
        if (0 == strcmp(name, isl_map_get_tuple_name(map, isl_dim_in)))
        {
            if (NULL == result)
            {
                result = map;
            }
            else
            {
                result = isl_map_union(result, map);
            }
        }
        else
        {
            isl_map_free(map);
        }
    }
    
    isl_map_list_free(maps);
    
    return result;
}

__isl_give isl_union_map* tc_get_union_map_for_input_tuple(__isl_keep isl_union_map* umap, const char* name)
{
    isl_union_map* result = NULL;
    
    isl_map_list* maps = tc_collect_maps(umap);
    
    for (int i = 0; i < isl_map_list_n_map(maps); ++i)
    {
        isl_map* map = isl_map_list_get_map(maps, i);
        
        if (0 == strcmp(name, isl_map_get_tuple_name(map, isl_dim_in)))
        {
            if (NULL == result)
            {
                result = isl_union_map_from_map(map);
            }
            else
            {
                result = isl_union_map_add_map(result, map);
            }
        }
        else
        {
            isl_map_free(map);
        }
    }
    
    isl_map_list_free(maps);
    
    return result;
}

__isl_give isl_union_map* tc_remove_map_with_tuple(__isl_take isl_union_map* umap, const char* name)
{
    isl_union_map* result = NULL;
    
    isl_map_list* maps = tc_collect_maps(umap);
    
    isl_union_map_free(umap);
    
    for (int i = 0; i < isl_map_list_n_map(maps); ++i)
    {
        isl_map* map = isl_map_list_get_map(maps, i);
        
        if ((isl_map_has_tuple_name(map, isl_dim_in) && 0 == strcmp(name, isl_map_get_tuple_name(map, isl_dim_in)))
            || (isl_map_has_tuple_name(map, isl_dim_out) && 0 == strcmp(name, isl_map_get_tuple_name(map, isl_dim_out))))
        {
            isl_map_free(map);
        }
        else
        {
            if (NULL == result)
            {
                result = isl_union_map_from_map(map);
            }
            else
            {
                result = isl_union_map_add_map(result, map);
            }
        }
    }
    
    isl_map_list_free(maps);
    
    return result;
}

int tc_get_statement_depth(const char* label, __isl_keep isl_union_map* umap)
{
    isl_map* map = tc_get_map_for_input_tuple(umap, label);
    
    int depth = isl_map_dim(map, isl_dim_in);
    
    isl_map_free(map);
    
    return depth;
}

__isl_give isl_set* tc_normalize_union_set(__isl_keep isl_union_set* uset, __isl_keep isl_union_map* S)
{
    if (isl_union_set_is_empty(uset))
    {
        isl_set* range_set = isl_set_from_union_set(isl_union_map_range(isl_union_map_copy(S)));
                
        isl_set* uset_normalized = isl_set_empty(isl_set_get_space(range_set));
        
        isl_set_free(range_set);
        
        return uset_normalized;
    }
    else
    {
        isl_union_set* uset_normalized = isl_union_set_apply(isl_union_set_copy(uset), isl_union_map_copy(S));
                
        return isl_set_from_union_set(uset_normalized);
    }
}

__isl_give isl_map* tc_normalize_union_map(__isl_keep isl_union_map* umap, __isl_keep isl_union_map* S)
{
    if (isl_union_map_is_empty(umap))
    {
        isl_set* range_set = isl_set_from_union_set(isl_union_map_range(isl_union_map_copy(S)));
                
        isl_map* range_map = isl_map_from_domain_and_range(isl_set_copy(range_set), isl_set_copy(range_set));
        
        isl_set_free(range_set);
                
        isl_map* umap_normalized = isl_map_empty(isl_map_get_space(range_map));
        
        isl_map_free(range_map);
        
        return umap_normalized;
    }
    else
    {    
        isl_union_map* umap_normalized = isl_union_map_apply_domain(isl_union_map_copy(umap), isl_union_map_copy(S));
    
        umap_normalized = isl_union_map_apply_range(umap_normalized, isl_union_map_copy(S));
        
        return isl_map_from_union_map(umap_normalized);
    }
}

__isl_give isl_set* tc_normalize_set(__isl_keep isl_set* set, __isl_keep isl_union_map* S)
{
    isl_union_set* uset = isl_union_set_from_set(isl_set_copy(set));
    
    isl_set* set_normalized = tc_normalize_union_set(uset, S);
    
    isl_union_set_free(uset);
    
    return set_normalized;
}

__isl_give isl_map* tc_normalize_map(__isl_keep isl_map* map, __isl_keep isl_union_map* S)
{
    isl_union_map* umap = isl_union_map_from_map(isl_map_copy(map));
    
    isl_map* map_normalized = tc_normalize_union_map(umap, S);
    
    isl_union_map_free(umap);
    
    return map_normalized;
}

__isl_give isl_union_set* tc_denormalize_set(__isl_keep isl_set* set, __isl_keep isl_union_map* S)
{
    isl_union_set* uset_denormalized = isl_union_set_apply(isl_union_set_from_set(isl_set_copy(set)), isl_union_map_reverse(isl_union_map_copy(S)));

    return uset_denormalized;
}

__isl_give isl_union_map* tc_denormalize_map(__isl_keep isl_map* map, __isl_keep isl_union_map* S)
{
    isl_union_map* umap_denormalized = isl_union_map_apply_domain(isl_union_map_from_map(isl_map_copy(map)), isl_union_map_reverse(isl_union_map_copy(S)));
    
    umap_denormalized = isl_union_map_apply_range(umap_denormalized, isl_union_map_reverse(isl_union_map_copy(S)));

    return umap_denormalized;
}

static isl_stat tc_collect_basic_sets_callback(__isl_take isl_basic_set* bset, void* user)
{
    isl_basic_set_list** bsets = (isl_basic_set_list**)user;
    
    *bsets = isl_basic_set_list_add(*bsets, bset);
    
    return isl_stat_ok;
}

__isl_give isl_basic_set_list* tc_collect_basic_sets(__isl_keep isl_set* set)
{
    isl_basic_set_list* bsets = isl_basic_set_list_alloc(isl_set_get_ctx(set), isl_set_n_basic_set(set));
    
    isl_set_foreach_basic_set(set, &tc_collect_basic_sets_callback, &bsets);
    
    return bsets;
}

static isl_stat tc_collect_sets_callback(__isl_take isl_set* set, void* user)
{
    isl_set_list** sets = (isl_set_list**)user;
    
    *sets = isl_set_list_add(*sets, set);
    
    return isl_stat_ok;
}

__isl_give isl_set_list* tc_collect_sets(__isl_keep isl_union_set* uset)
{
    isl_set_list* sets = isl_set_list_alloc(isl_union_set_get_ctx(uset), isl_union_set_n_set(uset));
    
    isl_union_set_foreach_set(uset, &tc_collect_sets_callback, &sets);
    
    return sets;
}

static isl_stat tc_collect_basic_maps_callback(__isl_take isl_basic_map* bmap, void* user)
{
    isl_basic_map_list** bmaps = (isl_basic_map_list**)user;
    
    *bmaps = isl_basic_map_list_add(*bmaps, bmap);
    
    return isl_stat_ok;
}

__isl_give isl_basic_map_list* tc_collect_maps(__isl_keep isl_map* map)
{
    isl_basic_map_list* bmaps = isl_basic_map_list_alloc(isl_map_get_ctx(map), isl_map_n_basic_map(map));
    
    isl_map_foreach_basic_map(map, &tc_collect_basic_maps_callback, &bmaps);
    
    return bmaps;
}

static isl_stat tc_collect_maps_callback(__isl_take isl_map* map, void* user)
{
    isl_map_list** maps = (isl_map_list**)user;
    
    *maps = isl_map_list_add(*maps, map);
    
    return isl_stat_ok;
}

__isl_give isl_map_list* tc_collect_maps(__isl_keep isl_union_map* umap)
{
    isl_map_list* maps = isl_map_list_alloc(isl_union_map_get_ctx(umap), isl_union_map_n_map(umap));
    
    isl_union_map_foreach_map(umap, &tc_collect_maps_callback, &maps);
    
    return maps;
}

__isl_give isl_union_map* tc_unwrap_range(__isl_take isl_union_map* umap)
{ 
    isl_map_list* maps = tc_collect_maps(umap);
        
    isl_union_map* unwrapped_umap = isl_union_map_empty(isl_union_map_get_space(umap));
    
    for (int i = 0; i < isl_map_list_n_map(maps); ++i)
    {
        isl_map* map = isl_map_list_get_map(maps, i);
        
        if (isl_map_range_is_wrapping(map))
        {            
            map = isl_set_unwrap(isl_map_domain(isl_map_uncurry(map)));
        }
        
        unwrapped_umap = isl_union_map_add_map(unwrapped_umap, map);
    }
    
    isl_map_list_free(maps);
    isl_union_map_free(umap);
    
    return unwrapped_umap;
}

static isl_stat tc_collect_qpolynomials_callback(__isl_take isl_set* domain, __isl_take isl_qpolynomial* poly, void* user)
{
    struct tc_qpolynomials* polys = (tc_qpolynomials*)user;
    
    if (NULL == polys->domains)
    {
        isl_ctx* ctx = isl_set_get_ctx(domain);
    
        polys->polys = (isl_qpolynomial**)malloc(sizeof(isl_qpolynomial*));
        polys->domains = isl_set_list_alloc(ctx, 1);
        
        polys->polys[0] = poly;
        polys->domains = isl_set_list_add(polys->domains, domain);
    }
    else
    {
        int size = isl_set_list_n_set(polys->domains);

        polys->polys = (isl_qpolynomial**)realloc(polys->polys, (size + 1) * sizeof(isl_qpolynomial*));
        
        polys->polys[size] = poly;
        polys->domains = isl_set_list_add(polys->domains, domain);
    }
    
    return isl_stat_ok;
}

struct tc_qpolynomials* tc_collect_qpolynomials(__isl_keep isl_pw_qpolynomial* pw)
{    
    struct tc_qpolynomials* polys = (struct tc_qpolynomials*)malloc(sizeof(struct tc_qpolynomials));
    
    polys->polys = NULL;
    polys->domains = NULL;
    
    isl_pw_qpolynomial_foreach_piece(pw, &tc_collect_qpolynomials_callback, polys);    
    
    return polys;
}

void tc_qpolynomials_free(struct tc_qpolynomials* polys)
{
    if (NULL == polys)
        return;
    
    for (int i = 0; i < isl_set_list_n_set(polys->domains); ++i)
    {
        isl_qpolynomial_free(polys->polys[i]);
    }
    free(polys->polys);
    
    isl_set_list_free(polys->domains);
    
    free(polys);
}

static isl_stat tc_set_to_point_callback(__isl_take isl_point* point, void* user)
{
    isl_point** data = (isl_point**)user;
    
    if (*data != NULL)
    {
        isl_point_free(point);
        
        return isl_stat_error;
    }
    
    *data = point;
    
    return isl_stat_ok;
}

__isl_give isl_point* tc_set_to_point(__isl_take isl_set* set)
{
    isl_point* point = NULL;
    
    isl_stat ret = isl_set_foreach_point(set, &tc_set_to_point_callback, &point);
    
    isl_set_free(set);
    
    if (isl_stat_ok != ret)
    {
        return NULL;
    }
    
    return point;
}

isl_bool tc_points_compare(__isl_keep isl_point* lhs, __isl_keep isl_point* rhs)
{
    isl_space* space = isl_point_get_space(lhs);
    int n = isl_space_dim(space, isl_dim_set);
    isl_space_free(space);

    for (int i = 0; i < n; ++i)
    {
        isl_val* lhs_val = isl_point_get_coordinate_val(lhs, isl_dim_set, i);
        isl_val* rhs_val = isl_point_get_coordinate_val(rhs, isl_dim_set, i);
        
        long lhs_value = isl_val_get_num_si(lhs_val);
        long rhs_value = isl_val_get_num_si(rhs_val);
        
        isl_val_free(lhs_val);
        isl_val_free(rhs_val);

        if (lhs_value != rhs_value)
        {
            return lhs_value < rhs_value ? isl_bool_true : isl_bool_false;
        }
    }

    return isl_bool_false;
}

__isl_give isl_map* tc_get_lex_forward(__isl_take isl_map* R)
{
    isl_ctx* ctx = isl_map_get_ctx(R);
    
    int n_in = isl_map_dim(R, isl_dim_in);
    
    isl_id_list* e_in = tc_ids_sequence(ctx, "in", n_in);
    isl_id_list* e_out = tc_ids_sequence(ctx, "out", n_in);
    
    isl_map* lex_forward_map = tc_make_map(ctx, NULL, e_in, e_out, tc_tuples_gt(e_out, e_in).c_str());
        
    isl_map* lex_forward_map_intersection = isl_map_intersect(R, lex_forward_map);
            
    isl_id_list_free(e_in);
    isl_id_list_free(e_out);
    
    return lex_forward_map_intersection;
}

__isl_give isl_map* tc_get_lex_backward(__isl_take isl_map* R)
{
    isl_ctx* ctx = isl_map_get_ctx(R);
    
    int n_in = isl_map_dim(R, isl_dim_in);
    
    isl_id_list* e_in = tc_ids_sequence(ctx, "in", n_in);
    isl_id_list* e_out = tc_ids_sequence(ctx, "out", n_in);
    
    isl_map* lex_backward_map = tc_make_map(ctx, NULL, e_in, e_out, tc_tuples_lt(e_out, e_in).c_str());
        
    isl_map* lex_backward_map_intersection = isl_map_intersect(R, lex_backward_map);
            
    isl_id_list_free(e_in);
    isl_id_list_free(e_out);
    
    return lex_backward_map_intersection;
}

isl_bool tc_is_lex_forward(__isl_keep isl_map* R)
{
    isl_map* R_forward = tc_get_lex_forward(isl_map_copy(R));
    
    isl_bool is_lex_forward = isl_map_is_equal(R, R_forward);
    
    isl_map_free(R_forward);
    
    return is_lex_forward;
}

__isl_give isl_union_map* tc_simplify_schedule(__isl_take isl_union_map* S)
{    
    isl_map_list* maps = tc_collect_maps(S);
    
    isl_union_map_free(S);
                
    isl_map* map_sample = isl_map_list_get_map(maps, 0);
    
    int n_out = isl_map_dim(map_sample, isl_dim_out);
    
    for (int i = n_out - 1; i >= 0; --i)
    {        
        isl_val* val = isl_map_plain_get_val_if_fixed(map_sample, isl_dim_out, i);

        if (!isl_val_is_nan(val))
        {            
            isl_bool all_equal = isl_bool_true;
            
            for (int j = 0; j < isl_map_list_n_map(maps); ++j)
            {
                isl_map* map = isl_map_list_get_map(maps, j);
                
                isl_val* val_other = isl_map_plain_get_val_if_fixed(map, isl_dim_out, i);
                
                if (!isl_val_eq(val, val_other))
                {
                    all_equal = isl_bool_false;
                }
                
                isl_val_free(val_other);                
                isl_map_free(map);
            }
            
            if (all_equal)
            {
                for (int j = 0; j < isl_map_list_n_map(maps); ++j)
                {
                    isl_map* map = isl_map_list_get_map(maps, j);
                    
                    map = isl_map_remove_dims(map, isl_dim_out, i, 1);
                    
                    maps = isl_map_list_set_map(maps, j, map);
                }
            }
        }
        
        isl_val_free(val);
    }
    
    isl_map_free(map_sample);
    
    isl_union_map* S_prim = NULL;
    
    for (int i = 0; i < isl_map_list_n_map(maps); ++i)
    {
        isl_map* map = isl_map_list_get_map(maps, i);
        
        if (NULL == S_prim)
        {
            S_prim = isl_union_map_from_map(map);
        }
        else
        {
            S_prim = isl_union_map_add_map(S_prim, map);
        }
    }
    
    isl_map_list_free(maps);
    
    return S_prim;
}

__isl_give isl_union_map* tc_extend_schedule(__isl_take isl_union_map* S, int n)
{
    isl_union_map* S_ext = NULL;
    
    isl_map_list* S_maps = tc_collect_maps(S);
    
    isl_union_map_free(S);
    
    for (int i = 0; i < isl_map_list_n_map(S_maps); ++i)
    {
        isl_map* map = isl_map_list_get_map(S_maps, i);
        
        map = isl_map_insert_dims(map, isl_dim_out, 0, n);
    
        if (NULL == S_ext)
        {
            S_ext = isl_union_map_from_map(map);
        }
        else
        {
            S_ext = isl_union_map_add_map(S_ext, map);
        }
    }
    
    isl_map_list_free(S_maps);
    
    return S_ext;
}

__isl_give isl_union_map* tc_extend_union_map(__isl_take isl_union_map* umap, int n)
{
    isl_union_map* umap_ext = NULL;
    
    isl_map_list* maps = tc_collect_maps(umap);
    
    isl_union_map_free(umap);
    
    for (int i = 0; i < isl_map_list_n_map(maps); ++i)
    {
        isl_map* map = isl_map_list_get_map(maps, i);
        
        map = isl_map_insert_dims(map, isl_dim_in, 0, n);
        map = isl_map_insert_dims(map, isl_dim_out, 0, n);
        
        for (int j = 0; j < n; ++j)
        {
            map = isl_map_equate(map, isl_dim_in, j, isl_dim_out, j);
        }
    
        if (NULL == umap_ext)
        {
            umap_ext = isl_union_map_from_map(map);
        }
        else
        {
            umap_ext = isl_union_map_add_map(umap_ext, map);
        }
    }
    
    isl_map_list_free(maps);
    
    return umap_ext;
}

__isl_give isl_set* tc_get_params_set(__isl_take isl_set* set, __isl_keep isl_id_list* params)
{
    isl_set* params_set = isl_set_params(set);
    
    params_set = tc_lift_up_set_params(params_set, params);
    
    return params_set;
}

__isl_give isl_id_list* tc_get_params_ids(__isl_take isl_set* set)
{
    isl_ctx* ctx = isl_set_get_ctx(set);

    int n_param = isl_set_n_param(set);

    isl_id_list* params = isl_id_list_alloc(ctx, n_param);

    for (int i = 0; i < n_param; ++i)
    {
        params = isl_id_list_add(params, isl_set_get_dim_id(set, isl_dim_param, i));
    }

    return params;
}

isl_bool tc_map_carries_dependences(__isl_keep isl_map* map, int pos)
{
    isl_bool carries_dependences = isl_bool_false;
    
    int n_dim = isl_map_dim(map, isl_dim_in);
    
    isl_map* map_eq = isl_map_universe(isl_map_get_space(map));
    
    for (int i = 0; i < n_dim; ++i)
    {
        if (i == pos)
        {
            continue;
        }
        
        map_eq = isl_map_equate(map_eq, isl_dim_in, i, isl_dim_out, i);
    }
    
    isl_map* map_neq = isl_map_universe(isl_map_get_space(map));
    map_neq = isl_map_equate(map_neq, isl_dim_in, pos, isl_dim_out, pos);    
    map_neq = isl_map_complement(map_neq);
                
    isl_map* R = isl_map_intersect(isl_map_copy(map), map_eq);
    R = isl_map_intersect(R, map_neq);
    
    carries_dependences = isl_map_is_empty(R) ? isl_bool_false : isl_bool_true;
    
    isl_map_free(R);
    
    return carries_dependences;
}

__isl_give isl_id_list* tc_id_list_remove_duplicates(__isl_take isl_id_list* list)
{
    for (int i = 0; i < isl_id_list_n_id(list); ++i)
    {
        isl_id* id_i = isl_id_list_get_id(list, i);
        
        for (int j = i + 1; j < isl_id_list_n_id(list); )
        {
            isl_id* id_j = isl_id_list_get_id(list, j);
            
            if (0 == strcmp(isl_id_get_name(id_i), isl_id_get_name(id_j)))
            {
                list = isl_id_list_drop(list, j, 1);
            }
            else
            {
                ++j;
            }
            
            isl_id_free(id_j);
        }
        
        isl_id_free(id_i);
    }
    
    return list;
}

long tc_set_card_value(__isl_take isl_set* set)
{
    isl_pw_qpolynomial* set_card = isl_set_card(set);
        
    isl_val* set_card_val = isl_pw_qpolynomial_max(set_card);

    long card = isl_val_get_num_si(set_card_val);

    isl_val_free(set_card_val);
    
    return card;
}

char* tc_qpolynomial_to_str(__isl_keep isl_qpolynomial* poly)
{
    isl_ctx* ctx = isl_qpolynomial_get_ctx(poly);

    isl_printer* printer = isl_printer_to_str(ctx);

    printer = isl_printer_set_output_format(printer, ISL_FORMAT_C);

    printer = isl_printer_print_qpolynomial(printer, poly);

    char* str = isl_printer_get_str(printer);

    isl_printer_free(printer);
    
    return str;
}

char* tc_qpolynomial_fold_to_str(__isl_keep isl_qpolynomial_fold* fold)
{
    isl_ctx* ctx = isl_qpolynomial_fold_get_ctx(fold);

    isl_printer* printer = isl_printer_to_str(ctx);

    printer = isl_printer_set_output_format(printer, ISL_FORMAT_C);

    printer = isl_printer_print_qpolynomial_fold(printer, fold);

    char* str = isl_printer_get_str(printer);

    isl_printer_free(printer);
    
    return str;
}

__isl_give isl_union_map* tc_loop_interchange(__isl_take isl_union_map* S, const char* name, const char* a, const char* b)
{
    isl_map* statement_schedule = tc_get_map_for_input_tuple(S, name);

    int n_in = isl_map_dim(statement_schedule, isl_dim_in);

    isl_set* statement_schedule_domain = isl_map_domain(isl_map_copy(statement_schedule));

    int pos_a = -1;
    int pos_b = -1;

    for (int i = 0; i < n_in; ++i)
    {
        if (strcmp(isl_set_get_dim_name(statement_schedule_domain, isl_dim_set, i), a) == 0)
        {
            pos_a = i;
        }
        if (strcmp(isl_set_get_dim_name(statement_schedule_domain, isl_dim_set, i), b) == 0)
        {
            pos_b = i;
        }
    }

    if (pos_a != -1 && pos_b != -1)
    {
        isl_map* interchange = isl_map_from_domain_and_range(isl_set_copy(statement_schedule_domain), isl_set_copy(statement_schedule_domain));
        
        for (int i = 0; i < n_in; ++i)
        {
            if (i != pos_a && i != pos_b)
            {
                interchange = isl_map_equate(interchange, isl_dim_in, i, isl_dim_out, i);
            }
        }

        interchange = isl_map_equate(interchange, isl_dim_in, pos_a, isl_dim_out, pos_b);
        interchange = isl_map_equate(interchange, isl_dim_in, pos_b, isl_dim_out, pos_a);
        
        statement_schedule = isl_map_apply_domain(statement_schedule, interchange);

        S = tc_remove_map_with_tuple(S, name);

        S = isl_union_map_add_map(S, isl_map_copy(statement_schedule));
    }

    isl_set_free(statement_schedule_domain);
    isl_map_free(statement_schedule);

    return S;
}
