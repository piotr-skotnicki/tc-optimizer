#include "utility.h"

#include <isl/ctx.h>
#include <isl/id.h>
#include <isl/set.h>
#include <isl/map.h>
#include <isl/union_set.h>
#include <isl/union_map.h>

#include <stdio.h>
#include <stddef.h>
#include <string.h>

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

__isl_give isl_set* tc_project_out_dim_names(__isl_take isl_set* set, isl_dim_type dim, __isl_keep isl_id_list* names)
{    
    for (int i = 0; i < isl_id_list_n_id(names); ++i)
    {
        isl_id* id = isl_id_list_get_id(names, i);
        
        int pos = isl_set_find_dim_by_id(set, dim, id);
        
        set = isl_set_project_out(set, dim, pos, 1);
        
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

__isl_give isl_set* tc_lift_down_all_set_vars(__isl_take isl_set* set, __isl_keep isl_id_list* vars)
{    
    for (int i = isl_id_list_n_id(vars) - 1; i >= 0; --i)
    {
        set = isl_set_set_dim_id(set, isl_dim_set, i, isl_id_list_get_id(vars, i));
        
        set = isl_set_move_dims(set, isl_dim_param, 0, isl_dim_set, i, 1);
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
    char set_str[1024];
    
    if (NULL == constraints || 0 == strlen(constraints))
    {
        sprintf(set_str, "[%s] -> { [%s] }", NULL != params ? tc_comma(params).c_str() : "", tc_comma(vars).c_str());
    }
    else
    {
        sprintf(set_str, "[%s] -> { [%s] : %s }", NULL != params ? tc_comma(params).c_str() : "", tc_comma(vars).c_str(), constraints);    
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
    char map_str[1024];
    
    if (NULL == constraints || 0 == strlen(constraints))
    {
        sprintf(map_str, "[%s] -> { [%s] -> [%s] }", NULL != params ? tc_comma(params).c_str() : "", tc_comma(in).c_str(), tc_comma(out).c_str());
    }
    else
    {
        sprintf(map_str, "[%s] -> { [%s] -> [%s] : %s }", NULL != params ? tc_comma(params).c_str() : "", tc_comma(in).c_str(), tc_comma(out).c_str(), constraints);
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
    char set_str[1024];
    
    sprintf(set_str, "[%s] -> { : %s }", tc_comma(params).c_str(), constraints);    
    
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

__isl_give isl_set* tc_fix_param_value(__isl_take isl_set* set, __isl_take isl_id* name, int value)
{
    int pos = isl_set_find_dim_by_id(set, isl_dim_param, name);
    
    set = isl_set_fix_si(set, isl_dim_param, pos, value);
    
    set = isl_set_project_out(set, isl_dim_param, pos, 1);
    
    isl_id_free(name);
    
    return set;
}

__isl_give isl_map* tc_make_identity(__isl_take isl_map* map)
{
    isl_map* result = isl_map_identity(isl_map_get_space(map));
    
    result = isl_map_intersect_domain(result, isl_map_domain(isl_map_copy(map)));
    
    result = isl_map_intersect_range(result, isl_map_domain(map));
    
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
        
        k_closure = isl_union_map_union(k_closure, isl_union_map_copy(k_power));
    }
    
    isl_union_map_free(k_power);
    isl_union_map_free(umap);
    
    return k_closure;
}

__isl_give isl_id_list* tc_ids_sequence(__isl_keep isl_ctx* ctx, const char* id, int size)
{
    char buff[50];
    
    isl_id_list* result = isl_id_list_alloc(ctx, size);
    
    for (int i = 0; i < size; ++i)
    {
        sprintf(buff, "%s%d", id, i);
        
        result = isl_id_list_add(result, isl_id_alloc(ctx, buff, NULL));
    }
    
    return result;
}

__isl_give isl_id_list* tc_ids_add_suffix(__isl_keep isl_id_list* list, const char* suffix)
{
    char buff[50];
    
    isl_ctx* ctx = isl_id_list_get_ctx(list);
    
    int n_len = isl_id_list_n_id(list);
    
    isl_id_list* result = isl_id_list_alloc(ctx, n_len);
    
    for (int i = 0; i < n_len; ++i)
    {
        isl_id* id = isl_id_list_get_id(list, i);
        
        sprintf(buff, "%s%s", isl_id_get_name(id), suffix);
        
        result = isl_id_list_add(result, isl_id_alloc(ctx, buff, NULL));
        
        isl_id_free(id);
    }    
    
    return result;
}

__isl_give isl_id_list* tc_ids_double(__isl_keep isl_id_list* list)
{
    char buff[50];
    
    isl_ctx* ctx = isl_id_list_get_ctx(list);
    
    int n_len = isl_id_list_n_id(list);
    
    isl_id_list* result = isl_id_list_alloc(ctx, n_len);
    
    for (int i = 0; i < n_len; ++i)
    {
        isl_id* id = isl_id_list_get_id(list, i);
        
        sprintf(buff, "%s%s", isl_id_get_name(id), isl_id_get_name(id));
        
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

struct tc_get_map_for_input_tuple_user
{
    const char* name;
    
    isl_map** result;
};

static isl_stat tc_get_map_for_input_tuple_callback(__isl_take isl_map* map, void* user)
{
    tc_get_map_for_input_tuple_user* data = (tc_get_map_for_input_tuple_user*)user;
    
    if (0 == strcmp(data->name, isl_map_get_tuple_name(map, isl_dim_in)))
    {
        if (NULL == *data->result)
        {
            *data->result = map;
        }
        else
        {
            *data->result = isl_map_union(*data->result, map);
        }
    }
    else
    {
        isl_map_free(map);
    }
    
    return isl_stat_ok;
}

__isl_give isl_map* tc_get_map_for_input_tuple(__isl_keep isl_union_map* umap, const char* name)
{
    isl_map* result = NULL;
    
    tc_get_map_for_input_tuple_user user;
    user.name = name;
    user.result = &result;
    
    isl_union_map_foreach_map(umap, &tc_get_map_for_input_tuple_callback, &user);
    
    return result;
}

int tc_get_statement_depth(const char* label, __isl_keep isl_union_map* umap)
{
    isl_map* map = tc_get_map_for_input_tuple(umap, label);
    
    int depth = isl_map_n_in(map);
    
    isl_map_free(map);
    
    return depth;
}

__isl_give isl_set* tc_normalize_union_set(__isl_keep isl_union_set* uset, __isl_keep isl_union_map* scattering)
{
    isl_union_set* uset_normalized = isl_union_set_apply(isl_union_set_copy(uset), isl_union_map_copy(scattering));
                
    return isl_set_from_union_set(uset_normalized);
}

__isl_give isl_map* tc_normalize_union_map(__isl_keep isl_union_map* umap, __isl_keep isl_union_map* scattering)
{
    isl_union_map* umap_normalized = isl_union_map_apply_domain(isl_union_map_copy(umap), isl_union_map_copy(scattering));
    
    umap_normalized = isl_union_map_apply_range(umap_normalized, isl_union_map_copy(scattering));
            
    return isl_map_from_union_map(umap_normalized);
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
