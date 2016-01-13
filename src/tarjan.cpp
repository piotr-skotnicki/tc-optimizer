#include "tarjan.h"

#include <isl/ctx.h>
#include <isl/space.h>
#include <isl/point.h>
#include <isl/val.h>
#include <isl/set.h>
#include <isl/map.h>

#include <stddef.h>
#include <limits.h>

#include <map>
#include <stack>

static isl_stat tc_tarjan_transitive_closure_init(__isl_take isl_point* u, void* user);

static isl_stat tc_tarjan_transitive_closure_dfs(__isl_take isl_point* u, void* user);

static isl_stat tc_tarjan_transitive_closure_visit(__isl_take isl_point* v, void* user);

static isl_stat tc_tarjan_components_init(__isl_take isl_point* u, void* user);

static isl_stat tc_tarjan_components_dfs(__isl_take isl_point* u, void* user);

static isl_stat tc_tarjan_components_visit(__isl_take isl_point* v, void* user);

struct points_compare
{
    bool operator()(__isl_keep isl_point* a, __isl_keep isl_point* b) const
    {    
        isl_space* space = isl_point_get_space(a);
        int n = isl_space_dim(space, isl_dim_set);
        isl_space_free(space);
        
        for (int i = 0; i < n; ++i)
        {
            isl_val* va = isl_point_get_coordinate_val(a, isl_dim_set, i);
            isl_val* vb = isl_point_get_coordinate_val(b, isl_dim_set, i);
            long asi = isl_val_get_num_si(va);
            long bsi = isl_val_get_num_si(vb);
            isl_val_free(va);
            isl_val_free(vb);
            
            if (asi != bsi)
            {
                return asi < bsi;
            }
        }
        
        return false;
    }
};

static isl_map* R;

static unsigned int order;

static std::map<isl_point*, unsigned int, points_compare> low, num;

static std::map<isl_point*, bool, points_compare> vis;

static std::map<isl_point*, isl_set*, points_compare> tc;

static std::map<isl_point*, isl_set*, points_compare> scc;

static std::stack<isl_point*> S;

static isl_stat tc_tarjan_transitive_closure_init(__isl_take isl_point* u, void* user)
{    
    low[u] = 0;
    vis[u] = false;
    num[u] = 0;
    tc[u] = NULL;
        
    return isl_stat_ok;
}

static isl_stat tc_tarjan_transitive_closure_visit(__isl_take isl_point* v, void* uu)
{
    isl_point* u = (isl_point*)uu;
    
    tc_tarjan_transitive_closure_dfs(isl_point_copy(v), NULL);
    
    if (low[u] > low[v])
    {
        low[u] = low[v];
    }
    
    tc[u] = isl_set_union(tc[u], isl_set_copy(tc[v]));
    tc[u] = isl_set_coalesce(tc[u]);
    
    isl_point_free(v);
    
    return isl_stat_ok;
}

static isl_stat tc_tarjan_transitive_closure_dfs(__isl_take isl_point* u, void* data)
{
    if (!vis[u])
    {        
        low[u] = num[u] = ++order;
        vis[u] = true;
        S.push(u);

        isl_set* v = isl_set_apply(isl_set_from_point(isl_point_copy(u)), isl_map_copy(R));
        tc[u] = v;
        isl_set_foreach_point(v, &tc_tarjan_transitive_closure_visit, u);

        if (low[u] == num[u])
        {
            isl_point* w = NULL;
            do
            {
                isl_point_free(w);
                w = S.top(); S.pop();
                low[w] = UINT_MAX;
                
                if (w != u)
                {
                    isl_set_free(tc[w]);
                    tc[w] = isl_set_copy(tc[u]);
                }
            }
            while (w != u);            
            isl_point_free(u);
        }
    }
    else
    {
        isl_point_free(u);
    }
    
    return isl_stat_ok;
}

__isl_give isl_map* tc_tarjan_transitive_closure(__isl_take isl_map* R)
{
    ::R = R;
    ::order = 0;
    ::low.clear();
    ::num.clear();
    ::vis.clear();
    ::tc.clear();
    while (!::S.empty()) S.pop();
    
    isl_map* R_plus = isl_map_empty(isl_map_get_space(R));    
    
    isl_set* points = isl_set_union(isl_map_domain(isl_map_copy(R)), isl_map_range(isl_map_copy(R)));
    
    isl_set_foreach_point(points, &tc_tarjan_transitive_closure_init, NULL);
    isl_set_foreach_point(points, &tc_tarjan_transitive_closure_dfs, NULL);
    
    for (std::map<isl_point*, isl_set*, points_compare>::iterator it = tc.begin(); it != tc.end(); ++it)
    {
        isl_point* u = it->first;
        isl_set* v = it->second;
        
        isl_map* relation = isl_map_from_domain_and_range(isl_set_from_point(u), v);
        
        R_plus = isl_map_union(R_plus, relation);
        
        //R_plus = isl_map_coalesce(R_plus);
    }
        
    isl_set_free(points);
    isl_map_free(R);
    
    return R_plus;
}

static isl_stat tc_tarjan_components_init(__isl_take isl_point* u, void* user)
{    
    low[u] = 0;
    vis[u] = false;
    num[u] = 0;
    scc[u] = NULL;
    
    return isl_stat_ok;
}

static isl_stat tc_tarjan_components_visit(__isl_take isl_point* v, void* uu)
{
    isl_point* u = (isl_point*)uu;
    
    tc_tarjan_components_dfs(isl_point_copy(v), NULL);
    
    if (low[u] > low[v])
    {
        low[u] = low[v];
    }
    
    isl_point_free(v);
    
    return isl_stat_ok;
}

static isl_stat tc_tarjan_components_dfs(__isl_take isl_point* u, void* data)
{
    if (!vis[u])
    {        
        low[u] = num[u] = ++order;
        vis[u] = true;
        S.push(u);

        isl_set* v = isl_set_apply(isl_set_from_point(isl_point_copy(u)), isl_map_copy(R));
        isl_set_foreach_point(v, &tc_tarjan_components_visit, u);
        isl_set_free(v);

        if (low[u] == num[u])
        {
            scc[u] = isl_set_empty(isl_point_get_space(u));
            while (true)
            {
                isl_point* w = S.top(); S.pop();
                low[w] = UINT_MAX;
                if (w == u)
                {
                    break;
                }
                scc[u] = isl_set_union(scc[u], isl_set_from_point(w));
            }
            isl_point_free(u);
        }
    }
    else
    {
        isl_point_free(u);
    }
    
    return isl_stat_ok;
}

__isl_give isl_map* tc_tarjan_components(__isl_take isl_map* R)
{
    ::R = R;
    ::order = 0;
    ::low.clear();
    ::num.clear();
    ::vis.clear();
    ::scc.clear();
    while (!::S.empty()) S.pop();
    
    isl_map* R_scc = isl_map_empty(isl_map_get_space(R));    
        
    isl_set* points = isl_set_union(isl_map_domain(isl_map_copy(R)), isl_map_range(isl_map_copy(R)));
    
    isl_set_foreach_point(points, &tc_tarjan_components_init, NULL);
    isl_set_foreach_point(points, &tc_tarjan_components_dfs, NULL);
        
    for (std::map<isl_point*, isl_set*, points_compare>::iterator it = scc.begin(); it != scc.end(); ++it)
    {    
        isl_point* u = it->first;
        isl_set* v = it->second;
        
        isl_map* relation = isl_map_from_domain_and_range(isl_set_from_point(u), v);
        
        if (NULL != relation)
        {         
            R_scc = isl_map_union(R_scc, relation);
        }
        
        //R_scc = isl_map_coalesce(R_scc);
    }
    
    isl_set_free(points);
    isl_map_free(R);
    
    return R_scc;
}
