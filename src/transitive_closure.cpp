#include "transitive_closure.h"
#include "tarjan_transitive_closure.h"
#include "floyd_warshall_transitive_closure.h"
#include "utility.h"
#include "debug.h"

#include <isl/map.h>
#include <isl/union_map.h>

#include <stddef.h>

__isl_give isl_map* (*tc_transitive_closure)(__isl_take isl_map* R, __isl_keep isl_union_map* S, int* exact);

__isl_give isl_map* (*tc_map_power)(__isl_take isl_map* R, __isl_keep isl_union_map* S, int* exact);

__isl_give isl_map* tc_transitive_closure_adapter_isl_map(__isl_take isl_map* R, __isl_keep isl_union_map* S, int* exact)
{
    return isl_map_transitive_closure(R, exact);
}

__isl_give isl_map* tc_transitive_closure_adapter_isl_union_map(__isl_take isl_map* R, __isl_keep isl_union_map* S, int* exact)
{
    isl_union_map* R_denorm = tc_denormalize_map(R, S);
    
    isl_map_free(R);
    
    isl_union_map* R_plus_denorm = isl_union_map_transitive_closure(R_denorm, exact);
    
    isl_map* R_plus = tc_normalize_union_map(R_plus_denorm, S);
    
    isl_union_map_free(R_plus_denorm);
    
    return R_plus;
}

__isl_give isl_map* tc_transitive_closure_adapter_floyd_warshall(__isl_take isl_map* R, __isl_keep isl_union_map* S, int* exact)
{
    if (NULL != exact)
    {
        *exact = -1;
    }
        
    isl_union_map* R_denorm = tc_denormalize_map(R, S);
    
    isl_map_free(R);
    
    isl_union_map* R_plus_denorm = tc_floyd_warshall_transitive_closure(R_denorm);
    
    isl_map* R_plus = tc_normalize_union_map(R_plus_denorm, S);
    
    isl_union_map_free(R_plus_denorm);
    
    return R_plus;
}

__isl_give isl_map* tc_transitive_closure_adapter_tarjan(__isl_take isl_map* R, __isl_keep isl_union_map* S, int* exact)
{
    if (NULL != exact)
    {
        *exact = 1;
    }
        
    return tc_tarjan_transitive_closure(R);
}

__isl_give isl_map* tc_map_power_adapter_isl_map(__isl_take isl_map* R, __isl_keep isl_union_map* S, int* exact)
{
    return isl_map_power(R, exact);
}

__isl_give isl_map* tc_map_power_adapter_isl_union_map(__isl_take isl_map* R, __isl_keep isl_union_map* S, int* exact)
{
    isl_union_map* R_denorm = tc_denormalize_map(R, S);
    
    isl_map_free(R);
    
    isl_union_map* R_k_denorm = isl_union_map_power(R_denorm, exact);
        
    isl_map* R_k = tc_normalize_union_map(R_k_denorm, S);
    
    isl_union_map_free(R_k_denorm);
    
    return R_k;
}
