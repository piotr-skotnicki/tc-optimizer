#include "transitive_closure.h"
#include "tarjan_transitive_closure.h"
#include "floyd_warshall_transitive_closure.h"
#include "iterative_transitive_closure.h"
#include "utility.h"
#include "debug.h"
#include "timer.h"

#include <isl/map.h>
#include <isl/union_map.h>

#include <stddef.h>

__isl_give isl_map* (*tc_transitive_closure)(__isl_take isl_map* R, __isl_keep isl_union_map* S, isl_bool* exact);

__isl_give isl_map* (*tc_map_power)(__isl_take isl_map* R, __isl_keep isl_union_map* S, isl_bool* exact);

__isl_give isl_map* tc_transitive_closure_adapter_isl_map(__isl_take isl_map* R, __isl_keep isl_union_map* S, isl_bool* exact)
{
    if (NULL != exact)
    {
        *exact = isl_bool_false;
    }

    struct tc_timer* timer = tc_timer_start();

    isl_map* R_plus = isl_map_transitive_closure(R, exact);

    tc_debug("R+ calculations time: %ld ms", tc_timer_stop(timer));

    return R_plus;
}

__isl_give isl_map* tc_transitive_closure_adapter_isl_union_map(__isl_take isl_map* R, __isl_keep isl_union_map* S, isl_bool* exact)
{
    if (NULL != exact)
    {
        *exact = isl_bool_false;
    }

    isl_union_map* R_denorm = tc_denormalize_map(R, S);

    isl_map_free(R);

    struct tc_timer* timer = tc_timer_start();

    isl_union_map* R_plus_denorm = isl_union_map_transitive_closure(R_denorm, exact);

    tc_debug("R+ calculations time: %ld ms", tc_timer_stop(timer));

    isl_map* R_plus = tc_normalize_union_map(R_plus_denorm, S);

    isl_union_map_free(R_plus_denorm);

    return R_plus;
}

__isl_give isl_map* tc_transitive_closure_adapter_floyd_warshall(__isl_take isl_map* R, __isl_keep isl_union_map* S, isl_bool* exact)
{
    if (NULL != exact)
    {
        *exact = isl_bool_false;
    }

    isl_union_map* R_denorm = tc_denormalize_map(R, S);

    isl_map_free(R);

    struct tc_timer* timer = tc_timer_start();

    isl_union_map* R_plus_denorm = tc_floyd_warshall_transitive_closure(R_denorm, exact);

    tc_debug("R+ calculations time: %ld ms", tc_timer_stop(timer));

    isl_map* R_plus = tc_normalize_union_map(R_plus_denorm, S);

    isl_union_map_free(R_plus_denorm);

    return R_plus;
}

__isl_give isl_map* tc_transitive_closure_adapter_iterative(__isl_take isl_map* R, __isl_keep isl_union_map* S, isl_bool* exact)
{
    if (NULL != exact)
    {
        *exact = isl_bool_false;
    }

    struct tc_timer* timer = tc_timer_start();

    isl_map* R_plus = tc_iterative_transitive_closure(R, 10, 100, exact);

    tc_debug("R+ calculations time: %ld ms", tc_timer_stop(timer));

    return R_plus;
}

__isl_give isl_map* tc_transitive_closure_adapter_tarjan(__isl_take isl_map* R, __isl_keep isl_union_map* S, isl_bool* exact)
{
    if (NULL != exact)
    {
        *exact = isl_bool_true;
    }

    struct tc_timer* timer = tc_timer_start();

    isl_map* R_plus = tc_tarjan_transitive_closure(R);

    tc_debug("R+ calculations time: %ld ms", tc_timer_stop(timer));

    return R_plus;
}

__isl_give isl_map* tc_map_power_adapter_isl_map(__isl_take isl_map* R, __isl_keep isl_union_map* S, isl_bool* exact)
{
    if (NULL != exact)
    {
        *exact = isl_bool_false;
    }

    struct tc_timer* timer = tc_timer_start();

    isl_map* R_k = isl_map_power(R, exact);

    tc_debug("R^k calculations time: %ld ms", tc_timer_stop(timer));

    return R_k;
}

__isl_give isl_map* tc_map_power_adapter_isl_union_map(__isl_take isl_map* R, __isl_keep isl_union_map* S, isl_bool* exact)
{
    if (NULL != exact)
    {
        *exact = isl_bool_false;
    }

    isl_union_map* R_denorm = tc_denormalize_map(R, S);

    isl_map_free(R);

    struct tc_timer* timer = tc_timer_start();

    isl_union_map* R_k_denorm = isl_union_map_power(R_denorm, exact);

    tc_debug("R^k calculations time: %ld ms", tc_timer_stop(timer));

    isl_map* R_k = tc_normalize_union_map(R_k_denorm, S);

    isl_union_map_free(R_k_denorm);

    return R_k;
}
