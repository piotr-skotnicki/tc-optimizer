#ifndef TC_TRANSITIVE_CLOSURE_H
#define TC_TRANSITIVE_CLOSURE_H

#include <isl/map.h>
#include <isl/union_map.h>

extern __isl_give isl_map* (*tc_transitive_closure)(__isl_take isl_map* R, __isl_keep isl_union_map* S, isl_bool* exact);

extern __isl_give isl_map* (*tc_map_power)(__isl_take isl_map* R, __isl_keep isl_union_map* S, isl_bool* exact);

__isl_give isl_map* tc_transitive_closure_adapter_isl_map(__isl_take isl_map* R, __isl_keep isl_union_map* S, isl_bool* exact);

__isl_give isl_map* tc_transitive_closure_adapter_isl_union_map(__isl_take isl_map* R, __isl_keep isl_union_map* S, isl_bool* exact);

__isl_give isl_map* tc_transitive_closure_adapter_floyd_warshall(__isl_take isl_map* R, __isl_keep isl_union_map* S, isl_bool* exact);

__isl_give isl_map* tc_transitive_closure_adapter_iterative(__isl_take isl_map* R, __isl_keep isl_union_map* S, isl_bool* exact);

__isl_give isl_map* tc_transitive_closure_adapter_tarjan(__isl_take isl_map* R, __isl_keep isl_union_map* S, isl_bool* exact);

__isl_give isl_map* tc_map_power_adapter_isl_map(__isl_take isl_map* R, __isl_keep isl_union_map* S, isl_bool* exact);

__isl_give isl_map* tc_map_power_adapter_isl_union_map(__isl_take isl_map* R, __isl_keep isl_union_map* S, isl_bool* exact);

#endif // TC_TRANSITIVE_CLOSURE_H
