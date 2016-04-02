#ifndef TC_TRANSITIVE_CLOSURE_H
#define	TC_TRANSITIVE_CLOSURE_H

#include <isl/map.h>
#include <isl/union_map.h>

__isl_give isl_map* tc_transitive_closure(__isl_take isl_map* R, int max);

__isl_give isl_union_map* tc_floyd_warshall_transitive_closure(__isl_take isl_union_map* R);

#endif // TC_TRANSITIVE_CLOSURE_H
