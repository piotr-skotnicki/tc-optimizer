#ifndef TC_FLOYD_WARSHALL_TRANSITIVE_CLOSURE_H
#define	TC_FLOYD_WARSHALL_TRANSITIVE_CLOSURE_H

#include <isl/union_map.h>

__isl_give isl_union_map* tc_floyd_warshall_transitive_closure(__isl_take isl_union_map* R, int* exact);

#endif // TC_FLOYD_WARSHALL_TRANSITIVE_CLOSURE_H
