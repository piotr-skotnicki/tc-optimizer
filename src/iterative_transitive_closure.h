#ifndef TC_ITERATIVE_TRANSITIVE_CLOSURE_H
#define TC_ITERATIVE_TRANSITIVE_CLOSURE_H

#include <isl/map.h>

__isl_give isl_map* tc_iterative_transitive_closure(__isl_take isl_map* map, int max_iterations, int max_disjunctions, isl_bool* exact);

#endif // TC_ITERATIVE_TRANSITIVE_CLOSURE_H
