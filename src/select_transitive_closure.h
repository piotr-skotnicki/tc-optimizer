#ifndef TC_SELECT_TRANSITIVE_CLOSURE_H
#define TC_SELECT_TRANSITIVE_CLOSURE_H

#include <isl/map.h>
#include <isl/union_map.h>

__isl_give isl_map* tc_select_transitive_closure(__isl_take isl_map* R, __isl_keep isl_union_map* S, isl_bool* exact);

#endif // TC_SELECT_TRANSITIVE_CLOSURE_H
