#ifndef TC_TRANSITIVE_CLOSURE_H
#define	TC_TRANSITIVE_CLOSURE_H

#include <isl/map.h>

__isl_give isl_map* tc_transitive_closure(__isl_take isl_map* R, int max);

#endif // TC_TRANSITIVE_CLOSURE_H
