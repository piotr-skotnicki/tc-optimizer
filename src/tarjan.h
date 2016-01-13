#ifndef TC_TARJAN_H
#define	TC_TARJAN_H

#include <isl/map.h>

__isl_give isl_map* tc_tarjan_transitive_closure(__isl_take isl_map* R);

__isl_give isl_map* tc_tarjan_components(__isl_take isl_map* R);

#endif // TC_TARJAN_H
