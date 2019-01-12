#ifndef TC_SCHEDULING_H
#define TC_SCHEDULING_H

#include "scop.h"
#include "options.h"

#include <isl/set.h>
#include <isl/map.h>
#include <isl/id.h>

extern void (*tc_scheduling)(struct tc_scop* scop, struct tc_options* options, __isl_take isl_union_set* LD, __isl_take isl_union_map* S, __isl_take isl_union_map* R, __isl_take isl_set* ii_set, __isl_take isl_set* tile, __isl_take isl_map* Rtile, __isl_take isl_id_list* II, __isl_take isl_id_list* I);

#endif // TC_SCHEDULING_H
