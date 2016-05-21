#ifndef TC_SERIAL_CODEGEN_H
#define TC_SERIAL_CODEGEN_H

#include "scop.h"
#include "options.h"

#include <isl/id.h>
#include <isl/set.h>
#include <isl/union_map.h>

void tc_codegen_serial(struct tc_scop* scop, struct tc_options* options, __isl_take isl_union_map* S, __isl_take isl_set* tile, __isl_keep isl_id_list* iterators);

#endif // TC_SERIAL_CODEGEN_H
