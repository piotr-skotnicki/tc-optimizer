#ifndef TC_OMP_GPU_CODEGEN_H
#define TC_OMP_GPU_CODEGEN_H

#include "scop.h"
#include "options.h"

#include <isl/id.h>
#include <isl/set.h>
#include <isl/union_map.h>

void tc_codegen_omp_gpu(struct tc_scop* scop, struct tc_options* options, __isl_take isl_union_map* S, __isl_take isl_set* tile, __isl_keep isl_id_list* iterators, __isl_keep isl_id_list* parallel_iterators);

#endif // TC_OMP_GPU_CODEGEN_H
