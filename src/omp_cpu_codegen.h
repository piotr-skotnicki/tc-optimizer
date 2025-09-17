#ifndef TC_OMP_CPU_CODEGEN_H
#define TC_OMP_CPU_CODEGEN_H

#include "scop.h"
#include "options.h"

#include <isl/id.h>
#include <isl/union_map.h>

void tc_codegen_omp_parallel_for(struct tc_scop* scop, struct tc_options* options, __isl_take isl_union_map* S, __isl_keep isl_id_list* iterators, __isl_keep isl_id_list* parallel_iterators, int nested);

void tc_codegen_omp_task_for(struct tc_scop* scop, struct tc_options* options, __isl_take isl_union_map* S, __isl_keep isl_id_list* iterators, __isl_keep isl_id_list* parallel_iterators, int nested);

#endif // TC_OMP_CPU_CODEGEN_H
