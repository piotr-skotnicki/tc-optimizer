#ifndef TC_ISL_SCHEDULING_H
#define TC_ISL_SCHEDULING_H

#include "scheduling.h"
#include "scop.h"
#include "options.h"

#include <isl/map.h>
#include <isl/union_map.h>
#include <isl/set.h>
#include <isl/id.h>
#include <isl/schedule.h>

void tc_scheduling_isl(struct tc_scop* scop, struct tc_options* options, __isl_take isl_union_set* LD, __isl_take isl_union_map* S, __isl_take isl_union_map* R, __isl_take isl_set* ii_set, __isl_take isl_set* tile, __isl_take isl_map* Rtile, __isl_take isl_id_list* II, __isl_take isl_id_list* I, enum tc_scheduling_enum scheduling);

void tc_scheduling_adapter_isl(struct tc_scop* scop, struct tc_options* options, __isl_take isl_union_set* LD, __isl_take isl_union_map* S, __isl_take isl_union_map* R, __isl_take isl_set* ii_set, __isl_take isl_set* tile, __isl_take isl_map* Rtile, __isl_take isl_id_list* II, __isl_take isl_id_list* I);

void tc_scheduling_adapter_feautrier(struct tc_scop* scop, struct tc_options* options, __isl_take isl_union_set* LD, __isl_take isl_union_map* S, __isl_take isl_union_map* R, __isl_take isl_set* ii_set, __isl_take isl_set* tile, __isl_take isl_map* Rtile, __isl_take isl_id_list* II, __isl_take isl_id_list* I);

void tc_scheduling_adapter_isl_wavefronting(struct tc_scop* scop, struct tc_options* options, __isl_take isl_union_set* LD, __isl_take isl_union_map* S, __isl_take isl_union_map* R, __isl_take isl_set* ii_set, __isl_take isl_set* tile, __isl_take isl_map* Rtile, __isl_take isl_id_list* II, __isl_take isl_id_list* I);

__isl_give isl_schedule* tc_schedule_compute(__isl_keep isl_union_set* IS, __isl_keep isl_union_map* R);

__isl_give isl_schedule* tc_schedule_compute_normalized(__isl_keep isl_set* IS, __isl_keep isl_map* R);

__isl_give isl_schedule* tc_schedule_compute_from_map(__isl_keep isl_map* R);

__isl_give isl_schedule* tc_schedule_compute_from_umap(__isl_keep isl_union_map* R);

__isl_give isl_map* tc_schedule_map(__isl_keep isl_schedule* S);

__isl_give isl_union_map* tc_schedule_umap(__isl_keep isl_schedule* schedule);

__isl_give isl_set* tc_schedule_validity(__isl_keep isl_map* m, __isl_keep isl_map* R);

__isl_give isl_set* tc_schedule_validity_bounded(__isl_keep isl_map* m, __isl_keep isl_map* R, __isl_keep isl_set* bounds);

__isl_give isl_map* tc_schedule_validity_map(__isl_keep isl_id_list* II, __isl_keep isl_map* R, __isl_keep isl_map* S);

isl_bool tc_schedule_validity_is_valid(__isl_take isl_set* validity);

isl_bool tc_schedule_validity_map_is_valid(__isl_take isl_map* validity);

#endif // TC_ISL_SCHEDULING_H
