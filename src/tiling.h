#ifndef TC_TILING_H
#define TC_TILING_H

#include <isl/id.h>
#include <isl/set.h>
#include <isl/map.h>
#include <isl/union_map.h>

#include <vector>
#include <string>
#include <map>

__isl_give isl_set* tc_tile_set(__isl_keep isl_id_list* II, __isl_keep isl_id_list* I, const std::vector<int>& BLOCK, __isl_keep isl_set* set);

__isl_give isl_set* tc_ii_set_set(__isl_keep isl_id_list* II, const std::vector<int>& BLOCK, __isl_keep isl_set* set);

__isl_give isl_set* tc_tile_set_of(__isl_keep isl_set* set, __isl_keep isl_set* ii_set, __isl_keep isl_id_list* II, std::string(*f)(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs));

__isl_give isl_set* tc_tile_lt_set(__isl_keep isl_set* tile, __isl_keep isl_set* ii_set, __isl_keep isl_id_list* II);

__isl_give isl_set* tc_tile_gt_set(__isl_keep isl_set* tile, __isl_keep isl_set* ii_set, __isl_keep isl_id_list* II);

__isl_give isl_map* tc_get_tile_schedule(const char* statement_label, /*const std::map<std::string, std::vector<std::string> >& groups,*/ __isl_keep isl_union_map* S);

__isl_give isl_set* tc_normalize_params(__isl_take isl_set* tile, __isl_keep isl_map* tile_schedule, __isl_keep isl_id_list* IIprim, __isl_keep isl_id_list* II);

__isl_give isl_map* tc_Rtile_map(__isl_keep isl_set* tile, __isl_keep isl_map* R);

void tc_tile_loop_nest(__isl_keep isl_union_set* LD, __isl_keep isl_union_map* S, const std::map<std::string, std::vector<int> >& blocks, __isl_keep isl_id_list* II, __isl_keep isl_id_list* I, __isl_give isl_set** tile, __isl_give isl_set** ii_set);

#endif // TC_TILING_H
