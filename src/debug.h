#ifndef TC_DEBUG_H
#define TC_DEBUG_H

#include <isl/space.h>
#include <isl/set.h>
#include <isl/map.h>
#include <isl/union_set.h>
#include <isl/union_map.h>
#include <isl/val.h>
#include <isl/polynomial.h>
#include <isl/schedule.h>
#include <isl/aff.h>
#include <isl/id.h>
#include <isl/ast.h>

#include <barvinok/isl.h>

extern int tc_debug_flag;

void tc_debug(const char* msg, ...);

void tc_warn(const char* msg, ...);

void tc_error(const char* msg, ...);

void tc_assert(int condition, const char* msg, ...);

void tc_die(int status);

void tc_debug_space(__isl_keep isl_space* set, const char* msg, ...);

void tc_debug_bset(__isl_keep isl_basic_set* bset, const char* msg, ...);

void tc_debug_bset_latex(__isl_keep isl_basic_set* bset, const char* msg, ...);

void tc_debug_set(__isl_keep isl_set* set, const char* msg, ...);

void tc_debug_set_latex(__isl_keep isl_set* set, const char* msg, ...);

void tc_debug_uset(__isl_keep isl_union_set* uset, const char* msg, ...);

void tc_debug_uset_latex(__isl_keep isl_union_set* uset, const char* msg, ...);

void tc_debug_bmap(__isl_keep isl_basic_map* bmap, const char* msg, ...);

void tc_debug_bmap_latex(__isl_keep isl_basic_map* bmap, const char* msg, ...);

void tc_debug_map(__isl_keep isl_map* map, const char* msg, ...);

void tc_debug_map_latex(__isl_keep isl_map* map, const char* msg, ...);

void tc_debug_umap(__isl_keep isl_union_map* umap, const char* msg, ...);

void tc_debug_umap_latex(__isl_keep isl_union_map* umap, const char* msg, ...);

void tc_debug_id(__isl_keep isl_id* id, const char* msg, ...);

void tc_debug_id_list(__isl_keep isl_id_list* list, const char* msg, ...);

void tc_debug_val(__isl_keep isl_val* val, const char* msg, ...);

void tc_debug_bool(isl_bool condition, const char* msg, ...);

void tc_debug_set_card(__isl_keep isl_set* set, const char* msg, ...);

void tc_debug_map_card(__isl_keep isl_map* map, const char* msg, ...);

void tc_debug_pw_qpolynomial(__isl_keep isl_pw_qpolynomial* qpoly, const char* msg, ...);

void tc_debug_qpolynomial(__isl_keep isl_qpolynomial* poly, const char* msg, ...);

void tc_debug_aff(__isl_keep isl_aff* aff, const char* msg, ...);

void tc_debug_multi_aff(__isl_keep isl_multi_aff* aff, const char* msg, ...);

void tc_debug_pw_aff(__isl_keep isl_pw_aff* aff, const char* msg, ...);

void tc_debug_pw_multi_aff(__isl_keep isl_pw_multi_aff* maff, const char* msg, ...);

void tc_debug_union_pw_multi_aff(__isl_keep isl_union_pw_multi_aff* umaff, const char* msg, ...);

void tc_debug_multi_union_pw_aff(__isl_keep isl_multi_union_pw_aff* muaff, const char* msg, ...);

void tc_debug_schedule(__isl_keep isl_schedule* schedule, const char* msg, ...);

void tc_debug_ast_expr(__isl_keep isl_ast_expr* expr, const char* msg, ...);

#endif // TC_DEBUG_H
