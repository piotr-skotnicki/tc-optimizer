#ifndef TC_DEBUG_H
#define TC_DEBUG_H

#include <isl/set.h>
#include <isl/map.h>
#include <isl/union_set.h>
#include <isl/union_map.h>
#include <isl/val.h>
#include <isl/polynomial.h>

#include <barvinok/isl.h>

extern int tc_debug_flag;

void tc_debug(const char* msg, ...);

void tc_debug_set(__isl_keep isl_set* set, const char* msg, ...);

void tc_debug_set_latex(__isl_keep isl_set* set, const char* msg, ...);

void tc_debug_uset(__isl_keep isl_union_set* uset, const char* msg, ...);

void tc_debug_uset_latex(__isl_keep isl_union_set* uset, const char* msg, ...);

void tc_debug_map(__isl_keep isl_map* map, const char* msg, ...);

void tc_debug_map_latex(__isl_keep isl_map* map, const char* msg, ...);

void tc_debug_umap(__isl_keep isl_union_map* umap, const char* msg, ...);

void tc_debug_umap_latex(__isl_keep isl_union_map* umap, const char* msg, ...);

void tc_debug_val(__isl_keep isl_val* val, const char* msg, ...);

void tc_debug_bool(isl_bool condition, const char* msg, ...);

void tc_debug_set_card(__isl_keep isl_set* set, const char* msg, ...);

void tc_debug_map_card(__isl_keep isl_map* map, const char* msg, ...);

void tc_debug_pw_qpolynomial(__isl_keep isl_pw_qpolynomial* qpoly, const char* msg, ...);

void tc_debug_qpolynomial(__isl_keep isl_qpolynomial* poly, const char* msg, ...);

#endif // TC_DEBUG_H
