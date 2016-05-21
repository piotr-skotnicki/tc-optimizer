#ifndef TC_UTILITY_H
#define TC_UTILITY_H

#include <isl/ctx.h>
#include <isl/id.h>
#include <isl/set.h>
#include <isl/map.h>
#include <isl/union_set.h>
#include <isl/union_map.h>
#include <isl/point.h>

#include <barvinok/isl.h>

#include <string>

__isl_give isl_set* tc_rename_dim(__isl_take isl_set* set, isl_dim_type dim, __isl_keep isl_id_list* from, __isl_keep isl_id_list* to);

__isl_give isl_set* tc_rename_params(__isl_take isl_set* set, __isl_keep isl_id_list* from, __isl_keep isl_id_list* to);

__isl_give isl_map* tc_map_rename_dim(__isl_take isl_map* map, isl_dim_type dim, __isl_keep isl_id_list* from, __isl_keep isl_id_list* to);

__isl_give isl_map* tc_map_rename_params(__isl_take isl_map* map, __isl_keep isl_id_list* from, __isl_keep isl_id_list* to);

__isl_give isl_set* tc_project_out_dim_names(__isl_take isl_set* set, isl_dim_type dim, __isl_keep isl_id_list* names);

__isl_give isl_set* tc_project_out_params(__isl_take isl_set* set, __isl_keep isl_id_list* names);

__isl_give isl_set* tc_project_out_dims_except_pos(__isl_take isl_set* set, isl_dim_type dim, int pos);

__isl_give isl_set* tc_project_out_params_except_pos(__isl_take isl_set* set, int pos);

__isl_give isl_set* tc_project_out_dims_except(__isl_take isl_set* set, isl_dim_type dim, __isl_take isl_id* id);

__isl_give isl_set* tc_project_out_params_except(__isl_take isl_set* set, __isl_take isl_id* id);

__isl_give isl_union_set* tc_union_set_project_out_dim_names(__isl_take isl_union_set* uset, isl_dim_type dim, __isl_keep isl_id_list* names);

__isl_give isl_union_set* tc_union_set_project_out_params(__isl_take isl_union_set* uset, __isl_keep isl_id_list* names);

__isl_give isl_map* tc_map_project_out_dim_names(__isl_take isl_map* map, isl_dim_type dim, __isl_keep isl_id_list* names);

__isl_give isl_map* tc_map_project_out_params(__isl_take isl_map* map, __isl_keep isl_id_list* names);

__isl_give isl_union_map* tc_union_map_project_out_dim_names(__isl_take isl_union_map* umap, isl_dim_type dim, __isl_keep isl_id_list* names);

__isl_give isl_union_map* tc_union_map_project_out_params(__isl_take isl_union_map* umap, __isl_keep isl_id_list* names);

__isl_give isl_map* tc_map_set_dim_names(__isl_take isl_map* map, isl_dim_type dim, __isl_keep isl_id_list* names);

__isl_give isl_set* tc_flatten_union_set(__isl_take isl_union_set* uset);

__isl_give isl_map* tc_flatten_union_map(__isl_take isl_union_map* umap);

__isl_give isl_set* tc_lift_up_set_params(__isl_take isl_set* set, __isl_keep isl_id_list* params);

__isl_give isl_set* tc_lift_down_all_set_vars(__isl_take isl_set* set, __isl_keep isl_id_list* vars);

__isl_give isl_map* tc_lift_down_all_map_vars(__isl_take isl_map* map
                                            , __isl_keep isl_id_list* in
                                            , __isl_keep isl_id_list* out);

std::string tc_comma(__isl_keep isl_id_list* list);

std::string tc_parens_comma(__isl_keep isl_id_list* list);

__isl_give isl_set* tc_make_set(__isl_keep isl_ctx* ctx
                              , __isl_keep isl_id_list* params
                              , __isl_keep isl_id_list* vars
                              , const char* constraints);

__isl_give isl_union_set* tc_make_union_set(__isl_keep isl_ctx* ctx
                                          , __isl_keep isl_id_list* params
                                          , __isl_keep isl_id_list* vars
                                          , const char* constraints);

__isl_give isl_map* tc_make_map(__isl_keep isl_ctx* ctx
                              , __isl_keep isl_id_list* params
                              , __isl_keep isl_id_list* in
                              , __isl_keep isl_id_list* out
                              , const char* constraints);

__isl_give isl_union_map* tc_make_union_map(__isl_keep isl_ctx* ctx
                                          , __isl_keep isl_id_list* params
                                          , __isl_keep isl_id_list* in
                                          , __isl_keep isl_id_list* out
                                          , const char* constraints);

__isl_give isl_set* tc_make_params(__isl_keep isl_ctx* ctx
                                 , __isl_keep isl_id_list* params
                                 , const char* constraints);

__isl_give isl_set* tc_parameterize(__isl_take isl_set* set, int pos, __isl_take isl_id* name);

__isl_give isl_set* tc_parameterize_all(__isl_take isl_set* set, __isl_keep isl_id_list* names);

__isl_give isl_map* tc_parameterize_map_all_in(__isl_take isl_map* set, __isl_keep isl_id_list* names);

__isl_give isl_id_list* tc_get_set_params_names(__isl_keep isl_set* set);

__isl_give isl_id_list* tc_get_union_set_params_names(__isl_keep isl_union_set* uset);

__isl_give isl_id_list* tc_get_map_params_names(__isl_keep isl_map* map);

__isl_give isl_id_list* tc_get_union_map_params_names(__isl_keep isl_union_map* umap);

__isl_give isl_set* tc_make_set_constraints(__isl_take isl_set* set, __isl_keep isl_id_list* names);

__isl_give isl_set* tc_make_set_constraints_for(__isl_take isl_set* set, int pos, __isl_take isl_id* name);

__isl_give isl_set* tc_make_map_constraints(__isl_take isl_map* map, __isl_keep isl_id_list* in, __isl_keep isl_id_list* out);

__isl_give isl_set* tc_set_fix_param_value(__isl_take isl_set* set, __isl_take isl_id* name, int value);

__isl_give isl_map* tc_map_fix_param_value(__isl_take isl_map* map, __isl_take isl_id* name, int value);

__isl_give isl_set* tc_set_fix_params_bounds(__isl_take isl_set* set, __isl_take isl_set* bounds);

__isl_give isl_map* tc_map_fix_params_bounds(__isl_take isl_map* map, __isl_take isl_set* bounds);

__isl_give isl_union_set* tc_union_set_fix_params_bounds(__isl_take isl_union_set* uset, __isl_take isl_set* bounds);

__isl_give isl_union_map* tc_union_map_fix_params_bounds(__isl_take isl_union_map* umap, __isl_take isl_set* bounds);

__isl_give isl_map* tc_make_identity(__isl_take isl_map* map);

void tc_scan_set(__isl_keep isl_set* set);

void tc_scan_map(__isl_keep isl_map* map);

void tc_scan_union_set(__isl_keep isl_union_set* uset);

void tc_scan_union_map(__isl_keep isl_union_map* umap);

__isl_give isl_map* tc_map_closure(__isl_take isl_map* map, int k);

__isl_give isl_union_map* tc_union_map_closure(__isl_take isl_union_map* umap, int k);

__isl_give isl_id_list* tc_ids_sequence(__isl_keep isl_ctx* ctx, const char* id, int size);

__isl_give isl_id_list* tc_ids_add_suffix(__isl_keep isl_id_list* list, const char* suffix);

__isl_give isl_id_list* tc_ids_double(__isl_keep isl_id_list* list);

__isl_give isl_id_list* tc_ids_prim(__isl_keep isl_id_list* list);

__isl_give isl_id_list* tc_ids_bis(__isl_keep isl_id_list* list);

__isl_give isl_id_list* tc_ids_ter(__isl_keep isl_id_list* list);

__isl_give isl_id_list* tc_ids_sub(__isl_keep isl_id_list* list, int begin, int end);

std::string tc_conjunction(const std::string& lhs, const std::string& rhs);

std::string tc_disjunction(const std::string& lhs, const std::string& rhs);

std::string tc_negation(const std::string& arg);

std::string tc_tuples_eq(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs);

std::string tc_tuples_neq(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs);

std::string tc_tuples_lt(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs);

std::string tc_tuples_gt(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs);

std::string tc_tuples_op(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs, const char* op);

__isl_give isl_set* tc_lexmin_set_pos(__isl_take isl_set* set, int pos);

__isl_give isl_set* tc_lexmax_set_pos(__isl_take isl_set* set, int pos);

__isl_give isl_set* tc_lexmin_set_var(__isl_take isl_set* set, __isl_take isl_id* id);

__isl_give isl_set* tc_lexmax_set_var(__isl_take isl_set* set, __isl_take isl_id* id);

long tc_lexmin_set_pos_value(__isl_keep isl_set* set, int pos);

long tc_lexmax_set_pos_value(__isl_keep isl_set* set, int pos);

long tc_lexmin_param_value(__isl_keep isl_set* set, __isl_take isl_id* param);

long tc_lexmax_param_value(__isl_keep isl_set* set, __isl_take isl_id* param);

__isl_give isl_set* tc_get_set_bounds(__isl_keep isl_set* set, __isl_keep isl_id_list* LB, __isl_keep isl_id_list* UB);

__isl_give isl_map* tc_get_map_for_input_tuple(__isl_keep isl_union_map* umap, const char* name);

__isl_give isl_union_map* tc_remove_map_with_tuple(__isl_take isl_union_map* umap, const char* name);

int tc_get_statement_depth(const char* label, __isl_keep isl_union_map* umap);

__isl_give isl_set* tc_normalize_union_set(__isl_keep isl_union_set* uset, __isl_keep isl_union_map* S);

__isl_give isl_map* tc_normalize_union_map(__isl_keep isl_union_map* umap, __isl_keep isl_union_map* S);

__isl_give isl_set* tc_normalize_set(__isl_keep isl_set* set, __isl_keep isl_union_map* S);

__isl_give isl_map* tc_normalize_map(__isl_keep isl_map* map, __isl_keep isl_union_map* S);

__isl_give isl_union_set* tc_denormalize_set(__isl_keep isl_set* set, __isl_keep isl_union_map* S);

__isl_give isl_union_map* tc_denormalize_map(__isl_keep isl_map* map, __isl_keep isl_union_map* S);

__isl_give isl_basic_set_list* tc_collect_basic_sets(__isl_keep isl_set* set);

__isl_give isl_set_list* tc_collect_sets(__isl_keep isl_union_set* uset);

__isl_give isl_basic_map_list* tc_collect_basic_maps(__isl_keep isl_map* map);

__isl_give isl_map_list* tc_collect_maps(__isl_keep isl_union_map* umap);

__isl_give isl_union_map* tc_unwrap_range(__isl_take isl_union_map* umap);

struct tc_qpolynomials
{
    isl_qpolynomial** polys;
    
    isl_set_list* domains;
};

struct tc_qpolynomials* tc_collect_qpolynomials(__isl_keep isl_pw_qpolynomial* pw);

void tc_qpolynomials_free(struct tc_qpolynomials* polys);

__isl_give isl_point* tc_set_to_point(__isl_take isl_set* set);

isl_bool tc_points_compare(__isl_keep isl_point* lhs, __isl_keep isl_point* rhs);

__isl_give isl_map* tc_get_lex_forward(__isl_keep isl_map* R);

__isl_give isl_map* tc_get_lex_backward(__isl_keep isl_map* R);

isl_bool tc_is_lex_forward(__isl_keep isl_map* R);

__isl_give isl_union_map* tc_simplify_schedule(__isl_take isl_union_map* S);

__isl_give isl_union_map* tc_extend_schedule(__isl_take isl_union_map* S, int n);

__isl_give isl_union_map* tc_extend_union_map(__isl_take isl_union_map* umap, int n);

__isl_give isl_set* tc_get_params_set(__isl_take isl_set* set, __isl_keep isl_id_list* params);

isl_bool tc_map_carries_dependences(__isl_keep isl_map* map, int pos);

#endif // TC_UTILITY_H
