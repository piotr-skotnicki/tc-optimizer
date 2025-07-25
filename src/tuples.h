#ifndef TC_TUPLES_H
#define TC_TUPLES_H

#include <isl/id.h>

#include <string>

std::string tc_tuples_eq(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs);

std::string tc_tuples_neq(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs);

std::string tc_tuples_lt(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs);

std::string tc_tuples_gt(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs);

std::string tc_tuples_tlt(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs);

std::string tc_tuples_tgt(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs);

/*
std::string tc_tuples_tslice(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs);

std::string tc_tuples_tslice_next(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs);
*/

std::string tc_tuples_op(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs, const char* op);

std::string tc_tuples_eq_param_0(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs);

std::string tc_tuples_next_param_0(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs);

std::string tc_tuples_prev_param_0(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs);

std::string tc_tuples_eq_params_except_1(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs);

std::string tc_tuples_next_param_0_eq_other_params_except_1(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs);

std::string tc_tuples_prev_param_0_eq_other_params_except_1(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs);

std::string tc_tuples_eq_t_k_ii(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs);

std::string tc_tuples_t_k_ii_lt(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs);

std::string tc_tuples_t_k_ii_gt(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs);

std::string tc_tuples_eq_t_ii(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs);

std::string tc_tuples_t_ii_lt(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs);

std::string tc_tuples_t_ii_gt(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs);

#endif // TC_TUPLES_H
