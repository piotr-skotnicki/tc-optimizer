#include "tuples.h"
#include "utility.h"

std::string tc_tuples_eq(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs)
{
    std::string str;

    for (int i = 0; i < isl_id_list_n_id(lhs); ++i)
    {
        isl_id* lhs_id = isl_id_list_get_id(lhs, i);
        isl_id* rhs_id = isl_id_list_get_id(rhs, i);

        std::string eq = std::string(isl_id_get_name(lhs_id)) + " = " + std::string(isl_id_get_name(rhs_id));

        if (str.size() > 0)
        {
            str = tc_conjunction(str, eq);
        }
        else
        {
            str = eq;
        }

        isl_id_free(lhs_id);
        isl_id_free(rhs_id);
    }

    return str;
}

std::string tc_tuples_neq(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs)
{
    return tc_negation(tc_tuples_eq(lhs, rhs));
}

std::string tc_tuples_lt(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs)
{
    std::string str;

    for (int i = 0; i < isl_id_list_n_id(lhs); ++i)
    {
        isl_id* lhs_id_i = isl_id_list_get_id(lhs, i);
        isl_id* rhs_id_i = isl_id_list_get_id(rhs, i);

        str += "(";

        for (int j = 0; j < i; ++j)
        {
            isl_id* lhs_id_j = isl_id_list_get_id(lhs, j);
            isl_id* rhs_id_j = isl_id_list_get_id(rhs, j);

            str += std::string(isl_id_get_name(lhs_id_j)) + " = " + std::string(isl_id_get_name(rhs_id_j)) + " and ";

            isl_id_free(lhs_id_j);
            isl_id_free(rhs_id_j);
        }

        str += std::string(isl_id_get_name(lhs_id_i)) + " < " + std::string(isl_id_get_name(rhs_id_i)) + ") or ";

        isl_id_free(lhs_id_i);
        isl_id_free(rhs_id_i);
    }

    str += "false";

    return str;
}

std::string tc_tuples_gt(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs)
{
    std::string str;

    for (int i = 0; i < isl_id_list_n_id(lhs); ++i)
    {
        isl_id* lhs_id_i = isl_id_list_get_id(lhs, i);
        isl_id* rhs_id_i = isl_id_list_get_id(rhs, i);

        str += "(";

        for (int j = 0; j < i; ++j)
        {
            isl_id* lhs_id_j = isl_id_list_get_id(lhs, j);
            isl_id* rhs_id_j = isl_id_list_get_id(rhs, j);

            str += std::string(isl_id_get_name(lhs_id_j)) + " = " + std::string(isl_id_get_name(rhs_id_j)) + " and ";

            isl_id_free(lhs_id_j);
            isl_id_free(rhs_id_j);
        }

        str += std::string(isl_id_get_name(lhs_id_i)) + " > " + std::string(isl_id_get_name(rhs_id_i)) + ") or ";

        isl_id_free(lhs_id_i);
        isl_id_free(rhs_id_i);
    }

    str += "false";

    return str;
}

std::string tc_tuples_tlt(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs)
{
    isl_id_list* ids_0_lhs = isl_id_list_from_id(isl_id_list_get_id(lhs, 0));
    isl_id_list* ids_0_rhs = isl_id_list_from_id(isl_id_list_get_id(rhs, 0));

    std::string str = tc_conjunction(tc_tuples_lt(lhs, rhs), tc_tuples_eq(ids_0_lhs, ids_0_rhs));

    isl_id_list_free(ids_0_lhs);
    isl_id_list_free(ids_0_rhs);

    return str;
}

std::string tc_tuples_tgt(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs)
{
    isl_id_list* ids_0_lhs = isl_id_list_from_id(isl_id_list_get_id(lhs, 0));
    isl_id_list* ids_0_rhs = isl_id_list_from_id(isl_id_list_get_id(rhs, 0));

    std::string str = tc_conjunction(tc_tuples_gt(lhs, rhs), tc_tuples_eq(ids_0_lhs, ids_0_rhs));

    isl_id_list_free(ids_0_lhs);
    isl_id_list_free(ids_0_rhs);

    return str;
}

/*
std::string tc_tuples_tslice(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs)
{
    isl_id_list* ids_0_lhs = isl_id_list_from_id(isl_id_list_get_id(lhs, 0));
    isl_id_list* ids_0_rhs = isl_id_list_from_id(isl_id_list_get_id(rhs, 0));

    std::string str = tc_tuples_eq(ids_0_lhs, ids_0_rhs);

    isl_id_list_free(ids_0_lhs);
    isl_id_list_free(ids_0_rhs);

    return str;
}

std::string tc_tuples_tslice_next(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs)
{
    isl_id_list* ids_0_lhs = isl_id_list_from_id(isl_id_list_get_id(lhs, 0));
    isl_id_list* ids_0_rhs = isl_id_list_from_id(isl_id_list_get_id(rhs, 0));

    std::string str = tc_conjunction(tc_tuples_gt(lhs, rhs), tc_tuples_eq(ids_0_lhs, ids_0_rhs));

    isl_id_list_free(ids_0_lhs);
    isl_id_list_free(ids_0_rhs);

    return str;
}
*/

std::string tc_tuples_op(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs, const char* op)
{
    std::string str;

    for (int i = 0; i < isl_id_list_n_id(lhs); ++i)
    {
        isl_id* lhs_id_i = isl_id_list_get_id(lhs, i);
        isl_id* rhs_id_i = isl_id_list_get_id(rhs, i);

        str += std::string(isl_id_get_name(lhs_id_i)) + " " + std::string(op) + " " + std::string(isl_id_get_name(rhs_id_i)) + " and ";

        isl_id_free(lhs_id_i);
        isl_id_free(rhs_id_i);
    }

    str += "true";

    return str;
}

std::string tc_tuples_eq_param_0(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs)
{
    isl_id* lhs_0_id = isl_id_list_get_id(lhs, 0);
    isl_id* rhs_0_id = isl_id_list_get_id(rhs, 0);

    return "(" + std::string(isl_id_get_name(lhs_0_id)) + " = " + std::string(isl_id_get_name(rhs_0_id)) + ")";
}

std::string tc_tuples_next_param_0(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs)
{
    isl_id* lhs_0_id = isl_id_list_get_id(lhs, 0);
    isl_id* rhs_0_id = isl_id_list_get_id(rhs, 0);

    return "(" + std::string(isl_id_get_name(lhs_0_id)) + " = (" + std::string(isl_id_get_name(rhs_0_id)) + " + 1))";
}

std::string tc_tuples_prev_param_0(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs)
{
    isl_id* lhs_0_id = isl_id_list_get_id(lhs, 0);
    isl_id* rhs_0_id = isl_id_list_get_id(rhs, 0);

    return "((" + std::string(isl_id_get_name(lhs_0_id)) + " + 1) = " + std::string(isl_id_get_name(rhs_0_id)) + ")";
}

std::string tc_tuples_eq_params_except_1(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs)
{
    std::string str;

    for (int i = 0; i < isl_id_list_n_id(lhs); ++i)
    {
        isl_id* lhs_id_i = isl_id_list_get_id(lhs, i);
        isl_id* rhs_id_i = isl_id_list_get_id(rhs, i);

        if (i != 1)
        {
            str += "(" + std::string(isl_id_get_name(lhs_id_i)) + " = " + std::string(isl_id_get_name(rhs_id_i)) + ") and ";
        }

        isl_id_free(lhs_id_i);
        isl_id_free(rhs_id_i);
    }

    str += "true";

    return str;
}

std::string tc_tuples_next_param_0_eq_other_params_except_1(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs)
{
    std::string str;

    for (int i = 0; i < isl_id_list_n_id(lhs); ++i)
    {
        isl_id* lhs_id_i = isl_id_list_get_id(lhs, i);
        isl_id* rhs_id_i = isl_id_list_get_id(rhs, i);

        if (i == 0)
        {
            str += "(" + std::string(isl_id_get_name(lhs_id_i)) + " = (" + std::string(isl_id_get_name(rhs_id_i)) + " + 1)) and ";
        }
        else if (i != 1)
        {
            str += "(" + std::string(isl_id_get_name(lhs_id_i)) + " = " + std::string(isl_id_get_name(rhs_id_i)) + ") and ";
        }

        isl_id_free(lhs_id_i);
        isl_id_free(rhs_id_i);
    }

    str += "true";

    return str;
}

std::string tc_tuples_prev_param_0_eq_other_params_except_1(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs)
{
    std::string str;

    for (int i = 0; i < isl_id_list_n_id(lhs); ++i)
    {
        isl_id* lhs_id_i = isl_id_list_get_id(lhs, i);
        isl_id* rhs_id_i = isl_id_list_get_id(rhs, i);

        if (i == 0)
        {
            str += "((" + std::string(isl_id_get_name(lhs_id_i)) + " + 1) = " + std::string(isl_id_get_name(rhs_id_i)) + ") and ";
        }
        else if (i != 1)
        {
            str += "(" + std::string(isl_id_get_name(lhs_id_i)) + " = " + std::string(isl_id_get_name(rhs_id_i)) + ") and ";
        }

        isl_id_free(lhs_id_i);
        isl_id_free(rhs_id_i);
    }

    str += "true";

    return str;
}

std::string tc_tuples_eq_t_k_ii(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs)
{
    isl_id_list* ids_t_k_ii_lhs = tc_ids_sub(lhs, 0, 3);
    isl_id_list* ids_t_k_ii_rhs = tc_ids_sub(rhs, 0, 3);

    std::string str = tc_tuples_eq(ids_t_k_ii_lhs, ids_t_k_ii_rhs);

    isl_id_list_free(ids_t_k_ii_lhs);
    isl_id_list_free(ids_t_k_ii_rhs);

    return str;
}

std::string tc_tuples_t_k_ii_lt(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs)
{
    isl_id_list* ids_t_k_ii_lhs = tc_ids_sub(lhs, 0, 3);
    isl_id_list* ids_t_k_ii_rhs = tc_ids_sub(rhs, 0, 3);

    isl_id_list* ids_lhs_jj_kk = tc_ids_sub(lhs, 3, isl_id_list_n_id(lhs));
    isl_id_list* ids_rhs_jj_kk = tc_ids_sub(rhs, 3, isl_id_list_n_id(rhs));

    std::string str = tc_conjunction(tc_tuples_eq(ids_t_k_ii_lhs, ids_t_k_ii_rhs), tc_tuples_lt(ids_lhs_jj_kk, ids_rhs_jj_kk));

    isl_id_list_free(ids_t_k_ii_lhs);
    isl_id_list_free(ids_t_k_ii_rhs);
    isl_id_list_free(ids_lhs_jj_kk);
    isl_id_list_free(ids_rhs_jj_kk);

    return str;
}

std::string tc_tuples_t_k_ii_gt(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs)
{
    isl_id_list* ids_t_k_ii_lhs = tc_ids_sub(lhs, 0, 3);
    isl_id_list* ids_t_k_ii_rhs = tc_ids_sub(rhs, 0, 3);

    isl_id_list* ids_lhs_jj_kk = tc_ids_sub(lhs, 3, isl_id_list_n_id(lhs));
    isl_id_list* ids_rhs_jj_kk = tc_ids_sub(rhs, 3, isl_id_list_n_id(rhs));

    std::string str = tc_conjunction(tc_tuples_eq(ids_t_k_ii_lhs, ids_t_k_ii_rhs), tc_tuples_gt(ids_lhs_jj_kk, ids_rhs_jj_kk));

    isl_id_list_free(ids_t_k_ii_lhs);
    isl_id_list_free(ids_t_k_ii_rhs);
    isl_id_list_free(ids_lhs_jj_kk);
    isl_id_list_free(ids_rhs_jj_kk);

    return str;
}

std::string tc_tuples_eq_t_ii(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs)
{
    isl_id_list* ids_t_ii_lhs = tc_ids_sub(lhs, 0, 2);
    isl_id_list* ids_t_ii_rhs = tc_ids_sub(rhs, 0, 2);

    std::string str = tc_tuples_eq(ids_t_ii_lhs, ids_t_ii_rhs);

    isl_id_list_free(ids_t_ii_lhs);
    isl_id_list_free(ids_t_ii_rhs);

    return str;
}

std::string tc_tuples_t_ii_lt(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs)
{
    isl_id_list* ids_t_ii_lhs = tc_ids_sub(lhs, 0, 2);
    isl_id_list* ids_t_ii_rhs = tc_ids_sub(rhs, 0, 2);

    isl_id_list* ids_lhs_jj_kk = tc_ids_sub(lhs, 2, isl_id_list_n_id(lhs));
    isl_id_list* ids_rhs_jj_kk = tc_ids_sub(rhs, 2, isl_id_list_n_id(rhs));

    std::string str = tc_conjunction(tc_tuples_eq(ids_t_ii_lhs, ids_t_ii_rhs), tc_tuples_lt(ids_lhs_jj_kk, ids_rhs_jj_kk));

    isl_id_list_free(ids_t_ii_lhs);
    isl_id_list_free(ids_t_ii_rhs);
    isl_id_list_free(ids_lhs_jj_kk);
    isl_id_list_free(ids_rhs_jj_kk);

    return str;
}

std::string tc_tuples_t_ii_gt(__isl_keep isl_id_list* lhs, __isl_keep isl_id_list* rhs)
{
    isl_id_list* ids_t_ii_lhs = tc_ids_sub(lhs, 0, 2);
    isl_id_list* ids_t_ii_rhs = tc_ids_sub(rhs, 0, 2);

    isl_id_list* ids_lhs_jj_kk = tc_ids_sub(lhs, 2, isl_id_list_n_id(lhs));
    isl_id_list* ids_rhs_jj_kk = tc_ids_sub(rhs, 2, isl_id_list_n_id(rhs));

    std::string str = tc_conjunction(tc_tuples_eq(ids_t_ii_lhs, ids_t_ii_rhs), tc_tuples_gt(ids_lhs_jj_kk, ids_rhs_jj_kk));

    isl_id_list_free(ids_t_ii_lhs);
    isl_id_list_free(ids_t_ii_rhs);
    isl_id_list_free(ids_lhs_jj_kk);
    isl_id_list_free(ids_rhs_jj_kk);

    return str;
}
