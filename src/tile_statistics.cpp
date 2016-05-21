#include "tile_statistics.h"
#include "utility.h"

#include <isl/ctx.h>
#include <isl/id.h>
#include <isl/set.h>
#include <isl/map.h>
#include <isl/union_map.h>
#include <isl/polynomial.h>
#include <isl/val.h>
#include <isl/ast_build.h>
#include <isl/printer.h>

#include <barvinok/isl.h>

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

static int tc_tile_statistics_block_compare(const void* a, const void* b)
{
    const struct tc_tile_statistics_block* const* lhs = (const struct tc_tile_statistics_block* const*)a;
    const struct tc_tile_statistics_block* const* rhs = (const struct tc_tile_statistics_block* const*)b;
        
    return strcmp(isl_id_get_name((*lhs)->label), isl_id_get_name((*rhs)->label));
}

static int tc_tile_statistics_group_compare(const void* a, const void* b)
{
    const struct tc_tile_statistics_group* const* lhs = (const struct tc_tile_statistics_group* const*)a;
    const struct tc_tile_statistics_group* const* rhs = (const struct tc_tile_statistics_group* const*)b;
        
    return ((*lhs)->card * (*lhs)->n_occurrences) < ((*rhs)->card * (*rhs)->n_occurrences);
}

struct tc_tile_statistics_user
{
    isl_set* tile;
    
    isl_set* ii_set;
    
    isl_id_list* II;
    
    isl_union_map* S;
    
    struct tc_tile_statistics* stats;
    
    enum tc_tile_statistics_category_enum category;
    
    isl_qpolynomial* card_poly;
};

static isl_stat tc_get_tile_statistics_callback(__isl_take isl_point* point, void* user)
{    
    struct tc_tile_statistics_user* data = (struct tc_tile_statistics_user*)user;
    
    isl_id_list* II = data->II;
    isl_set* tile = data->tile;
    isl_union_map* S = data->S;
    struct tc_tile_statistics* stats = data->stats;
        
    isl_ctx* ctx = isl_point_get_ctx(point);
    
    isl_set* point_set = isl_set_from_point(point);
    
    isl_set* point_constraints = tc_make_set_constraints(point_set, II);
        
    isl_set* tile_ii = isl_set_intersect_params(isl_set_copy(tile), point_constraints);
    
    if (isl_set_is_empty(tile_ii))
    {
        isl_set_free(tile_ii);
        
        return isl_stat_ok;
    }
        
    struct tc_tile_statistics_group* group = (struct tc_tile_statistics_group*)malloc(sizeof(struct tc_tile_statistics_group));
    group->blocks = NULL;
    group->n_blocks = 0;
    group->n_blocks_capacity = 0;
    group->card = 0;
    group->n_occurrences = 1;
    group->tile_sample = isl_set_copy(tile_ii);
    group->category = data->category;
    group->card_qpolynomial = isl_qpolynomial_copy(data->card_poly);
                        
    isl_union_set* tile_ii_denormalized = tc_denormalize_set(tile_ii, S);
       
    isl_set_list* tile_ii_sets = tc_collect_sets(tile_ii_denormalized);
    
    for (int i = 0; i < isl_set_list_n_set(tile_ii_sets); ++i)
    {
        isl_set* tile_set = isl_set_list_get_set(tile_ii_sets, i);
                
        struct tc_tile_statistics_block* block = (struct tc_tile_statistics_block*)malloc(sizeof(struct tc_tile_statistics_block));
        block->card = 0;
        block->label = isl_set_get_tuple_id(tile_set);
        block->dimensions = NULL;
        
        isl_pw_qpolynomial* tile_set_card = isl_set_card(isl_set_copy(tile_set));
        
        isl_val* tile_set_card_val = isl_pw_qpolynomial_max(tile_set_card);

        block->card = isl_val_get_num_si(tile_set_card_val);
        
        isl_val_free(tile_set_card_val);

        block->dimensions = isl_vec_alloc(ctx, isl_set_n_dim(tile_set));

        for (int j = 0; j < isl_set_n_dim(tile_set); ++j)
        {
            int dim = tc_lexmax_set_pos_value(tile_set, j) - tc_lexmin_set_pos_value(tile_set, j) + 1;

            block->dimensions = isl_vec_set_element_si(block->dimensions, j, dim);
        }

        /*
        if (NULL == group->blocks)
        {
            group->blocks = (struct tc_tile_statistics_block**)malloc(sizeof(struct tc_tile_statistics_block*));
            group->n_blocks = 1;
        }
        else
        {
            group->blocks = (struct tc_tile_statistics_block**)realloc(group->blocks, (group->n_blocks + 1) * sizeof(struct tc_tile_statistics_block*));                        
            group->n_blocks = group->n_blocks + 1;
        }
        group->blocks[group->n_blocks - 1] = block;
        */
        
        if (NULL == group->blocks)
        {
            group->n_blocks_capacity = 1;
            group->blocks = (struct tc_tile_statistics_block**)malloc(sizeof(struct tc_tile_statistics_block*));
        }
        else
        {
            if (group->n_blocks + 1 > group->n_blocks_capacity)
            {
                group->n_blocks_capacity = group->n_blocks_capacity * 2;
                group->blocks = (struct tc_tile_statistics_block**)realloc(group->blocks, (group->n_blocks_capacity) * sizeof(struct tc_tile_statistics_block*));
            }
        }
        group->n_blocks = group->n_blocks + 1;
        group->blocks[group->n_blocks - 1] = block;
        
        
        group->card += block->card;
        
        isl_set_free(tile_set);
    }
    
    qsort(group->blocks, group->n_blocks, sizeof(group->blocks[0]), &tc_tile_statistics_block_compare);
    
    isl_bool group_found = isl_bool_false;
    
    for (int i = 0; i < stats->n_groups; ++i)
    {
        struct tc_tile_statistics_group* existing_group = stats->groups[i];
        
        isl_bool groups_equal = isl_bool_true;
        
        if (group->n_blocks == existing_group->n_blocks
            && group->card == existing_group->card
            && group->category == existing_group->category)
        {
            for (int j = 0; j < group->n_blocks; ++j)
            {
                struct tc_tile_statistics_block* block = group->blocks[j];                
                struct tc_tile_statistics_block* existing_block = existing_group->blocks[j];
                
                if (block->card != existing_block->card
                    || 0 != strcmp(isl_id_get_name(block->label), isl_id_get_name(existing_block->label))
                    || !isl_vec_is_equal(block->dimensions, existing_block->dimensions))
                {
                    groups_equal = isl_bool_false;
                    break;
                }
            }
        }
        else
        {
            groups_equal = isl_bool_false;
        }
        
        if (groups_equal)
        {
            existing_group->n_occurrences = existing_group->n_occurrences + 1;
            group_found = isl_bool_true;
        }
    }

    if (!group_found)
    {
        if (NULL == stats->groups)
        {
            stats->n_groups_capacity = 1;
            stats->groups = (struct tc_tile_statistics_group**)malloc(sizeof(struct tc_tile_statistics_group*));
        }
        else
        {
            if (stats->n_groups + 1 > stats->n_groups_capacity)
            {
                stats->n_groups_capacity = stats->n_groups_capacity * 2;
                stats->groups = (struct tc_tile_statistics_group**)realloc(stats->groups, (stats->n_groups_capacity) * sizeof(struct tc_tile_statistics_group*));
            }
        }
        stats->n_groups = stats->n_groups + 1;
        stats->groups[stats->n_groups - 1] = group;
    }
    else
    {
        tc_tile_statistics_group_free(group);
    }
    
    isl_set_list_free(tile_ii_sets);    
    isl_union_set_free(tile_ii_denormalized);
    isl_set_free(tile_ii);
    
    return isl_stat_ok;
}

struct tc_tile_statistics* tc_compute_tile_statistics(__isl_keep isl_set* tile, __isl_keep isl_set* ii_set, __isl_keep isl_id_list* II, __isl_keep isl_set* bounds, __isl_keep isl_union_set* LD, __isl_keep isl_union_map* S, struct tc_scop* scop, const std::map<std::string, std::vector<int> >& blocks)
{    
    struct tc_tile_statistics* stats = (struct tc_tile_statistics*)malloc(sizeof(struct tc_tile_statistics));
    stats->scop = scop;
    stats->LD = isl_union_set_copy(LD);
    stats->S = isl_union_map_copy(S);
    stats->bounds = bounds;
    stats->n_statement_instances = 0;
    stats->groups = NULL;
    stats->n_groups = 0;
    stats->n_groups_capacity = 0;
    stats->statements = NULL;
    stats->n_statements = 0;
    stats->blocks = &blocks;
    
    isl_ctx* ctx = isl_set_get_ctx(tile);
    
    isl_set* tile_bounded = tc_set_fix_params_bounds(isl_set_copy(tile), isl_set_copy(bounds));
        
    isl_set* ii_set_bounded = tc_set_fix_params_bounds(isl_set_copy(ii_set), isl_set_copy(bounds));
    
    isl_set* LD_normalized = tc_normalize_union_set(LD, S);
    
    isl_set* LD_normalized_bounded = tc_set_fix_params_bounds(LD_normalized, isl_set_copy(bounds));
    
    isl_pw_qpolynomial* LD_card = isl_set_card(LD_normalized_bounded);
    
    isl_val* LD_card_val = isl_pw_qpolynomial_max(LD_card);
    
    stats->n_statement_instances = isl_val_get_num_si(LD_card_val);
    
    isl_val_free(LD_card_val);
    
    isl_set_list* LD_sets = tc_collect_sets(LD);
    
    stats->n_statements = isl_set_list_n_set(LD_sets);
    stats->statements = (struct tc_tile_statistics_block**)malloc(stats->n_statements * sizeof(struct tc_tile_statistics_block*));    
    
    for (int i = 0; i < isl_set_list_n_set(LD_sets); ++i)
    {
        isl_set* LD_set = isl_set_list_get_set(LD_sets, i);
        
        LD_set = tc_set_fix_params_bounds(LD_set, isl_set_copy(bounds));
        
        struct tc_tile_statistics_block* statement = (struct tc_tile_statistics_block*)malloc(sizeof(struct tc_tile_statistics_block));
        
        statement->label = isl_set_get_tuple_id(LD_set);
        
        isl_pw_qpolynomial* LD_set_card = isl_set_card(isl_set_copy(LD_set));
    
        isl_val* LD_set_card_val = isl_pw_qpolynomial_max(LD_set_card);

        statement->card = isl_val_get_num_si(LD_set_card_val);

        isl_val_free(LD_set_card_val);
        
        statement->dimensions = isl_vec_alloc(ctx, isl_set_n_dim(LD_set));

        for (int j = 0; j < isl_set_n_dim(LD_set); ++j)
        {
            int dim = tc_lexmax_set_pos_value(LD_set, j) - tc_lexmin_set_pos_value(LD_set, j) + 1;

            statement->dimensions = isl_vec_set_element_si(statement->dimensions, j, dim);
        }
        
        stats->statements[i] = statement;
        
        isl_set_free(LD_set);
    }
    
    isl_set_list_free(LD_sets);
    
    qsort(stats->statements, stats->n_statements, sizeof(stats->statements[0]), &tc_tile_statistics_block_compare);
   
    struct tc_tile_statistics_user data;
    data.tile = tile_bounded;
    data.ii_set = ii_set_bounded;
    data.II = II;
    data.S = S;
    data.stats = stats;
        
    isl_pw_qpolynomial* tile_card = isl_set_card(isl_set_copy(tile));
    
    isl_id_list* params = tc_get_union_set_params_names(LD);
    
    struct tc_qpolynomials* qpolys = tc_collect_qpolynomials(tile_card);
    
    for (int i = 0; i < isl_set_list_n_set(qpolys->domains); ++i)
    {        
        isl_qpolynomial* poly = qpolys->polys[i];
        isl_set* domain = isl_set_list_get_set(qpolys->domains, i);
        
        isl_set* tile_bounded_repr = isl_set_intersect_params(isl_set_copy(tile_bounded), isl_set_copy(domain));
        tile_bounded_repr = tc_set_fix_params_bounds(tile_bounded_repr, isl_set_copy(bounds));

        isl_set* ii_set_bounded_repr = tc_lift_up_set_params(isl_set_copy(domain), II);
        ii_set_bounded_repr = tc_set_fix_params_bounds(ii_set_bounded_repr, isl_set_copy(bounds));
        ii_set_bounded_repr = isl_set_intersect(ii_set_bounded_repr, isl_set_copy(ii_set_bounded));
        
        if (!isl_set_is_empty(tile_bounded_repr))
        {
            isl_bool is_parametric = isl_bool_false;
            isl_bool is_varied = isl_bool_false;

            for (int j = 0; j < isl_id_list_n_id(params); ++j)
            {
                isl_id* id = isl_id_list_get_id(params, j);

                int pos = isl_set_find_dim_by_id(domain, isl_dim_param, id);

                if (isl_bool_true == isl_qpolynomial_involves_dims(poly, isl_dim_param, pos, 1))
                {
                    is_parametric = isl_bool_true;
                }

                isl_id_free(id);
            }

            for (int j = 0; j < isl_id_list_n_id(II); ++j)
            {
                isl_id* id = isl_id_list_get_id(II, j);

                int pos = isl_set_find_dim_by_id(domain, isl_dim_param, id);

                if (isl_bool_true == isl_qpolynomial_involves_dims(poly, isl_dim_param, pos, 1))
                {
                    is_varied = isl_bool_true;
                }

                isl_id_free(id);
            }
            
            if (is_varied || is_parametric)
            {
                isl_printer* printer = isl_printer_to_str(ctx);
                printer = isl_printer_set_output_format(printer, ISL_FORMAT_C);
                printer = isl_printer_print_qpolynomial(printer, poly);
                char* str = isl_printer_get_str(printer);
                
                if (is_varied)
                {
                    is_varied = isl_bool_false;
                    for (int j = 0; j < isl_id_list_n_id(II); ++j)
                    {
                        isl_id* id = isl_id_list_get_id(II, j);
                        if (strstr(str, isl_id_get_name(id)))
                        {
                            is_varied = isl_bool_true;
                            isl_id_free(id);
                            break;
                        }
                        isl_id_free(id);
                    }
                }
                if (is_parametric)
                {
                    is_parametric = isl_bool_false;
                    for (int j = 0; j < isl_id_list_n_id(params); ++j)
                    {
                        isl_id* id = isl_id_list_get_id(params, j);
                        if (strstr(str, isl_id_get_name(id)))
                        {
                            is_parametric = isl_bool_true;
                            isl_id_free(id);
                            break;
                        }
                        isl_id_free(id);
                    }
                }                
                
                free(str);
                isl_printer_free(printer);
            }

            data.card_poly = poly;
            if (is_varied && is_parametric)
            {
                data.category = tc_tile_statistics_category_parametric_varied;
            }
            else if (is_varied)
            {
                data.category = tc_tile_statistics_category_varied;
            }
            else if (is_parametric)
            {
                data.category = tc_tile_statistics_category_parametric;
            }
            else
            {
                data.category = tc_tile_statistics_category_fixed;
            }

            isl_set_foreach_point(ii_set_bounded_repr, &tc_get_tile_statistics_callback, &data);
        }

        isl_set_free(domain);
        isl_set_free(ii_set_bounded_repr);
        isl_set_free(tile_bounded_repr);
    }
    
    tc_qpolynomials_free(qpolys);
    
    qsort(stats->groups, stats->n_groups, sizeof(stats->groups[0]), &tc_tile_statistics_group_compare);

    isl_id_list_free(params);    
    isl_pw_qpolynomial_free(tile_card);    
    isl_set_free(tile_bounded);
    isl_set_free(ii_set_bounded);
    
    return stats;
}

void tc_compute_tile_statistics_trend(FILE* out, __isl_keep isl_set* tile, __isl_keep isl_set* ii_set, __isl_keep isl_id_list* II, __isl_keep isl_set_list* bounds_list, __isl_keep isl_union_set* LD, __isl_keep isl_union_map* S, struct tc_scop* scop, const std::map<std::string, std::vector<int> >& blocks)
{
    for (int i = 0; i < isl_set_list_n_set(bounds_list); ++i)
    {
        isl_set* bounds = isl_set_list_get_set(bounds_list, i);
        
        struct tc_tile_statistics* stats = tc_compute_tile_statistics(tile, ii_set, II, bounds, LD, S, scop, blocks);
        
        long n_tiles = 0;
        long n_fixed_tiles = 0;
        long n_fixed_tiles_statements = 0;
        long n_parametric_tiles = 0;
        long n_parametric_tiles_statements = 0;
        long n_varied_tiles = 0;
        long n_varied_tiles_statements = 0;
        long n_parametric_varied_tiles = 0;
        long n_parametric_varied_tiles_statements = 0;

        for (int i = 0; i < stats->n_groups; ++i)
        {
            struct tc_tile_statistics_group* group = stats->groups[i];

            if (tc_tile_statistics_category_parametric_varied == group->category)
            {
                n_parametric_varied_tiles = n_parametric_varied_tiles + group->n_occurrences;
                n_parametric_varied_tiles_statements = n_parametric_varied_tiles_statements + (group->n_occurrences * group->card);
            }
            else if (tc_tile_statistics_category_parametric == group->category)
            {
                n_parametric_tiles = n_parametric_tiles + group->n_occurrences;
                n_parametric_tiles_statements = n_parametric_tiles_statements + (group->n_occurrences * group->card);
            }
            else if (tc_tile_statistics_category_varied == group->category)
            {
                n_varied_tiles = n_varied_tiles + group->n_occurrences;
                n_varied_tiles_statements = n_varied_tiles_statements + (group->n_occurrences * group->card);
            }
            else if (tc_tile_statistics_category_fixed == group->category)
            {
                n_fixed_tiles = n_fixed_tiles + group->n_occurrences;
                n_fixed_tiles_statements = n_fixed_tiles_statements + (group->n_occurrences * group->card);
            }

            n_tiles = n_tiles + group->n_occurrences;
        }
                
        if (isl_union_set_dim(stats->LD, isl_dim_param) > 0)
        {
            isl_id_list* params = tc_get_union_set_params_names(stats->LD);

            for (int j = 0; j < isl_id_list_n_id(params); ++j)
            {
                isl_id* param = isl_id_list_get_id(params, j);

                int pos = isl_set_find_dim_by_id(stats->bounds, isl_dim_param, param);

                isl_val* val = isl_set_plain_get_val_if_fixed(stats->bounds, isl_dim_param, pos);

                long value = isl_val_get_num_si(val);

                isl_val_free(val);

                fprintf(out, "%s = %ld%s", isl_id_get_name(param), value, j == isl_id_list_n_id(params)-1 ? "" : ", ");

                isl_id_free(param);
            }
            fprintf(out, "\n");

            isl_id_list_free(params);
        }
        
        fprintf(out, "tile type             tile percentage      %% of subspace occupied by this tile type\n");
        fprintf(out, "1. fixed                 %12.8f                                  %12.8f\n", (100.0 * n_fixed_tiles) / n_tiles, (100.0 * n_fixed_tiles_statements) / stats->n_statement_instances);
        fprintf(out, "2. varied                %12.8f                                  %12.8f\n", (100.0 * n_varied_tiles) / n_tiles, (100.0 * n_varied_tiles_statements) / stats->n_statement_instances);
        fprintf(out, "3. parametric            %12.8f                                  %12.8f\n", (100.0 * n_parametric_tiles) / n_tiles, (100.0 * n_parametric_tiles_statements) / stats->n_statement_instances);
        fprintf(out, "4. parametric/varied     %12.8f                                  %12.8f\n", (100.0 * n_parametric_varied_tiles) / n_tiles, (100.0 * n_parametric_varied_tiles_statements) / stats->n_statement_instances);
        fprintf(out, "\n");        
        
        tc_tile_statistics_free(stats);
        
        isl_set_free(bounds);
    }
}

void tc_tile_statistics_print(FILE* out, struct tc_tile_statistics* stats)
{
    long n_tiles = 0;
    long n_fixed_tiles = 0;
    long n_fixed_tiles_statements = 0;
    long n_parametric_tiles = 0;
    long n_parametric_tiles_statements = 0;
    long n_varied_tiles = 0;
    long n_varied_tiles_statements = 0;
    long n_parametric_varied_tiles = 0;
    long n_parametric_varied_tiles_statements = 0;
    
    for (int i = 0; i < stats->n_groups; ++i)
    {
        struct tc_tile_statistics_group* group = stats->groups[i];
        
        if (tc_tile_statistics_category_parametric_varied == group->category)
        {
            n_parametric_varied_tiles = n_parametric_varied_tiles + group->n_occurrences;
            n_parametric_varied_tiles_statements = n_parametric_varied_tiles_statements + (group->n_occurrences * group->card);
        }
        else if (tc_tile_statistics_category_parametric == group->category)
        {
            n_parametric_tiles = n_parametric_tiles + group->n_occurrences;
            n_parametric_tiles_statements = n_parametric_tiles_statements + (group->n_occurrences * group->card);
        }
        else if (tc_tile_statistics_category_varied == group->category)
        {
            n_varied_tiles = n_varied_tiles + group->n_occurrences;
            n_varied_tiles_statements = n_varied_tiles_statements + (group->n_occurrences * group->card);
        }
        else if (tc_tile_statistics_category_fixed == group->category)
        {
            n_fixed_tiles = n_fixed_tiles + group->n_occurrences;
            n_fixed_tiles_statements = n_fixed_tiles_statements + (group->n_occurrences * group->card);
        }
        
        n_tiles = n_tiles + group->n_occurrences;
    }
    
    fprintf(out, "Total statement instances: %ld\n", stats->n_statement_instances);
    fprintf(out, "Total tiles: %ld\n\n", n_tiles);
    
    if (isl_union_set_dim(stats->LD, isl_dim_param) > 0)
    {
        fprintf(out, "Parameters values:\n");
    
        isl_id_list* params = tc_get_union_set_params_names(stats->LD);

        for (int i = 0; i < isl_id_list_n_id(params); ++i)
        {
            isl_id* param = isl_id_list_get_id(params, i);

            int pos = isl_set_find_dim_by_id(stats->bounds, isl_dim_param, param);

            isl_val* val = isl_set_plain_get_val_if_fixed(stats->bounds, isl_dim_param, pos);

            long value = isl_val_get_num_si(val);

            isl_val_free(val);

            fprintf(out, "%s = %ld\n", isl_id_get_name(param), value);

            isl_id_free(param);
        }    
        
        isl_id_list_free(params);
        
        fprintf(out, "\n");
    }
    
    for (int i = 0; i < stats->n_statements; ++i)
    {
        struct tc_tile_statistics_block* statement = stats->statements[i];
        
        fprintf(out, "Statement %s has %ld instances in space ", isl_id_get_name(statement->label), statement->card);
        
        int size = isl_vec_size(statement->dimensions);
        
        for (int j = 0; j < size; ++j)
        {
            isl_val* val = isl_vec_get_element_val(statement->dimensions, j);

            long dim = isl_val_get_num_si(val);

            fprintf(out, "%ld%s", dim, j == size - 1 ? "" : " x ");

            isl_val_free(val);
        }
        
        fprintf(out, ", tiles are of size ");
        const char* statement_label = isl_id_get_name(statement->label);
        int statement_depth = tc_get_statement_depth(statement_label, stats->S);
        std::vector<int> block;
        if ((*stats->blocks).count(statement_label) > 0)
        {
            block = (*stats->blocks).at(statement_label);
        }
        else
        {
            block = (*stats->blocks).at("__DEFAULT__");
        }
        
        for (int j = 0; j < statement_depth; ++j)
        {
            fprintf(out, "%d %s", block[j], j == statement_depth - 1 ? "" : "x ");
        }        
        
        fprintf(out, "\n");
    }
    fprintf(out, "\n");
    
    fprintf(out, "Fixed tiles: %ld (%.8f %%) with total of %ld statement instances (%.8f %%)\n", n_fixed_tiles, (100.0 * n_fixed_tiles) / n_tiles, n_fixed_tiles_statements, (100.0 * n_fixed_tiles_statements) / stats->n_statement_instances);
    fprintf(out, "Varied tiles: %ld (%.8f %%) with total of %ld statement instances (%.8f %%)\n", n_varied_tiles, (100.0 * n_varied_tiles) / n_tiles, n_varied_tiles_statements, (100.0 * n_varied_tiles_statements) / stats->n_statement_instances);
    fprintf(out, "Parametric tiles: %ld (%.8f %%) with total of %ld statement instances (%.8f %%)\n", n_parametric_tiles, (100.0 * n_parametric_tiles) / n_tiles, n_parametric_tiles_statements, (100.0 * n_parametric_tiles_statements) / stats->n_statement_instances);    
    fprintf(out, "Parametric/varied tiles: %ld (%.8f %%) with total of %ld statement instances (%.8f %%)\n\n", n_parametric_varied_tiles, (100.0 * n_parametric_varied_tiles) / n_tiles, n_parametric_varied_tiles_statements, (100.0 * n_parametric_varied_tiles_statements) / stats->n_statement_instances);
        
    fprintf(out, "--------------------------------------------------------\n\n");
    
    for (int i = 0; i < stats->n_groups; ++i)
    {
        struct tc_tile_statistics_group* group = stats->groups[i];
        
        const char* category = "";
        switch (group->category)
        {
            case tc_tile_statistics_category_parametric: category = "parametric"; break;
            case tc_tile_statistics_category_varied: category = "varied"; break;
            case tc_tile_statistics_category_parametric_varied: category = "parametric/varied"; break;
            default:
            case tc_tile_statistics_category_fixed: category = "fixed"; break;
        }
        
        int plural = (group->n_occurrences > 1);
        const char* verb = (plural ? "" : "s");
        const char* noun = (plural ? "s" : "");
        fprintf(out, "%ld %s tile%s (%.8f %% of all tiles)%s including %ld statement instances (%.8f %% of all statement instances)\nTile%s contain%s %.8f %% of all statement instances\n", group->n_occurrences, category, noun, (100.0 * group->n_occurrences) / n_tiles, plural ? " each" : "", group->card, (100.0 * group->card) / stats->n_statement_instances, noun, verb, (100.0 * group->card * group->n_occurrences) / stats->n_statement_instances);
        
        int max_depth = 0;
        for (int j = 0; j < group->n_blocks; ++j)
        {
            struct tc_tile_statistics_block* block = group->blocks[j];
            int statement_depth = tc_get_statement_depth(isl_id_get_name(block->label), stats->S);
            max_depth = (statement_depth > max_depth ? statement_depth : max_depth);
        }
        
        for (int j = 0; j < group->n_blocks; ++j)
        {
            struct tc_tile_statistics_block* block = group->blocks[j];
            
            fprintf(out, "\t%ld statement instances of %s in a tile of size ", block->card, isl_id_get_name(block->label));
            
            int size = isl_vec_size(block->dimensions);
            
            long hblock_volume = 0;
            
            for (int k = 0; k < size; ++k)
            {
                isl_val* val = isl_vec_get_element_val(block->dimensions, k);

                long dim = isl_val_get_num_si(val);

                if (0 == hblock_volume)
                {
                    hblock_volume = dim;
                }
                else
                {
                    hblock_volume *= dim;
                }

                fprintf(out, "%ld%s", dim, k == size - 1 ? "" : " x ");

                isl_val_free(val);
            }
            
            for (int k = size; k < max_depth; ++k)
            {
                fprintf(out, " x 1");
            }

            fprintf(out, " (tile coverage with statement instances = %.8f %%)\n", 0 == hblock_volume ? 0 : (100.0 * block->card) / hblock_volume);
        }
        
        fprintf(out, "\n");
        
        ///*
        isl_ctx* ctx = isl_set_get_ctx(group->tile_sample);
                
        isl_ast_build* ast_build = isl_ast_build_from_context(isl_set_copy(stats->scop->pet->context));

        isl_union_map* S_prim = isl_union_map_intersect_range(isl_union_map_copy(stats->S), isl_union_set_from_set(isl_set_copy(group->tile_sample)));

        isl_ast_node* ast_tile = isl_ast_build_ast_from_schedule(ast_build, S_prim);

        isl_printer* printer = isl_printer_to_str(ctx);    

        printer = isl_printer_set_output_format(printer, ISL_FORMAT_C);

        isl_ast_print_options* ast_options = isl_ast_print_options_alloc(ctx);
        
        printer = isl_ast_node_print(ast_tile, printer, ast_options);
        
        char* code = isl_printer_get_str(printer);
        fprintf(out, "%s\n", code);
        free(code);
        
        printer = isl_printer_flush(printer);
        
        printer = isl_printer_print_qpolynomial(printer, group->card_qpolynomial);
        
        char* expression = isl_printer_get_str(printer);        
        fprintf(out, "The number of statement instances = %s\n\n", expression);        
        free(expression);
        
        isl_printer_free(printer);
        isl_ast_build_free(ast_build);
        isl_ast_node_free(ast_tile);
        //*/ 
        
        fprintf(out, "--------------------------------------------------------\n\n");
    }
}

void tc_tile_statistics_block_free(struct tc_tile_statistics_block* block)
{
    if (NULL == block)
        return;
    
    isl_id_free(block->label);
    isl_vec_free(block->dimensions);

    free(block);
}

void tc_tile_statistics_group_free(struct tc_tile_statistics_group* group)
{
    if (NULL == group)
        return;
    
    for (int i = 0; i < group->n_blocks; ++i)
    {
        tc_tile_statistics_block_free(group->blocks[i]);
    }
    free(group->blocks);
    
    isl_set_free(group->tile_sample);
    isl_qpolynomial_free(group->card_qpolynomial);
    
    free(group);
}

void tc_tile_statistics_free(struct tc_tile_statistics* stats)
{
    if (NULL == stats)
        return;
    
    for (int i = 0; i < stats->n_groups; ++i)
    {
        tc_tile_statistics_group_free(stats->groups[i]);
    }
    free(stats->groups);
    
    for (int i = 0; i < stats->n_statements; ++i)
    {
        tc_tile_statistics_block_free(stats->statements[i]);
    }
    free(stats->statements);
    
    isl_union_set_free(stats->LD);
    isl_union_map_free(stats->S);
    
    free(stats);
}

/*
void tc_compute_tile_statistics_2(__isl_keep isl_set* tile, __isl_keep isl_set* ii_set, __isl_keep isl_id_list* II, __isl_keep isl_set* bounds, __isl_keep isl_union_set* LD, __isl_keep isl_union_map* S, struct tc_scop* scop)
{
    isl_pw_qpolynomial* card = isl_set_card(isl_set_copy(tile));
    
    isl_id_list* params = tc_get_union_set_params_names(LD);    
    
    struct tc_qpolynomials* qpolys = tc_collect_qpolynomials(card);
    
    for (int i = 0; i < isl_set_list_n_set(qpolys->domains); ++i)
    {
        isl_qpolynomial* poly = qpolys->polys[i];
        isl_set* domain = isl_set_list_get_set(qpolys->domains, i);
        
        isl_bool is_parametric = isl_bool_false;
        isl_bool is_varied = isl_bool_false;

        for (int j = 0; j < isl_id_list_n_id(params); ++j)
        {
            isl_id* id = isl_id_list_get_id(params, j);

            int pos = isl_set_find_dim_by_id(domain, isl_dim_param, id);

            if (isl_qpolynomial_involves_dims(poly, isl_dim_param, pos, 1))
            {
                is_parametric = isl_bool_true;
            }

            isl_id_free(id);
        }

        for (int j = 0; j < isl_id_list_n_id(II); ++j)
        {
            isl_id* id = isl_id_list_get_id(II, j);

            int pos = isl_set_find_dim_by_id(domain, isl_dim_param, id);

            if (isl_qpolynomial_involves_dims(poly, isl_dim_param, pos, 1))
            {
                is_varied = isl_bool_true;
            }

            isl_id_free(id);
        }
        
        const char* category = "";
        if (is_parametric && is_varied) category = "parametric/varied";
        else if (is_parametric) category = "parametric";
        else if (is_varied) category = "varied";
        else category = "fixed";
        
        fprintf(out, stderr, "%s\n", category);
        
        fprintf(out, stderr, "card\n");
        isl_qpolynomial_dump(poly);

        isl_set* domain_set = tc_lift_up_set_params(isl_set_copy(domain), II);

        isl_set* S_i = isl_set_intersect(isl_set_copy(ii_set), isl_set_copy(domain_set));
        fprintf(out, stderr, "Si\n");
        isl_set_dump(S_i);
        
        isl_pw_qpolynomial* S_i_card = isl_set_card(isl_set_copy(S_i));
        fprintf(out, stderr, "card(Si)=\n");
        isl_pw_qpolynomial_dump(S_i_card);
        
        isl_set* S_i_bounded = isl_set_intersect_params(isl_set_copy(S_i), isl_set_copy(bounds));
        
        isl_pw_qpolynomial* S_i_bounded_card = isl_set_card(isl_set_copy(S_i_bounded));
        fprintf(out, stderr, "card(Si_bounded)=\n");
        isl_pw_qpolynomial_dump(S_i_bounded_card);

        //isl_pw_qpolynomial* S_i_card = isl_set_card(isl_set_copy(S_i));
        //fprintf(out, stderr, "card(Si)=\n");
        //isl_pw_qpolynomial_dump(S_i_card);

        isl_set* T_i = isl_set_intersect_params(isl_set_copy(tile), isl_set_copy(domain));
        fprintf(out, stderr, "Ti\n");
        isl_set_dump(T_i);
        
        fprintf(out, stderr, "---------------------------------------\n");

        isl_set_free(domain);
        isl_set_free(domain_set);
        isl_qpolynomial_free(poly);
        isl_set_free(T_i);
        isl_set_free(S_i);
        isl_set_free(S_i_bounded);
        isl_pw_qpolynomial_free(S_i_card);
        isl_pw_qpolynomial_free(S_i_bounded_card);
    }
    
    isl_pw_qpolynomial_free(card);
}
*/
