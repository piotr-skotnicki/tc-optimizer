#include "scop.h"
#include "utility.h"
#include "options.h"
#include "debug.h"

#include <isl/ctx.h>
#include <isl/space.h>
#include <isl/set.h>
#include <isl/map.h>
#include <isl/union_set.h>
#include <isl/union_map.h>
#include <isl/schedule.h>
#include <isl/flow.h>

#include <pet.h>

#include <stddef.h>
#include <string.h>
#include <stdlib.h>

struct tc_scop* tc_scop_extract(__isl_keep isl_ctx* ctx, const char* filename)
{    
    struct pet_scop* pet = pet_scop_extract_from_C_source(ctx, filename, NULL);

    if (NULL == pet)
        return NULL;

    struct tc_scop* scop = (struct tc_scop*)malloc(sizeof(struct tc_scop));
            
    scop->domain = isl_union_set_empty(isl_set_get_space(pet->context));
            
    scop->schedule = isl_schedule_get_map(pet->schedule);
    
    scop->relation = tc_dependence_analysis(pet);
    
    scop->reads = pet_scop_get_may_reads(pet);
    
    scop->writes = pet_scop_get_may_writes(pet);
    
    scop->ctx = isl_union_set_get_ctx(scop->domain);
    
    scop->pet = pet;

    for (int i = 0; i < pet->n_stmt; ++i) 
    {
        struct pet_stmt* stmt = pet->stmts[i];

        if (!pet_stmt_is_kill(stmt))
        {
            isl_set* statement_domain = isl_set_copy(stmt->domain);
            
            if (pet->stmts[i]->n_arg > 0)
            {
                statement_domain = isl_map_domain(isl_set_unwrap(statement_domain));
            }
            
            scop->domain = isl_union_set_add_set(scop->domain, statement_domain);
        }
        else
        {
            const char* statement_label = isl_set_get_tuple_name(stmt->domain);
            
            scop->schedule = tc_remove_map_with_tuple(scop->schedule, statement_label);
            
            scop->relation = tc_remove_map_with_tuple(scop->relation, statement_label);
        }
    }
    
    scop->schedule = tc_simplify_schedule(scop->schedule);
    
    return scop;
}

__isl_give isl_union_map* tc_dependence_analysis(struct pet_scop* pet)
{        
    isl_space* space = isl_set_get_space(pet->context);
    
    isl_union_map* empty = isl_union_map_empty(space);
    
    isl_union_map* schedule = isl_schedule_get_map(pet->schedule);
    isl_union_map* reads = pet_scop_get_may_reads(pet);
    isl_union_map* writes = pet_scop_get_may_writes(pet);
            
    isl_union_map* dep_raw;
    isl_union_map* dep_war;
    isl_union_map* dep_waw;
    
    /*
    isl_union_map_compute_flow(isl_union_map_copy(reads),
            isl_union_map_copy(writes),
            isl_union_map_copy(empty),
            isl_union_map_copy(schedule),
            &dep_raw, NULL, NULL, NULL);
            
    isl_union_map_compute_flow(isl_union_map_copy(writes),
            isl_union_map_copy(writes),
            isl_union_map_copy(reads),
            isl_union_map_copy(schedule),
            &dep_waw, &dep_war, NULL, NULL);
    */
    
    isl_union_map_compute_flow(isl_union_map_copy(reads),
            isl_union_map_copy(empty),
            isl_union_map_copy(writes),
            isl_union_map_copy(schedule),
            NULL, &dep_raw, NULL, NULL);
    
    isl_union_map_compute_flow(isl_union_map_copy(writes),
            isl_union_map_copy(empty),
            isl_union_map_copy(reads),
            isl_union_map_copy(schedule),
            NULL, &dep_war, NULL, NULL);
    
    isl_union_map_compute_flow(isl_union_map_copy(writes),
            isl_union_map_copy(empty),
            isl_union_map_copy(writes),
            isl_union_map_copy(schedule),
            NULL, &dep_waw, NULL, NULL);
    
    dep_raw = isl_union_map_coalesce(dep_raw);
    dep_war = isl_union_map_coalesce(dep_war);
    dep_waw = isl_union_map_coalesce(dep_waw);
    
    isl_union_map* relation = dep_raw;
    relation = isl_union_map_union(relation, dep_war);
    relation = isl_union_map_union(relation, dep_waw);
    
    relation = isl_union_map_coalesce(relation);

    isl_union_map_free(empty);
    isl_union_map_free(reads);
    isl_union_map_free(writes);
    isl_union_map_free(schedule);
        
    relation = tc_unwrap_range(relation);
    
    return relation;
}

__isl_give isl_union_map* tc_scop_data_to_cache_lines(struct tc_scop* scop, struct tc_options* options, __isl_keep isl_set* bounds)
{
    isl_ctx* ctx = isl_set_get_ctx(bounds);
    
    isl_set* defines = tc_options_get_defines(options, ctx);
        
    isl_set* params_bounds = isl_set_intersect_params(isl_set_copy(bounds), defines);
    
    isl_union_map* access_to_address = NULL;
        
    for (int i = 0; i < scop->pet->n_array; ++i)
    {
        struct pet_array* array = scop->pet->arrays[i];
        
        isl_set* extent = isl_set_copy(array->extent);
        
        extent = tc_set_fix_params_bounds(extent, isl_set_copy(params_bounds));
                    
        isl_map* address = isl_map_from_domain(isl_set_copy(extent));

        address = isl_map_add_dims(address, isl_dim_out, 2);

        address = isl_map_fix_si(address, isl_dim_out, 0, i); // [x = i, offset]
        
        // N[] = [x, 0]
        // A[i,j,k] = [x, i * j_size * k_size * data_size + j * k_size * data_size + k * data_size]        
        
        isl_constraint* constraint = isl_constraint_alloc_equality(isl_local_space_from_space(isl_map_get_space(address)));

        constraint = isl_constraint_set_coefficient_si(constraint, isl_dim_out, 1, -1);
        
        const int n_dim_extent = isl_set_n_dim(extent);

        for (int j = 0; j < n_dim_extent; ++j)
        {
            int multiplier = 1;

            for (int k = j + 1; k < n_dim_extent; ++k)
            {
                multiplier *= tc_lexmax_set_pos_value(extent, k) - tc_lexmin_set_pos_value(extent, k) + 1;
            }

            multiplier *= array->element_size;

            constraint = isl_constraint_set_coefficient_si(constraint, isl_dim_in, j, multiplier);
        }

        address = isl_map_add_constraint(address, constraint);

        if (NULL == access_to_address)
        {
            access_to_address = isl_union_map_from_map(address);
        }
        else
        {
            access_to_address = isl_union_map_add_map(access_to_address, address);
        }
        
        isl_set_free(extent);
    }

    char cache_line_str[256];
    
    snprintf(cache_line_str, sizeof(cache_line_str), "{ [x, i] -> [x, j] : j = floor(i/%d) }", tc_options_cache_line(options));
    
    isl_union_map* address_to_cache_line = isl_union_map_read_from_str(ctx, cache_line_str);
    
    isl_union_map* C = isl_union_map_apply_range(access_to_address, address_to_cache_line);
    
    isl_set_free(params_bounds);
    
    return C;
}

struct pet_stmt* tc_scop_get_pet_stmt(struct tc_scop* scop, const char* label)
{
    struct pet_stmt* value = NULL;
    
    struct pet_scop* pet = scop->pet;
    
    for (int i = 0; i < pet->n_stmt; ++i) 
    {
        struct pet_stmt* stmt = pet->stmts[i];
    
        const char* stmt_label = isl_set_get_tuple_name(stmt->domain);
        
        if (NULL != stmt_label && 0 == strcmp(label, stmt_label))
        {
            value = stmt;
            break;
        }
    }
    
    return value;
}

void tc_scop_free(struct tc_scop* scop)
{
    if (NULL == scop)
        return;
    
    pet_scop_free(scop->pet);

    isl_union_set_free(scop->domain);
    isl_union_map_free(scop->schedule);
    isl_union_map_free(scop->relation);
    isl_union_map_free(scop->reads);
    isl_union_map_free(scop->writes);

    free(scop);
}
