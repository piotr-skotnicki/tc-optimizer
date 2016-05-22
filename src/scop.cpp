#include "scop.h"
#include "utility.h"
#include "options.h"

#include <isl/ctx.h>
#include <isl/space.h>
#include <isl/set.h>
#include <isl/map.h>
#include <isl/union_set.h>
#include <isl/union_map.h>
#include <isl/printer.h>
#include <isl/ast_build.h>
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

__isl_give isl_printer* tc_print_statements_macros(struct tc_scop* scop, __isl_take isl_printer* printer)
{
    struct pet_scop* pet = scop->pet;
    
    char buff[2048];
    
    for (int i = 0; i < pet->n_stmt; ++i) 
    {
        struct pet_stmt* stmt = pet->stmts[i];

        if (!pet_stmt_is_kill(stmt))
        {
            isl_set* statement_domain = isl_set_copy(stmt->domain);
            
            isl_ctx* ctx = isl_set_get_ctx(statement_domain);
            
            if (pet->stmts[i]->n_arg > 0)
            {
                statement_domain = isl_map_domain(isl_set_unwrap(statement_domain));
            }
            
            const char* label = isl_set_get_tuple_name(statement_domain);
            
            int n_set_dim = isl_set_dim(statement_domain, isl_dim_set);
            
            isl_id_list* iterators = isl_id_list_alloc(ctx, n_set_dim);
            
            for (int j = 0; j < n_set_dim; ++j)
            {            
                iterators = isl_id_list_add(iterators, isl_id_alloc(ctx, isl_set_get_dim_name(statement_domain, isl_dim_set, j), NULL));
            }
            
            sprintf(buff, "#define %s_I(%s) ", label, tc_comma(iterators).c_str());            
            
            printer = isl_printer_print_str(printer, buff);
            
            isl_ast_build* ast_stmt_build = isl_ast_build_from_context(isl_set_copy(statement_domain));
            
            isl_id_to_ast_expr* id2ast = pet_stmt_build_ast_exprs(stmt, ast_stmt_build, NULL, NULL, NULL, NULL);
            
            isl_ast_build_free(ast_stmt_build);
                        
            printer = pet_stmt_print_body(stmt, printer, id2ast);
            
            isl_id_to_ast_expr_free(id2ast);
            
            sprintf(buff, "#define %s(%s) %s_I(%s)\n", label, tc_comma(iterators).c_str(), label, tc_parens_comma(iterators).c_str());
            
            printer = isl_printer_print_str(printer, buff);
            
            isl_id_list_free(iterators);            
            isl_set_free(statement_domain);
        }
    }
    
    return printer;
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

static __isl_give isl_multi_pw_aff* tc_scop_build_ast_exprs_index_callback(__isl_take isl_multi_pw_aff* index, __isl_keep isl_id* id, void* user)
{
    isl_pw_multi_aff* iterator_map = (isl_pw_multi_aff*)user;

    iterator_map = isl_pw_multi_aff_copy(iterator_map);
    
    return isl_multi_pw_aff_pullback_pw_multi_aff(index, iterator_map);
}

__isl_give isl_ast_node* tc_scop_at_each_domain(__isl_take isl_ast_node* node, __isl_keep isl_ast_build* build, void* user)
{
    struct tc_scop* scop = (struct tc_scop*)user;

    struct tc_stmt* stmt = tc_stmt_alloc();
    
    isl_ctx* ctx = isl_ast_node_get_ctx(node);

    isl_ast_expr* expr = isl_ast_node_user_get_expr(node);

    isl_ast_expr* arg = isl_ast_expr_get_op_arg(expr, 0);

    isl_id* id = isl_ast_expr_get_id(arg);
    
    stmt->stmt = tc_scop_get_pet_stmt(scop, isl_id_get_name(id));

    isl_map* map = isl_map_from_union_map(isl_ast_build_get_schedule(build));
    
    map = isl_map_reverse(map);
    
    isl_pw_multi_aff* iterator_map = isl_pw_multi_aff_from_map(map);
    
    stmt->id2expr = pet_stmt_build_ast_exprs(stmt->stmt, build, &tc_scop_build_ast_exprs_index_callback, iterator_map, NULL, NULL);
    
    isl_pw_multi_aff_free(iterator_map);

    isl_id* annotation = isl_id_alloc(ctx, NULL, stmt);
    
    annotation = isl_id_set_free_user(annotation, &tc_stmt_free);
    
    isl_id_free(id);
    isl_ast_expr_free(expr);    
    isl_ast_expr_free(arg);
    
    node = isl_ast_node_set_annotation(node, annotation);
    
    return node;
}

struct tc_stmt* tc_stmt_alloc()
{
    struct tc_stmt* stmt = (struct tc_stmt*)malloc(sizeof(struct tc_stmt));
    
    stmt->stmt = NULL;    
    stmt->id2expr = NULL;
    
    return stmt;
}

void tc_stmt_free(void* user)
{    
    if (NULL == user)
        return;
    
    struct tc_stmt* stmt = (struct tc_stmt*)user;
    
    isl_id_to_ast_expr_free(stmt->id2expr);
    
    free(stmt);
}

__isl_give isl_printer* tc_scop_print_user(__isl_take isl_printer* printer, __isl_take isl_ast_print_options* options, __isl_keep isl_ast_node* node, void* user)
{
    isl_id* annotation = isl_ast_node_get_annotation(node);
    
    struct tc_stmt* stmt = (struct tc_stmt*)isl_id_get_user(annotation);
    
    printer = pet_stmt_print_body(stmt->stmt, printer, stmt->id2expr);

    isl_ast_print_options_free(options);
    isl_id_free(annotation);

    return printer;
}

__isl_give isl_printer* tc_print_prologue(struct tc_scop* scop, struct tc_options* options, __isl_take isl_printer* printer)
{
    char* command_line = tc_options_get_command_line(options);
    
    printer = isl_printer_print_str(printer, "/* ");
    printer = isl_printer_print_str(printer, command_line);
    printer = isl_printer_print_str(printer, " */\n");
    
    free(command_line);
    
    return printer;
}

__isl_give isl_printer* tc_print_epilogue(struct tc_scop* scop, struct tc_options* options, __isl_take isl_printer* printer)
{
    //printer = isl_printer_print_str(printer, "/**/\n");
    
    return printer;
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
