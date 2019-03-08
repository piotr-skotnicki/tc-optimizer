#include "ast.h"
#include "scop.h"
#include "debug.h"

#include <isl/ctx.h>
#include <isl/id.h>
#include <isl/ast.h>

#include <pet.h>

#include <stdlib.h>
#include <stddef.h>
#include <string.h>

struct tc_ast_visitor_context* tc_ast_visitor_context_alloc(__isl_keep isl_ctx* ctx)
{
    struct tc_ast_visitor_context* context = (struct tc_ast_visitor_context*)malloc(sizeof(struct tc_ast_visitor_context));
    
    context->scop = NULL;
    context->parallel_iterators = NULL;
    
    return context;
}

void tc_ast_visitor_context_free(struct tc_ast_visitor_context* context)
{
    if (NULL == context)
        return;
    
    free(context);
}

struct tc_ast_for_annotation* tc_ast_for_annotation_alloc(__isl_keep isl_ctx* ctx)
{
    struct tc_ast_for_annotation* annotation = (struct tc_ast_for_annotation*)malloc(sizeof(struct tc_ast_for_annotation));
    
    annotation->nested_statements = isl_id_list_alloc(ctx, 1);
    annotation->is_parallel = isl_bool_false;
    
    return annotation;
}

void tc_ast_for_annotation_free(void* user)
{
    if (NULL == user)
        return;
    
    struct tc_ast_for_annotation* annotation = (struct tc_ast_for_annotation*)user;
    
    isl_id_list_free(annotation->nested_statements);
    
    free(annotation);
}

struct tc_ast_stmt_annotation* tc_ast_stmt_annotation_alloc()
{
    struct tc_ast_stmt_annotation* annotation = (struct tc_ast_stmt_annotation*)malloc(sizeof(struct tc_ast_stmt_annotation));
    
    annotation->stmt = NULL;    
    annotation->id2expr = NULL;
    
    return annotation;
}

void tc_ast_stmt_annotation_free(void* user)
{    
    if (NULL == user)
        return;
    
    struct tc_ast_stmt_annotation* annotation = (struct tc_ast_stmt_annotation*)user;
    
    isl_id_to_ast_expr_free(annotation->id2expr);
    
    free(annotation);
}

__isl_give isl_id* tc_ast_visitor_before_for(__isl_keep isl_ast_build* build, void* user)
{
    isl_ctx* ctx = isl_ast_build_get_ctx(build);
    
    struct tc_ast_visitor_context* context = (struct tc_ast_visitor_context*)user;
    
    struct tc_ast_for_annotation* annotation = tc_ast_for_annotation_alloc(ctx);
    
    isl_id* annotation_id = isl_id_alloc(ctx, NULL, annotation);
    
    annotation_id = isl_id_set_free_user(annotation_id, &tc_ast_for_annotation_free);
    
    context->annotations_stack->push_back(annotation);
        
    return annotation_id;
}

__isl_give isl_ast_node* tc_ast_visitor_after_for(__isl_take isl_ast_node* node, __isl_keep isl_ast_build* build, void* user)
{
    struct tc_ast_visitor_context* context = (struct tc_ast_visitor_context*)user;
    
    isl_ast_expr* expr = isl_ast_node_for_get_iterator(node);
    
    isl_id* annotation_id = isl_ast_node_get_annotation(node);
    
    struct tc_ast_for_annotation* annotation = (struct tc_ast_for_annotation*)isl_id_get_user(annotation_id);
    
    isl_id* iterator = isl_ast_expr_get_id(expr);
    
    if (context->parallel_iterators != NULL)
    {
        for (int i = 0; i < isl_id_list_n_id(context->parallel_iterators); ++i)
        {
            isl_id* parallel_iterator = isl_id_list_get_id(context->parallel_iterators, i);

            if (0 == strcmp(isl_id_get_name(iterator), isl_id_get_name(parallel_iterator)))
            {
                annotation->is_parallel = isl_bool_true;
            }

            isl_id_free(parallel_iterator);
        }
    }
    
    context->annotations_stack->pop_back();
    
    isl_id_free(annotation_id);
    isl_id_free(iterator);
    isl_ast_expr_free(expr);
    
    return node;
}

static __isl_give isl_multi_pw_aff* tc_ast_visitor_build_ast_exprs_index_callback(__isl_take isl_multi_pw_aff* index, __isl_keep isl_id* id, void* user)
{
    isl_pw_multi_aff* iterator_map = (isl_pw_multi_aff*)user;

    iterator_map = isl_pw_multi_aff_copy(iterator_map);
    
    return isl_multi_pw_aff_pullback_pw_multi_aff(index, iterator_map);
}

__isl_give isl_ast_node* tc_ast_visitor_at_each_domain(__isl_take isl_ast_node* node, __isl_keep isl_ast_build* build, void* user)
{
    struct tc_ast_visitor_context* context = (struct tc_ast_visitor_context*)user;

    struct tc_scop* scop = context->scop;

    struct tc_ast_stmt_annotation* annotation = tc_ast_stmt_annotation_alloc();
    
    isl_ctx* ctx = isl_ast_node_get_ctx(node);

    isl_ast_expr* expr = isl_ast_node_user_get_expr(node);

    isl_ast_expr* arg = isl_ast_expr_get_op_arg(expr, 0);

    isl_id* id = isl_ast_expr_get_id(arg);

    for (int i = 0; i < (*context->annotations_stack).size(); ++i)
    {
        isl_id_list* nested_statements = (*context->annotations_stack)[i]->nested_statements;

        isl_bool has_nested_statement = isl_bool_false;
        
        for (int j = 0; j < isl_id_list_n_id(nested_statements); ++j)
        {
            isl_id* nested_statement = isl_id_list_get_id(nested_statements, j);
            
            if (0 == strcmp(isl_id_get_name(id), isl_id_get_name(nested_statement)))
            {
                has_nested_statement = isl_bool_true;
            }

            isl_id_free(nested_statement);
        }

        if (has_nested_statement == isl_bool_false)
        {
            nested_statements = isl_id_list_add(nested_statements, isl_id_copy(id));

            (*context->annotations_stack)[i]->nested_statements = nested_statements;
        }
    }
    
    annotation->stmt = tc_scop_get_pet_stmt(scop, isl_id_get_name(id));

    isl_map* map = isl_map_from_union_map(isl_ast_build_get_schedule(build));
    
    map = isl_map_reverse(map);
    
    isl_pw_multi_aff* iterator_map = isl_pw_multi_aff_from_map(map);
    
    annotation->id2expr = pet_stmt_build_ast_exprs(annotation->stmt, build, &tc_ast_visitor_build_ast_exprs_index_callback, iterator_map, NULL, NULL);
    
    isl_pw_multi_aff_free(iterator_map);

    isl_id* annotation_id = isl_id_alloc(ctx, NULL, annotation);
    
    annotation_id = isl_id_set_free_user(annotation_id, &tc_ast_stmt_annotation_free);
    
    isl_id_free(id);
    isl_ast_expr_free(expr);
    isl_ast_expr_free(arg);
    
    node = isl_ast_node_set_annotation(node, annotation_id);
    
    return node;
}
