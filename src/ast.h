#ifndef TC_AST_H
#define TC_AST_H

#include <isl/ctx.h>
#include <isl/id.h>
#include <isl/ast.h>

#include <pet.h>

#include <vector>

struct tc_ast_stmt_annotation
{
    struct pet_stmt* stmt;
        
    isl_id_to_ast_expr* id2expr;
};

struct tc_ast_for_annotation
{
    isl_id_list* nested_statements;
    
    isl_bool is_parallel;
};

struct tc_ast_visitor_context
{
    struct tc_scop* scop;

    std::vector<struct tc_ast_for_annotation*>* annotations_stack;
    
    isl_id_list* parallel_iterators;
};

struct tc_ast_visitor_context* tc_ast_visitor_context_alloc(__isl_keep isl_ctx* ctx);

void tc_ast_visitor_context_free(struct tc_ast_visitor_context* context);

__isl_give isl_ast_node* tc_ast_visitor_at_each_domain(__isl_take isl_ast_node* node, __isl_keep isl_ast_build* build, void* user);

__isl_give isl_id* tc_ast_visitor_before_for(__isl_keep isl_ast_build* build, void* user);

__isl_give isl_ast_node* tc_ast_visitor_after_for(__isl_take isl_ast_node* node, __isl_keep isl_ast_build* build, void* user);

struct tc_ast_stmt_annotation* tc_ast_stmt_annotation_alloc();

void tc_ast_stmt_annotation_free(void* user);

struct tc_ast_for_annotation* tc_ast_for_annotation_alloc(__isl_keep isl_ctx* ctx);

void tc_ast_for_annotation_free(void* user);

#endif // TC_AST_H
