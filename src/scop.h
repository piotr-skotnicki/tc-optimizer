#ifndef TC_SCOP_H
#define TC_SCOP_H

#include <isl/ctx.h>
#include <isl/union_set.h>
#include <isl/union_map.h>
#include <isl/printer.h>

#include <pet.h>

struct tc_scop
{
    isl_ctx* ctx;
    
    isl_union_set* domain;
    
    isl_union_map* schedule;
    
    isl_union_map* relation;
    
    isl_union_map* reads;
    
    isl_union_map* writes;
    
    struct pet_scop* pet;
};

struct tc_stmt
{
    struct pet_stmt* stmt;
        
    isl_id_to_ast_expr* id2expr;
};

struct tc_scop* tc_scop_extract(__isl_keep isl_ctx* ctx, const char* filename);

__isl_give isl_union_map* tc_dependence_analysis(struct pet_scop* scop);

__isl_give isl_printer* tc_print_statements_macros(struct tc_scop* scop, __isl_take isl_printer* printer);

__isl_give isl_ast_node* tc_scop_at_each_domain(__isl_take isl_ast_node* node, __isl_keep isl_ast_build* build, void* user);

__isl_give isl_printer* tc_scop_print_user(__isl_take isl_printer* printer, __isl_take isl_ast_print_options* options, __isl_keep isl_ast_node* node, void* user);

__isl_give isl_printer* tc_print_prologue(struct tc_scop* scop, struct tc_options* options, __isl_take isl_printer* printer);

__isl_give isl_printer* tc_print_epilogue(struct tc_scop* scop, struct tc_options* options, __isl_take isl_printer* printer);

void tc_scop_free(struct tc_scop* scop);

struct tc_stmt* tc_stmt_alloc();

void tc_stmt_free(void* user);

struct pet_stmt* tc_scop_get_pet_stmt(struct tc_scop* scop, const char* label);

#endif // TC_SCOP_H
