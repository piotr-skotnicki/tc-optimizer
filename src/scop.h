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
    
    struct pet_scop* pet;
};

struct tc_scop* tc_scop_extract(__isl_keep isl_ctx* ctx, const char* filename);

__isl_give isl_union_map* tc_dependence_analysis(struct pet_scop* scop);

__isl_give isl_printer* tc_print_statements_macros(struct tc_scop* scop, __isl_take isl_printer* printer, __isl_keep isl_ast_build* ast_build);

void tc_scop_free(struct tc_scop* scop);

#endif // TC_SCOP_H
