#ifndef TC_SCOP_H
#define TC_SCOP_H

#include <isl/ctx.h>

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

struct tc_scop* tc_scop_extract(__isl_keep isl_ctx* ctx, const char* filename);

__isl_give isl_union_map* tc_dependence_analysis(struct pet_scop* scop);

void tc_scop_free(struct tc_scop* scop);

struct pet_stmt* tc_scop_get_pet_stmt(struct tc_scop* scop, const char* label);

#endif // TC_SCOP_H
