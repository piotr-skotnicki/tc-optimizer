#ifndef TC_TILE_STATISTICS_H
#define TC_TILE_STATISTICS_H

#include "scop.h"

#include <isl/id.h>
#include <isl/set.h>
#include <isl/union_set.h>
#include <isl/union_map.h>
#include <isl/vec.h>

#include <barvinok/isl.h>

#include <stdio.h>

#include <string>
#include <vector>
#include <map>

enum tc_tile_statistics_category_enum
{
    tc_tile_statistics_category_fixed,
    tc_tile_statistics_category_parametric,
    tc_tile_statistics_category_varied,
    tc_tile_statistics_category_parametric_varied,
    tc_tile_statistics_category_unknown
};

struct tc_tile_statistics_block
{
    isl_id* label;
    
    isl_vec* dimensions;
    
    long card;
};

struct tc_tile_statistics_group
{    
    struct tc_tile_statistics_block** blocks;
    
    int n_blocks;
    
    int n_blocks_capacity;
    
    int dimensionality;
    
    long card;
    
    long memory;
    
    long n_occurrences;
    
    isl_set* tile_sample;
    
    isl_qpolynomial* card_qpolynomial;
    
    enum tc_tile_statistics_category_enum category;
};

struct tc_tile_statistics
{    
    struct tc_scop* scop;
    
    struct tc_options* options;
    
    isl_union_set* LD;
    
    isl_union_map* S;
    
    isl_set* bounds;
    
    struct tc_tile_statistics_group** groups;
    
    int n_groups;
    
    int n_groups_capacity;
    
    struct tc_tile_statistics_block** statements;
    
    int n_statements;
    
    long n_statement_instances;
    
    long n_dimensionality_tiles[10];
    
    const std::map<std::string, std::vector<int> >* blocks;
};

struct tc_tile_statistics* tc_compute_tile_statistics(__isl_keep isl_set* tile
                                                    , __isl_keep isl_set* ii_set
                                                    , __isl_keep isl_id_list* II
                                                    , __isl_keep isl_set* bounds
                                                    , __isl_keep isl_union_set* LD
                                                    , __isl_keep isl_union_map* S
                                                    , __isl_keep isl_union_map* RA
                                                    , __isl_keep isl_union_map* WA
                                                    , struct tc_scop* scop
                                                    , struct tc_options* options
                                                    , const std::map<std::string, std::vector<int> >& blocks);

void tc_compute_tile_statistics_trend(FILE* out
                                    , __isl_keep isl_set* tile
                                    , __isl_keep isl_set* ii_set    
                                    , __isl_keep isl_id_list* II
                                    , __isl_keep isl_set_list* bounds_list
                                    , __isl_keep isl_union_set* LD
                                    , __isl_keep isl_union_map* S
                                    , __isl_keep isl_union_map* RA
                                    , __isl_keep isl_union_map* WA
                                    , struct tc_scop* scop
                                    , struct tc_options* options
                                    , const std::map<std::string, std::vector<int> >& blocks);

void tc_tile_statistics_print(FILE* out, struct tc_tile_statistics* stats);

void tc_tile_statistics_block_free(struct tc_tile_statistics_block* block);

void tc_tile_statistics_group_free(struct tc_tile_statistics_group* group);

void tc_tile_statistics_free(struct tc_tile_statistics* stats);

//void tc_compute_tile_statistics_2(__isl_keep isl_set* tile, __isl_keep isl_set* ii_set, __isl_keep isl_id_list* II, __isl_keep isl_set* bounds, __isl_keep isl_union_set* LD, __isl_keep isl_union_map* S, struct tc_scop* scop);

#endif // TC_TILE_STATISTICS_H
