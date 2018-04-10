#include "floyd_warshall_transitive_closure.h"
#include "iterative_transitive_closure.h"
#include "utility.h"
#include "debug.h"

#include <isl/set.h>
#include <isl/map.h>
#include <isl/union_map.h>

#include <stddef.h>

#include <string>
#include <vector>
#include <map>

__isl_give isl_union_map* tc_floyd_warshall_transitive_closure(__isl_take isl_union_map* R, int* exact)
{
    isl_map* tcRR = NULL;
    isl_map* tcPR = NULL;
    isl_map* tcRQ = NULL;
        
    int size = 0;    
    
    std::map<std::string, int> id2index;
    std::map<int, std::string> index2id;
    
    isl_map_list* maps = tc_collect_maps(R);
    
    for (int i = 0; i < isl_map_list_n_map(maps); ++i)
    {
        isl_map* map = isl_map_list_get_map(maps, i);
        
        const char* in_name = isl_map_get_tuple_name(map, isl_dim_in);
        const char* out_name = isl_map_get_tuple_name(map, isl_dim_out);
        
        if (id2index.count(in_name) == 0)
        {
            id2index[in_name] = size;
            index2id[size] = in_name;            
            ++size;
        }
        
        if (id2index.count(out_name) == 0)
        {
            id2index[out_name] = size;
            index2id[size] = out_name;
            ++size;
        }
        
        isl_map_free(map);
    }
    
    isl_map*** relations = (isl_map***)calloc(size, sizeof(isl_map**));    
    for (int i = 0; i < size; ++i)
    {
        relations[i] = (isl_map**)calloc(size, sizeof(isl_map*));
        for (int j = 0; j < size; ++j)
        {
            relations[i][j] = NULL;
        }
    }
    
    for (int i = 0; i < isl_map_list_n_map(maps); ++i)
    {
        isl_map* map = isl_map_list_get_map(maps, i);
        
        const char* in_name = isl_map_get_tuple_name(map, isl_dim_in);
        const char* out_name = isl_map_get_tuple_name(map, isl_dim_out);

        map = isl_map_reset_tuple_id(map, isl_dim_in);
        map = isl_map_reset_tuple_id(map, isl_dim_out);
        
        relations[id2index[in_name]][id2index[out_name]] = map;
    }
    
    isl_map_list_free(maps);

    for (int r = 0; r < size; ++r)
    {
        if (relations[r][r] != NULL)
        {        
            //int exact;            
            
            tc_debug_map(relations[r][r], "RR");
            //tcRR = isl_map_transitive_closure(isl_map_copy(relations[r][r]), exact);
            tcRR = tc_iterative_transitive_closure(isl_map_copy(relations[r][r]), 10, 100, exact);
            tc_debug_map(tcRR, "RR+");

            isl_space* space = isl_map_get_space(tcRR);
            tcRR = isl_map_union(tcRR, isl_map_identity(space));
            //tcRR = isl_map_compute_divs(tcRR);
            tcRR = isl_map_coalesce(tcRR);
            //tcRR = isl_map_from_basic_map(isl_map_simple_hull(tcRR));
        }

        for (int p = 0; p < size; ++p)
        {
            if (relations[p][r] == NULL)
            {
                continue;
            }
            
            if (p != r)
            {
                if (relations[r][r] != NULL)
                {
                    //printf("\nCompose %d %d -> %d %d\n", r, r, p, r);
                    
                    tcPR = isl_map_apply_range(isl_map_copy(relations[p][r]), isl_map_copy(tcRR));
                    tcPR = isl_map_compute_divs(tcPR);
                    tcPR = isl_map_coalesce(tcPR);
                    //tcPR = isl_map_from_basic_map(isl_map_simple_hull(tcPR));
                }
                else
                {
                    tcPR = isl_map_copy(relations[p][r]);
                }
            }
            else
            {
                tcPR = isl_map_copy(tcRR);
            }

            for (int q = 0; q < size; ++q)
            {
                if (relations[r][q] == NULL)
                {
                    continue;
                }
                
                if (relations[p][q] == NULL)
                {
                    if (q != r)
                    {
                        //printf("\nCompose %d %d -> %d %d -> %d %d\n", r, q, r, r, p, r);

                        relations[p][q] = isl_map_apply_range(isl_map_copy(tcPR), isl_map_copy(relations[r][q]));
                        relations[p][q] = isl_map_compute_divs(relations[p][q]);
                        relations[p][q] = isl_map_coalesce(relations[p][q]);
                        //relations[p][q] = isl_map_from_basic_map(isl_map_simple_hull(relations[p][q]));
                    }
                    else
                    {
                        relations[p][q] = isl_map_copy(tcPR);
                    }
                }
                else
                {
                    if (q != r)
                    {
                        //printf("\nCompose %d %d -> %d %d -> %d %d\n", r, q, r, r, p, r);

                        tcRQ = isl_map_apply_range(isl_map_copy(tcPR), isl_map_copy(relations[r][q]));
                        tcRQ = isl_map_compute_divs(tcRQ);
                        tcRQ = isl_map_coalesce(tcRQ);
                        //tcRQ = isl_map_from_basic_map(isl_map_simple_hull(tcRQ));

                        relations[p][q] = isl_map_union(relations[p][q], tcRQ);

                        tcRQ = NULL;

                        //relations[p][q] = isl_map_union(relations[p][q], isl_map_apply_range(isl_map_copy(tcPR), isl_map_copy(relations[r][q])));
                        //relations[p][q] = isl_map_coalesce(relations[p][q]);

                        //isl_map_apply_range(isl_map_copy(tcPR), isl_map_copy(relations[r][q])));
                        //relations[p][q] = isl_map_coalesce(relations[p][q]);
                    }
                    else
                    {
                        relations[p][q] = isl_map_union(relations[p][q], isl_map_copy(tcPR));
                    }

                    relations[p][q] = isl_map_compute_divs(relations[p][q]);
                    relations[p][q] = isl_map_coalesce(relations[p][q]);
                    //relations[p][q] = isl_map_from_basic_map(isl_map_polyhedral_hull(relations[p][q]));                        
                }
            }
            isl_map_free(tcPR);
            tcPR = NULL;
        } 

        isl_map_free(tcRR); 
        tcRR = NULL;
    }
    
    isl_union_map* R_plus = NULL;
    
    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < size; ++j)
        {
            isl_map* map = relations[i][j];
            
            if (map != NULL)
            {
                map = isl_map_set_tuple_name(map, isl_dim_in, index2id[i].c_str());
                map = isl_map_set_tuple_name(map, isl_dim_out, index2id[j].c_str());

                if (R_plus == NULL)
                {
                    R_plus = isl_union_map_from_map(map);
                }
                else
                {
                    R_plus = isl_union_map_add_map(R_plus, map);
                }
            }
        }
        
        free(relations[i]);
    }
    
    free(relations);
    
    isl_union_map_free(R);
    
    return R_plus;
}
