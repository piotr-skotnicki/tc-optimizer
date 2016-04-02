#include "transitive_closure.h"
#include "utility.h"

#include <isl/set.h>
#include <isl/map.h>
#include <isl/union_map.h>

#include <stddef.h>

#include <string>
#include <vector>
#include <map>

__isl_give isl_map* tc_transitive_closure(__isl_take isl_map* R, int max)
{
    isl_ctx* ctx = isl_map_get_ctx(R);
    
    isl_id_list* I = tc_ids_sequence(ctx, "i", isl_map_n_in(R));
    isl_id_list* J = tc_ids_sequence(ctx, "j", isl_map_n_out(R));
        
    isl_map* RP = tc_make_map(ctx, NULL, I, J, tc_tuples_op(I, J, "<=").c_str());
    
    isl_map* R1 = isl_map_intersect(isl_map_copy(R), RP);
    isl_map* R2 = isl_map_subtract(isl_map_copy(R), isl_map_copy(R1));

    R1 = isl_map_project_out(R1, isl_dim_param, 0, isl_map_n_param(R1));
    R2 = isl_map_project_out(R2, isl_dim_param, 0, isl_map_n_param(R2));
    
    int exact;
    isl_map* R1_plus = isl_map_transitive_closure(isl_map_copy(R1), &exact);
    isl_map* R2_plus = isl_map_transitive_closure(isl_map_copy(R2), &exact);
    
    isl_map* R_plus_old = isl_map_union(isl_map_union(R1, R2), isl_map_union(R1_plus, R2_plus));
        
    int k = 1;
    while (1)
    {
        isl_map* R_plus_new = isl_map_union(isl_map_copy(R_plus_old), isl_map_apply_range(isl_map_copy(R_plus_old), isl_map_copy(R_plus_old)));
        
        isl_map* delta = isl_map_subtract(isl_map_copy(R_plus_new), R_plus_old);
        
        R_plus_old = R_plus_new;
        
        isl_bool empty = isl_map_is_empty(delta);
        isl_map_free(delta);
        
        if (empty || k == max)
        {
            break;
        }
        k = k + 1;
    }
            
    isl_set* IS = isl_set_union(isl_map_domain(isl_map_copy(R)), isl_map_range(isl_map_copy(R)));
    
    isl_map* R_IS = tc_make_map(ctx, NULL, I, J, tc_tuples_gt(J, I).c_str());
        
    R_IS = isl_map_intersect_domain(R_IS, isl_set_copy(IS));
    R_IS = isl_map_intersect_range(R_IS, IS);
    
    isl_map* R_plus = isl_map_intersect(R_plus_old, R_IS);
    R_plus = isl_map_coalesce(R_plus);
    
    isl_map_free(R);
    
    isl_id_list_free(I);
    isl_id_list_free(J);
    
    return R_plus;
}

__isl_give isl_union_map* tc_floyd_warshall_transitive_closure(__isl_take isl_union_map* R)
{
    isl_map* tcRR = NULL;
    isl_map* tcPR = NULL;
    isl_map* tcRQ = NULL;
        
    int size = 0;    
    
    std::map<std::string, int> id2index;
    //std::map<int, std::string> index2id;
    
    isl_map_list* maps = tc_collect_maps(R);
    
    for (int i = 0; i < isl_map_list_n_map(maps); ++i)
    {
        isl_map* map = isl_map_list_get_map(maps, i);
        
        const char* in_name = isl_map_get_tuple_name(map, isl_dim_in);
        const char* out_name = isl_map_get_tuple_name(map, isl_dim_out);
        
        if (id2index.count(in_name) == 0)
        {
            id2index[in_name] = size;
            //index2id[size] = in_name;            
            ++size;
        }
        
        if (id2index.count(out_name) == 0)
        {
            id2index[out_name] = size;
            //index2id[size] = out_name;
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
        
        relations[id2index[in_name]][id2index[out_name]] = map;
    }
    
    isl_map_list_free(maps);

    for (int r = 0; r < size; ++r)
    {
	if (relations[r][r] != NULL)
        {	    
            int exact;
	    tcRR = isl_map_transitive_closure(isl_map_copy(relations[r][r]), &exact);
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
                //map = isl_map_set_tuple_name(map, isl_dim_in, index2id[i].c_str());
                //map = isl_map_set_tuple_name(map, isl_dim_out, index2id[j].c_str());

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
    
    return R_plus;
}
