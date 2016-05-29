#include "stencil_tiling.h"
#include "regular_tiling.h"
#include "merge_tiling.h"
#include "correction_tiling.h"
#include "options.h"
#include "transitive_closure.h"
#include "scop.h"
#include "debug.h"
#include "utility.h"

#include <isl/ctx.h>
#include <isl/options.h>

#include <pet.h>

#include <stddef.h>
#include <stdio.h>

#include <sys/time.h>

int main(int argc, char* argv[])
{
    struct tc_options* options = tc_options_alloc(argc, argv);
    
    tc_options_check_spelling(options);

    options->output = stdout;
    
    if (tc_options_is_set(options, "-v", "--version"))
    {
        tc_options_credits();
    }
    else if (argc < 3 || tc_options_is_set(options, "-h", "--help"))
    {
        tc_options_help();
    }
    else
    {        
        tc_debug_flag = tc_options_is_verbose(options);
        
        const char* file = tc_options_source_file(options);

        struct pet_options* poptions = pet_options_new_with_defaults();

        isl_ctx* ctx = isl_ctx_alloc_with_options(&pet_options_args, poptions);    

        pet_options_set_autodetect(ctx, 0);
        
        if (tc_options_is_set(options, NULL, "--braces"))
        {
            isl_options_set_ast_always_print_block(ctx, 1);
        }

        struct tc_scop* scop = tc_scop_extract(ctx, file);
        
        if (NULL == scop)
        {
            tc_options_error("No SCoP was found in file `%s'.", file);
        }
        else
        {
            isl_set* defines = tc_options_get_defines(options, ctx);
        
            scop->domain = tc_union_set_fix_params_bounds(scop->domain, isl_set_copy(defines));

            scop->relation = tc_union_map_fix_params_bounds(scop->relation, isl_set_copy(defines));
            
            scop->reads = tc_union_map_fix_params_bounds(scop->reads, isl_set_copy(defines));
            
            scop->writes = tc_union_map_fix_params_bounds(scop->writes, isl_set_copy(defines));

            scop->schedule = tc_union_map_fix_params_bounds(scop->schedule, defines);
            
            if (isl_union_set_is_empty(scop->domain))
            {
                tc_options_error("SCoP is empty.");
            }
            
            enum tc_transitive_closure_enum transitive_closure = tc_options_transitive_closure(options);

            switch (transitive_closure)
            {
                case tc_transitive_closure_enum_isl_map:
                case tc_transitive_closure_enum_unknown:
                default:
                {
                    tc_transitive_closure = &tc_transitive_closure_adapter_isl_map;
                    tc_map_power = &tc_map_power_adapter_isl_map;
                }
                break;

                case tc_transitive_closure_enum_isl_union_map:
                {
                    tc_transitive_closure = &tc_transitive_closure_adapter_isl_union_map;
                    tc_map_power = &tc_map_power_adapter_isl_map;
                }
                break;

                case tc_transitive_closure_enum_floyd_warshall:
                {
                    tc_transitive_closure = &tc_transitive_closure_adapter_floyd_warshall;
                    tc_map_power = &tc_map_power_adapter_isl_map;
                }
                break;

                case tc_transitive_closure_enum_tarjan:
                {
                    tc_transitive_closure = &tc_transitive_closure_adapter_tarjan;
                    tc_map_power = &tc_map_power_adapter_isl_map;
                }
                break;
            }

            enum tc_algorithm_enum algorithm = tc_options_algorithm(options);
            
            struct timeval start_time, end_time;
            gettimeofday(&start_time, NULL);
            
            switch (algorithm)
            {
                case tc_algorithm_enum_stencil_tiling:
                {
                    tc_algorithm_stencil_tiling(scop, options);
                }
                break;

                case tc_algorithm_enum_regular_tiling:
                {
                    tc_algorithm_regular_tiling(scop, options);
                }
                break;

                case tc_algorithm_enum_correction_tiling:
                {
                    tc_algorithm_correction_tiling(scop, options);
                }
                break;

                case tc_algorithm_enum_merge_tiling:
                {
                    tc_algorithm_merge_tiling(scop, options);
                }
                break;
                
                default:
                {
                    tc_options_error("Invalid algorithm");
                }
                break;
            }
            
            gettimeofday(&end_time, NULL);
            
            if (tc_options_is_set(options, NULL, "--time"))
            {
                unsigned long elapsed = 1000 * (end_time.tv_sec - start_time.tv_sec) + (end_time.tv_usec - start_time.tv_usec) / 1000;
                
                fprintf(stderr, "Total calculations time: %ld ms\n", elapsed);
            }
        }
                
        tc_scop_free(scop);     

        isl_ctx_free(ctx);
    }
        
    fclose(options->output);
    
    tc_options_free(options);
}

