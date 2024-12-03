#include "stencil_tiling.h"
#include "regular_tiling.h"
#include "merge_tiling.h"
#include "split_tiling.h"
#include "correction_tiling.h"
#include "correction_inv_tiling.h"
#include "mod_correction_tiling.h"
#include "options.h"
#include "transitive_closure.h"
#include "scop.h"
#include "debug.h"
#include "utility.h"
#include "scheduling.h"
#include "lex_scheduling.h"
#include "isl_scheduling.h"
#include "sfs_scheduling.h"
#include "free_scheduling.h"
#include "dynamic_free_scheduling.h"
#include "timer.h"

#include <isl/ctx.h>
#include <isl/options.h>
#include <isl/id.h>
#include <isl/isl_options_private.h>

#include <pet.h>

#include <stddef.h>
#include <stdio.h>

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

        struct isl_options* isloptions = isl_ctx_options(ctx);
        isloptions->closure = ISL_CLOSURE_ISL;

        pet_options_set_autodetect(ctx, 0);

        isl_options_set_ast_always_print_block(ctx, 1);

        struct tc_scop* scop = tc_scop_extract(ctx, file);

        if (NULL == scop)
        {
            tc_error("No SCoP was found in file `%s'.", file);
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
                tc_error("SCoP is empty.");
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

                case tc_transitive_closure_enum_omega_map:
                {
                    tc_transitive_closure = &tc_transitive_closure_adapter_omega_map;
                    tc_map_power = &tc_map_power_adapter_isl_map;
                }
                break;

                case tc_transitive_closure_enum_omega_union_map:
                {
                    tc_transitive_closure = &tc_transitive_closure_adapter_omega_union_map;
                    tc_map_power = &tc_map_power_adapter_isl_map;
                }
                break;

                case tc_transitive_closure_enum_floyd_warshall:
                {
                    tc_transitive_closure = &tc_transitive_closure_adapter_floyd_warshall;
                    tc_map_power = &tc_map_power_adapter_isl_map;
                }
                break;

                case tc_transitive_closure_enum_iterative:
                {
                    tc_transitive_closure = &tc_transitive_closure_adapter_iterative;
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

            enum tc_scheduling_enum scheduling = tc_options_scheduling(options);

            switch (scheduling)
            {
                case tc_scheduling_enum_lex:
                {
                    tc_scheduling = &tc_scheduling_lex;
                }
                break;

                case tc_scheduling_enum_sfs_tile:
                {
                    tc_scheduling = &tc_scheduling_sfs_tiles;
                }
                break;

                case tc_scheduling_enum_sfs_single:
                {
                    tc_scheduling = &tc_scheduling_sfs_single;
                }
                break;

                case tc_scheduling_enum_sfs_multiple:
                {
                    tc_scheduling = &tc_scheduling_sfs_multiple;
                }
                break;

                case tc_scheduling_enum_free_rk:
                {
                    tc_scheduling = &tc_scheduling_free_schedule_rk;
                }
                break;

                case tc_scheduling_enum_free_karl:
                {
                    tc_scheduling = &tc_scheduling_free_schedule_karl;
                }
                break;

                case tc_scheduling_enum_free_dynamic:
                {
                    tc_scheduling = &tc_scheduling_dynamic_free_schedule;
                }
                break;

                case tc_scheduling_enum_isl:
                {
                    tc_scheduling = &tc_scheduling_adapter_isl;
                }
                break;

                case tc_scheduling_enum_isl_wavefronting:
                {
                    tc_scheduling = &tc_scheduling_adapter_isl_wavefronting;
                }
                break;

                case tc_scheduling_enum_feautrier:
                {
                    tc_scheduling = &tc_scheduling_adapter_feautrier;
                }
                break;
            }

            enum tc_algorithm_enum algorithm = tc_options_algorithm(options);

            struct tc_timer* timer = tc_timer_start();

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

                case tc_algorithm_enum_correction_inv_tiling:
                {
                    tc_algorithm_correction_inv_tiling(scop, options);
                }
                break;

                case tc_algorithm_enum_merge_tiling:
                {
                    tc_algorithm_merge_tiling(scop, options);
                }
                break;

                case tc_algorithm_enum_split_tiling:
                {
                    tc_algorithm_split_tiling(scop, options);
                }
                break;

                case tc_algorithm_enum_mod_correction_tiling:
                {
                    tc_algorithm_mod_correction_tiling(scop, options);
                }
                break;
                
                default:
                {
                    tc_error("Invalid algorithm");
                }
                break;
            }

            // if (tc_options_is_set(options, NULL, "--time"))
            {
                tc_debug("Total calculations time: %ld ms", tc_timer_stop(timer));
            }
        }

        tc_scop_free(scop);

        isl_ctx_free(ctx);
    }

    fclose(options->output);

    tc_options_free(options);
}

