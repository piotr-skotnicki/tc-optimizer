#include <isl/ctx.h>
#include <isl/options.h>

#include <pet.h>

#include <stddef.h>

#include "stencil.h"
#include "free_schedule_tiling.h"
#include "static_tile_correction.h"
#include "options.h"

int main(int argc, char* argv[])
{
    if (argc < 3)
    {
        tc_options_help();
    }
    else
    {    
        const char* file = tc_options_source_file(argc, argv);

        enum tc_algorithm_enum algorithm = tc_options_algorithm(argc, argv);

        struct pet_options* options = pet_options_new_with_defaults();

        isl_ctx* ctx = isl_ctx_alloc_with_options(&pet_options_args, options);        

        pet_options_set_autodetect(ctx, 0);

        isl_options_set_coalesce_bounded_wrapping(ctx, 0);

        struct tc_scop* scop = tc_scop_extract(ctx, file);
        
        if (NULL == scop)
        {
            tc_options_error("Invalid input file `%s'", file);
        }
        else
        {
            switch (algorithm)
            {
                case tc_algorithm_enum_stencil_tiling:
                {
                    tc_algorithm_stencil(argc - 3, argv + 3, scop);
                }
                break;

                case tc_algorithm_enum_free_schedule_tiling:
                {
                    tc_algorithm_free_schedule_tiling(argc - 3, argv + 3, scop);
                }
                break;

                case tc_algorithm_enum_static_correction_tiling:
                {
                    tc_algorithm_static_tile_correction(argc - 3, argv + 3, scop);
                }
                break;
            }
        }
        
        tc_scop_free(scop);

        isl_ctx_free(ctx);
    }
}
