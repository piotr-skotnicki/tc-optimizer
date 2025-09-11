#include "select_transitive_closure.h"
#include "utility.h"
#include "debug.h"
#include "input_output.h"

#include "transitive_closure.h"

#include <isl/ctx.h>
#include <isl/map.h>
#include <isl/union_map.h>

#include <stdio.h>

__isl_give isl_map* tc_select_transitive_closure(__isl_take isl_map* R, __isl_keep isl_union_map* S, isl_bool* exact)
{
    const char* names[] = {
        "--isl-map-tc           ISL normalized map transitive closure",
        "--isl-union-map-tc     ISL union map transitive closure",
        "--iterative-tc         Iterative algorithm",
        "--floyd-warshall-tc    Floyd-Warshall algorithm",
        "--omega-map-tc         Omega normalized map transitive closure",
        "--omega-union-map-tc   Omega union map transitive closure",
        "--tarjan-tc            Tarjan algorithm for finite graphs",
    };

    isl_map* (*algorithms[])(__isl_take isl_map* R, __isl_keep isl_union_map* S, isl_bool* exact) = {
        &tc_transitive_closure_adapter_isl_map,
        &tc_transitive_closure_adapter_isl_union_map,
        &tc_transitive_closure_adapter_iterative,
        &tc_transitive_closure_adapter_floyd_warshall,
        &tc_transitive_closure_adapter_omega_map,
        &tc_transitive_closure_adapter_omega_union_map,
        &tc_transitive_closure_adapter_tarjan,
    };

    const int size = sizeof(names) / sizeof(*names);

    int opt = 0;
    do
    {
        for (int i = 0; i < size; ++i)
        {
            printf("%d. %s\n", i+1, names[i]);
        }
        printf("\n");
        printf("Select tc algorithm: ");
        scanf("%d", &opt);
    } while (!(0 < opt && opt <= size));

    return algorithms[opt-1](R, S, exact);
}
