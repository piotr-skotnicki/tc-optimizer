#include "slicing.h"

#include <isl/set.h>
#include <isl/map.h>

__isl_give isl_set* tc_uds_set(__isl_keep isl_map* R)
{
    isl_set* uds = isl_set_subtract(isl_map_domain(isl_map_copy(R)), isl_map_range(isl_map_copy(R)));
    return uds;
}

__isl_give isl_set* tc_udd_set(__isl_keep isl_map* R)
{
    isl_set* udd = isl_set_subtract(isl_map_range(isl_map_copy(R)), isl_map_domain(isl_map_copy(R)));
    return udd;
}

__isl_give isl_set* tc_ind_set(__isl_keep isl_set* LD, __isl_keep isl_map* R)
{
    isl_set* ind = isl_set_subtract(isl_set_subtract(isl_set_copy(LD), isl_map_domain(isl_map_copy(R))), isl_map_range(isl_map_copy(R)));
    return ind;
}
