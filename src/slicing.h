#ifndef TC_SLICING_H
#define	TC_SLICING_H

#include <isl/set.h>
#include <isl/map.h>

__isl_give isl_set* tc_uds_set(__isl_keep isl_map* R);

__isl_give isl_set* tc_udd_set(__isl_keep isl_map* R);

__isl_give isl_set* tc_ind_set(__isl_keep isl_set* LD, __isl_keep isl_map* R);

__isl_give isl_map* tc_Rusc_map(__isl_keep isl_map* R);

#endif // TC_SLICING_H
