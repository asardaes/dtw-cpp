#ifndef _DTW_H
#define _DTW_H

#include <vector>
#include "ts.h"

namespace TSdist {

/** Simple DTW distance and optionally a slanted band constraint

    Parameter window_size is for the global constraint. <= 0 means no constraint
    Parameter p is for the Lp norm
    Parameter diag_weight is the weight of the diagonal in the step pattern
 */
double computeDTW(const TimeSeriesBase& x, const TimeSeriesBase& y,
                  int window_size, int p, int diag_weight);


/** Normalized DTW distance and optionally a slanted band constraint

    Parameter window_size is for the global constraint. <= 0 means no constraint
    Parameter p is for the Lp norm
 */
double computeNormalizedDTW(const TimeSeriesBase& x, const TimeSeriesBase& y,
                            int window_size, int p);


/** Simple DTW distance with backtracking and optionally a slanted band constraint

    Parameter window_size is for the global constraint. <= 0 means no constraint
    Parameter p is for the Lp norm
    Parameter diag_weight is the weight of the diagonal in the step pattern
    The indices of the correspondence between x and y are returned in idx and idy
 */
double backtrackDTW(const TimeSeriesBase& x, const TimeSeriesBase& y,
                    int window_size, int p, int diag_weight,
                    std::vector<int>& idx, std::vector<int>& idy);


/** Normalized DTW distance with backtracking and optionally a slanted band constraint

    Parameter window_size is for the global constraint. <= 0 means no constraint
    Parameter p is for the Lp norm
    The indices of the correspondence between x and y are returned in idx and idy
 */
double backtrackNormalizedDTW(const TimeSeriesBase& x, const TimeSeriesBase& y,
                              int window_size, int p,
                              std::vector<int>& idx, std::vector<int>& idy);

}

#endif // _DTW_H
