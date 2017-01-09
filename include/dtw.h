#ifndef _DTW_H
#define _DTW_H

#include <vector>
#include "ts.h"

namespace TSdist {

/** Simple DTW distance and optionally a slanted band constraint
 *
 * Parameter p is for the Lp norm
 * Parameter diag_step is the weight of the diagonal in the step pattern
 */
double computeDTW(TimeSeriesBase const &x, TimeSeriesBase const &y,
                  int window_size, int p, int diag_step);


/** Normalized DTW distance and optionally a slanted band constraint
 *
 * Parameter p is for the Lp norm
 */
double computeNormalizedDTW(TimeSeriesBase const &x, TimeSeriesBase const &y,
                            int window_size, int p);


/** Simple DTW distance with backtracking and optionally a slanted band constraint
 *
 * Parameter p is for the Lp norm
 * Parameter diag_step is the weight of the diagonal in the step pattern
 * The indices of the correspondance between x and y are returned in idx and idy
 */
double computeBacktrackDTW(TimeSeriesBase const &x, TimeSeriesBase const &y,
                           int window_size, int p, int diag_step,
                           std::vector<int> &idx, std::vector<int> &idy);


/** Normalized DTW distance with backtracking and optionally a slanted band constraint
 *
 * Parameter p is for the Lp norm
 * Parameter diag_step is the weight of the diagonal in the step pattern
 * The indices of the correspondance between x and y are returned in idx and idy
 */
double computeBacktrackNormalizedDTW(TimeSeriesBase const &x, TimeSeriesBase const &y,
                                     int window_size, int p,
                                     std::vector<int> &idx, std::vector<int> &idy);

}

#endif // _DTW_H
