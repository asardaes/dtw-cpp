#ifndef _DTW_H
#define _DTW_H

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

}

#endif // _DTW_H
