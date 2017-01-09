#ifndef _DTW_H
#define _DTW_H

#include "ts.h"

namespace TSdist {

/** Simple DTW distance with L1 or L2 norm and optionally a slanted band constraint
 *
 * Parameter p is for the Lp norm
 */
double computeDTW(TimeSeriesBase const &x, TimeSeriesBase const &y,
                  int window_size, int p, int diag_step);

}

#endif // _DTW_H
