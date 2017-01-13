#ifndef _LB_H
#define _LB_H

namespace TSdist {

#include "ts.h"

/** Compute warping envelop based on a Sakoe-Chiba window

    Only univariate series supported.
    All series must have the same length.

    Parameter 'window_size' must be greater than zero.
    The objects 'lower_envelop' and 'upper_envelop' are updated with the envelops corresponding to
    'x'.

    (c) Daniel Lemire, 2008
    Taken from https://github.com/lemire/lbimproved/blob/master/dtw.h
 */
void computeEnvelop(const TimeSeriesBase& x, int window_size,
                    TimeSeriesBase& lower_envelop, TimeSeriesBase& upper_envelop);

/** DTW lower bound: LB_Keogh

    Both 'x' and 'y' must have the same length
    Only univariate series supported.
    This version assumes that envelops are already available. See function computeEnvelop above.

    Parameter x is the reference
    Parameter y is the query
    Parameter p is for the Lp norm
    Envelops must correspond to 'y'
 */
double lbKeogh(const TimeSeriesBase& x, const TimeSeriesBase& y, int p,
               const TimeSeriesBase& lower_envelop, const TimeSeriesBase& upper_envelop);

}

#endif // _LB_H
