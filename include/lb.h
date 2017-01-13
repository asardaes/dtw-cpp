#ifndef _LB_H
#define _LB_H

namespace TSdist {

#include "ts.h"

/** Compute warping envelop based on a Sakoe-Chiba window

    Parameter 'window_size' must be greater than zero.
    The objects 'lower_envelop' and 'upper_envelop' are updated with the envelops corresponding to
    'x'. All must have the same length.
    Only univariate series supported.

    (c) Daniel Lemire, 2008
    Taken from https://github.com/lemire/lbimproved/blob/master/dtw.h
 */
void computeEnvelop(const TimeSeriesBase& x, int window_size,
                    TimeSeriesBase& lower_envelop, TimeSeriesBase& upper_envelop);

}

#endif // _LB_H
