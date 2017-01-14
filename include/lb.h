#ifndef _LB_H
#define _LB_H

#include "ts.h"

namespace TSdist {

/** Compute warping envelop based on a Sakoe-Chiba window

    (c) Daniel Lemire, 2008
    Taken from https://github.com/lemire/lbimproved/blob/master/dtw.h

    All series must have the same length.
    Only univariate series supported.

    Parameter 'window_size' must be greater than zero.
    The objects 'lower_envelop' and 'upper_envelop' are updated with the envelops corresponding to
    'x'.
 */
void computeEnvelop(const TimeSeriesBase& x, int window_size,
                    TimeSeriesBase& lower_envelop, TimeSeriesBase& upper_envelop);

/** DTW lower bound: LB_Keogh

    All series must have the same length.
    Only univariate series supported.
    This version assumes that envelops are already available. See function computeEnvelop above.

    Parameter x is the reference
    Parameter y is the query
    Parameter p is for the Lp norm
    Envelops must correspond to 'y'
 */
double lbKeogh(const TimeSeriesBase& x, const TimeSeriesBase& y, int p,
               const TimeSeriesBase& lower_envelop, const TimeSeriesBase& upper_envelop);

/** DTW lower bound: LB_Improved

    All series must have the same length.
    Only univariate series supported.
    This version computes both sets of required envelops and saves them in lower/uper_envelop.

    Parameter x is the reference
    Parameter y is the query
    Parameter p is for the Lp norm
    Parameter window_size is for the window constraint
    Parameter H is a helper series that is needed for the calculation
 */
double lbImproved(const TimeSeriesBase& x, const TimeSeriesBase& y,
                  int window_size, int p,
                  TimeSeriesBase& lower_envelop, TimeSeriesBase& upper_envelop,
                  TimeSeriesBase& H);

}

#endif // _LB_H
