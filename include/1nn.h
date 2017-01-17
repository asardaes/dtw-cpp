#ifndef _1NN_H
#define _1NN_H

#include <cmath>
#include <limits>
#include "ts.h"

namespace TSdist {

/** 1-Nearest-Neighbor in DTW space exploiting its lower bounds

    All series in the database should have the same length as 'query'

    This assumes the time-series database (TSDB) supports iterators that reference/point to
    TimeSeriesBase derivatives.
    Parameters 'L*', 'U*' and 'H' are helpers that should have the same length as 'query'
 */
template<typename TSDB>
const TimeSeriesBase& NearestNeighborDTW(const TSDB& tsdb, const TimeSeriesBase& query,
                                   int window_size, int p, int diag_weight,
                                   TimeSeriesBase& L1, TimeSeriesBase& U1,
                                   TimeSeriesBase& L2, TimeSeriesBase& U2,
                                   TimeSeriesBase& H)
{
    // Window size and length checked here
    computeEnvelop(query, window_size, L1, U1);

    // Initial DTW distance
    double d = std::numeric_limits<double>::max();

    // To return
    const TimeSeriesBase* NN = nullptr;

    for (const TimeSeriesBase& REF : tsdb)
    {
        double lb = 0;

        // LB_Keogh
        for (int i = 0; i < REF.length(); i++)
        {
            if (REF[i][0] > U1[i][0]) {
                H[i][0] = U1[i][0];
                lb += std::pow(REF[i][0] - U1[i][0], p);

            } else if (REF[i][0] < L1[i][0]) {
                H[i][0] = L1[i][0];
                lb += std::pow(L1[i][0] - REF[i][0], p);

            } else {
                H[i][0] = REF[i][0];
            }
        }

        if (lb < d) {
            // LB_Improved
            computeEnvelop(H, window_size, L2, U2);

            for (int i = 0; i < query.length(); i++)
            {
                if (query[i][0] > U2[i][0])
                    lb += std::pow(query[i][0] - U2[i][0], p);
                else if (query[i][0] < L2[i][0])
                    lb += std::pow(L2[i][0] - query[i][0], p);
            }

            if (lb < d) {
                // DTW distance
                double dtw = computeDTW(REF, query, window_size, p, diag_weight);
                dtw = std::pow(dtw, p);

                if (dtw < d) {
                    NN = &REF;
                    d = dtw;
                }
            }
        }
    }

    return reinterpret_cast<const TimeSeriesBase&>(*NN);
}

}

#endif // _1NN_H
