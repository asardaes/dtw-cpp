#include <cmath>
#include <deque>
#include "ts.h"
#include "lb.h"

namespace TSdist {

// ================================================================================================
/* Warping envelop */
// ================================================================================================
void computeEnvelop(const TimeSeriesBase& x, int window_size,
                    TimeSeriesBase& lower_envelop, TimeSeriesBase& upper_envelop)
{
    if (window_size < 1)
        throw("Window size must be positive.");

    if (x.numVars() != 1)
        throw("Only univariate series are supported.");

    if (x.length() != lower_envelop.length() || x.length() != upper_envelop.length())
        throw("Length mismatch between x and the envelops.");

    window_size = window_size * 2 + 1;
    int constraint = (window_size - 1) / 2;
    int array_size = static_cast<int>(x.length());
    std::deque<int> maxfifo, minfifo;

    maxfifo.push_back(0);
    minfifo.push_back(0);

    for(int i = 1; i < array_size; ++i) {
        if(i >= constraint + 1) {
            upper_envelop[i - constraint - 1][0] = x[maxfifo.front()][0];
            lower_envelop[i - constraint - 1][0] = x[minfifo.front()][0];
        }

        if(x[i][0] > x[i - 1][0]) { //overshoot
            maxfifo.pop_back();

            while(maxfifo.size() > 0) {
                if(x[i][0] <= x[maxfifo.back()][0]) break;
                maxfifo.pop_back();
            }

        } else {
            minfifo.pop_back();

            while(minfifo.size() > 0) {
                if(x[i][0] >= x[minfifo.back()][0]) break;
                minfifo.pop_back();
            }
        }

        maxfifo.push_back(i);
        minfifo.push_back(i);

        if(i == window_size + maxfifo.front())
            maxfifo.pop_front();
        else if(i == window_size + minfifo.front())
            minfifo.pop_front();
    }

    for(int i = x.length(); i <= array_size + constraint; ++i) {
        upper_envelop[i - constraint - 1][0] = x[maxfifo.front()][0];
        lower_envelop[i - constraint - 1][0] = x[minfifo.front()][0];

        if(i - maxfifo.front() >= window_size) maxfifo.pop_front();
        if(i - minfifo.front() >= window_size) minfifo.pop_front();
    }
}

// ================================================================================================
/* LB_Keogh */
// ================================================================================================
double lbKeogh(const TimeSeriesBase& x, const TimeSeriesBase& y, int p,
               const TimeSeriesBase& lower_envelop, const TimeSeriesBase& upper_envelop)
{
    if (p < 1)
        throw("Parameter p must be positive.");

    if (x.numVars() != 1 || y.numVars() != 1)
        throw("Only univariate series are supported.");

    if (x.length() != y.length())
        throw("Length mismatch between x and y.");

    if (y.length() != lower_envelop.length() || y.length() != upper_envelop.length())
        throw("Length mismatch between y and the envelops.");

    double lb = 0;

    for (int i = 0; i < x.length(); i++)
    {
        if (x[i][0] > upper_envelop[i][0])
            lb += std::pow(x[i][0] - upper_envelop[i][0], p);
        else if (x[i][0] < lower_envelop[i][0])
            lb += std::pow(lower_envelop[i][0] - x[i][0], p);
    }

    return std::pow(lb, 1.0 / p);
}

// ================================================================================================
/* LB_Improved */
// ================================================================================================
double lbImproved(const TimeSeriesBase& x, const TimeSeriesBase& y,
                  int window_size, int p,
                  TimeSeriesBase& lower_envelop, TimeSeriesBase& upper_envelop,
                  TimeSeriesBase& H)
{
    if (p < 1)
        throw("Parameter p must be positive.");

    if (x.numVars() != 1 || y.numVars() != 1)
        throw("Only univariate series are supported.");

    if (x.length() != y.length())
        throw("Length mismatch between x and y.");

    double lb = 0;

    // window size and length checked here
    computeEnvelop(y, window_size, lower_envelop, upper_envelop);

    for (int i = 0; i < x.length(); i++)
    {
        if (x[i][0] > upper_envelop[i][0]) {
            H[i][0] = upper_envelop[i][0];
            lb += std::pow(x[i][0] - upper_envelop[i][0], p);

        } else if (x[i][0] < lower_envelop[i][0]) {
            H[i][0] = lower_envelop[i][0];
            lb += std::pow(lower_envelop[i][0] - x[i][0], p);

        } else {
            H[i][0] = x[i][0];
        }
    }

    // window size and length checked here
    computeEnvelop(H, window_size, lower_envelop, upper_envelop);

    for (int i = 0; i < y.length(); i++)
    {
        if (y[i][0] > upper_envelop[i][0])
            lb += std::pow(y[i][0] - upper_envelop[i][0], p);
        else if (y[i][0] < lower_envelop[i][0])
            lb += std::pow(lower_envelop[i][0] - y[i][0], p);
    }

    return std::pow(lb, 1.0 / p);
}

}
