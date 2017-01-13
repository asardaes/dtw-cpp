#include <deque>
#include "ts.h"
#include "lb.h"

namespace TSdist {

void computeEnvelop(const TimeSeriesBase& x, int window_size,
                    TimeSeriesBase& lower_envelop, TimeSeriesBase& upper_envelop)
{
    if (window_size < 1)
        throw("Window size must be positive.");

    if (x.numVars() > 1)
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

}