#include <iostream>
#include <vector>
#include "TSdist.h"

using namespace std;

// Univariate time series class
class UnivariateTimeSeries: public TSdist::TimeSeriesBase
{
public:
    UnivariateTimeSeries(vector<double> &singleSeries) {
        if (singleSeries.empty()) throw("Series cannot be empty.");
        series = singleSeries;
    }

    int length() const {
        return series.size();
    }

    int numVars() const {
        return 1;
    }

    double indexSeries(int time_index, int var_index) const {
        return series[time_index];
    }

private:
    vector<double> series;
};

// Example
int main()
{
    std::vector<double> series1 = {1.0, 2.0, 3.0, 4.0};
    std::vector<double> series2 = {0.0, -1.0, 5.0, 2.5};

    UnivariateTimeSeries ts1(series1);
    UnivariateTimeSeries ts2(series2);

    cout << "DTW distance is: " << TSdist::computeDTW(ts1, ts2, 1, 2, 2) << endl;
    cout << "nDTW distance is: " << TSdist::computeNormalizedDTW(ts1, ts2, 1, 2) << endl;
    return 0;
}
