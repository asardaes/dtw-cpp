#include <iostream>
#include <vector>
#include "TSdist.h"

using namespace std;

// Univariate time series class
class UnivariateTimeSeries: public TSdist::TimeSeriesBase
{
public:
    UnivariateTimeSeries(std::vector<double> &singleSeries) {
        if (singleSeries.empty()) throw("Series cannot be empty.");
        series = singleSeries;
    }

    int length() const override {
        return series.size();
    }

    int numVars() const override {
        return 1;
    }

    const double& indexSeries(int time_index, int var_index) const override {
        return series.at(time_index);
    }

    double& indexSeries(int time_index, int var_index) override {
        return series.at(time_index);
    }

private:
    std::vector<double> series;
};

// Example
int main()
{
    std::vector<double> series1 = {1.0, 2.0, 3.0, 4.0};
    std::vector<double> series2 = {0.0, -1.0, 5.0, 2.5};

    std::vector<int> idx;
    std::vector<int> idy;

    UnivariateTimeSeries ts1(series1);
    UnivariateTimeSeries ts2(series2);

    // check that assignment works
    ts1[0][0] = 1;

    cout << "DTW distance is: " << TSdist::computeDTW(ts1, ts2, 0, 2, 2) << endl;
    cout << "nDTW distance is: " << TSdist::computeNormalizedDTW(ts1, ts2, 1, 2) << endl;

    cout << "DTW distance with backtrack is: " <<
        TSdist::backtrackDTW(ts1, ts2, 0, 2, 2, idx, idy) << endl;
    cout << "X indices are: ";
    for (auto i : idx) cout << i << ", ";
    cout << endl;
    cout << "Y indices are: ";
    for (auto j : idy) cout << j << ", ";
    cout << endl;

    cout << "nDTW distance with backtrack is: " <<
        TSdist::backtrackNormalizedDTW(ts1, ts2, 1, 2, idx, idy) << endl;
    cout << "X indices are: ";
    for (auto i : idx) cout << i << ", ";
    cout << endl;
    cout << "Y indices are: ";
    for (auto j : idy) cout << j << ", ";

    return 0;
}
