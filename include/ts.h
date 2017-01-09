#ifndef _TS_H
#define _TS_H

namespace TSdist {

// Policy enforcement (methods required to compute distances)
class TimeSeriesBase
{
public:
    // Move semantics?
    TimeSeriesBase(const TimeSeriesBase&) = default;
    TimeSeriesBase(TimeSeriesBase&&) = default;
    TimeSeriesBase& operator=(const TimeSeriesBase&) = default;
    TimeSeriesBase& operator=(TimeSeriesBase&&) = default;
    virtual ~TimeSeriesBase() = default;

    // For multivariate series
    virtual int numVars() const = 0;

    // Length of a single variable of the series
    virtual int length() const = 0;

    // This assumes index starts at 0
    virtual double indexSeries(int time_index, int var_index) const = 0;

protected:
    // Abstract classes cannot be instantiated
    TimeSeriesBase() = default;
};

}

#endif // _TS_H
