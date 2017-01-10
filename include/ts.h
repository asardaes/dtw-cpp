#ifndef _TS_H
#define _TS_H

namespace TSdist {

// Policy enforcement (methods required to compute distances)
class TimeSeriesBase
{
public:
    // Move semantics?
    TimeSeriesBase(TimeSeriesBase const &) = default;
    TimeSeriesBase(TimeSeriesBase&&) = default;
    TimeSeriesBase& operator=(TimeSeriesBase const &) = default;
    TimeSeriesBase& operator=(TimeSeriesBase&&) = default;
    virtual ~TimeSeriesBase() = default;

    // ============================================================================================
    /* Pure virtual methods */
    // ============================================================================================

    // For multivariate series
    virtual int numVars() const = 0;

    // Length of a single variable of the series
    virtual int length() const = 0;

    // This assumes indices start at 0
    virtual const double& indexSeries(int time_index, int var_index) const = 0;
    virtual double& indexSeries(int time_index, int var_index) = 0;

    // ============================================================================================
    /* Virtual methods with defaults */
    // ============================================================================================

    // --------------------------------------------------------------------------------------------
    /* Double subscript operator will call indexSeries */
    // --------------------------------------------------------------------------------------------

    // Proxy class for double subscript
    class SubscriptProxy
    {
    public:
        SubscriptProxy(TimeSeriesBase const &outer, int time_index) :
            _outer(outer) ,
            _time_index(time_index)
        { };

        // If Time Series is const, return const value
        const double& operator[](int var_index) const {
            return _outer.indexSeries(_time_index, var_index);
        }

        // If Time Series is not const (e.g. when assigning using double brackets), return non-const
        double& operator[](int var_index) {
            // I couldn't find a way to do this without using const_cast
            return const_cast<TimeSeriesBase&>(_outer).indexSeries(_time_index, var_index);
        }

    private:
        const TimeSeriesBase& _outer;
        int _time_index;
    };

    // Const version
    virtual const SubscriptProxy operator[](int time_index) const {
        return SubscriptProxy(*this, time_index);
    }

    // Non-const version
    virtual SubscriptProxy operator[](int time_index) {
        return SubscriptProxy(*this, time_index);
    }

protected:
    // Abstract classes cannot be instantiated
    TimeSeriesBase() = default;
};

}

#endif // _TS_H
