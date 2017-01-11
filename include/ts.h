#ifndef _TS_H
#define _TS_H

namespace TSdist {

// Policy enforcement (methods required to compute distances)
class TimeSeriesBase
{
public:

    // Move semantics
    TimeSeriesBase(const TimeSeriesBase&) = default;
    TimeSeriesBase(TimeSeriesBase&&) = default;
    virtual TimeSeriesBase& operator=(const TimeSeriesBase&) = default;
    virtual TimeSeriesBase& operator=(TimeSeriesBase&&) = default;
    virtual ~TimeSeriesBase() = default;

    // ============================================================================================
    /* Pure virtual methods */
    // ============================================================================================

    // For multivariate series. Number of variables.
    virtual int numVars() const = 0;

    // Length of a single variable of the series
    virtual int length() const = 0;

    // This assumes indices start at 0
    virtual const double& indexSeries(int time_index, int var_index) const = 0;
    virtual double& indexSeries(int time_index, int var_index) = 0;

private: // More public methods below

    // Proxy class for double subscript
    class SubscriptProxy
    {
    public:
        SubscriptProxy(const TimeSeriesBase* outer_this, int time_index) :
            _outer_this(outer_this) ,
            _time_index(time_index)
        { };

        const double& operator[](int var_index) const {
            return _outer_this->indexSeries(_time_index, var_index);
        }

        double& operator[](int var_index) {
            return const_cast<TimeSeriesBase*>(_outer_this)->indexSeries(_time_index, var_index);
        }

    private:
        const TimeSeriesBase* _outer_this;
        int _time_index;
    };

public:

    // ============================================================================================
    /* Non-virtual methods */
    // ============================================================================================

    // --------------------------------------------------------------------------------------------
    /* Double subscript operator will call indexSeries through SubscriptProxy struct */
    // --------------------------------------------------------------------------------------------

    // Const version
    const SubscriptProxy operator[](int time_index) const {
        return SubscriptProxy{this, time_index};
    }

    // Non-const version
    SubscriptProxy operator[](int time_index) {
        return SubscriptProxy{this, time_index};
    }

protected:

    // Abstract classes cannot be instantiated
    TimeSeriesBase() = default;

};

}

#endif // _TS_H
