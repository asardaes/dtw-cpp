#include <cmath>
#include <limits>
#include "ts.h"
#include "dtw.h"

namespace TSdist {

static const double NOT_VISITED = -1;
static const int STEP_DIAG = 0;
static const int STEP_LEFT = 1;
static const int STEP_UP = 2;
static const double DBL_MAX = std::numeric_limits<double>::max();

double lnorm(TimeSeriesBase const &x, TimeSeriesBase const &y,
             int p, int time_x, int time_y)
{
    double result = 0;

    for (int k = 0; k < x.numVars(); k++) {
        result += std::abs(std::pow(x.indexSeries(time_x, k) - y.indexSeries(time_y, k), p));
    }

    return std::pow(result, 1.0 / p);
}

double computeDTW(TimeSeriesBase const &x, TimeSeriesBase const &y,
                  int window_size, int p, int diag_step)
{
    if (x.numVars() != y.numVars())
        throw("Series must have the same number of variables.");

    if (p < 1)
        throw("Parameter p must be positive.");

    if (diag_step != 1 && diag_step != 2)
        throw("Diagonal step can only be 1 or 2.");

    int nx = x.length();
    int ny = y.length();

    // local/global cost matrix, only 2 rows to implement memory-saving version
    double D[2][ny + 1];
    // possible directions to take when traversing D
    double tuple_direction[3];

    int i, j, direction;
    double local_cost;

    // initialization of first row and first column
    for (j = 0; j <= ny; j++) D[0][j] = NOT_VISITED;
    D[1][0] = NOT_VISITED;

    // first value, must set here to avoid multiplying by step
    D[1][1] = std::pow(lnorm(x, y, p, 0, 0), p);

    // dynamic programming
    for (i = 1; i <= nx; i++)
    {
        int j1, j2;

        // adjust limits depending on window
        if (window_size < 1) {
            j1 = 1;
            j2 = ny;

        } else {
            j1 = std::ceil((double) i * ny / nx - window_size);
            j2 = std::floor((double) i * ny / nx + window_size);

            j1 = j1 > 1 ? j1 : 1;
            j2 = j2 < ny ? j2 : ny;
        }

        for (j = 1; j <= ny; j++)
        {
            // very first value already set above
            if (i == 1 && j == 1) continue;

            if (j < j1 || j > j2) {
                // cell outside of window
                D[i % 2][j] = NOT_VISITED;
                continue;
            }

            // l-norm for single time-point, p only affects multivariate series
            local_cost = std::pow(lnorm(x, y, p, i-1, j-1), p);

            // first assign the current values of possible directions
            tuple_direction[STEP_DIAG] = D[(i - 1) % 2][j - 1];
            tuple_direction[STEP_LEFT] = D[i % 2][j - 1];
            tuple_direction[STEP_UP] = D[(i - 1) % 2][j];

            // then update values appropriately
            tuple_direction[STEP_DIAG] = (tuple_direction[STEP_DIAG] == NOT_VISITED) ? DBL_MAX : tuple_direction[STEP_DIAG] + diag_step * local_cost;
            tuple_direction[STEP_LEFT] = (tuple_direction[STEP_LEFT] == NOT_VISITED) ? DBL_MAX : tuple_direction[STEP_LEFT] + local_cost;
            tuple_direction[STEP_UP] = (tuple_direction[STEP_UP] == NOT_VISITED) ? DBL_MAX : tuple_direction[STEP_UP] + local_cost;

            // which direction has the least associated cost?
            direction = (tuple_direction[STEP_LEFT] < tuple_direction[STEP_DIAG]) ? STEP_LEFT : STEP_DIAG;
            direction = (tuple_direction[STEP_UP] < tuple_direction[direction]) ? STEP_UP : direction;

            D[i % 2][j] = tuple_direction[direction];
        }
    }

    // calculate p-root on the very last value
    return std::pow(D[nx % 2][ny], 1.0 / p);
}

double computeNormalizedDTW(TimeSeriesBase const &x, TimeSeriesBase const &y,
                            int window_size, int p)
{
    return computeDTW(x, y, window_size, p, 2) / (x.length() + y.length());
}

}
