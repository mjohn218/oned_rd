#include <cstdlib>
#include <cmath>

double GaussV()
{
    double R=2.0;
    double V1;

    while (R >= 1.0) {
	  V1 = 2.0 *rand()/(1.0*RAND_MAX) - 1.0;
	  double V2 = 2.0 * rand()/(1.0*RAND_MAX) - 1.0;
        R = (V1 * V1) + (V2 * V2);
    }
    return (V1 * sqrt(-2.0 * log(R) / R));
}
