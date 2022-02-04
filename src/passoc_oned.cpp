#include <iostream>
#include <cmath>

using namespace std;


/*evaluate the association probability (1-survival probability) between a pair of particles separated
  by r1, in 1-dimensions.
 */
double passoc_oned(double r1, double ka, double D, double dt, double sigma)
{
  double dist=r1-sigma;
  double arg1=dist/sqrt(4.0*D*dt);
  double pAssoc=erfc(arg1)-exp(ka/D*dist+ka*ka*dt/D)*erfc(arg1+ka*sqrt(dt/D));

  return pAssoc;
}
