#include <iostream>
#include <cmath>

using namespace std;


/*evaluate the association probability (1-survival probability) between a pair of particles separated
  by r1, in 1-dimensions.
 */
double reweight_ratio_ps_oned(double r1, double prevR, double ka, double dt, double D, double sigma, double psurvivePrev, double rTol)
{
  /*Calculate the irreversible probability given prevR and currR=r1*/
  double dist=r1-prevR;
  double fourDdt=4.0*D*dt;
  double arg1=r1+prevR-2.0*sigma;
  double cof=1.0/(sqrt(4.0*M_PI*D*dt));
  double term1=exp(-dist*dist/fourDdt)+exp(-arg1*arg1/fourDdt);
  term1=term1*cof;
  double term2=-ka/D*exp(ka*ka*dt/D+ka/D*arg1)*erfc(arg1/sqrt(fourDdt)+ka*sqrt(dt/D));
  double pirr=term1+term2;
  
  
  /*Calculate the free probability given prevR and r1*/
  double adist=r1+prevR;
  double pfree=exp(-dist*dist/fourDdt)-exp(-adist*adist/fourDdt);
  pfree=pfree*cof;
  
  //  double pnorm=1+0.5*(erf((prevR-sigma)/sqrt(4.0*D*dt))-erf((prevR+sigma)/sqrt(4.0*D*dt)));//Normalization of the free propagator (which is restricted to 0->infinity), from sigma->infinity)
  double pnorm=0.5*(erf((prevR-sigma)/sqrt(4.0*D*dt))+erf((prevR+sigma)/sqrt(4.0*D*dt)));
  pfree=pfree/pnorm;
  /*Reweighting ratio is pirr/pfree*, where pfree* has been rescaled by the normalization from sigma->inf range and from the previous survival probability.
   */
  
  double ratio;
  //  if(pirr-pfree * psurvivePrev <rTol)
  //ratio=1.0;//ignore precision issues between them
  //else
  ratio=pirr/(pfree*psurvivePrev);
  //cout <<r1<<' '<<pirr<<' '<<pfree<<' '<<ratio<<endl;
  return ratio;
}
