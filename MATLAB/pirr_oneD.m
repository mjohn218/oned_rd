function[pirr, pAssoc]=pirr_oneD(ka, D, sigma, dt, r, r0)

dist=r-r0;
fourDdt=4.0*D*dt;
arg1=r+r0-2.0*sigma;
cof=1.0/(sqrt(4.0*pi*D*dt));
 term1=exp(-dist.*dist/fourDdt)+exp(-arg1.*arg1/fourDdt);
  term1=term1*cof;
 term2=-ka/D*exp(ka*ka*dt/D+ka/D*arg1).*erfc(arg1/sqrt(fourDdt)+ka*sqrt(dt/D));
 pirr=term1+term2;
 
distA=r0-sigma;
 argA=distA/sqrt(4.0*D*dt);
pAssoc=erfc(argA)-exp(ka/D*distA+ka*ka*dt/D)*erfc(argA+ka*sqrt(dt/D));
