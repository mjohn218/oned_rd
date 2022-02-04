function[pAssoc]=preact_oneD(ka, D, sigma, dt, r0)

distA=r0-sigma;
 argA=distA/sqrt(4.0*D*dt);
pAssoc=erfc(argA)-exp(ka/D*distA+ka*ka*dt/D).*erfc(argA+ka*sqrt(dt/D));
