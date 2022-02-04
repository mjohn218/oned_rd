%length is in nm
%time is in us

function[pfree]=pfree_1d(D, sigma, dt, r, r0)



cof=1/sqrt(4.0*pi*D*dt);

dist=r-r0;
adist=r+r0;
fDt=4.0*D*dt;
%terms=exp(-dist.*dist/fDt)+exp(-adist.*adist/fDt);
terms=exp(-dist.*dist/fDt)-exp(-adist.*adist/fDt);

pfree=cof*terms;

%norm=1+0.5*( erf((r0-sigma)/sqrt(4*D*dt))-erf((r0+sigma)/sqrt(4*D*dt)))
norm=0.5*( erf((r0-sigma)/sqrt(4*D*dt))+erf((r0+sigma)/sqrt(4*D*dt)))

pfree=pfree/norm;
vec=find(r<=sigma);
rID=vec(end)
pfree(1:rID)=0;