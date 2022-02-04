/*
  initialize a single pair at a separation of x0.
Propagate the pair until they react. 
Repeat many times to determine the survival probability as a function of time, and initial separation
x0.

MODEL PARAMETERS: All lengths in units of nm, all time in units of us.
ka (nm/us)
Dtot=DA+DA (nm^2/us)
sigma (nm)

SIM PARAMETERS
X0: initial separations (nm)
dt: timestep (us)
N: total iterations to simulate (TotalTime=dt*Nsteps) 
 
OUTPUT:
Passoc=1-psurvive
store in 2 columnes
time  passoc(t | x0)



*/
#include <fstream>
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <sys/time.h>
#include "oneD_functions.h"

using namespace std;

double GaussV();

int main(int argc, char *argv[])
{
  
  int i, j, n, k;
  timeval tim;
  gettimeofday(&tim, 0);
  double t1=tim.tv_sec+tim.tv_usec;
  
  int seed=int(t1);
  
  //  seed=1364043256;
  cout <<"seed: "<<seed<<endl;
  srand(seed);
  //double randmax=pow(2.0, 32);
  //double irandmax=1.0/randmax;
  
  double x0=atof(argv[1]);//nm
  cout <<"Read in x0: "<<x0<<endl;
  double sigma=1;//nm
  double ka=0.1;//nm/us
  double deltat=0.1;//us
  int Nit=500;
  int Nrep=100000;
  int rep;
  double D1=0.1;
  double D2=0.1;
  double Dtot=D1+D2;//sum of diffusion constants for both particles.
  
  char tname[100];
  ofstream probfile;
  sprintf(tname, "p_reactVsTime_dt%g_x0_%g_ka%g_D%g_sig%g.dat",deltat, x0, ka,Dtot, sigma);
 
  int Nit1=Nit+1;
  

  int *pSurviveHist=new int[Nit1];
  for(j=0;j<Nit1;j++){
    pSurviveHist[j]=0;
  }
  
  double p1=0; //p1 will always be on the left of p2, so p1<p2. If they swap, reject move. 
  double p2=p1+x0;
  double prevIter=-2;//This is needed for reweighting
  double prevSep=0;//for reweighting
  double psurvivePrev=0;//for reweighting
  double prevNorm=1;//for reweighting

  int it=1;
  double flagstop=0;//terminate trajectories that end in association.
  double rTol=1e-12;
  double *pTheory=new double[Nit1];
  double currSep, ncross;
  double currNorm;
  double reweightRatio, probNoRewgt, probAssoc;
  double dx1, dx2;
  double currp1, currp2;
  double rnum;
  bool flagContinue;
  double scalar=4;//scale the expected displacement in 1D to capture variance beyond the mean
  double Rmax=sigma+scalar*sqrt(2*Dtot*deltat);//max displacement beyond which pAssoc -> 0

  /*theoretical association probability*/
  double dx=0.001;
  double Nx=10000;
  //  double *rRatio=new double[Nit];
  pTheory[0]=0;
  for(i=1;i<Nit+1;i++)
	pTheory[i]=passoc_oned(x0, ka, Dtot, i*deltat, sigma);	//calculate the theoretical association probability from the starting separation as a function of time.
  
  double ptemp=passoc_oned(x0, ka, Dtot, deltat, sigma);
  ofstream tempFile("reweightVsR.dat");
  psurvivePrev=1.0-ptemp;
  for(i=0;i<Nx;i++){
	reweightRatio=reweight_ratio_ps_oned(sigma+i*dx, x0, ka, deltat, Dtot, sigma, psurvivePrev, rTol);
	tempFile<<sigma+i*dx<<' '<<reweightRatio<<endl;
  }
  
  /*PERFORM MULTIPLE REPEATS TO GET STATISTICS ON A SINGLE PAIR*/
  for(rep=0;rep<Nrep;rep++){

	//    cout <<"Repeat: "<<rep<<endl;
    
    /*generate initial coordinates, separated by x0*/
	p1=0; //p1 will always be on the left of p2, so p1<p2. If they swap, reject move. 
	p2=p1+x0;
    prevIter=-2;//This is needed for reweighting
	prevSep=0;//for reweighting
	psurvivePrev=0;//for reweighting
	prevNorm=1;//for reweighting
	Dtot=D1+D2;//sum of diffusion constants for both particles.
	it=0;
	flagContinue=true;//terminate trajectories that end in association.
	
	/*BEGIN ITERATIONS, TERMINATE AT MAX TIME OR IF ASSSOCIATION OCCURS*/
	while(it<Nit && flagContinue==true){
	  it++;
	  currSep=p2-p1;//this will always be positive, because we will reject moves that allow them to swap positions, otherwise we would need to take abs value.
	  ncross=0;
	  /*Calc probability of associating*/
	  
	  if(currSep<Rmax){
	      /*If A is free, calculate association probability, otherwise calculate trajectory probability*/
		ncross=1;//evaluate binding probability!
		if(currSep<sigma){
		  cout <<"WARNING: CURRENT SEPARATION IS LESS THAN SIGMA: "<<p1<<' '<<p2<<' '<<" iter: "<<it<<" rep: "<<rep<<endl;
		  currSep=sigma;
		}
	    currNorm=1.0;//accumulating reweighting ratio
		reweightRatio=1.0;
		
		/*Determine if you were already in reaction zone, so the reweighting ratio is accumulated over 
		  all previous steps*/
		
		if(prevIter==(it-1)){
		  /*you were already in reaction zone in previous step, so update the reweighting ratio as product
			of current ratio with previous steps*/
		  reweightRatio=reweight_ratio_ps_oned(currSep, prevSep, ka, deltat, Dtot, sigma, psurvivePrev, rTol);
		  currNorm=prevNorm*reweightRatio;//Current reweighting ratio is for this step, multiplied by all previous products stored in prevnoem
		  //cout<<" rep: "<<rep<<" it: "<<it<<' '<<currSep<<" reweightRatio: "<<reweightRatio<<' '<<currNorm<<' '<<psurvivePrev<<endl;
		}
	    probNoRewgt=passoc_oned(currSep, ka, Dtot, deltat, sigma);	

		probAssoc=probNoRewgt*currNorm;
	    
		/*Store your current separations for next step evaluation of reweighting ratios*/
		prevSep=currSep;
		prevIter=it;
		prevNorm=currNorm;
		psurvivePrev=1.0-probAssoc;//Survival probability for this step
		
	    
	  }else{
		if(prevIter==(it-1))
		  cout <<rep<<' '<<it<<' '<<currSep<<' '<<currNorm<<endl;
		//just undergo free diffusion
		dx1=sqrt(2.0*deltat*D1)*GaussV();
		dx2=sqrt(2.0*deltat*D2)*GaussV();
		currp1=p1+dx1;
		currp2=p2+dx2;
		while(currp2-currp1<sigma){
		  //reject move due to overlap
		  dx1=sqrt(2.0*deltat*D1)*GaussV();
		  dx2=sqrt(2.0*deltat*D2)*GaussV();
		  currp1=p1+dx1;
		  currp2=p2+dx2;
		}
		//complete the position update
		p1=p1+dx1;
		p2=p2+dx2;
	  }

	  /*Evaluate if the reaction will occur*/
      if(ncross>0){
		
		/*might perform this reaction, depending on k_associate*/
		rnum=1.0*rand()/(1.0*RAND_MAX);
		//cout <<rnum<<endl;
		if(probAssoc>rnum){
		  /*PERFORM ASSOCIATION*/
		  flagContinue=false;
		}
		if(flagContinue==true){
		  //Update using diffusion, reject moves that cause overlap.
		  //just undergo free diffusion
	      dx1=sqrt(2.0*deltat*D1)*GaussV();
	      dx2=sqrt(2.0*deltat*D2)*GaussV();
		  currp1=p1+dx1;
		  currp2=p2+dx2;
		  while(currp2-currp1<sigma){
			//reject move due to overlap
			dx1=sqrt(2.0*deltat*D1)*GaussV();
			dx2=sqrt(2.0*deltat*D2)*GaussV();
			currp1=p1+dx1;
			currp2=p2+dx2;
		  }
		  //complete the position update
		  p1=p1+dx1;
		  p2=p2+dx2;
		  	    
		}
		
	  }//done evaluating reaction probability
	  
      if(flagContinue==true){
		pSurviveHist[it]+=1;//another pair survived to this time-step.
	  }
      
    }//end iterating over time steps
    
	//
    if(rep%10==0){
      probfile.open(tname);
      probfile<<0<<'\t'<<0<<'\t'<<0<<endl;//at time zero, association probability is zero 
      for(i=1;i<Nit+1;i++){
		
		probfile <<i*deltat<<'\t'<<1.0-pSurviveHist[i]/(1.0*rep)<<'\t'<<pTheory[i]<<endl;
      }
      probfile.close();
      
    }
  
    
  }//end looping over repetitions
  
  
  probfile.open(tname);
  probfile<<0<<'\t'<<0<<'\t'<<0<<endl;//at time zero, association probability is zero 
  for(i=1;i<Nit+1;i++){
	probfile <<i*deltat<<'\t'<<1.0-pSurviveHist[i]/(1.0*Nrep)<<'\t'<<pTheory[i]<<endl;
  }
  probfile.close();
    

  cout <<"End Main, complete run "<<endl;
  
}//end main
