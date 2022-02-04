/*
  initialize a single pair A at center.
Reversible interactions with N B particles, B's don't interact with one another.


MODEL PARAMETERS: All lengths in units of nm, all time in units of us.
ka (nm/us)
Dtot=DA+DA (nm^2/us)
sigma (nm)
kb (1/s)

SIM PARAMETERS
NB: number of B particles (nm)
dt: timestep (us)
N: total iterations to simulate (TotalTime=dt*Nsteps) 
 
OUTPUT:
pBound vs time
store in 2 columnes
time  pBound



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
  
  int NB=atoi(argv[1]);//nm
  cout <<"Read in NB: "<<NB<<endl;
  double sigma=1;//nm
  double ka=0.1;//nm/us
  double kb=500;///s
  double deltat=0.1;//us
  double boxLength=100;//nm
  int Nit=10000000;
  int Nrep=5000;
  int rep;
  double D1=0;
  double D2=0.1;
  double Dtot=D1+D2;//sum of diffusion constants for both particles.
  
  char tname[100];
  ofstream probfile;
  sprintf(tname, "pBound_VsTime_dt%g_NB%d_ka%g_D%g_kb%g_sig%g.dat",deltat, NB, ka, Dtot,kb, sigma);

  char mname[100];
  ofstream meansfile;
  sprintf(mname, "meanpBound_VsRepeat_dt%g_NB%d_ka%g_D%g_kb%g_sig%g.dat",deltat, NB, ka, Dtot,kb, sigma);
  meansfile.open(mname);
  
  double dur1=Nit*deltat*0.75;
  double dur2=Nit*deltat*0.5;
  meansfile <<" Repeat. "<<"avg over last "<<dur1<<" us. Avg over last "<<dur2<<" us. "<<endl;
  int Nit1=Nit+1;
  

  double *pBoundSq=new double[Nit1];
  double *pBoundOneTraj=new double[Nit1];
  double *pBoundHist=new double[Nit1];
  for(j=0;j<Nit1;j++){
    pBoundHist[j]=0;
	pBoundOneTraj[j]=0;
	pBoundSq[j]=0;
  }
  
  double p1=0; //coordinates of the A particles
  double *p2=new double[NB+1];
  int *prevIter=new int[NB+1];//-2;//This is needed for reweighting
  double *prevSep=new double[NB+1];//=0;//for reweighting
  double *psurvivePrev=new double[NB+1];//=0;//for reweighting
  double *prevNorm=new double[NB+1];//1;//for reweighting
  bool * boundStatus=new bool[NB+1];
  bool * stopMoving=new bool[NB+1];
  
  int maxPartner=50;//just to be safe!
  int *crossPartner=new int[maxPartner];//max number of partners
  double *probAssoc=new double[maxPartner];//max number of partners

  int it=1;
  //  double dit=1;
  double mean, var;
  double flagstop=0;//terminate trajectories that end in association.
  double rTol=1e-12;
  double *pTheory=new double[Nit1];
  double currSep, ncross;
  double currNorm;
  double reweightRatio, probNoRewgt;
  double dx1, dx2;
  double currp1, currp2;
  double rnum;
  bool flagContinue;
  double scalar=4;//scale the expected displacement in 1D to capture variance beyond the mean
  double Rmax=sigma+scalar*sqrt(2*Dtot*deltat);//max displacement beyond which pAssoc -> 0
  double us_to_s=1E-6;
  /*theoretical association probability*/
  double dx=0.001;
  double Nx=10000;
  //  double *rRatio=new double[Nit];
  // pTheory[0]=0;
//   for(i=1;i<Nit+1;i++)
// 	pTheory[i]=passoc_oned(x0, ka, Dtot, i*deltat, sigma);	//calculate the theoretical association probability from the starting separation as a function of time.
  
  int ncCurr;
  double dist;
  double x0=1;
  double ptemp=passoc_oned(x0, ka, Dtot, deltat, sigma);
  ofstream tempFile("reweightVsR.dat");
  double psPrev=1.0-ptemp;
  for(i=0;i<Nx;i++){
	reweightRatio=reweight_ratio_ps_oned(sigma+i*dx, x0, ka, deltat, Dtot, sigma, psPrev, rTol);
	tempFile<<sigma+i*dx<<' '<<reweightRatio<<endl;
  }
  int currPartnerA, aPartner;
  double aBound;
  int nevent;
  double pDissoc;
  double Keq=ka/kb*1E6;//units of nm
  double B0=NB*1.0/boxLength;//particles/nm
  double currBound=0;
  double probBoundTheory=Keq*B0/(1+Keq*B0);
  int MaxEvent=5000;
  cout<<" Expected Bound probability : "<<probBoundTheory<<" total repetitions: "<<Nrep<<endl;
  /*PERFORM MULTIPLE REPEATS TO GET STATISTICS ON A SINGLE PAIR*/
  for(rep=0;rep<Nrep;rep++){

	cout <<"Repeat: "<<rep<<endl;
    currBound=0;
    for(it=0;it<Nit1;it++)
	  pBoundOneTraj[it]=0;//reset to zero each repeat
	p1=0; //p1 will always be on the left of p2, so p1<p2. If they swap, reject move. 
	/*generate initial coordinates for all B molecules*/
    for(i=0;i<NB+1;i++){
	  prevIter[i]=-2;//This is needed for reweighting
	  prevSep[i]=0;//for reweighting
	  psurvivePrev[i]=0;//for reweighting
	  prevNorm[i]=1;//for reweighting
	  boundStatus[i]=false;//not bound
	  stopMoving[i]=false;//did not just dissociate
	}
	aBound=0;//switch this to 1 if A becomes bound
	currPartnerA=-1;//not bound to anyone.
	Dtot=D1+D2;//sum of diffusion constants for both particles.
	//dit=0;
	it=0;
	//select random coordinates, do not allow overlap with sigma
	p2[0]=0;//this is the A particle.
	for(i=1;i<NB+1;i++){
	  p2[i]=rand()/(1.0*RAND_MAX)*boxLength;//position is between 0 and boxlength
	  while(p2[i]<sigma)
		p2[i]=rand()/(1.0*RAND_MAX)*boxLength;//position is between 0 and boxlength
	}

	//flagContinue=true;//terminate trajectories that end in association.
	nevent=0;
	/*BEGIN ITERATIONS, TERMINATE AT MAX TIME OR IF ASSSOCIATION OCCURS*/
	while(it<Nit1){
	  //dit=it+1;
	  it++;
	  ncross=0;
	  /*First evaluate dissociation*/
	  if(aBound==1){
		
		pDissoc=1.0-exp(-kb*deltat*us_to_s);
		rnum=1.0*rand()/(1.0*RAND_MAX);
		if(pDissoc>rnum){
		  //cout <<" DISSOCIATION: "<<it<<' '<<rep<<endl;
		  aBound=0;//Change status of A from bound to free.
		  ncross=-1;//do not try to bind this step.
		  aPartner=currPartnerA;
		  boundStatus[aPartner]=false;
		  stopMoving[aPartner]=true;
		  currPartnerA=-1;//no more partner
		  //p2[Apartner];
		  prevNorm[aPartner]=1.0;
		  
		}
		
	  }else{
		/*Calc probability of associating for each molecule of B with A, if A is unbound*/
		//A is free to bind
		if(ncross>-1){
		  //don't evaluate binding if A just dissociated.
		  for(i=1;i<NB+1;i++){
			currSep=p2[i]-p1;//p1 should never change--always at zero
			if(currSep<Rmax){
			  /*If A is free, calculate association probability, otherwise calculate trajectory probability*/
			  ncCurr=ncross;
			  
			  if(currSep<sigma){
				cout <<"WARNING: CURRENT SEPARATION IS LESS THAN SIGMA: "<<p1<<' '<<p2<<' '<<" iter: "<<it<<" rep: "<<rep<<endl;
				currSep=sigma;
			  }
			  currNorm=1.0;//accumulating reweighting ratio
			  reweightRatio=1.0;
			  
			  /*Determine if you were already in reaction zone, so the reweighting ratio is accumulated over 
				all previous steps*/
			  
			  if(prevIter[i]==(it-1)){
				/*you were already in reaction zone in previous step, so update the reweighting ratio as product
				  of current ratio with previous steps*/
				reweightRatio=reweight_ratio_ps_oned(currSep, prevSep[i], ka, deltat, Dtot, sigma, psurvivePrev[i], rTol);
				currNorm=prevNorm[i]*reweightRatio;//Current reweighting ratio is for this step, multiplied by all previous products stored in prevnoem
				//cout<<" rep: "<<rep<<" it: "<<it<<' '<<currSep<<" reweightRatio: "<<reweightRatio<<' '<<currNorm<<' '<<psurvivePrev<<endl;
			  }
			  probNoRewgt=passoc_oned(currSep, ka, Dtot, deltat, sigma);	
			  
			  probAssoc[ncCurr]=probNoRewgt*currNorm;
			  crossPartner[ncCurr]=i;
			  /*Store your current separations for next step evaluation of reweighting ratios*/
			  prevSep[i]=currSep;
			  prevIter[i]=it;
			  prevNorm[i]=currNorm;
			  psurvivePrev[i]=1.0-probAssoc[ncCurr];//Survival probability for this step
			  ncross++;
			  
			}//done updating reaction status for molecule i
			
			
		  }//done looping over all molecules NB
		}//only check association if did not just dissociate
	  }//end performing bimolecular reaction evaluations, skip this if the A molecule started bound.
	
	  /*Now Evaluate if the reaction will occur, if we just dissociated or B's are all far away, then no attempt is made to react*/
	  if(ncross>0){
		
		/*Loop over all partners within the reaction zone*/
		for(i=0;i<ncross;i++){
		  
		  /*might perform this reaction, depending on k_associate*/
		  rnum=1.0*rand()/(1.0*RAND_MAX);
		  //cout <<rnum<<endl;
		  if(probAssoc[i]>rnum){
			/*PERFORM ASSOCIATION*/
			cout <<" ASSOCIATION ! "<<aPartner<<" iter: "<<it<<' '<<rep<<" currbound: "<<currBound/(1.0*it)<<endl;
			nevent++;
			aBound=1;
			aPartner=crossPartner[i];
			currPartnerA=aPartner;//this is who A is bound to.
			//move the position of aPartner to sigma
			p2[aPartner]=sigma+1e-12;//small shift for precision
			boundStatus[aPartner]=true;
			stopMoving[aPartner]=true;
			i=ncross;//break out of this loop so no more association events are attempted.
		  }
		}
	  }//done attempting to perform a bimolecular reaction

	  /*Now we need to move all particles, whether a reaction occured or not, and avoid overlap*/
	  for(i=1;i<NB+1;i++){
		if(boundStatus[i]==false && stopMoving[i]==false){
		  /*just undergo free diffusion*/
		  dx2=sqrt(2.0*deltat*D2)*GaussV();
		  //currp1=p1+dx1;
		  currp2=p2[i]+dx2;
		  //put back within the box if it hops out.
		  if(currp2>boxLength){
			dist=currp2-boxLength;
			currp2=currp2-2.0*dist;//factor of two, one dist puts you at wall, two dists puts you reflected back inside
			dx2=dx2-2.0*dist;
		  }
		  while(currp2<sigma){
			//reject move due to overlap
			//				dx1=sqrt(2.0*deltat*D1)*GaussV();
			dx2=sqrt(2.0*deltat*D2)*GaussV();
			//currp1=p1+dx1;
			currp2=p2[i]+dx2;
			//do not need to check boxLength here --the box should not be so small you could both hop within sigma and bounce off the wall in one step.
		  }
		  //complete the position update
		  //p1=p1+dx1;
		  p2[i]=p2[i]+dx2;
		  
		}//do not move a bound B molecule.
	  }//done looping over all B molecules.
      currBound+=aBound;//adds in either a 1 or a 0. in the end divide by total iterations.
	  pBoundOneTraj[it]+=aBound;//reset to zero each time.
	  pBoundSq[it]+=aBound*aBound;
	  pBoundHist[it]+=aBound;//add in 1 if it is bound, 0 if it is unbound
	  /*reset stopMoving*/
	  for(i=0;i<NB+1;i++){
		stopMoving[i]=false;//this only gets changed for molecules that dissociate.
	  }
      
    }//end iterating over time steps
    
	//calculate means over the later parts of each traj.
	double avg75=0;
	double avg50=0;
	int start75=Nit1*0.25;//includes 75% of the trajectory
	int start50=Nit1*0.5;
	for(i=start75;i<Nit1;i++){
	  avg75+=pBoundOneTraj[i];
	  if(i>start50)
		avg50+=pBoundOneTraj[i];
	}
	double dura1=Nit1-start75;
	double dura2=Nit1-start50;
	meansfile<<rep<<'\t'<<avg75/(1.0*dura1)<<'\t'<<avg50/(1.0*dura2)<<endl;

    if(rep%10==0){
      probfile.open(tname);
      probfile<<"time (us) "<<" pBound "<<" Stdev "<<endl; 
	  probfile<<0<<'\t'<<0<<endl;//at time zero, association probability is zero 
      for(i=1;i<Nit+1;i+=10){
		mean=pBoundHist[i]/(1.0*(rep+1));
		var=pBoundSq[i]/(1.0*(rep+1))-mean*mean;
		probfile <<i*deltat<<'\t'<<mean<<'\t'<<sqrt(var)<<endl;
      }
      probfile.close();
      
    }
  
    
  }//end looping over repetitions
  
  
  probfile.open(tname);
  probfile<<"time (us) "<<" pBound "<<" Stdev "<<endl; 
  probfile<<0<<'\t'<<0<<'\t'<<0<<endl;//at time zero, association probability is zero 
  for(i=1;i<Nit+1;i+=10){
	mean=pBoundHist[i]/(1.0*(rep+1));
	var=pBoundSq[i]/(1.0*(rep+1))-mean*mean;
	probfile <<i*deltat<<'\t'<<mean<<'\t'<<sqrt(var)<<endl;
  }
  probfile.close();
    

  cout <<"End Main, complete run "<<endl;
  
}//end main
