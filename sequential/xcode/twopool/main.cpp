#include <math.h>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "nr.h"	
using namespace std;

//------------------------------------------------------------------------------------------------------------------------
#define importFolderName    "/Users/stefanhallermann/Desktop/github/in/"
#define exportFolderName    "/Users/stefanhallermann/Desktop/github/sequential/out/"
//-------------------------------------------------------------------------------------------------------------------------


//number of free parameters
//if this is changed, some stuff in main() and setrates() has to be changed
#define ndim 5

#define howManyFiles 1
#define howManyExperiments 6


//caWave
double caGlobalDt = 1000e-6;//1ms
#define caGlobalLength 5100    //5.1s at caWaveDt=1ms => 5 / 0.001 = 5.100
float caGlobal_300[caGlobalLength+1];
//global ca dynamics
//double cagl;
double timegl;

//global ca dynamics
double MicroCaPeak;
double MicroCaTau;

//pools and rates
double p_rel;             //initial pr of vesicles in pool1
double Ntot;
double k1, k2;
double b1, b2, b3;
double k10, k20;
double schluck1, schluck2;
double K_1, K_2;
double TSLIfrac;

#define statesN 3       //number of pools
Vec_DP states(statesN);

//define arrays for data
float minisize[100];    //size in miniature EPSC in pA for each experiment; is here set constant to 20 pA
int apNumber300;        //number of stimuli. Is set in load().
int recNumber;        //number of stimuli. Is set in load().
float apTimes300[1000]; //time of stimuli in s
//float «s[1000]; //time of stimuli in s
float apSamp300[1000];  //sampling (weighting) of stimuli for chi2 calculation
float apPostSyn300[1000];   //factor describing postsynaptic depression for each stimulus
float apAmp300[1000];       //amplitude of EPSC in units of released vesicles (= measured phasic EPSC amplitude/mini size of corresponding file)
//CAVE this value is contaminated by postsynaptic depression
//in export(), this value is multiplied by the mini size of the corresponding file resulting in units pA
float apAmpSim300[1000];    //amplitude of simulated EPSC in units of released vesicles (taking postsyn. depression into account)
//in output(), this value is multiplied by mini size of the corresponding file resulting in units pA
float apAmpStore300[1000][howManyExperiments];      //corresponding values stored for each experiment
float apAmpSimStore300[1000][howManyExperiments];   //corresponding values stored for each experiment
float export300[1000][5];                           //size of pool0, pool1, pool2, and release probability of pool1 and pool2 for each stimulus
float exportStore300[1000][howManyExperiments][5];  //corresponding values stored for each experiment

double chi2Sum1,chi2Sum2,chi2Sum3;

//define global arrays used in simplex and chi2
const int mpts = ndim+1;		
double p[mpts][ndim];	//that are e.g. 4 3-dimensional vectors
double psum[ndim];
double y[mpts];			//contains the chi values of the points of the simplex
double ptry[ndim];		//used in chi2()
double pBest[ndim];		//set in simplex()
double pBestStore[ndim][howManyExperiments];
double start[ndim];
double startMerk[ndim];
double startDelta[ndim];
double startDeltaMerk[ndim];
double chi2SumBest;     //set in simplex()
double chi2Sum1BestStore[howManyExperiments];
double chi2Sum2BestStore[howManyExperiments];
double chi2Sum3BestStore[howManyExperiments];
double chi2SumBestStore[howManyExperiments];

//parameters for repeated simplex and odeint
double reduceFactor = 1.5;
double ftol = 1e-5;
const int KMAX=5000;
DP eps=0.0001;          //precision
DP h1=500e-6;           //inital step guess within odeint
DP hmin=0;              //step should be larger than //1e-10;
int nrhs;               //counts function evaluations
int nbad,nok;           //count bad and good steps
DP dxsav;               //defining declarations
int kmax,kount;

// define general function
double sqr(double x){ 	return x*x;
}

/////////////////////////////////////////////////////////////////////////////////////////////
// load data into the global variables //////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
void load(){
	int i,ii,iii;
	FILE *fd;
	char filename[300];
	char buf[5];
	
    // aptimes + sampling (weighting) + postsyn. depression
	strcpy(filename,importFolderName);
    strcat(filename,"apTimesSamp300.txt");
	fd=fopen(filename,"r");
	i=0;
	while(fscanf(fd, "%f", &apTimes300[i]) != EOF)	{
		fscanf(fd, "%f", &apSamp300[i]);
        fscanf(fd, "%f", &apPostSyn300[i]);
        //fscanf(fd, "%f", &apFacil300[i]);
        //apTimes300[i] /= 1000.;//convert ms in s
		i++;
	}
	fclose(fd);
	apNumber300=i;

    // mini size
    for(ii=0;ii<howManyExperiments;ii++) {
        minisize[ii]=15.2; //assumung 1 nA mini size
    }

	// currents
	// -----------   300 Hz  ---------------
    for(i=1;i<=howManyFiles;i++) {
        strcpy(filename,importFolderName);
        sprintf(buf,"%d",i);// convert ii to string [buf]
        strcat(filename,buf);
        strcat(filename,"_300Hz.txt");
        cout << "I loaded " << filename << endl;
        fd=fopen(filename,"r");
        for(ii=0;ii<apNumber300;ii++){
            for(iii=0;iii<howManyExperiments;iii++) {
                fscanf(fd, "%f", &apAmpStore300[ii][iii]);
                apAmpStore300[ii][iii] /= -1.*minisize[iii];
            }
        }
        fclose(fd);
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////
// make Ca /////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
void makeCaWave()
{
    int i,caGlobalIndex;
    double time,expTime;

    caGlobalIndex = 0;
    while(caGlobalIndex <= caGlobalLength) {
        //caGlobal_300[caGlobalIndex] = (float)car;
        caGlobal_300[caGlobalIndex] = (float)0;
        caGlobalIndex++;
    }

//50Hz
    for(i=0;i<apNumber300;i++)    {
        time = apTimes300[i];
        expTime=0;
        caGlobalIndex = (int)(0.5 + time/caGlobalDt);
        //caGlobalIndex = 1+(int)(1+time/caGlobalDt);
        while(caGlobalIndex < caGlobalLength && expTime < 200*MicroCaTau) {
            caGlobal_300[caGlobalIndex] += (float)(MicroCaPeak * exp(-expTime/MicroCaTau));
            caGlobalIndex++;
            expTime += caGlobalDt;
        }
    }
    printf("the ca global waves were calculated!");
}

/////////////////////////////////////////////////////////////////////////////////////////////
// output in txt files    ///////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
void output()
{
	int i,ii;
	FILE *fd,*fdSim;
	FILE *fdrates;
	char filename[300];

	//export EPSC amplitude
	strcpy(filename,exportFolderName);
	strcat(filename,"apAmp300.txt");
	fd=fopen(filename,"w");
	strcpy(filename,exportFolderName);
	strcat(filename,"apAmp300Sim.txt");
	fdSim=fopen(filename,"w");
	//experimental data
	for(i = 0; i < apNumber300; i += 1){
		for(ii=0;ii<howManyExperiments;ii++)	fprintf(fd,"%f	",minisize[ii]*apAmpStore300[i][ii]);
		fprintf(fd,"\n");
	}
	//simulated data
	for(i = 0; i < apNumber300; i += 1){
		for(ii=0;ii<howManyExperiments;ii++)	fprintf(fdSim,"%f	",minisize[ii]*apAmpSimStore300[i][ii]);
		fprintf(fdSim,"\n");
	}
	fclose(fd);
	fclose(fdSim);
   
    //export pool sizes and release probability
	strcpy(filename,exportFolderName);	strcat(filename,"nR_300.txt");	fd=fopen(filename,"w");
	for(i = 0; i < apNumber300; i += 1){		
		for(ii=0;ii<howManyExperiments;ii++)	fprintf(fd,"%e	",exportStore300[i][ii][0]);
		fprintf(fd,"\n");
	}
	fclose(fd);
	strcpy(filename,exportFolderName);	strcat(filename,"n1_300.txt");	fd=fopen(filename,"w");
	for(i = 0; i < apNumber300; i += 1){		
		for(ii=0;ii<howManyExperiments;ii++)	fprintf(fd,"%e	",exportStore300[i][ii][1]);
		fprintf(fd,"\n");
	}
	fclose(fd);
	strcpy(filename,exportFolderName);	strcat(filename,"n2_300.txt");	fd=fopen(filename,"w");
	for(i = 0; i < apNumber300; i += 1){		
		for(ii=0;ii<howManyExperiments;ii++)	fprintf(fd,"%e	",exportStore300[i][ii][2]);
		fprintf(fd,"\n");
	}
	fclose(fd);
	strcpy(filename,exportFolderName);	strcat(filename,"p1_300.txt");	fd=fopen(filename,"w");
	for(i = 0; i < apNumber300; i += 1){		
		for(ii=0;ii<howManyExperiments;ii++)	fprintf(fd,"%e	",exportStore300[i][ii][3]);
		fprintf(fd,"\n");
	}
	fclose(fd);
	strcpy(filename,exportFolderName);	strcat(filename,"p2_300.txt");	fd=fopen(filename,"w");
	for(i = 0; i < apNumber300; i += 1){		
		for(ii=0;ii<howManyExperiments;ii++)	fprintf(fd,"%e	",exportStore300[i][ii][4]);
		fprintf(fd,"\n");
	}
	fclose(fd);

    
    //export timebase
    strcpy(filename,exportFolderName);
    strcat(filename,"time300.txt");
    fd=fopen(filename,"w");
    for(i = 0; i < apNumber300; i += 1){
        fprintf(fd,"%f    ",apTimes300[i]);
        fprintf(fd,"\n");
    }
    fclose(fd);

    
    //export caWave
    strcpy(filename,exportFolderName);
    strcat(filename,"caWave.txt");
    fd=fopen(filename,"w");
    for(i = 0; i < caGlobalLength; i += 1){
        fprintf(fd,"%f    ",1e6*caGlobal_300[i]);
        fprintf(fd,"\n");
    }
    fclose(fd);

    
    //export best rates and chi2
	strcpy(filename,exportFolderName);
	strcat(filename,"rates.txt");
	fdrates=fopen(filename,"w");
	for(ii=0;ii<howManyExperiments;ii++){
		for(i=0;i<ndim;i++)	fprintf(fdrates, "%e	", pBestStore[i][ii]);
        //fprintf(fdrates,"%e\t", chi2Sum1BestStore[ii]);
        //fprintf(fdrates,"%e\t", chi2Sum2BestStore[ii]);   //=0
        //fprintf(fdrates,"%e\t", chi2Sum3BestStore[ii]);   //=0
        fprintf(fdrates,"%e\n", chi2SumBestStore[ii]);
	}
	fclose(fdrates);

	printf("output into files done!");
}

// differential equation for pools
void derivs(const DP x, Vec_I_DP &yy, Vec_O_DP &dydx)        //use yy and not y to prevent confusion with global y
{
    double timeTmp;
    double microCalcium;
    double k1, k2;
    
    timeTmp = x + timegl;
    microCalcium = caGlobal_300[(int)(0.5 + timeTmp/caGlobalDt)]; //is global Ca2+ minus resting Ca2+

    //Michaelis-Menten according to eq. 42 of Lin et al. PNAS 2022
    k1 = (k10 + schluck1*microCalcium/(MicroCaPeak*MicroCaTau)) / (1 + microCalcium/K_1);
    k2 = (k20 + schluck2*microCalcium/(MicroCaPeak*MicroCaTau)) / (1 + microCalcium/K_2);


    //LS
    dydx[0] = -(b1+k1+k2)*yy[0] + (b2-k1)*yy[1] + (b3-k1)*yy[2] + Ntot*k1;
    //TS
    dydx[1] = k2*yy[0] - b2*yy[1];
	//TSL
    dydx[2] = -b3*yy[2];
}

void initPools(){
    //LS
    states[0] = Ntot*((k10*k20) / (k10*b2 + b1*b2 + k10*k20)) * b2/k20;
    //TS
    states[1] = Ntot*(k10*k20) / (k10*b2 + b1*b2 + k10*k20);
    //TSL
    states[2] = 0;
}

void fadvance(double timestep){
    kmax=KMAX;
    NR::odeint(states,0.0,timestep,eps,h1,hmin,nok,nbad,derivs,NR::rkqs);
}

/////////////////////////////////////////////////////////////////////////////////////////////
//  setRates      ///////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
void setRates() 
{
    //varried values
    k10 = ptry[0];
    b1 = ptry[1];
    k20 = ptry[2];
    b2 = ptry[3];
    p_rel = ptry[4];
    if(p_rel>1) {p_rel=1.;ptry[4]=1.;};//Pr cannot be >1
}

/////////////////////////////////////////////////////////////////////////////////////////////
//  chi2 ////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
double chi2()
{
	int i;
	double nextApTime;
    double chi2Sum;

    setRates();
	//define start conditions
	initPools();
	i=0;
	while(1)	{
		timegl = apTimes300[i];

        //number of released vesciles multipied by postsynaptic depression
        apAmpSim300[i] = (float)(states[1]*p_rel + states[2]*p_rel)*apPostSyn300[i];
        
        
        //for export: pool sizes before AP
        export300[i][0]=(float)states[0];
        export300[i][1]=(float)states[1];
        export300[i][2]=(float)states[2];
        export300[i][3]=(float)p_rel;
        export300[i][4]=(float)p_rel;

        states[1] = (1 - p_rel)*states[1];     //subtract released vesicles
        states[2] = (1 - p_rel)*states[2];     //subtract released vesicles

        //per AP TSL is increased by TSLIfrac*LS
        states[2] += TSLIfrac*states[0];
        states[0] -= TSLIfrac*states[0];
        
		if(i>=apNumber300-1) break;         //last AP
        nextApTime = apTimes300[i+1];
        //calculate the development of the number of vesicles in the pools until next AP
		fadvance(nextApTime - timegl);
		i++;
	}
	chi2Sum=0;
    chi2Sum1=0;
    chi2Sum2=0;
    chi2Sum3=0;
    for(i=0;i<apNumber300;i++)
		chi2Sum1 += apSamp300[i]*sqr(apAmp300[i] - apAmpSim300[i]);
    chi2Sum = chi2Sum1;
    return(chi2Sum);
}

void start_vectors()
{
	int i,j;
	for(i=0;i<mpts;i++)	{
		for(j=0;j<ndim;j++)
			p[i][j] = (start[j]);
	}
	for(j=0;j<ndim;j++)
		p[j+1][j] = (start[j] + startDelta[j]);
}

void get_psum()
{
	int i,j;
	double sum; 
	//p and psum are global

	for(j=0;j<ndim;j++)	{
		for (sum=0,i=0;i<mpts;i++)
			sum += p[i][j];
		psum[j] = sum;
	}
}

double amotry(int ihi, double fac)
{
	int j;
	double fac1,fac2,ytry;
	fac1 = (1.0-fac)/ndim;
	fac2 = fac1-fac;
	for(j=0;j<ndim;j++)
		ptry[j] = psum[j]*fac1 - p[ihi][j]*fac2;
	//check for broders at the moment only >0 
	for(j=0;j<ndim;j++)
		if(ptry[j] <= 0) return y[ihi];

	ytry = chi2();	//uses ptry
	if(ytry < y[ihi])	{
		y[ihi] = ytry;
		for(j=0;j<ndim;j++)		{
			psum[j] += ptry[j]-p[ihi][j];
			p[ihi][j] = ptry[j];
		}
	}
	return ytry;
}

void simplex()
{
	//see Press et al., Numerical recipes in C++ 2nd ed. p. 415
	const int NMAX = 5000;
	const double TINY=1.0e-10;
	int i,ihi,ilo,inhi,j,k;
	double rtol,ysave,ytry;
	int nfunk = 0;
	for(i=0;i<mpts;i++)
	{
		for(k=0;k<ndim;k++)
			ptry[k] = p[i][k];
		y[i] = chi2();	//uses ptry
	}
	get_psum();
	for(;;)
	{
		//find highest next-highest and lowest
		ilo=0;
		//ihi = y[0]>y[1] ? (inhi=1,0) : (inhi=0,1);
        if (y[0]>y[1])
        {
            ihi=0;
            inhi=1;
        }
        else
        {
            ihi=1;
            inhi=0;
        }
        for(i=0;i<mpts;i++)
		{
			if(y[i] <= y[ilo]) ilo=i;
			if(y[i] > y[ihi])
			{
				inhi=ihi;
				ihi=i;
			}
			else if(y[i] > y[inhi] && i!=ihi) inhi=i;
		}
		rtol = 2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo])+TINY);

		if(rtol < ftol)         //finish and store best values and best chi2 in global variables
		{
			for(k=0;k<ndim;k++)
				pBest[k] = p[ilo][k];	//used in output()
			chi2SumBest=y[ilo];
			break;
		}

		if(nfunk >= NMAX)
		{
			printf("\n\nNMAX exceeded!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1\n\n");
			for(k=0;k<ndim;k++)
				pBest[k] = p[ilo][k];	//used in output()
			break;
		}
		nfunk += 2;
		//new iteration
		ytry = amotry(ihi,-1.0);
		if(ytry <= y[ilo])	//better than best point, so try aditional extrapol.
		{				
			ytry = amotry(ihi,2.0);
		}
		else
			if(ytry  >= y[inhi])		//worse than 2nd highest so look for intermediate lower point
			{
				ysave = y[ihi];
				ytry = amotry(ihi,0.5);
				if(ytry >= ysave)		//cant get rid of highest point, so contract around the best (lowest) point
				{
					for(i=0;i<mpts;i++)
					{
						if(i != ilo)
						{
							for(j=0;j<ndim;j++)
								p[i][j] = ptry[j] = 0.5*(p[i][j]+p[ilo][j]);
							y[i] = chi2();	//uses ptry
						}
					}
					nfunk += ndim;
					get_psum();
				}
			}
			else --nfunk;
	}
	printf("\n\ndone!\n");
}

/////////////////////////////////////////////////////////////////////////////////////////////
// main /////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
	int k,m,ii,l,ll;
    int indexTmp;
	double chi2Merk,chi2RelChange;
    double    p_rel_start, k10_start, b1_start, k20_start, b2_start;
    
    //values BEFORE
    TSLIfrac = 0.160891; //per AP TSL is increased by TSLIfrac*LS
    p_rel_start = 0.544;
    k10_start = 2.965;
    b1_start = 2.618;
    k20_start = 0.398;
    b2_start = 1.090;
    
    b3 = 1/0.0775412;//cf. Tau_TSL in Lin et al. PNAS 2022
    MicroCaPeak = 2.09e-07;
    MicroCaTau = 0.1395482;

    schluck1 = 0.07320108;
    K_1 = 6.12E-07;
    schluck2 = 0.1051034;
    K_2 = 10;


    //define Ntot according to eq. 4 and multiply with the average amplitude of the before ("experiment" 1) and the after traces ("experiment" 22)
    Ntot = (((-255.0042286 + -162.2618905)/2)/-15.2)/(p_rel_start*((k10_start*k20_start) / (k10_start*b2_start + b1_start*b2_start + k10_start*k20_start)));
     
    //define start values that should be optimized
    indexTmp = 0;
    startMerk[indexTmp] = k10_start;
    startDeltaMerk[indexTmp] = 0.1 * startMerk[indexTmp];

    indexTmp += 1;
    startMerk[indexTmp] = b1_start;
    startDeltaMerk[indexTmp] = 0.1 * startMerk[indexTmp];

    indexTmp += 1;
    startMerk[indexTmp] = k20_start;
    startDeltaMerk[indexTmp] = 0.1 * startMerk[indexTmp];

    indexTmp += 1;
    startMerk[indexTmp] = b2_start;
    startDeltaMerk[indexTmp] = 0.1 * startMerk[indexTmp];

    indexTmp += 1;
    startMerk[indexTmp] = p_rel_start;
    startDeltaMerk[indexTmp] = 0.1 * startMerk[indexTmp];
    
	load();
    makeCaWave();

    for(ii=0;ii<howManyExperiments;ii++){
		printf("\n\nNow file: %d\n",ii);
		for(k=0;k<ndim;k++){
            pBest[k] = startMerk[k];
			startDelta[k] = startDeltaMerk[k];
		}
		//load Experiment into currently used array
		for(l=0;l<apNumber300;l++) apAmp300[l] = apAmpStore300[l][ii];
        // minimazation
		chi2Merk=1e100;
		chi2RelChange = 1;
		m=0;
        while(m<10 && chi2RelChange > 0.0001){          //restart minimazation with best values + deltas and reduce deltas -> chance of jumping out of local minimum
            for(k=0;k<ndim;k++)    {
                start[k] = pBest[k];
                if(m>0) pBest[k] += startDelta[k];      //do not add deltas when done for the first time
               startDelta[k] /= reduceFactor;
            }
            start_vectors();
            
            simplex();   //will set pBest and chi2SumBest
            
			printf("simplex %d. th time: chi2 = %f\n",m+1,chi2SumBest);
			for(k=0;k<ndim;k++)	{
				printf("%e	",pBest[k]);
			}
			printf("\n");
			chi2RelChange = fabs(chi2Merk-chi2SumBest)/chi2SumBest;
			chi2Merk = chi2SumBest;
			m++;
		}
		//do it again with best paramters
		for(k=0;k<ndim;k++)
			ptry[k] = pBest[k];
		chi2SumBestStore[ii] = chi2();      //chi2() also generate apAmpSim with the best parameters
        chi2Sum1BestStore[ii] = chi2Sum1;
        chi2Sum2BestStore[ii] = chi2Sum2;
        chi2Sum3BestStore[ii] = chi2Sum3;
        //store current for exort in outStartAndEnd()
		for(l=0;l<apNumber300;l++) apAmpSimStore300[l][ii] = apAmpSim300[l];
        //store pool size and pr for exort in outStartAndEnd()
		for(l=0;l<apNumber300;l++) for(ll=0;ll<5;ll++) exportStore300[l][ii][ll] = export300[l][ll];
        //store parameters for exort in outStartAndEnd()
		for(k=0;k<ndim;k++)	pBestStore[k][ii] = pBest[k];
	}
	output();
}
