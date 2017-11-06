//
//  PedMixTest.cpp
//
//
//  Created by Yufeng Wu on 12/5/14.
//
//

#include "PedMixTest.h"
#include <cmath>
#include <iomanip>
#include "Utils4.h"
#include <vector>

using namespace std;

struct SValue{
	vector< vector<PedigreeMixHaplotype> > listHaps;
	PopMixingModel* modelPedMix;
	const PedigreeMixLikelihood* pedMixCalc;
	int geno_index;
	double DEF_BFGS_STEP_SIZE;
	double pho;
	int ll;
	double probswitch;
};


static lbfgsfloatval_t evaluate(void *instance,const lbfgsfloatval_t *x,lbfgsfloatval_t *g,const int n,const lbfgsfloatval_t step)
{
	SValue *sa;
	sa=(SValue *)instance;
	//cout<<"haplotype size is "<<sa->listHaps.size()<<endl;
	//cout<<"pedMixCalc :"<<sa->pedMixCalc->GetMaxNumPedigreeGen()<<endl;
	//cout<<"model"<<sa->modelPedMix->GetNumSites(0)<<endl;

	//calculate fx
	lbfgsfloatval_t fx;
	double phopara=sa->pho*sa->ll;//**********************param to change
	//int ll=sa->ll;//**********************param to change
	//double probswitch=sa->probswitch;
	//cout<<"phopara "<<phopara<<endl;
	//cout<<ll;
	int hap_index1=2*sa->geno_index;
	int hap_index2=2*sa->geno_index+1;
	PopMixingModel mp=*(sa->modelPedMix);
	vector<double> llh(sa->listHaps.size());

	vector<double> PP;
	for (int i=0;i<n;i++)
	{
		PP.push_back(exp(x[i])/(1+exp(x[i])));
	}

	cout<<"PP vector----[";
	for (int i=0;i<n;i++)
	{
		cout<<PP[i]<<" ";
	}
	cout<<"]"<<endl;
#pragma omp parallel num_threads(22)
    {
#pragma omp for ordered schedule(dynamic)
	for (int i=0;i<(int)sa->listHaps.size();i++) //enumerate all the locus
	{
                long tstart2 = GetCurrentTimeTick();
		llh[i]=sa->pedMixCalc->Compute( mp, i, sa->listHaps[i][hap_index1], sa->listHaps[i][hap_index2], PP,phopara,sa->ll,sa->probswitch );
                cout<<"haplotype "<<(i+1)<<" takes "<<GetElapseTime( tstart2 )<<" seconds"<<endl;
	}
    }
	fx=-GetLogSumOfLogs( llh );
	cout<<"computed likelihood is : "<<fx<<endl;


	vector<double> PP_tmp;
	//calculate gradient
    for(int i=0; i<n; ++i)
    {
        double deltaInc = GetDeltaIncBFGS( x[i], sa->DEF_BFGS_STEP_SIZE );
        //cout<<"increase of x["<<i<<"] is:"<<deltaInc<<endl;
        PP_tmp.clear();
        PP_tmp=PP;
        PP_tmp[i]=exp(x[i]+deltaInc)/(1+exp(x[i]+deltaInc));
        //print new PP
/*
        cout<<"test new PP is :[";
        for (int k=0;k<n;k++)
        {
        	cout<<PP_tmp[k]<<" ";
        }
        cout<<"]"<<endl;
*/
        vector<double> llh2(sa->listHaps.size());
#pragma omp parallel num_threads(22)
        {
#pragma omp for ordered schedule(dynamic)
        for (int j=0;j<(int)sa->listHaps.size();j++)
        {
                long tstart2 = GetCurrentTimeTick();
        	llh2[j]=sa->pedMixCalc->Compute( mp, j, sa->listHaps[j][hap_index1], sa->listHaps[j][hap_index2], PP_tmp,phopara,sa->ll,sa->probswitch );
                //cout<<"haplotype "<<(j+1)<<" takes "<<GetElapseTime( tstart2 )<<" seconds"<<endl;
        }
        }
        double fx2=-GetLogSumOfLogs( llh2 );
        //cout<<"the new likelihood is :"<<fx2<<endl;
        double gradientVal = ( fx2-fx )/deltaInc;
        g[i] = gradientVal;
    }
	return fx;
}
/*
int test1(void* par)
{
	SValue* sa=(SValue*)par;
	for(int i=0;i<sa->listHaps.size();i++)
	{
		cout<<sa->listHaps[i].size()<<endl;
		for(int j=0;j<sa->listHaps[i].size();j++)
		{
			PedigreeMixHaplotype pl=sa->listHaps[i][j];
			cout<<pl.Dump()<<" ";
		}
		cout<<endl;
	}
	cout<<sa->listHaps.size()<<endl;

	vector<int> a;
	a.push_back(1);
	vector<int> b;
	b=a;

}
*/

//
int TestLikelihoodForHaps( const vector<vector<PedigreeMixHaplotype> > &listHaps, PopMixingModel &modelPedMix, int numMixGens, double mixratio, vector<double> Par  )
{

	PedigreeMixLikelihood pedMixCalc;
    modelPedMix.SetMixRatio( mixratio );

    if( pedMixCalc.GetMaxNumPedigreeGen() < GetNumofPerfectPedGens()  )
    {
        cout << "The number of generations in pedigree is too large. Stop." << endl;
        exit(1);
    }

    modelPedMix.SetNumMixGens( numMixGens );
    //pedMixCalc.UseRecursive(true);
    //bool fRecursive = false;
    bool fRecursive = true;
    pedMixCalc.UseRecursive(fRecursive);
    bool fLogMode = true;
    pedMixCalc.SetLogMode( fLogMode );
    //cout << "The model to be used: ";
    //modelPedMix.Dump();

    for (int geno_index=0;geno_index<listHaps[0].size()/2;geno_index++)
    {
    //start BFGS main function
    int i, ret = 0;
    int num_anc=pow(2,GetNumofPerfectPedGens());
    int N=2*num_anc;
    lbfgsfloatval_t fx;
    lbfgsfloatval_t *x = lbfgs_malloc(N);
    lbfgs_parameter_t param;
    if (x == NULL) {
        printf("ERROR: Failed to allocate a memory block for variables.\n");
        return 1;
    }
    //initialize PP variables
    double basev=5/Par[2];
    //cout<<"base value is"<<basev<<endl;
    for (i=0;i<N;++i)
    {
    	x[i]=log((i*0.5+1)*basev/(1-(i*0.5+1)*basev));     //*********************param to change
    }
    //initialize parameters for L-BFGS optimization
    lbfgs_parameter_init(&param);
    //start L-BFGS optimization progress

    //pass data****************************************************************
    SValue* sa=new SValue;
    sa->listHaps=listHaps;
    sa->modelPedMix=&modelPedMix;
    sa->pedMixCalc= &pedMixCalc;
    sa->geno_index=geno_index;
    sa->probswitch=Par[0];
    sa->pho=Par[1];
    sa->ll=(int)Par[2];
    sa->DEF_BFGS_STEP_SIZE=Par[3];
    ret = lbfgs(N, x, &fx, evaluate, progress, (void *)sa, &param);
    //pass data****************************************************************

    //Report result
    cout<<"--*important info*--the maximum likelihood is  [";
    cout<<setprecision(10)<<fx<<"]"<<endl;
    //####### remap PP to [0, 1]
    vector<double> PP;
    cout<<"--*important info*--PP vector is   [";
    for (int idx=0;idx<N;idx++)
    {
    	PP.push_back(exp(x[idx])/(1+exp(x[idx])));
    	cout<<PP[idx]<<" ";
    }
    cout<<"]"<<endl;
    double ratiov=0.0;
    cout<<"--*important info*--the admixture proportions of ancestral population A is  [";
    for (int anc=0;anc<num_anc;anc++)
    {
    	ratiov=100*PP[2*anc+1]/(PP[2*anc]+PP[2*anc+1]);
    	cout<<ratiov<<"%  ";
    }
    cout<<"]"<<endl;
    //clear L-BFGS
    lbfgs_free(x);
    }
    return 0;

}

static double GetDeltaIncBFGS(double valCurr, double DEF_BFGS_STEP_SIZE)
{
    //cout<<"step size is "<<DEF_BFGS_STEP_SIZE<<endl;
    //const double DEF_BFGS_STEP_SIZE = 0.00000001;    //*********************************param to change
    double trueval=exp(valCurr)/(1+exp(valCurr));
    return log((trueval+DEF_BFGS_STEP_SIZE)/(1-trueval-DEF_BFGS_STEP_SIZE)) - valCurr;
}

static int progress(
    void *instance,
    const lbfgsfloatval_t *x,
    const lbfgsfloatval_t *g,
    const lbfgsfloatval_t fx,
    const lbfgsfloatval_t xnorm,
    const lbfgsfloatval_t gnorm,
    const lbfgsfloatval_t step,
    int n,
    int k,
    int ls
    )
{
    printf("Iteration %d:\n", k);
    printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
    printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
    printf("\n");
    return 0;
}


