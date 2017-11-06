#include "UtilsNumerical.h"
#include "Utils3.h"
#include <iostream>
#include <iomanip>

// Some matrix utilities
// YW: seem to be some risk of memory issue: not freeing???

#if 0		// move implementaiton to header file
inline int* dec2binarr(long n, int dim)
{
    // note: res[dim] will save the sum res[0]+...+res[dim-1]
    int* res = (int*)calloc(dim + 1, sizeof(int));
    int pos = dim - 1;

    // note: this will crash if dim < log_2(n)...
    while (n > 0)
    {
        res[pos] = n % 2;
        res[dim] += res[pos];
        n = n / 2; // integer division
        pos--;
    }

    return res;
}

template<class T>
T MatrixPermanent(const vector<T>& A, int n)
{
	// expects n by n matrix encoded as vector
    T sum = 0;
    T rowsumprod, rowsum;
    int* chi = new int[n + 1];
    double C = (double)pow((double)2, n);

    // loop all 2^n submatrices of A
    for (int k = 1; k < C; k++)
    {
        rowsumprod = 1;
        chi = dec2binarr(k, n); // characteristic vector

        // loop columns of submatrix #k
        for (int m = 0; m < n; m++)
        {
            rowsum = 0;

            // loop rows and compute rowsum
            for (int p = 0; p < n; p++)
                rowsum += chi[p] * A[m * n + p];

            // update product of rowsums
            rowsumprod *= rowsum;

            // (optional -- use for sparse matrices)
            // if (rowsumprod == 0) break;
        }

        sum += (T)pow((double)-1, n - chi[n]) * rowsumprod;
    }

	//delete [] chi;

    return sum;
}

#endif


//////////////////////////////////////////////////////////////////////////////////////////////////

double NumericalAlgoUtils :: Func1DMinBrent( double ax, double bx, double cx, double tol, double *xmin )
{
//cout << "Func1DMinBrent: " << ", [" << ax << ", " << bx << ", " << cx  << ", tol  " << tol << "], \n";
	// YW: this function is based Numerical Receipe in C book.
	// search for best 1 D function (in this case, the likelihood) using Brent's method
	//Given a function f, and given a bracketing triplet of abscissas ax, bx, cx (such that bx is
	//between ax and cx, and f(bx) is less than both f(ax) and f(cx)), this routine isolates
	//the minimum to a fractional precision of about tol using Brent’s method. The abscissa of
	//the minimum is returned as xmin, and the minimum function value is returned as brent, the
	//returned function value.
	#define ITMAX 100
	#define CGOLD 0.3819660
	#define ZEPS 1.0e-10
	//Here ITMAX is the maximum allowed number of iterations; CGOLD is the golden ratio; ZEPS is
	//a small number that protects against trying to achieve fractional accuracy for a minimum that
	//happens to be exactly zero.
	#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
	#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

//cout << "Func1DMinBrent: ax=" << ax << ", bx=" << bx << ", cx=" << cx << ", tol=" << tol << endl;
	int iter;
	double a,b,d=0.0,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	double e=0.0;   // This will be the distance moved on the step before last.
	a=(ax < cx ? ax : cx);  //a and b must be in ascending order,
	b=(ax > cx ? ax : cx); // but input abscissas need not be.
	x=w=v=bx;				//Initializations...
	fw=fv=fx= 0; //EvaluateAt( x, NULL ); change later

#if 0
	// in case f(a) < f(b) < f(c), stop
	double fa1 = EvaluateAt( a, NULL );
	double fb1 = EvaluateAt( b, NULL );
cout << "fa1 = " << fa1 << " for a = " << a << ", fb1= " << fb1 << " for b = " << b << ", fx = " << fx << ", x = " << x << endl;
	if( IsSignificantlyLarge(fa1, fx) == false || IsSignificantlyLarge( fb1, fx) == false )
	{
		//  take the minimum
		*xmin = a;
		double res1 = fa1;
		if( fa1 > fb1)
		{
			*xmin = b;
			res1 = fb1;
		}
		return res1;
	}
#endif

	for (iter=1;iter<=ITMAX;iter++)
	{ //Main program loop.
//cout << "iteration " << iter << endl;
		xm=0.5*(a+b);
		tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
		if (fabs(x-xm) <= (tol2-0.5*(b-a)))
		{ //Test for done here.
			*xmin=x;
//cout << "x = " << x << ", xm = " << xm << ", tol2 = " << tol2 << ", b = " << b << ", a = " << a << endl;
//cout << "Here: STOP EARLY\n";
			return fx;
		}
		if (fabs(e) > tol1)
		{ // Construct a trial parabolic fit.
//cout << "here...\n";
			r=(x-w)*(fx-fv);
			q=(x-v)*(fx-fw);
			p=(x-v)*q-(x-w)*r;
			q=2.0*(q-r);
			if (q > 0.0) p = -p;
			q=fabs(q);
			etemp=e;
			e=d;
			if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
				d=CGOLD*(e=(x >= xm ? a-x : b-x));
			//The above conditions determine the acceptability of the parabolic fit. Here we
			//take the golden section step into the larger of the two segments.
			else
			{
				d=p/q; //Take the parabolic step.
				u=x+d;
				if (u-a < tol2 || b-u < tol2)
					d=SIGN(tol1,xm-x);
			}
		}
		else
		{
//cout << "here2\n";
			d=CGOLD*(e=(x >= xm ? a-x : b-x));
		}
		u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
		fu=0; //EvaluateAt(u, NULL ); change later
		//This is the one function evaluation per iteration.
		if (fu <= fx) { //Now decide what to do with our func
			if(u >= x) a=x; else b=x; //tion evaluation.
			SHFT(v,w,x,u)                       //Housekeeping follows:
			SHFT(fv,fw,fx,fu)
		}
		else
		{
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x)
			{
				v=w;
				w=u;
				fv=fw;
				fw=fu;
			}
			else if (fu <= fv || v == x || v == w)
			{
				v=u;
				fv=fu;
			}
		} //Done with housekeeping. Back for
//cout << "** -fx = " << -1.0*fx << endl;
	} //another iteration.
	//YW_ASSERT_INFO(false, "Too many iterations in brent");
    cout << "WARNING: Too many iterations in brent.\n";
	*xmin=x; //Never get here.
	return fx;
}

/*
double NumericalAlgoUtils :: EvaluateAt()
{

}
*/
bool NumericalAlgoUtils :: IsSignificantlyLarge(double v1, double v2) const
{
	// is v1 significantly larger than v2 (i.e. larger by some threshold)?
	// by default, the computed values are in log-space, and thus we ask then to differ by at least 5%
	const double thresDef = log(1.05);
	return v1 >= v2+ thresDef;
}


double BrentMethod :: Func1DMinBrent( const PedigreeMixLikelihood &pedMixCalc, const PopMixingModel &modelPedMix, const vector<vector<PedigreeMixHaplotype> > &haplotype, int geno_index, vector<double> PP, double ax, int anc_ind, double cx, double tol, double *xmin,double phopara,int ll ,double probswitch)
{
//cout << "Func1DMinBrent: " << ", [" << ax << ", " << bx << ", " << cx  << ", tol  " << tol << "], \n";
	// YW: this function is based Numerical Receipe in C book.
	// search for best 1 D function (in this case, the likelihood) using Brent's method
	//Given a function f, and given a bracketing triplet of abscissas ax, bx, cx (such that bx is
	//between ax and cx, and f(bx) is less than both f(ax) and f(cx)), this routine isolates
	//the minimum to a fractional precision of about tol using Brent’s method. The abscissa of
	//the minimum is returned as xmin, and the minimum function value is returned as brent, the
	//returned function value.
	#define ITMAX 100
	#define CGOLD 0.3819660
	#define ZEPS 1.0e-10
	//Here ITMAX is the maximum allowed number of iterations; CGOLD is the golden ratio; ZEPS is
	//a small number that protects against trying to achieve fractional accuracy for a minimum that
	//happens to be exactly zero.
	#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
	#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

//cout << "Func1DMinBrent: ax=" << ax << ", bx=" << bx << ", cx=" << cx << ", tol=" << tol << endl;
	double bx=PP[anc_ind];
	int iter;
	double a,b,d=0.0,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	double e=0.0;   // This will be the distance moved on the step before last.
	a=(ax < cx ? ax : cx);  //a and b must be in ascending order,
	b=(ax > cx ? ax : cx); // but input abscissas need not be.
	x=w=v=bx;				//Initializations...
	fw=fv=fx= EvaluateAt(x, pedMixCalc, modelPedMix, haplotype, geno_index, PP, anc_ind, phopara,ll,probswitch);


	for (iter=1;iter<=ITMAX;iter++)
	{ //Main program loop.
//cout << "iteration " << iter << endl;
		xm=0.5*(a+b);
		tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
		if (fabs(x-xm) <= (tol2-0.5*(b-a)))
		{ //Test for done here.
			//cout << "x = " << x << ", xm = " << xm << ", tol2 = " << tol2 << ", b = " << b << ", a = " << a << endl;
			//cout<<xmin;
			*xmin=x;
//cout << "x = " << x << ", xm = " << xm << ", tol2 = " << tol2 << ", b = " << b << ", a = " << a << endl;
//cout << "Here: STOP EARLY\n";
			return fx;
		}
		if (fabs(e) > tol1)
		{ // Construct a trial parabolic fit.
//cout << "here...\n";
			r=(x-w)*(fx-fv);
			q=(x-v)*(fx-fw);
			p=(x-v)*q-(x-w)*r;
			q=2.0*(q-r);
			if (q > 0.0) p = -p;
			q=fabs(q);
			etemp=e;
			e=d;
			if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
				d=CGOLD*(e=(x >= xm ? a-x : b-x));
			//The above conditions determine the acceptability of the parabolic fit. Here we
			//take the golden section step into the larger of the two segments.
			else
			{
				d=p/q; //Take the parabolic step.
				u=x+d;
				if (u-a < tol2 || b-u < tol2)
					d=SIGN(tol1,xm-x);
			}
		}
		else
		{
//cout << "here2\n";
			d=CGOLD*(e=(x >= xm ? a-x : b-x));
		}
		u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
		fu=EvaluateAt(u, pedMixCalc, modelPedMix, haplotype, geno_index, PP, anc_ind, phopara,ll,probswitch);
		//This is the one function evaluation per iteration.
		if (fu <= fx) { //Now decide what to do with our func
			if(u >= x) a=x; else b=x; //tion evaluation.
			SHFT(v,w,x,u)                       //Housekeeping follows:
			SHFT(fv,fw,fx,fu)
		}
		else
		{
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x)
			{
				v=w;
				w=u;
				fv=fw;
				fw=fu;
			}
			else if (fu <= fv || v == x || v == w)
			{
				v=u;
				fv=fu;
			}
		} //Done with housekeeping. Back for
//cout << "** -fx = " << -1.0*fx << endl;
	} //another iteration.
	//YW_ASSERT_INFO(false, "Too many iterations in brent");
    cout << "WARNING: Too many iterations in brent.\n";
	*xmin=x; //Never get here.
	return fx;
}

double BrentMethod :: EvaluateAt(double pt, const PedigreeMixLikelihood &pedMixCalc, const PopMixingModel &modelPedMix,const vector<vector<PedigreeMixHaplotype> > &haplotype, int geno_index, vector<double> PP, int anc_ind,double phopara,int ll,double probswitch)
{
	//long tstart = GetCurrentTimeTick();
	PP[anc_ind]=pt;

	vector<double> llh(haplotype.size());
	//cout<<"size of llh is: "<<haplotype.size()<<endl;
#pragma omp parallel num_threads(20)
	{
#pragma omp for
	for (int i=0;i<(int)haplotype.size();i++) //enumerate all the locus
	{
		//cout<<haplotype.size()<<"----"<<i<<endl;
		/*
		cout<<"Now consider haplotype: ";
		for (int j=0;j<4;j++)
		{
			cout<<haplotype[i][hap_index].GetAlleleAt(j)<<" ";
		}
		cout<<endl;
		*/
		int hap_index1=2*geno_index;
		int hap_index2=2*geno_index+1;
		llh[i]=pedMixCalc.Compute( modelPedMix, i, haplotype[i][hap_index1], haplotype[i][hap_index2], PP,phopara,ll,probswitch );
		//llh.push_back(pedMixCalc.Compute( modelPedMix, i, haplotype[i][hap_index1], haplotype[i][hap_index2], PP,phopara,ll ));//change later; specify locus site
	}
	}
/*
	cout<<"The list of probs is :";
	for (int j=0;j<llh.size();j++)
	{
		cout<<llh[j]<<" ";
	}
	cout<<endl<<"Sum of log prob is:";
	cout<<GetLogSumOfLogs( llh )<<endl;
*/
/*
		cout<<"For this setting of PP: ";
		for (int k=0;k<PP.size();k++)
		{
			cout<<PP[k]<<"-";
		}

		cout<<endl;

	cout<<setprecision(10)<<-GetLogSumOfLogs( llh )<<endl;
*/

	//cout<<-GetLogSumOfLogs( llh )<<endl;
	//cout << "Elapsed time = " << GetElapseTime( tstart ) << " seconds." << endl;
	return -GetLogSumOfLogs( llh );
	//return -ll because BrentMethod search for the smallest value

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double RoundDoubleValTo(double val, int numFractionDigits)
{
    // numFractiondigits: how many digits after . we want to keep
    YW_ASSERT_INFO(numFractionDigits >= 0, "numFracDigits:; must be positive");
    double ratioInc = pow(10.0, numFractionDigits);
    return round(val*ratioInc)/ratioInc;
}

