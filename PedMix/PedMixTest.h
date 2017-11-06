//
//  PedMixTest.h
//
//
//  Created by Yufeng Wu on 12/5/14.
//  Test pedigree likelihood functions
//

#ifndef ____PedMixTest__
#define ____PedMixTest__

#include "PedigreeMixLikelihood.h"
#include "UtilsNumerical.h"
#include <stdio.h>
#include "lbfgs.h"

int TestLikelihoodForHaps( const vector<vector<PedigreeMixHaplotype> > &listHaps, PopMixingModel &modelPedMix, int numMixGens, double mixratio, vector<double> Par );
static int progress(void *instance, const lbfgsfloatval_t *x, const lbfgsfloatval_t *g, const lbfgsfloatval_t fx, const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step, int n, int k, int ls);
static lbfgsfloatval_t evaluate(void *instance,const lbfgsfloatval_t *x,lbfgsfloatval_t *g,const int n,const lbfgsfloatval_t step);
static double GetDeltaIncBFGS(double valCurr, double DEF_BFGS_STEP_SIZE);

#endif /* defined(____PedMixTest__) */
