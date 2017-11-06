//
//  PedigreeMixIO.h
//
//
//  Created by Yufeng Wu on 12/7/14.
//  Handle I/O of pedigree-based admixture
//
//

#ifndef ____PedigreeMixIO__
#define ____PedigreeMixIO__

#include "PedigreeMixLikelihood.h"

//************************************************************************
// Input processing

void ReadPedMixInputFromFile(const char *fileName, vector<vector<PedigreeMixHaplotype> > &listHaps, PopMixingModel &modelPedMix);
vector<double> ReadParFromFile(const char *ParFile);

#endif /* defined(____PedigreeMixIO__) */
