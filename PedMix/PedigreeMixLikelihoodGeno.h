//
//  PedigreeMixLikelihoodGeno.h
//
//
//  Created by Yufeng Wu on 3/25/15.
//
//

#ifndef ____PedigreeMixLikelihoodGeno__
#define ____PedigreeMixLikelihoodGeno__

#include "PedigreeMixLikelihood.h"

//************************************************************************
//  Pedigree likelihood computing for genotypes
//  here, genotype uses the exactly the same interface
//  the main difference is that

class PedigreeMixLikelihoodGeno
{
public:
    virtual int GetMaxNumPedigreeGen() const { return 2; };
};



#endif /* defined(____PedigreeMixLikelihoodGeno__) */
