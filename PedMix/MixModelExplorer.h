//
//  MixModelExplorer.h
//
//
//  Created by Yufeng Wu on 12/8/14.
//
//

#ifndef ____MixModelExplorer__
#define ____MixModelExplorer__

#include "PedigreeMixLikelihood.h"
#include "UtilsNumerical.h"

//************************************************************************
// Search for the best model given haplotypes and the initial model

class PedigreeMixModelExplorer : public NumericalAlgoUtils
{
public:
    PedigreeMixModelExplorer( PopMixingModel &modelToExp, const vector<vector<PedigreeMixHaplotype> > &listHapsInParam );
    virtual bool Explore();
    virtual double EvaluateAt(double pt, void *pParam);
    double GetLogLikeli() const { return logLikeliCur; }

protected:
    double CalcHapLikelihood() const;
    void SetLogLikeli(double logl) { logLikeliCur = logl; }
    PopMixingModel &GetModel() { return modelCur; }
    const PopMixingModel &GetModel() const { return modelCur; }

private:
    PopMixingModel &modelCur;
    const vector<vector<PedigreeMixHaplotype> > &listHapsIn;
    double logLikeliCur;
    bool fLogMode;
};

//************************************************************************
// Opt over mixing ratio

class PedigreeMixModelExplorerMixRatio : public PedigreeMixModelExplorer
{
public:
    PedigreeMixModelExplorerMixRatio(  PopMixingModel &modelToExp, const vector<vector<PedigreeMixHaplotype> > &listHapsInParam  );
    virtual bool Explore();
    virtual double EvaluateAt(double pt, void *pParam);

};

//************************************************************************
// Opt over mixing time

class PedigreeMixModelExplorerMixGen : public PedigreeMixModelExplorer
{
public:
    PedigreeMixModelExplorerMixGen(  PopMixingModel &modelToExp, const vector<vector<PedigreeMixHaplotype> > &listHapsInParam  );
    virtual bool Explore();
    virtual double EvaluateAt(double pt, void *pParam);

};



#endif /* defined(____MixModelExplorer__) */
