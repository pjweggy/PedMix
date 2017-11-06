//
//  MixModelExplorer.cpp
//
//
//  Created by Yufeng Wu on 12/8/14.
//
//

#include "MixModelExplorer.h"
#include <cmath>
using namespace std;

//************************************************************************
// Search for the best model given haplotypes and the initial model
const int MIX_GEN_MIN = 4;
const int MIX_GEN_MAX = 40;

PedigreeMixModelExplorer :: PedigreeMixModelExplorer( PopMixingModel &modelToExp, const vector<vector<PedigreeMixHaplotype> > &listHapsInParam ) : modelCur(modelToExp), listHapsIn(listHapsInParam), logLikeliCur(MAX_NEG_DOUBLE_VAL), fLogMode(true)
{
    //
}

bool PedigreeMixModelExplorer :: Explore()
{
    // return true if a significantly better model is found, false otherwise
    // here we scale by the max length
    // iterate over optimization of generaiton number and mix ratio

    bool fGenOptStarted = false;
    //bool fMROptStarted = false;

    while(true)
    {
cout << "Before optimizing, loglikelihood: " << this->GetLogLikeli() << endl;
        // terminate whenever optimization is
        PedigreeMixModelExplorerMixRatio mrOpt( modelCur, listHapsIn );
        mrOpt.SetLogLikeli( this->GetLogLikeli() );
        bool fUpdatedMR = mrOpt.Explore();
        //fMROptStarted = true;

        if( fUpdatedMR == false && fGenOptStarted == true )
        {
            // that is it
cout << "MR: not updated\n";
            break;
        }
cout << "mrOpt.GetLogLikeli(): " << mrOpt.GetLogLikeli() << endl;
        this->SetLogLikeli( mrOpt.GetLogLikeli() );

cout << "After optimizing mixing ration, likelihood is " << this->GetLogLikeli() << ", model: ";
modelCur.DumpModelOnly();

        // otherwise opt wrt generation
        PedigreeMixModelExplorerMixGen mgenOpt( modelCur, listHapsIn );
        mgenOpt.SetLogLikeli( this->GetLogLikeli() );
        bool fUpdatedGen = mgenOpt.Explore();
        fGenOptStarted = true;

        if( fUpdatedGen == false )
        {
            // that is it
cout << "Gen: not updated\n";
            break;
        }
cout << "mgenOpt.GetLogLikeli(): " << mgenOpt.GetLogLikeli() << endl;
        this->SetLogLikeli( mgenOpt.GetLogLikeli() );
cout << "After optimizing mixing generation, likelihood is " << this->GetLogLikeli() << ", model: ";
modelCur.DumpModelOnly();
    }

    return true;


#if 0
    // repeat until no change is possible
    //while(true)
    //{
        // explore the
        int ngenCurModel = modelCur.GetNumMixGens();
        //double mrCur = modelCur.GetMixRatio(1);
        double likeliCur = GetLogLikeli();

        double ngOpt;
        double likeliMax = -1.0*Func1DMinBrent((1.0*MIX_GEN_MIN)/MIX_GEN_MAX, ((double)ngenCurModel)/MIX_GEN_MAX, 1.0, 1.0/MIX_GEN_MAX, &ngOpt );

//cout << "In PedigreeMixModelExplorerMixRatio :: Explore: maximum likelihood = " << likeliMax << ", original likelihood = " << GetLogLikeli() << endl;

        //YW_ASSERT_INFO( likeliMax+tolnum >= loglikeliCurBest, "Wrong: Brent's result does not make sense" );
        if( IsSignificantlyLarge(likeliMax, likeliCur ) )
        {
            //
            SetLogLikeli( likeliMax );
            //speciesTreeBest.SetBranchLen(branch, szCurr);
            int ngOptInt = (int)(ngOpt*MIX_GEN_MAX);
            modelCur.SetNumMixGens(ngOptInt);
            // continue
            return true;
        }
        else
        {
            // restore to original settings, do not change branch length since no increase is found
            SetLogLikeli( likeliCur );
            modelCur.SetNumMixGens(ngenCurModel);
            //modelCur.SetMixRatio(mrCur);
            //break;
            return false;
        }
    //}
    //return true;
#endif
}

double PedigreeMixModelExplorer :: EvaluateAt(double pt, void *pParam)
{
//cout << "PedigreeMixModelExplorer: Evaluatiing at: " << pt << endl;
    // the pt refers to the current generation num
    YW_ASSERT_INFO( pt >= 1.0/MIX_GEN_MAX, "Cannot be too small" );
    int numgen = (int) (pt*MIX_GEN_MAX);
    modelCur.SetNumMixGens(numgen);

#if 0  // for now fix the mixing ratio
    // also optimize over mixture ratio
    PedigreeMixModelExplorerMixRatio mrOpt( modelCur, listHapsIn );
    mrOpt.SetLogLikeli( this->GetLogLikeli() );
    mrOpt.Explore();
//cout << "During PedigreeMixModelExplorer: Evaluate: model is ";
//GetModel().Dump();
//cout << "Current likelihood: " << mrOpt.GetLogLikeli() << endl;

    // update the likelihood
    SetLogLikeli( mrOpt.GetLogLikeli() );

    double ll = -1.0*mrOpt.GetLogLikeli();
#endif

    double ll = -1.0*CalcHapLikelihood();
    return ll;
}

double PedigreeMixModelExplorer :: CalcHapLikelihood() const
{
    //const int numPerfGens = 2;
    PedigreeMixLikelihood pedMixCalc;
    pedMixCalc.UseRecursive(true);
    pedMixCalc.SetLogMode(fLogMode);

    vector<double> Ancratio;//waiting for further check

    // note: compute the sum of log-likelihood
    double resSumLogLikeli = 0.0;
    // just compute the probability of each haplotype under the model
    for(int i=0; i<(int)listHapsIn.size(); ++i)
    {
        for(int j=0; j<(int)listHapsIn[i].size(); ++j )
        {
            //double prob2 = pedMixCalc.Compute(GetModel(), i, listHapsIn[i][j] , Ancratio); //waiting for further check
        	double prob2 =0.0;
//cout << "The likelihood (by recursive mode) of the haplotype: " << prob2 << endl;

            if( fLogMode == false )
            {
                if( prob2 >= MIN_POS_VAL)
                {
                    resSumLogLikeli += log(prob2);
                }
            }
            else
            {
                resSumLogLikeli += prob2;
            }
        }
    }
cout << "Overall haplotype likelihood: " << resSumLogLikeli << endl;
    return resSumLogLikeli;
}


//************************************************************************
// Opt over mixing ratio


PedigreeMixModelExplorerMixRatio :: PedigreeMixModelExplorerMixRatio(  PopMixingModel &modelToExp, const vector<vector<PedigreeMixHaplotype> > &listHapsInParam  ) : PedigreeMixModelExplorer( modelToExp, listHapsInParam )
{
    //
}

bool PedigreeMixModelExplorerMixRatio :: Explore()
{
    //
    const double MR_MIN = 0.05;
    const double MR_MAX = 0.95;

    //while(true)
    //{

    double mrOpt;
    double mrCurModel = GetModel().GetMixRatio(1);
    double likeliCur = GetLogLikeli();
    double likeliMax = -1.0*Func1DMinBrent(MR_MIN, mrCurModel, MR_MAX, 0.05, &mrOpt );

    //cout << "Branch: " << branch << ", cur len: " << regMid << ", curBestProb " << loglikeliCurBest <<  ", loglikMax = " << likeliMax << endl;;
    //YW_ASSERT_INFO( likeliMax+tolnum >= loglikeliCurBest, "Wrong: Brent's result does not make sense" );

cout << "In PedigreeMixModelExplorerMixRatio :: Explore: maximum likelihood = " << likeliMax << ", original likelihood = " << GetLogLikeli() << endl;

    if( IsSignificantlyLarge(likeliMax, likeliCur) )
    {
        //
        SetLogLikeli( likeliMax );
cout << "Set Max: likeli: " << GetLogLikeli() << endl;
        //speciesTreeBest.SetBranchLen(branch, szCurr);
        GetModel().SetMixRatio(mrOpt);
        // just continue
        return true;
    }
    else
    {
        // restore to original settings, do not change branch length since no increase is found
cout << "MR: not changed\n";
        SetLogLikeli( likeliCur );
        GetModel().SetMixRatio(mrCurModel);
        // terminate the search
        return false;
    }

}

double PedigreeMixModelExplorerMixRatio :: EvaluateAt(double pt, void *pParam)
{
cout << "PedigreeMixModelExplorerMixRatio: Evaluatiing at: " << pt << endl;
    // set the mix ratio to the value
    GetModel().SetMixRatio(pt);
    return -1.0*CalcHapLikelihood();
}

//************************************************************************
// Opt over mixing generation


PedigreeMixModelExplorerMixGen :: PedigreeMixModelExplorerMixGen(  PopMixingModel &modelToExp, const vector<vector<PedigreeMixHaplotype> > &listHapsInParam  ) : PedigreeMixModelExplorer( modelToExp, listHapsInParam )
{
    //
}

bool PedigreeMixModelExplorerMixGen :: Explore()
{
    //
    //const double MR_MIN = 0.05;
    //const double MR_MAX = 0.95;

    //while(true)
    //{

    double ngOpt;
    int ngenCurModel = GetModel().GetNumMixGens();
    double likeliCur = GetLogLikeli();
    double likeliMax = -1.0*Func1DMinBrent((1.0*MIX_GEN_MIN)/MIX_GEN_MAX, ((double)ngenCurModel)/MIX_GEN_MAX, 1.0, 1.0/MIX_GEN_MAX, &ngOpt);

    //cout << "Branch: " << branch << ", cur len: " << regMid << ", curBestProb " << loglikeliCurBest <<  ", loglikMax = " << likeliMax << endl;;
    //YW_ASSERT_INFO( likeliMax+tolnum >= loglikeliCurBest, "Wrong: Brent's result does not make sense" );

cout << "In PedigreeMixModelExplorerMixRatio :: Explore: maximum likelihood = " << likeliMax << ", original likelihood = " << likeliCur <<": ngOpt: " << ngOpt << endl;

    if( IsSignificantlyLarge(likeliMax, likeliCur ) )
    {
        //
        SetLogLikeli( likeliMax );
        //speciesTreeBest.SetBranchLen(branch, szCurr);
        int ngOptInt = (int)(ngOpt*MIX_GEN_MAX);
        GetModel().SetNumMixGens(ngOptInt);
        // just continue
        return true;
    }
    else
    {
        // restore to original settings, do not change branch length since no increase is found
        SetLogLikeli( likeliCur );
        GetModel().SetNumMixGens(ngenCurModel);
        // terminate the search
        return false;
    }

}

double PedigreeMixModelExplorerMixGen :: EvaluateAt(double pt, void *pParam)
{
    // set the mix ratio to the value
    YW_ASSERT_INFO( pt >= 1.0/MIX_GEN_MAX, "Cannot be too small" );
    int numgen = (int) (pt*MIX_GEN_MAX);
    GetModel().SetNumMixGens(numgen);
cout << "PedigreeMixModelExplorerMixRatio: Evaluatiing at: " << pt <<", ng: " << numgen << endl;
    return -1.0*CalcHapLikelihood();
}



