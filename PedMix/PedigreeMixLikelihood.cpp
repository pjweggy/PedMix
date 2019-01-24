//
//  PedigreeMixLikelihood.cpp
//
//
//  Created by Yufeng Wu on 12/5/14.
//
//

#include "PedigreeMixLikelihood.h"
#include "FFTZ2n.h"
#include "Utils4.h"
#include <cmath>
using namespace std;


//************************************************************************
//  Settings

int numPerfectPedigreeGensSettings = 2;

void SetNumofPerfectPedGens(int numPerfGen)
{
    numPerfectPedigreeGensSettings = numPerfGen;
}
int GetNumofPerfectPedGens()
{
    return numPerfectPedigreeGensSettings;
}

//************************************************************************
//  Pedigree likelihood haplotype


PedigreeMixHaplotype :: PedigreeMixHaplotype( const vector<int> &vecAlleles ) : vecHapAllele( vecAlleles )
{
    //
}

int PedigreeMixHaplotype :: GetAlleleAt(int site) const
{
    //
    return vecHapAllele[site];
}

void PedigreeMixHaplotype :: Dump() const
{
    //
    DumpSequence( vecHapAllele );
}

//************************************************************************
//  Pedigree likelihood model

PopMixingModel :: PopMixingModel(int nPops, int numMixGensInit, double ratioMixInit) : numSrcPops(nPops), numMixGens(numMixGensInit), ratioMix(ratioMixInit)
{
    //
}

int PopMixingModel :: GetNumSites(int locus) const
{
    //
    YW_ASSERT_INFO(locus<(int)listAlleleFreqsPops.size() && listAlleleFreqsPops[locus].size() > 0, "Not initialized");
    return listAlleleFreqsPops[locus][0].size();
}

double PopMixingModel :: GetPopulationAlleleFreq(int locus, int pop, int site, int allele) const
{
    //
//cout << "For pop: " << pop << ", site = " << site << ", allele = " << allele << endl;
    YW_ASSERT_INFO(locus <(int)listAlleleFreqsPops.size(), "Overflow");
    YW_ASSERT_INFO(pop<=1, "Can have only two populations");
    YW_ASSERT_INFO( site<(int) listAlleleFreqsPops[locus][pop].size(), "Overflow" );
    if( allele == 0 )
    {
        // the stored is the allele freq of 1s
    	//cout<<1.0-listAlleleFreqsPops[locus][pop][site]<<endl;
        return 1.0-listAlleleFreqsPops[locus][pop][site];
    }
    else
    {
        return listAlleleFreqsPops[locus][pop][site];
    }
}

double PopMixingModel :: GetMixRatio(int pop) const
{
    // for now assume two pops
    YW_ASSERT_INFO(pop<=1, "Only two populations are supported right now");
    if( pop == 0 )
    {
        // assume ratioMix is the ratio of the population 1 (not 0!)
        return 1.0-ratioMix;
    }
    else
    {
        return ratioMix;
    }
}

double PopMixingModel :: CalcProbRecNoSwitchOneGen(int locus, int site) const
{
    //
    double recFrac = GetRecombFracAt(locus, site);
    return 0.5+0.5*exp(-2.0*recFrac);
}

double PopMixingModel :: CalcProbRecNoSwitch(int locus, int site, int gen) const
{
    //
    return pow( CalcProbRecNoSwitchOneGen(locus, site), gen );
}

double PopMixingModel :: CalcAncesSwitchProb(int locus, int site, int numGensBackwards, int popCur, int popNew) const
{
    cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1\n";
    // prob of switching from one src pop to another at a site within
    // some generations
    double probNoSwitch = CalcProbRecNoSwitch(locus, site, numGensBackwards);
//cout << "Locus: " << locus << ", site: " << site << ", numGensBackwards: " << numGensBackwards << ", popCur: " << popCur << ", popNew: " << popNew << ", probNoSwitch: " << probNoSwitch << endl;
    if( popCur == popNew )
    {
        //
        return probNoSwitch + (1.0-probNoSwitch)* GetMixRatio( popCur );
    }
    else
    {
        // need to first switch and then to the desired one
//cout << "mix ratio for new pop: " << GetMixRatio( popNew ) << endl;
        return (1.0-probNoSwitch) * GetMixRatio( popNew );
    }
}

void PopMixingModel :: SetRecFractions(int locus, const vector<double> &listRecs)
{
    // is the locus already in? if so, just update
    if( locus <(int) listRecFractions.size() )
    {
        listRecFractions[locus] = listRecs;
    }
    else
    {
        // othrwise must be the latest to add
        YW_ASSERT_INFO( locus == (int)listRecFractions.size(), "Cannot set recombination fractions");
        listRecFractions.push_back(listRecs);
    }
}

void PopMixingModel :: SetAlleleFreqForPop( int locus, int pop, const vector<double> &vecPopAlleleFreq )
{
    //
    YW_ASSERT_INFO(locus <=(int)listAlleleFreqsPops.size(), "Wrong1"  );
    if( locus == (int)listAlleleFreqsPops.size() )
    {
        vector<vector<double> > listpp;
        listAlleleFreqsPops.push_back(listpp);
    }
    if( listAlleleFreqsPops[locus].size() == 0 )
    {
        listAlleleFreqsPops[locus].resize( GetNumSrcPops() );
        // also set up the number of sites
        //this->numSites = (int)vecPopAlleleFreq.size();
    }
    YW_ASSERT_INFO(pop<(int)listAlleleFreqsPops[locus].size(), "Only two populations are allowed");
    listAlleleFreqsPops[locus][pop] = vecPopAlleleFreq;
}

void PopMixingModel :: AddAlleleFreqs(int locus, double af1, double af2)
{
//cout << "locus = " << locus << ", af1 = " << af1 << ", af2 = " << af2 << endl;
    YW_ASSERT_INFO(locus <=(int)listAlleleFreqsPops.size(), "Wrong1"  );
    if( locus == (int)listAlleleFreqsPops.size() )
    {
        vector<vector<double> > listpp;
        listAlleleFreqsPops.push_back(listpp);
    }
    //
    if( listAlleleFreqsPops[locus].size() == 0 )
    {
        listAlleleFreqsPops[locus].resize( GetNumSrcPops() );
        // also set up the number of sites
        //this->numSites = 0;
    }
    //vector<double> plist;
    //plist.push_back(af1);
    //plist.push_back(af2);
    listAlleleFreqsPops[locus][0].push_back(af1);
    listAlleleFreqsPops[locus][1].push_back(af2);
    //++numSites;
}

void PopMixingModel :: AddRecFrac(int locus, double rf)
{
    YW_ASSERT_INFO(locus <=(int)listRecFractions.size(), "Wrong1"  );
    if( locus == (int)listRecFractions.size() )
    {
        vector<double>  listpp;
        listRecFractions.push_back(listpp);
    }

    listRecFractions[locus].push_back(rf);
}


void PopMixingModel :: Dump() const
{
    //
    DumpModelOnly();
    cout << "Number of loci: " << GetNumLoci() << endl;
    for(int ll=0; ll<GetNumLoci(); ++ll)
    {
        cout << "At locus " << ll << "   number of sites: " << GetNumSites(ll) << "  ";
        cout << "Allele frequency of sites: ";
        for(int i=0; i<(int)listAlleleFreqsPops[ll][0].size(); ++i)
        {
            cout << "[" << listAlleleFreqsPops[ll][0][i] << "," << listAlleleFreqsPops[ll][1][i] << "] ";
        }
        cout << endl;
        cout <<"Recombinaiton fraction: ";
        DumpDoubleVec(listRecFractions[ll]);
    }
}


void PopMixingModel :: DumpModelOnly() const
{
    //
    cout << "Number of source populaitons: " << numSrcPops << ", number of mixing generations: " << numMixGens << ", ratio of mixing: " << ratioMix << endl;
}

//************************************************************************
//  Pedigree likelihood computing interface

PedigreeMixLikelihood :: PedigreeMixLikelihood() : numPerfectPedigreeGens(GetNumofPerfectPedGens() ), fFFTMode(false), fRecursive(false), fLogMode(false)
{
    //
}
/*
void PedigreeMixLikelihood :: AncestryInfer( const PopMixingModel &model, int locus, const PedigreeMixHaplotype &haplotype, vector<int> &listMPStates, vector<vector<int> > &listFounderPops ) const
{
    // find the list of most likely states using posterier decoding; this will give for each position, what is the most likely
    // ancestry
    vector<double> vecProbsCur;
    vector<double> vecProbsNext;

    // this is a Lander-Green style approach
    vector<vector<double> > listVecProbsCurForward;        // at each site, there are a number of ACs, this stores probs for each AC at this site

    for( int i=0; i<GetNumACs(); ++i )
    {
        //
        double probInit = CalcProbObserveCurSite(i, 0, model, locus, haplotype.GetAlleleAt(0) );
        double probPrior = GetACPrior();
        //cout << "ProbInit: " << probInit << ", probPrior = " << probPrior << endl;
        double probInitCombo = probInit*probPrior;
        if( fLogMode == true )
        {
            probInitCombo = log( probInitCombo );
        }
        vecProbsCur.push_back(probInitCombo);
    }
    listVecProbsCurForward.push_back( vecProbsCur );
    // now move on to next site
    for(int i=1; i<model.GetNumSites(locus); ++i)
    {
        // now update the current probs
        if( fFFTMode == false )
        {
            if(fRecursive == true)
            {
                //
                CalcProbsFromPrevRecur( vecProbsCur, i, model, locus, haplotype.GetAlleleAt(i), vecProbsNext );
            }
            else
            {
                CalcProbsFromPrev( vecProbsCur, i, model, locus, haplotype.GetAlleleAt(i), vecProbsNext );
            }
        }
        else
        {
            YW_ASSERT_INFO(fLogMode == false, "Fatal error: log-mode of FFT is not yet done");
            // test FFT
            //vector<double> vecProbsNext2;
            CalcProbsFromPrevFFT( vecProbsCur, i, model, locus, haplotype.GetAlleleAt(i), vecProbsNext );
        }
        //cout << "***** FFT-based entries: ";
        //DumpDoubleVec(vecProbsNext2);
        // update the current
        vecProbsCur = vecProbsNext;
        //cout << "At site " << i << ": current prob list is ";
        //DumpDoubleVec(vecProbsCur);

        listVecProbsCurForward.push_back(vecProbsCur);
    }

    // now do the same for backward: NOTE BACKEARD is stored backwards
    vector<vector<double> > listVecProbsCurBackward;        // at each site, there are a number of ACs, this stores probs for each AC at this site
    // init at the first position
    vecProbsCur.clear();
    vecProbsNext.clear();
    for(int i=0; i<GetNumACs(); ++i)
    {
        // init to prob w/o transition prob
        //double probInit = CalcProbObserveCurSite(i, 0, model, locus, haplotype.GetAlleleAt(0) );
        double probPrior = 1.0;         // always stop here
        //cout << "ProbInit: " << probInit << ", probPrior = " << probPrior << endl;
        double probInitCombo = probPrior;
        if( fLogMode == true )
        {
            probInitCombo = log( probInitCombo );
        }
        vecProbsCur.push_back(probInitCombo);
    }
    listVecProbsCurBackward.push_back(vecProbsCur);
    //cout << "Prob at the first column: ";
    //DumpDoubleVec(vecProbsCur);
    // now move on to next site
    for(int i=model.GetNumSites(locus)-2; i>=0; --i)
    {
        // now update the current probs
        if( fFFTMode == false )
        {
            if( fRecursive == true )
            {
                CalcProbsFromPrevRecurBackwards( vecProbsCur, i, model, locus, haplotype.GetAlleleAt(i+1), vecProbsNext );
            }
            else
            {
                CalcProbsFromPrevBackwards( vecProbsCur, i, model, locus, haplotype.GetAlleleAt(i+1), vecProbsNext );
            }
        }
        else
        {
            YW_ASSERT_INFO(false, "Not implemented yet.\n");
            //YW_ASSERT_INFO(fLogMode == false, "Fatal error: log-mode of FFT is not yet done");
            // test FFT
            //vector<double> vecProbsNext2;
            //CalcProbsFromPrevFFT( vecProbsCur, i, model, locus, haplotype.GetAlleleAt(i), vecProbsNext );
        }
        //cout << "***** FFT-based entries: ";
        //DumpDoubleVec(vecProbsNext2);
        // update the current
        vecProbsCur = vecProbsNext;
//cout << "At site " << i << ": current prob list is ";
//DumpDoubleVec(vecProbsCur);
        listVecProbsCurBackward.push_back( vecProbsCur );
    }
    ReverseVec( listVecProbsCurBackward );

#if 0
// dump out probs
cout << "listVecProbsCurForward: \n";
for(int i=0; i<(int)listVecProbsCurForward.size(); ++i)
{
DumpDoubleVec( listVecProbsCurForward[i] );
}
cout << "ListVecProbsCurBackward: \n";
for(int i=0; i<(int)listVecProbsCurBackward.size(); ++i)
{
DumpDoubleVec( listVecProbsCurBackward[i] );
}
#endif
    // now check for each site to see which AC has the highest posterior decoding value
    listMPStates.clear();
    for(int i=0; i<model.GetNumSites(locus); ++i)
    {
        vector<double> listProbsACPosterior;
        for(int j=0;j<(int)listVecProbsCurForward[i].size(); ++j )
        {
            //
            double probCombo;
            if( fLogMode == true)
            {
                probCombo = listVecProbsCurForward[i][j] + listVecProbsCurBackward[i][j];
            }
            else
            {
                probCombo = listVecProbsCurForward[i][j]*listVecProbsCurBackward[i][j];
            }
            listProbsACPosterior.push_back( probCombo );
        }
        int ACMaxIndex = GetLargestIndiceInDoubleVec( listProbsACPosterior );
        listMPStates.push_back(ACMaxIndex);
    }
    //cout << "RESULT of POSTERIOR DECODING: ";
    //DumpIntVec( listMPStates );

    // also output the founder populaiton assignment
    listFounderPops.clear();
    for(int i=0; i<(int)listMPStates.size(); ++i)
    {
        vector<int> vv;
        GetFounderPopsForAC( listMPStates[i], vv );
        listFounderPops.push_back(vv);
    }
}
*/
/*
double PedigreeMixLikelihood :: Compute(const PopMixingModel &model, int locus, const PedigreeMixHaplotype &haplotype) const
{
    vector<double> vecProbsCur;//A list of prob for all ACs
    vector<double> vecProbsNext;

    // init at the first position
    for(int i=0; i<GetNumACs(); ++i)
    {
        // init to prob w/o transition prob
        double probInit = CalcProbObserveCurSite(i, 0, model, locus, haplotype.GetAlleleAt(0) );
        double probPrior = GetACPrior();
//cout << "ProbInit: " << probInit << ", probPrior = " << probPrior << endl;
        double probInitCombo = probInit*probPrior;
        if( fLogMode == true )
        {
            probInitCombo = log( probInitCombo );
        }
        vecProbsCur.push_back(probInitCombo);
    }
//cout << "Prob at the first column: ";
//DumpDoubleVec(vecProbsCur);
    // now move on to next site
    for(int i=1; i<model.GetNumSites(locus); ++i)
    {
        // now update the current probs
        if( fFFTMode == false )
        {
            if(fRecursive == true)
            {
                //
                CalcProbsFromPrevRecur( vecProbsCur, i, model, locus, haplotype.GetAlleleAt(i), vecProbsNext );
            }
            else
            {
                CalcProbsFromPrev( vecProbsCur, i, model, locus, haplotype.GetAlleleAt(i), vecProbsNext );
            }
        }
        else
        {
            YW_ASSERT_INFO(fLogMode == false, "Fatal error: log-mode of FFT is not yet done");
            // test FFT
            //vector<double> vecProbsNext2;
            CalcProbsFromPrevFFT( vecProbsCur, i, model, locus, haplotype.GetAlleleAt(i), vecProbsNext );
        }
//cout << "***** FFT-based entries: ";
//DumpDoubleVec(vecProbsNext2);
        // update the current
        vecProbsCur = vecProbsNext;
//cout << "At site " << i << ": current prob list is ";
//DumpDoubleVec(vecProbsCur);
    }

    // this is the final prob: add them up
    if( fLogMode == true )
    {
        return GetLogSumOfLogs( vecProbsCur );
    }
    else
    {
        return GetOverallProb( vecProbsCur );
    }
}
*/

double PedigreeMixLikelihood :: Compute(const PopMixingModel &model, int locus, const PedigreeMixHaplotype &haplotype1, const PedigreeMixHaplotype &haplotype2, vector<double> PP,double phopara,int ll,double probswitch) const
{
    vector<double> vecProbsCur;//A list of prob for all ACs
    vector<double> vecProbsNext;

    //cout<<"the number of ACs is "<<GetNumACs()<<endl;
    //cout<<"the number of internal nodes is "<<GetNumInternalNodes()<<endl;
    //cout<<"the number of leaves is "<<GetNumPerfectPedLeaves()<<endl;

    // init at the first position

    for(int i=0; i<GetNumACs(); ++i)
    {
        // init to prob w/o transition prob
        double probInit = CalcProbObserveCurSite(i,0, model, locus, haplotype1.GetAlleleAt(0), 0)*CalcProbObserveCurSite(i, 0, model, locus, haplotype2.GetAlleleAt(0), 1);//source pop allele frequency
        //cout<<"indexAC----"<<i<<"----emission prob----"<<CalcProbObserveCurSite(i,0, model, locus, haplotype1.GetAlleleAt(0), 0)<<"x"<<CalcProbObserveCurSite(i, 0, model, locus, haplotype2.GetAlleleAt(0), 1)<<endl;
        double probInitCombo = 0.5*probInit;
        if( fLogMode == true )
        {
            probInitCombo = log( probInitCombo );
        }
        vecProbsCur.push_back(probInitCombo);
    }
//********* Init All recom indexAC with prob 0
   if (locus == 0)
   {
    int numAC=GetNumPerfectPedLeaves()+GetNumInternalNodes();
    int numRec=GetNumInternalNodes();
    //cout<<numAC<<"      "<<numRec<<endl; 
    if ( numRec >0 )
    {
         for(int j=0;j<vecProbsCur.size()/2; ++j)
         {
              for(int nb=0;nb<numRec;nb++)
              {
                   int mask = (0x1 << nb );
                   if ( (j & mask ) != 0  && vecProbsCur[j] != log(10e-50))
                   {
                         vecProbsCur[j]=log(10e-50);
                         //cout<<"index "<<j<<" now becomes "<<vecProbsCur[j]<<endl;
                   }
                   mask = (0x1 << (nb+numAC) );
                   if ( (j & mask ) != 0  && vecProbsCur[j] != log(10e-50))
                   {
                         vecProbsCur[j]=log(10e-50);
                         //cout<<"index "<<j<<" now becomes "<<vecProbsCur[j]<<endl;
                   }
              }
         }
    }
    }
//*********
//********* Init All switch indexAC with prob 0
    for(int j=0;j<vecProbsCur.size()/2; ++j)
    {
        //vecProbsCur[j+vecProbsCur.size()/2]=vecProbsCur[j];
        vecProbsCur[j+vecProbsCur.size()/2]=log(10e-50);
    }

    //printvec(vecProbsCur);


    for(int i=1; i<model.GetNumSites(locus); ++i)
    {
        // now update the current probs
        if( fFFTMode == false )
        {
            if(fRecursive == true)
            {
                //
                CalcProbsFromPrevRecur( vecProbsCur, i, model, locus, haplotype1.GetAlleleAt(i), haplotype2.GetAlleleAt(i), vecProbsNext,PP, phopara,ll,probswitch );
            }
            else
            {
                CalcProbsFromPrev( vecProbsCur, i, model, locus, haplotype1.GetAlleleAt(i), haplotype2.GetAlleleAt(i), vecProbsNext,PP, phopara, ll  );
            }
        }
        else
        {
            YW_ASSERT_INFO(fLogMode == false, "Fatal error: log-mode of FFT is not yet done");
            // test FFT
            //vector<double> vecProbsNext2;
            CalcProbsFromPrevFFT( vecProbsCur, i, model, locus, haplotype1.GetAlleleAt(i), vecProbsNext );
        }
//cout << "***** FFT-based entries: ";
//DumpDoubleVec(vecProbsNext2);
        // update the current
        vecProbsCur = vecProbsNext;
        //printvec(vecProbsCur);
        /*
        cout<<endl<<"At site---"<<i<<"---prob list---";
        for (int k=0;k<vecProbsCur.size();k++)
        {
        	cout<<vecProbsCur[k]<<",";
        }
        */
//cout << "At site " << i << ": current prob list is ";
//DumpDoubleVec(vecProbsCur);
    }
    // this is the final prob: add them up
    if( fLogMode == true )
    {
        return GetLogSumOfLogs( vecProbsCur );
    }
    else
    {
        return GetOverallProb( vecProbsCur );
    }
}

void PedigreeMixLikelihood :: CalcProbsFromPrev( const vector<double> &vecProbsPrev, int siteCur, const PopMixingModel &model, int locus, int alleleCur1, int alleleCur2, vector<double> &vecProbsCur,vector<double> PP, double phopara, int ll ) const
{
    // compute the prob of current position based on the prob of previous site
    vecProbsCur.clear();
    for(int i=0; i<(int)vecProbsPrev.size(); ++i  )
    {
        //
    	//cout<<"now compute prob at AC "<<i<<" :"<<endl;
        double probCurStep = CalcProbsFromPrevForIndexSingle( vecProbsPrev, siteCur, model, locus, i, PP, phopara, ll );
        double probObsCurSite = CalcProbObserveCurSite( i, siteCur, model, locus, alleleCur1,0)*CalcProbObserveCurSite( i, siteCur, model, locus, alleleCur2,1);
        //cout<<"emission prob"<<CalcProbObserveCurSite( i, siteCur, model, locus, alleleCur1,0)<<"x"<<CalcProbObserveCurSite( i, siteCur, model, locus, alleleCur2,1)<<endl;
        if( fLogMode == false )
        {
            vecProbsCur.push_back( probCurStep*probObsCurSite );
        }
        else
        {
            vecProbsCur.push_back( probCurStep + log(probObsCurSite ) );
        }
    }
}

void PedigreeMixLikelihood :: CalcProbsFromPrevBackwards( const vector<double> &vecProbsLater, int siteCur, const PopMixingModel &model, int locus, int alleleNext, vector<double> &vecProbsCur ) const
{
    // same as above, but this time compute the backward prob as in standard HMM
    // compute the prob of current position based on the prob of previous site
    vecProbsCur.clear();
    for(int i=0; i<(int)vecProbsLater.size(); ++i  )
    {
        //
        double probCurStep = CalcProbsFromPrevForIndexSingleBackwards( vecProbsLater, siteCur, model, locus, alleleNext, i );
//cout << "In CalcProbsFromPrev: probCurStep = " << probCurStep << endl;
        vecProbsCur.push_back( probCurStep );
    }

}

double PedigreeMixLikelihood :: CalcProbsFromPrevForIndexSingle( const vector<double> &vecProbsPrev, int siteCur, const PopMixingModel &model, int locus, int indexAC, vector<double> PP, double phopara, int ll) const
{
    // directly compute transition prob and explicitly evaluate
    // prob of a single AC (indexed by the specific position)
    double res = 0.0;
    YW_ASSERT_INFO( siteCur > 0, "Site index: cannot be zero" );
    //double recFrac = model.GetRecombFracAt( siteCur-1 );
    // here we have
//cout << "In CalcProbsFromPrevForIndexSingle: indexAC = " << indexAC << ", ancSwitchProb0to1 = " << ancSwitchProb0to1 << ", ancSwitchProb1to0 = " << ancSwitchProb1to0 << ", vecProbsPrev = ";
//DumpDoubleVec(vecProbsPrev);
    for(int i=0; i<(int)vecProbsPrev.size(); ++i)
    {
        //
    	//cout<<"AC--"<<i<<"--transite--"<<indexAC<<"--";
        double probTrans = CalcTransProb( i, indexAC, siteCur, model, locus, PP, phopara, ll);
        if( fLogMode == false )
        {
            res += probTrans*vecProbsPrev[i];
        }
        else
        {
            if( i == 0 )
            {
                res = log(probTrans) + vecProbsPrev[i];
            }
            else
            {
            	if (probTrans == 0)
            	{
            		res=res;
            	}
            	else
            	{
            		res = GetLogSumOfTwo(res, log(probTrans)+vecProbsPrev[i]);
            	}

            }
        }
//cout << "i=" << i << ", probtrans = " << probTrans << ", res = " << res << endl;
    }

    return res;
}

double PedigreeMixLikelihood :: CalcProbsFromPrevForIndexSingleBackwards( const vector<double> &vecProbsPrev, int siteCur, const PopMixingModel &model, int locus, int alleleNext, int indexAC) const
{
    // directly compute transition prob and explicitly evaluate
    // prob of a single AC (indexed by the specific position)
    double res = 0.0;
    vector<double> PP;//change later
    YW_ASSERT_INFO( siteCur >= 0, "Site index: cannot be zero" );
    //double recFrac = model.GetRecombFracAt( siteCur-1 );
    // here we have
    double ancSwitchProb0to1 = model.CalcAncesSwitchProb( locus, siteCur, model.GetNumMixGens()-numPerfectPedigreeGens, 0, 1 );
    double ancSwitchProb1to0 = model.CalcAncesSwitchProb( locus, siteCur, model.GetNumMixGens()-numPerfectPedigreeGens, 1, 0 );
//cout << "In CalcProbsFromPrevForIndexSingle: indexAC = " << indexAC << ", ancSwitchProb0to1 = " << ancSwitchProb0to1 << ", ancSwitchProb1to0 = " << ancSwitchProb1to0 << ", vecProbsPrev = ";
//DumpDoubleVec(vecProbsPrev);
    for(int i=0; i<(int)vecProbsPrev.size(); ++i)
    {
        //
        //double probTrans = CalcTransProb( indexAC, i, siteCur+1, model, locus, ancSwitchProb0to1, ancSwitchProb1to0  );
    	double probTrans = 0.0;
        //double probObsCurSite = CalcProbObserveCurSite( i, siteCur+1, model, locus, alleleNext);
        double probObsCurSite =0.0;
//cout << "At step " << i << ": probObsCurSite = " << probObsCurSite << endl;
        if( fLogMode == false )
        {
            res += probTrans*probObsCurSite*vecProbsPrev[i];
        }
        else
        {
            if( i == model.GetNumSites(locus)-1 )
            {
                res = log(probTrans) + log(probObsCurSite) + vecProbsPrev[i];
            }
            else
            {
                res = GetLogSumOfTwo(res, log(probObsCurSite)+ log(probTrans)+vecProbsPrev[i]);
            }
        }
//cout << "i=" << i << ", probtrans = " << probTrans << ", res = " << res << endl;
    }

    return res;
}


double PedigreeMixLikelihood :: CalcProbObserveCurSite( int indexAC, int siteCur, const PopMixingModel &model, int locus, int alleleCur, int LorR) const
{
    // prob of observing the current allele at this site
    int leafReachedInAC = GetPedLeafReachedIn( indexAC, LorR );
    int srcPopInAC = GetPopSrcFromACIndex( indexAC, leafReachedInAC, LorR );
    //cout<<"Hap leaf reached at"<<leafReachedInAC<<endl;
    //cout<<"Hap src Population InAC is "<<srcPopInAC<<endl;
    double prob = model.GetPopulationAlleleFreq( locus, srcPopInAC, siteCur, alleleCur );
    //cout<<"get ancester population allele freq is-----"<<prob<<"-------prob*"<<endl;

    //cout<<endl<<"prob is "<<prob<<endl;
//cout << "CalcProbObserveCurSite: indexAC = " << indexAC << ", alleleCur = " << alleleCur <<", leafIndex = " << leafReachedInAC << ", srcPopInAC: " << srcPopInAC << ", prob = " << prob << endl;
    return prob;
}

double PedigreeMixLikelihood :: CalcTransProb( int indexPrev, int indexCur, int siteCur, const PopMixingModel &model, int locus,vector<double> PP, double phopara, int ll) const
{
    // calc prob of transiting from one AC to the other
    double res = 1.0;
    double probNoRec = model.CalcProbRecNoSwitchOneGen(locus, siteCur-1);
    int numAC=GetNumInternalNodes()+GetNumPerfectPedLeaves();
    // first the internal nodes
    for(int i=0; i<GetNumInternalNodes(); ++i)
    {
    	//for parent 1
        int rcPrev = GetRecChoiceofACIndex(indexPrev, i);
        int rcCurr = GetRecChoiceofACIndex(indexCur, i);
        double probStep1 = 0.0;
        if( rcPrev == rcCurr)
        {
            probStep1 = probNoRec;
        }
        else
        {
            probStep1 = 1.0-probNoRec;
        }
        res *= probStep1;
        //cout<<probStep1<<"x";
        //for parent 2
        rcPrev = GetRecChoiceofACIndex((indexPrev >> numAC), i);
        rcCurr = GetRecChoiceofACIndex((indexCur >> numAC), i);
        if( rcPrev == rcCurr)
        {
            probStep1 = probNoRec;
        }
        else
        {
            probStep1 = 1.0-probNoRec;
        }
        res *= probStep1;
        //cout<<probStep1<<"x";
    }
    // then AC swap

    double recfrag = model.GetRecombFracAt(locus,siteCur-1)/phopara;
    //    cout<<"recfrag"<<recfrag<<endl;
    double uni = 1.0/ll;
    if (recfrag == 0)
    {
    	recfrag = uni;
    }
    double reclength = recfrag/uni;
    //for parent 2
    for(int i=0; i<GetNumPerfectPedLeaves(); ++i)
    {
        int acPrev = GetPopSrcACIndex(indexPrev, i+GetNumInternalNodes());
        int acCurr = GetPopSrcACIndex(indexCur, i+GetNumInternalNodes());
        double probStep2 = 0.0;
        int pos=PP.size()/2+2*i;
        if(( acPrev == 0 ) && ( acCurr == 0 ))
        {
        	probStep2 = 1-reclength*PP[pos];
        }
        else if (( acPrev == 0 ) && ( acCurr == 1 ))
        {
            probStep2 = reclength*PP[pos];
        }
        else if (( acPrev == 1 ) && ( acCurr == 1 ))
        {
        	probStep2 = 1-reclength*PP[pos+1];
        }
        else if (( acPrev == 1 ) && ( acCurr == 0 ))
        {
        	probStep2 = reclength*PP[pos+1];
        }
        res *= probStep2;
        //cout<<probStep2<<"x";
    }
    //for parent 1
    for(int i=0; i<GetNumPerfectPedLeaves(); ++i)
    {
        int acPrev = GetPopSrcACIndex((indexPrev >> numAC), i+GetNumInternalNodes());
        int acCurr = GetPopSrcACIndex((indexCur >> numAC), i+GetNumInternalNodes());
        double probStep2 = 0.0;
        int pos=2*i;
        if(( acPrev == 0 ) && ( acCurr == 0 ))
        {
        	probStep2 = 1-reclength*PP[pos];
        }
        else if (( acPrev == 0 ) && ( acCurr == 1 ))
        {
            probStep2 = reclength*PP[pos];
        }
        else if (( acPrev == 1 ) && ( acCurr == 1 ))
        {
        	probStep2 = 1-reclength*PP[pos+1];
        }
        else if (( acPrev == 1 ) && ( acCurr == 0 ))
        {
        	probStep2 = reclength*PP[pos+1];
        }
        res *= probStep2;
        //cout<<probStep2<<"x";
    }
    int switchsign1 = (indexPrev >> 2*numAC);
    int switchsign2 = (indexCur >> 2*numAC);
    double switchprob = 0.00002;
    double probStep3;
    if (switchsign1 == switchsign2)
    {
    	probStep3 = 1 - reclength*switchprob;
    }
    else
    {
    	probStep3 = reclength*switchprob;
    }
    res *= probStep3;
    //cout<<probStep3<<endl;
    return res;
}

int PedigreeMixLikelihood :: GetPedLeafReachedIn(int indexAC,int LorR) const
{
    //  given a AC index, find the leaf that is turned on
	// LorR ==0: choose leaf parent no switch no phasing err; LorR==1 choose right parent switch phasing err
    //int mask = (0x1 <<  GetNumInternalNodes() ) - 1;
    //return indexAC & mask;
	int switchsign=(indexAC >> (2*(GetNumPerfectPedLeaves()+GetNumInternalNodes())));
	//cout<<"switchsign move to right "<<(2*(GetNumPerfectPedLeaves()+GetNumInternalNodes()))<<endl;
	//cout<<"indexAC--"<<indexAC<<"--switchsign--"<<switchsign<<endl;

	//
	if ((LorR == 0 && switchsign == 0) || (LorR ==1 && switchsign == 1))
	{
		indexAC = (indexAC >> (GetNumPerfectPedLeaves()+GetNumInternalNodes()));
		//cout<<"indexAC move to right "<<(GetNumPerfectPedLeaves()+GetNumInternalNodes())<<endl;
	}
    int res = 0;
    int nodecur = 0;
    while(nodecur < GetNumInternalNodes() )
    {
        //
        int rcChoice=GetRecChoiceofACIndex(indexAC, nodecur);
        res = (res << 1 )+ rcChoice;
        // move to next
        nodecur = 2*nodecur+1+rcChoice;
    }
    //cout<<"find leaf "<<res<<endl;
    return res;
}

int PedigreeMixLikelihood :: GetRecChoiceofACIndex(int indexAC, int indexNode) const
{
    // get what choice of recombinaiton is made at a node
    // for an index of AC, which src populaiton is a leave
    int mask = (0x1 << indexNode );
    //cout<<(0x1)<<"-----"<<mask<<endl;
    // assume only two populaitons for now
    int res = 0;
    //
    //cout<<( indexAC & mask )<<endl;
    if( ( indexAC & mask ) != 0 )
    {
        res = 1;
    }
    //return 1 if recombination return 0 if not.
    //cout<<"recombination event "<<res<<endl;
    return res;

}

int PedigreeMixLikelihood :: GetPopSrcACIndex(int indexAC, int indexNode) const
{
	int mask = (0x1 << indexNode );
	int res = 0;
	if( ( indexAC & mask ) != 0 )
	{
	    res = 1;
	}
	return res;
}

int PedigreeMixLikelihood :: GetPopSrcFromACIndex(int indexAC, int pedLeafIndex,int LorR) const
{
    // for an index of AC, which src populaiton is a leave
	int switchsign=(indexAC >> (2*(GetNumPerfectPedLeaves()+GetNumInternalNodes())));
	//cout<<"switchsign move right to "<<(2*(GetNumPerfectPedLeaves()+GetNumInternalNodes()))<<endl;
	//cout<<"indexAC--"<<indexAC<<"--switchsign--"<<switchsign<<endl;
	if ((LorR == 0 && switchsign == 0) || (LorR == 1 && switchsign == 1))
	{
		indexAC = (indexAC >>  (GetNumPerfectPedLeaves()+GetNumInternalNodes()));
	}
    int mask = (0x1 << ( GetNumInternalNodes() + GetNumPerfectPedLeaves()- 1 - pedLeafIndex ) );
    //cout<<"find pop src pos  "<<( GetNumInternalNodes() + pedLeafIndex )<<endl;
    // assume only two populaitons for now
    int res = 0;
    if( ( indexAC & mask ) != 0 )
    {
        res = 1;
    }
    //cout<<"Get population source "<<res<<endl;
    return res;
}

int PedigreeMixLikelihood :: GetNumACs() const
{
    //
    return 0x1 << (2*( GetNumInternalNodes() + GetNumPerfectPedLeaves() )+1);
}

int PedigreeMixLikelihood :: GetNumPerfectPedLeaves() const
{
    //
    return 0x1 << (numPerfectPedigreeGens-1);
}

double PedigreeMixLikelihood :: GetOverallProb(const vector<double> &listProbs) const
{
    //
    return GetSumOfElements(listProbs);
}

double PedigreeMixLikelihood :: GetACPrior() const
{
    return 1.0/GetNumACs();
}

void PedigreeMixLikelihood :: CalcProbsFromPrevFFT(const vector<double> &vecProbsPrev, int siteCur, const PopMixingModel &model, int locus, int alleleCur, vector<double> &vecProbsCur) const
{
    // use FFT to compute the probability of the next site
    FFTZ2nUtils fftUtils( GetLengthIndexAC() );

    // first do fft for the current prob
    vector<double> vecProbsPrevFFT;
    fftUtils.ForwardFFT( vecProbsPrev, vecProbsPrevFFT );

    // now the transition prob part; this part can be done by a closed form formula
    vector<double> vecTransFFT;
    FormTransProbFFT( vecTransFFT, siteCur, model, locus );
//cout << "vecTransFFT: ";
//DumpDoubleVec(vecTransFFT);

//vector<double> vecTransInit;
//ConsTransProbVec( vecTransInit, siteCur, model );
//vector<double> vecTransFFT2;
//fftUtils.ForwardFFT( vecTransInit, vecTransFFT2 );
//cout << "vecTransFFT2: ";
//DumpDoubleVec(vecTransFFT2);
    vector<double> PP;//change later
    // also, include the current allele in the transition prob
    for(int i=0; i<(int)vecProbsCur.size(); ++i)
    {
        //
        //double probCurSite = CalcProbObserveCurSite( 0, siteCur, model, locus, alleleCur);
    	double probCurSite = 0.0;
        vecProbsCur[i] *= probCurSite;
    }


    // pointwise multiply
    fftUtils.PointwiseMulti( vecProbsPrevFFT, vecTransFFT );

    // do backwards FFT
    fftUtils.BackwardFFT( vecProbsPrevFFT,  vecProbsCur);

//#if 0
    // at last, include the current allele for each
    for(int i=0; i<(int)vecProbsCur.size(); ++i)
    {
        //
        //double probCurSite = CalcProbObserveCurSite( i, siteCur, model, locus, alleleCur);
        double probCurSite = 0.0;
        vecProbsCur[i] *= probCurSite;
    }
//#endif
}

void PedigreeMixLikelihood :: FormTransProbFFT( vector<double> &listTransFFT, int siteCur, const PopMixingModel &model, int locus ) const
{
    // YW: this routine assume symmetric in transitiing from either of the two
    // populations. MAY NOT HOLD (i.e. when mixing ratio is not 0.5). TBD
    // setup the transition prob in FFT domain
    double probNoRec = model.CalcProbRecNoSwitchOneGen(locus, siteCur-1);
    double ancSwitchProb0to1 = model.CalcAncesSwitchProb( locus, siteCur-1, model.GetNumMixGens()-numPerfectPedigreeGens, 0, 1 );
    //double ancSwitchProb1to0 = model.CalcAncesSwitchProb( siteCur-1, model.GetNumMixGens()-numPerfectPedigreeGens, 1, 0 );
    listTransFFT.clear();
//cout << "probNoRec: " << probNoRec << ", ancSwitchProb0to1 = " << ancSwitchProb0to1 << endl;
    for(int i=0; i<GetNumACs(); ++i)
    {
        double probtransFFT = 1.0;
        for(int j=0; j<GetNumInternalNodes(); ++j)
        {
            if( GetRecChoiceofACIndex(i, j) == 1 )
            {
                probtransFFT *= (2*probNoRec-1.0);
            }
        }
        for(int j=0; j<GetNumPerfectPedLeaves(); ++j)
        {
            if( GetPopSrcFromACIndex(i,j,0) == 1 )//changed
            {
                probtransFFT *= 1-2*ancSwitchProb0to1;
            }
        }
        listTransFFT.push_back(probtransFFT);
    }
}

void PedigreeMixLikelihood :: ConsTransProbVec(vector<double> &listTransInit, int siteCur, const PopMixingModel &model, int locus ) const
{
    //
    double probNoRec = model.CalcProbRecNoSwitchOneGen(locus, siteCur-1);
    double ancSwitchProb0to1 = model.CalcAncesSwitchProb( locus, siteCur-1, model.GetNumMixGens()-numPerfectPedigreeGens, 0, 1 );
    //double ancSwitchProb1to0 = model.CalcAncesSwitchProb( siteCur-1, model.GetNumMixGens()-numPerfectPedigreeGens, 1, 0 );
    listTransInit.clear();
//cout << "probNoRec: " << probNoRec << ", ancSwitchProb0to1 = " << ancSwitchProb0to1 << endl;
    for(int i=0; i<GetNumACs(); ++i)
    {
        double probtrans = 1.0;
        for(int j=0; j<GetNumInternalNodes(); ++j)
        {
            if( GetRecChoiceofACIndex(i, j) == 1 )
            {
                probtrans *= (1-probNoRec);
            }
            else
            {
                probtrans *= probNoRec;
            }
        }
        for(int j=0; j<GetNumPerfectPedLeaves(); ++j)
        {
            if( GetPopSrcFromACIndex(i,j,0) == 1 )//changed
            {
                probtrans *= ancSwitchProb0to1;
            }
            else
            {
                probtrans *= 1.0-ancSwitchProb0to1;
            }
        }
        listTransInit.push_back(probtrans);
    }

}

void PedigreeMixLikelihood :: CalcProbsFromPrevRecur(const vector<double> &vecProbsPrev, int siteCur, const PopMixingModel &model, int locus, int alleleCur1,int alleleCur2, vector<double> &vecProbsCur,vector<double> PP,double phopara,int ll,double probswitch) const
{
    // use recurrence to compute the transition probability
    double probNoRec = model.CalcProbRecNoSwitchOneGen(locus, siteCur-1);
    //cout<<"prob with NO RECOMB is---"<<probNoRec<<endl;
    //double ancSwitchProb0to1 = model.CalcAncesSwitchProb( locus, siteCur-1, model.GetNumMixGens()-numPerfectPedigreeGens, 0, 1 );
    //double ancSwitchProb1to0 = model.CalcAncesSwitchProb( locus, siteCur-1, model.GetNumMixGens()-numPerfectPedigreeGens, 1, 0 );
//cout << "In CalcProbsFromPrevRecur:  " << " ancSwitchProb0to1 = " << ancSwitchProb0to1 << ", ancSwitchProb1to0 = " << ancSwitchProb1to0 <<  endl;

    double recfrag = model.GetRecombFracAt(locus,siteCur-1)/phopara;
//    cout<<"recfrag"<<recfrag<<endl;
    double uni = 1.0/ll;//change later
    if (recfrag == 0)
    {
     	recfrag = uni;
    }
    double reclength = recfrag*ll;
    CalcProbsFromPrevRecurRoutine( true, GetLengthIndexAC()-1, vecProbsPrev, vecProbsCur, probNoRec, PP, reclength, probswitch);


    // finally add the current site
    for(int i=0; i<(int)vecProbsCur.size(); ++i)
    {
        //
        double probCurSite = CalcProbObserveCurSite( i, siteCur, model, locus, alleleCur1,0)*CalcProbObserveCurSite( i, siteCur, model, locus, alleleCur2,1);
        if ( probCurSite == 0 )
        {
        	cout<<"warning!!!!!!!!!!!!!"<<endl;
        }

        if( fLogMode == false )
        {
            vecProbsCur[i] *= probCurSite;
        }
        else
        {
            vecProbsCur[i] += log(probCurSite );
        }
    }

}

/*
void PedigreeMixLikelihood :: CalcProbsFromPrevRecurBackwards(const vector<double> &vecProbsPrev, int siteCur, const PopMixingModel &model, int locus, int alleleNext, vector<double> &vecProbsCur) const
{
//cout << "^^CalcProbsFromPrevRecurBackwards: alleleNext: " << alleleNext << ", sitecur: " << siteCur << endl;
    // use recurrence to compute the transition probability
    double probNoRec = model.CalcProbRecNoSwitchOneGen(locus, siteCur);
    double ancSwitchProb0to1 = model.CalcAncesSwitchProb( locus, siteCur, model.GetNumMixGens()-numPerfectPedigreeGens, 0, 1 );
    double ancSwitchProb1to0 = model.CalcAncesSwitchProb( locus, siteCur, model.GetNumMixGens()-numPerfectPedigreeGens, 1, 0 );
//cout << "In CalcProbsFromPrevRecur:  " << " ancSwitchProb0to1 = " << ancSwitchProb0to1 << ", ancSwitchProb1to0 = " << ancSwitchProb1to0 <<  endl;
//cout << "CalcProbsFromPrevRecurBackwards: vecProbsPrev = ";
//for(int i=0; i<(int)vecProbsPrev.size(); ++i)
//{
//cout << vecProbsPrev[i] << "  ";
//}
//cout << endl;

    // we multiple the prob of observing the current allele to the backward prob pointwise
    // YW: this is due to the nature of backward algorithm
    vector<double> vecProbsCurObsMulti = vecProbsPrev;
    for(int i=0; i<(int)vecProbsCurObsMulti.size(); ++i)
    {
        //
        double probCurSite = CalcProbObserveCurSite( i, siteCur+1, model, locus, alleleNext );
        if( fLogMode == false )
        {
            vecProbsCurObsMulti[i] *= probCurSite;
        }
        else
        {
            vecProbsCurObsMulti[i] += log(probCurSite );
        }
    }
//cout << "vecProbsCurObsMulti: ";
//for(int i=0; i<(int)vecProbsCurObsMulti.size(); ++i)
//{
//cout << vecProbsCurObsMulti[i] << "  ";
//}
//cout << endl;

    // that is it
    CalcProbsFromPrevRecurRoutine( false, GetLengthIndexAC()-1, vecProbsCurObsMulti, ancSwitchProb0to1, ancSwitchProb1to0, probNoRec, vecProbsCur );

}
*/
void PedigreeMixLikelihood :: CalcProbsFromPrevRecurRoutine(bool fForward, int bitpos, const vector<double> &vecProbsPrev, vector<double> &vecProbsCur,double probNoRecSingleGen, vector<double> PP, double recfrag,double probswitch) const
{
    // fForward: determine whether it is for forward or backward (important for how to use transition prob)
//cout << "CalcProbsFromPrevRecurRoutine: bitpos=" << bitpos << ", probAncesSwitch0to1: " << probAncesSwitch0to1 << ", probAnc1to0: " << probAncesSwitch1to0 << ", probNorec: " << probNoRecSingleGen << ", vecProbPrev: ";
//DumpDoubleVec(vecProbsPrev);

    //YW_ASSERT_INFO( (int)vecProbsPrev.size() == (0x1 << (bitpos+1)), "CalcProbsFromPrevRecurRoutine: size is wrong" );
    vecProbsCur.clear();
    vector<vector<double> > vecParams;
    int num_gen=GetNumofPerfectPedGens();
    //cout<<"get generations to proceed is ----"<<num_gen<<endl;
    SetupRecParams( fForward, bitpos, PP, probNoRecSingleGen, vecParams,recfrag, probswitch);
    //cout<<"Now we reach to bitpos ------"<<bitpos<<"!"<<endl;
//cout << "vecParams: " << vecParams[0][0] << " " << vecParams[0][1] << " " << vecParams[1][0] << " " << vecParams[1][1] << endl;
    // setup the recurrence for computing at bit position
    if( bitpos ==0 )
    {
//cout << "bitpos is 0\n";
        // direct computation
        YW_ASSERT_INFO(vecProbsPrev.size() == 2, "vecProbsPrev: size wrong");
        if( fLogMode == false )
        {
            vecProbsCur.push_back( vecParams[0][0]*vecProbsPrev[0] + vecParams[0][1]*vecProbsPrev[1]  );
            vecProbsCur.push_back( vecParams[1][0]*vecProbsPrev[0] + vecParams[1][1]*vecProbsPrev[1]  );
        }
        else
        {
            vecProbsCur.push_back( GetLogSumOfTwo( log(vecParams[0][0])+vecProbsPrev[0], log(vecParams[0][1])+vecProbsPrev[1] )  );
            vecProbsCur.push_back( GetLogSumOfTwo( log(vecParams[1][0])+vecProbsPrev[0], log(vecParams[1][1])+vecProbsPrev[1] ) );
        }
    }
    else
    {
        // apply recurrence
        vector<double> vecProbsPrevPart1, vecProbsPrevPart2;
        SplitItemsofVecIntoTwoParts( vecProbsPrev, vecProbsPrevPart1, vecProbsPrevPart2, vecProbsPrev.size()/2 );
        // apply recurrence
        vector<double> vecSubRes1, vecSubRes2;
        CalcProbsFromPrevRecurRoutine(fForward, bitpos-1, vecProbsPrevPart1, vecSubRes1,probNoRecSingleGen, PP,recfrag,probswitch);
        CalcProbsFromPrevRecurRoutine(fForward, bitpos-1, vecProbsPrevPart2, vecSubRes2,probNoRecSingleGen, PP,recfrag,probswitch);
        // now combine
        vector<double> vecSubRes1Comb1 = vecSubRes1;
        if( fLogMode == false )
        {
            ScaleVectorValBy( vecSubRes1Comb1, vecParams[0][0] );
        }
        else
        {
            //
            OffsetVectorValBy( vecSubRes1Comb1, log(vecParams[0][0]) );
        }
        vector<double> vecSubRes1Comb2 = vecSubRes2;
        if( fLogMode == false )
        {
            ScaleVectorValBy( vecSubRes1Comb2, vecParams[0][1] );
        }
        else
        {
            OffsetVectorValBy( vecSubRes1Comb2, log(vecParams[0][1]) );
        }
        if( fLogMode == false )
        {
            PointwiseAddVectorBy(vecSubRes1Comb1, vecSubRes1Comb2);
        }
        else
        {
            SumofLogVecs(vecSubRes1Comb1, vecSubRes1Comb2);
        }
        vector<double> vecSubRes2Comb1 = vecSubRes1;
        if(fLogMode == false )
        {
            ScaleVectorValBy( vecSubRes2Comb1, vecParams[1][0] );
        }
        else
        {
            OffsetVectorValBy( vecSubRes2Comb1, log(vecParams[1][0]) );
        }
        vector<double> vecSubRes2Comb2 = vecSubRes2;
        if( fLogMode == false )
        {
            ScaleVectorValBy( vecSubRes2Comb2, vecParams[1][1] );
        }
        else
        {
            OffsetVectorValBy( vecSubRes2Comb2, log(vecParams[1][1])  );
        }
        if( fLogMode == false )
        {
            PointwiseAddVectorBy(vecSubRes2Comb1, vecSubRes2Comb2);
        }
        else
        {
            SumofLogVecs(vecSubRes2Comb1, vecSubRes2Comb2);
        }
        // now merge into the final result
        MergeTwoVectorsInto( vecProbsCur, vecSubRes1Comb1, vecSubRes2Comb1 );
    }
/*
cout << "**vecProbsCur: ";
for (int k=0;k<vecProbsCur.size();k++)
{
	cout<<vecProbsCur[k]<<",";
}
cout<<endl;
*/
//DumpDoubleVec(vecProbsCur);
}


/*This is the original CalcProbsFromPrevRecurRoutine
void PedigreeMixLikelihood :: CalcProbsFromPrevRecurRoutine(bool fForward, int bitpos, const vector<double> &vecProbsPrev, double probAncesSwitch0to1, double probAncesSwitch1to0, double probNoRecSingleGen, vector<double> &vecProbsCur) const
{
    // fForward: determine whether it is for forward or backward (important for how to use transition prob)
//cout << "CalcProbsFromPrevRecurRoutine: bitpos=" << bitpos << ", probAncesSwitch0to1: " << probAncesSwitch0to1 << ", probAnc1to0: " << probAncesSwitch1to0 << ", probNorec: " << probNoRecSingleGen << ", vecProbPrev: ";
//DumpDoubleVec(vecProbsPrev);
    YW_ASSERT_INFO( (int)vecProbsPrev.size() == (0x1 << (bitpos+1)), "CalcProbsFromPrevRecurRoutine: size is wrong" );
    vecProbsCur.clear();
    vector<vector<double> > vecParams;
    SetupRecParams( fForward, bitpos, probAncesSwitch0to1, probAncesSwitch1to0, probNoRecSingleGen, vecParams );
//cout << "vecParams: " << vecParams[0][0] << " " << vecParams[0][1] << " " << vecParams[1][0] << " " << vecParams[1][1] << endl;
    // setup the recurrence for computing at bit position
    if( bitpos == 0 )
    {
//cout << "bitpos is 0\n";
        // direct computation
        YW_ASSERT_INFO(vecProbsPrev.size() == 2, "vecProbsPrev: size wrong");
        if( fLogMode == false )
        {
            vecProbsCur.push_back( vecParams[0][0]*vecProbsPrev[0] + vecParams[0][1]*vecProbsPrev[1]  );
            vecProbsCur.push_back( vecParams[1][0]*vecProbsPrev[0] + vecParams[1][1]*vecProbsPrev[1]  );
        }
        else
        {
            vecProbsCur.push_back( GetLogSumOfTwo( log(vecParams[0][0])+vecProbsPrev[0], log(vecParams[0][1])+vecProbsPrev[1] )  );
            vecProbsCur.push_back( GetLogSumOfTwo( log(vecParams[1][0])+vecProbsPrev[0], log(vecParams[1][1])+vecProbsPrev[1] ) );
        }
    }
    else
    {
        // apply recurrence
        vector<double> vecProbsPrevPart1, vecProbsPrevPart2;
        SplitItemsofVecIntoTwoParts( vecProbsPrev, vecProbsPrevPart1, vecProbsPrevPart2, vecProbsPrev.size()/2 );
        // apply recurrence
        vector<double> vecSubRes1, vecSubRes2;
        CalcProbsFromPrevRecurRoutine(fForward, bitpos-1, vecProbsPrevPart1, probAncesSwitch0to1, probAncesSwitch1to0, probNoRecSingleGen, vecSubRes1);
        CalcProbsFromPrevRecurRoutine(fForward, bitpos-1, vecProbsPrevPart2, probAncesSwitch0to1, probAncesSwitch1to0, probNoRecSingleGen, vecSubRes2);
        // now combine
        vector<double> vecSubRes1Comb1 = vecSubRes1;
        if( fLogMode == false )
        {
            ScaleVectorValBy( vecSubRes1Comb1, vecParams[0][0] );
        }
        else
        {
            //
            OffsetVectorValBy( vecSubRes1Comb1, log(vecParams[0][0]) );
        }
        vector<double> vecSubRes1Comb2 = vecSubRes2;
        if( fLogMode == false )
        {
            ScaleVectorValBy( vecSubRes1Comb2, vecParams[0][1] );
        }
        else
        {
            OffsetVectorValBy( vecSubRes1Comb2, log(vecParams[0][1]) );
        }
        if( fLogMode == false )
        {
            PointwiseAddVectorBy(vecSubRes1Comb1, vecSubRes1Comb2);
        }
        else
        {
            SumofLogVecs(vecSubRes1Comb1, vecSubRes1Comb2);
        }
        vector<double> vecSubRes2Comb1 = vecSubRes1;
        if(fLogMode == false )
        {
            ScaleVectorValBy( vecSubRes2Comb1, vecParams[1][0] );
        }
        else
        {
            OffsetVectorValBy( vecSubRes2Comb1, log(vecParams[1][0]) );
        }
        vector<double> vecSubRes2Comb2 = vecSubRes2;
        if( fLogMode == false )
        {
            ScaleVectorValBy( vecSubRes2Comb2, vecParams[1][1] );
        }
        else
        {
            OffsetVectorValBy( vecSubRes2Comb2, log(vecParams[1][1])  );
        }
        if( fLogMode == false )
        {
            PointwiseAddVectorBy(vecSubRes2Comb1, vecSubRes2Comb2);
        }
        else
        {
            SumofLogVecs(vecSubRes2Comb1, vecSubRes2Comb2);
        }
        // now merge into the final result
        MergeTwoVectorsInto( vecProbsCur, vecSubRes1Comb1, vecSubRes2Comb1 );
    }
//cout << "**vecProbsCur: ";
//DumpDoubleVec(vecProbsCur);
}
*/

void PedigreeMixLikelihood :: SetupRecParams( bool fForward, int bitpos, vector<double> PP, double probNoRecSingleGen, vector<vector<double> > &vecParams,double recfrag,double probswitch) const
{
    // setup parameters of the recurrence (2x2)
    vecParams.clear();
    vecParams.resize(2);
    //cout<<"prob switch is"<<probswitch<<endl;
    //double probswitch=0.0000000000001;
    int numAC=GetNumPerfectPedLeaves()+GetNumInternalNodes();
    int numRec=GetNumInternalNodes();
    int len_PP=PP.size();
    //cout<<"numAC="<<numAC<<"  numRec="<<numRec<<"   len_PP="<<len_PP<<endl;
    if( fForward == true )
    {
        if (bitpos == 2*numAC) //switch bit
        {
                //cout<<"Phasing error bit"<<endl;
        	if (probswitch*recfrag > 1)
        	{
        		vecParams[0].push_back( 0.0 );
        		vecParams[0].push_back(1.0);
        		vecParams[1].push_back(1.0);
        		vecParams[1].push_back(0.0);
        	}
        	else
        	{
        		vecParams[0].push_back( 1.0-probswitch*recfrag );
        		vecParams[0].push_back(probswitch*recfrag);
        		vecParams[1].push_back(probswitch*recfrag);
        		vecParams[1].push_back(1.0-probswitch*recfrag);
        	}

        }
        else if((bitpos < numRec) || ((bitpos < numRec+numAC) && (bitpos >= numAC))) //recombination bits
        {
        	//cout<<"if bitpos is smaller than-----"<<GetNumInternalNodes()<<"--------get recombination parameters"<<endl;
            // this is for recombination only
            //cout<<"Recombination setting bit"<<endl;
            vecParams[0].push_back( probNoRecSingleGen );
            vecParams[0].push_back(1.0-probNoRecSingleGen);
            vecParams[1].push_back(1.0-probNoRecSingleGen);
            vecParams[1].push_back(probNoRecSingleGen);
        }
        else if (bitpos>=numRec && bitpos < numAC) //right parent
        {
        	//int pos=2*(len_PP/2-1);
        	int pos=len_PP/2+2*(numAC-bitpos-1);
                //cout<<"Now bit is "<<bitpos<<"Right Ancestral setting bit parent "<<pos<<endl;
//                double Ma=PP[pos+1]/(PP[pos]+PP[pos+1]);
//                double Mb=PP[pos]/(PP[pos]+PP[pos+1]);
        	if (recfrag*PP[pos]>1)
        	{
        		vecParams[0].push_back(0.0 );
        		vecParams[1].push_back(1.0);
        	}
        	else
        	{
//        		vecParams[0].push_back((1.0-recfrag*PP[pos])+(1.0-probNoRecSingleGen)*(Ma-(1.0-recfrag*PP[pos])) );
//        		vecParams[1].push_back(recfrag*PP[pos]+(1.0-probNoRecSingleGen)*(Mb-recfrag*PP[pos]));
                        vecParams[0].push_back(1.0-recfrag*PP[pos] );
                        vecParams[1].push_back(recfrag*PP[pos]);

        	}
        	if (recfrag*PP[pos+1]>1)
        	{
        		vecParams[0].push_back(1.0);
        		vecParams[1].push_back(0.0);
        	}
            else
            {
//        		vecParams[0].push_back(recfrag*PP[pos+1]+(1.0-probNoRecSingleGen)*(Ma-recfrag*PP[pos+1]));
//        		vecParams[1].push_back((1.0-recfrag*PP[pos+1])+(1.0-probNoRecSingleGen)*(Mb-(1.0-recfrag*PP[pos+1])));
                        vecParams[0].push_back(recfrag*PP[pos+1]);
                        vecParams[1].push_back(1.0-recfrag*PP[pos+1]);
        	}
        }
/*
        else if (bitpos==numRec+1)
        {
        	int pos=2*(len_PP/2-2);
        	if (recfrag*PP[pos]>1)
        	{
        		vecParams[0].push_back(0.0 );
        		vecParams[1].push_back(1.0);
        	}
        	else
        	{
        		vecParams[0].push_back(1.0-recfrag*PP[pos] );
        		vecParams[1].push_back(recfrag*PP[pos]);
        	}
        	if (recfrag*PP[pos+1]>1)
        	{
        		vecParams[0].push_back(1.0);
        		vecParams[1].push_back(0.0);
        	}
            else
            {
        		vecParams[0].push_back(recfrag*PP[pos+1]);
        		vecParams[1].push_back(1.0-recfrag*PP[pos+1]);
        	}
        }
*/
        else if (bitpos >= numAC+numRec && bitpos < 2*numAC)
        {
        	//int pos=2*(len_PP/2-3);
        	int pos=2*(2*numAC-bitpos-1);
                //cout<<"Now bit is "<<bitpos<<" Left Ancestral setting bit parent "<<pos<<endl;
//                double Ma=PP[pos+1]/(PP[pos]+PP[pos+1]);
//                double Mb=PP[pos]/(PP[pos]+PP[pos+1]);
        	if (recfrag*PP[pos]>1)
        	{
        		vecParams[0].push_back(0.0 );
        		vecParams[1].push_back(1.0);
        	}
        	else
        	{
//                        vecParams[0].push_back((1.0-recfrag*PP[pos])+(1.0-probNoRecSingleGen)*(Ma-(1.0-recfrag*PP[pos])) );
//                        vecParams[1].push_back(recfrag*PP[pos]+(1.0-probNoRecSingleGen)*(Mb-recfrag*PP[pos]));
        		vecParams[0].push_back(1.0-recfrag*PP[pos] );
        		vecParams[1].push_back(recfrag*PP[pos]);
        	}
        	if (recfrag*PP[pos+1]>1)
        	{
        		vecParams[0].push_back(1.0);
        		vecParams[1].push_back(0.0);
        	}
            else
            {
//                        vecParams[0].push_back(recfrag*PP[pos+1]+(1.0-probNoRecSingleGen)*(Ma-recfrag*PP[pos+1]));
//                        vecParams[1].push_back((1.0-recfrag*PP[pos+1])+(1.0-probNoRecSingleGen)*(Mb-(1.0-recfrag*PP[pos+1])));
        		vecParams[0].push_back(recfrag*PP[pos+1]);
        		vecParams[1].push_back(1.0-recfrag*PP[pos+1]);
        	}
        }
/*
        else if (bitpos==numAC+numRec+1)
        {
        	int pos=2*(len_PP/2-4);
        	if (recfrag*PP[pos]>1)
        	{
        		vecParams[0].push_back(0.0 );
        		vecParams[1].push_back(1.0);
        	}
        	else
        	{
        		vecParams[0].push_back(1.0-recfrag*PP[pos] );
        		vecParams[1].push_back(recfrag*PP[pos]);
        	}
        	if (recfrag*PP[pos+1]>1)
        	{
        		vecParams[0].push_back(1.0);
        		vecParams[1].push_back(0.0);
        	}
            else
            {
        		vecParams[0].push_back(recfrag*PP[pos+1]);
        		vecParams[1].push_back(1.0-recfrag*PP[pos+1]);
        	}
        }
*/
/*
        else if((bitpos < numAC) &&(bitpos >=numRec))
        {
        	int num_anc=PP.size()/2;
        	for (int k=0;k<num_anc/2;k++)
        	{
        		if ((bitpos>=numRec+2*k) && (bitpos<numRec+2*(k+1)))
        		{
        			int pos=2*(num_anc-k-1);
        			if (recfrag*PP[pos]>1)
        			{
        				vecParams[0].push_back(0.0 );
        				vecParams[1].push_back(1.0);
        			}
        			else
        			{
        				vecParams[0].push_back(1.0-recfrag*PP[pos] );
        				vecParams[1].push_back(recfrag*PP[pos]);
        			}
        			if (recfrag*PP[pos+1]>1)
        			{
        				vecParams[0].push_back(1.0);
        				vecParams[1].push_back(0.0);
        			}
        			else
        			{
        				vecParams[0].push_back(recfrag*PP[pos+1]);
        				vecParams[1].push_back(1.0-recfrag*PP[pos+1]);
        			}
        		}
        	}
        }
        else if((bitpos < 2*numAC) &&(bitpos >=(numRec+numAC)))
        {
        	int num_anc=PP.size()/2;
        	for (int k=0;k<num_anc/2;k++)
        	{
        		if ((bitpos>=numRec+numAC+2*k) && (bitpos<numRec+numAC+2*(k+1)))
        		{
        			int pos=2*(num_anc/2-k-1);
        			if (recfrag*PP[pos]>1)
        			{
        				vecParams[0].push_back(0.0 );
        				vecParams[1].push_back(1.0);
        			}
        			else
        			{
        				vecParams[0].push_back(1.0-recfrag*PP[pos] );
        				vecParams[1].push_back(recfrag*PP[pos]);
        			}
        			if (recfrag*PP[pos+1]>1)
        			{
        				vecParams[0].push_back(1.0);
        				vecParams[1].push_back(0.0);
        			}
        			else
        			{
        				vecParams[0].push_back(recfrag*PP[pos+1]);
        				vecParams[1].push_back(1.0-recfrag*PP[pos+1]);
        			}
        		}
        	}
        }
*/
/*
        //cout<<probNoRecSingleGen<<endl;
        for (int k=0;k<2;++k)
        {
        	cout<<vecParams[k][0]<<"--"<<vecParams[k][1]<<endl;
        }
*/
    }
    else
    {///haven't changed this part
#if 0
        if( bitpos < GetNumInternalNodes() )
        {
            // this is for recombination only
            vecParams[0].push_back( probNoRecSingleGen );
            vecParams[1].push_back(1.0-probNoRecSingleGen);
            vecParams[0].push_back(1.0-probNoRecSingleGen);
            vecParams[1].push_back(probNoRecSingleGen);
        }
        else
        {
            // for ancestry change
            vecParams[0].push_back( 1.0-probAncesSwitch0to1 );
            vecParams[1].push_back(probAncesSwitch1to0);
            vecParams[0].push_back(probAncesSwitch0to1);
            vecParams[1].push_back(1.0-probAncesSwitch1to0);
        }
#endif
    }
}
/* This is the original code
 void PedigreeMixLikelihood :: SetupRecParams( bool fForward, int bitpos, double probAncesSwitch0to1, double probAncesSwitch1to0, double probNoRecSingleGen, vector<vector<double> > &vecParams) const
{
    // setup parameters of the recurrence (2x2)
    vecParams.clear();
    vecParams.resize(2);

    if( fForward == true )
    {
        if( bitpos < GetNumInternalNodes() )
        {
            // this is for recombination only
            vecParams[0].push_back( probNoRecSingleGen );
            vecParams[0].push_back(1.0-probNoRecSingleGen);
            vecParams[1].push_back(1.0-probNoRecSingleGen);
            vecParams[1].push_back(probNoRecSingleGen);
        }
        else
        {
            // for ancestry change
            vecParams[0].push_back( 1.0-probAncesSwitch0to1 );
            vecParams[0].push_back(probAncesSwitch1to0);
            vecParams[1].push_back(probAncesSwitch0to1);
            vecParams[1].push_back(1.0-probAncesSwitch1to0);
        }
    }
    else
    {
        if( bitpos < GetNumInternalNodes() )
        {
            // this is for recombination only
            vecParams[0].push_back( probNoRecSingleGen );
            vecParams[1].push_back(1.0-probNoRecSingleGen);
            vecParams[0].push_back(1.0-probNoRecSingleGen);
            vecParams[1].push_back(probNoRecSingleGen);
        }
        else
        {
            // for ancestry change
            vecParams[0].push_back( 1.0-probAncesSwitch0to1 );
            vecParams[1].push_back(probAncesSwitch1to0);
            vecParams[0].push_back(probAncesSwitch0to1);
            vecParams[1].push_back(1.0-probAncesSwitch1to0);
        }
    }
}

 */
void PedigreeMixLikelihood :: GetFounderPopsForAC(int indexAC, vector<int> &listFounderPops) const
{
    // get the founders (great grand parents) populations coded in the AC index
    listFounderPops.clear();
    for(int i=0; i<GetNumPerfectPedLeaves(); ++i)
    {
        //
        int mask = 0x1 << GetNumInternalNodes();
        mask = mask << i;
        int pop = 0;
        if( (indexAC & mask) == 0 )
        {
            pop = 0;
        }
        else
        {
            pop = 1;
        }
        listFounderPops.push_back(pop);
    }
}

