//
//  PedigreeMixLikelihood.h
//
//
//  Created by Yufeng Wu on 12/5/14.
//  Compute the likelihood of one haplotype under the
//  perfect tree model
//

#ifndef ____PedigreeMixLikelihood__
#define ____PedigreeMixLikelihood__

#include "Utils3.h"

//************************************************************************
//  Settings

void SetNumofPerfectPedGens(int numPerfGen);
int GetNumofPerfectPedGens();

//************************************************************************
//  Pedigree likelihood haplotype

class PedigreeMixHaplotype
{
public:
    PedigreeMixHaplotype() {}
    PedigreeMixHaplotype(const PedigreeMixHaplotype &rhs) : vecHapAllele(rhs.vecHapAllele) {}
    PedigreeMixHaplotype( const vector<int> &vecAlleles );
    void AddAllele(int allele) { vecHapAllele.push_back(allele); }
    int GetAlleleAt(int site) const;
    void Dump() const;

private:
    vector<int> vecHapAllele;
};

//************************************************************************
//  Pedigree likelihood model

class PopMixingModel
{
public:
    PopMixingModel(int nPops, int numMixGensInit, double ratioMixInit);
    int GetNumLoci() const { return listAlleleFreqsPops.size(); }
    int GetNumSites(int loci) const;
    double GetPopulationAlleleFreq(int locus, int pop, int site, int allele) const;
    int GetNumSrcPops() const { return numSrcPops; }
    int GetNumMixGens() const { return numMixGens; }
    void SetNumMixGens(int ng) { numMixGens = ng; }
    double GetMixRatio(int pop) const;
    void SetMixRatio(double mr) { ratioMix = mr; }
    double GetRecombFracAt(int loci, int site) const { return listRecFractions[loci][site]; }
    double CalcAncesSwitchProb(int loci, int site, int numGensBackwards, int popCur, int popNew) const;
    double CalcProbRecNoSwitchOneGen(int loci, int site) const;
    double CalcProbRecNoSwitch(int loci, int site, int gen) const;
    void SetRecFractions(int loci, const vector<double> &listRecs);
    void SetAlleleFreqForPop( int loci, int pop,  const vector<double> &vecPopAlleleFreq );
    void AddAlleleFreqs(int loci, double af1, double af2);
    void AddRecFrac(int loci, double rf);
    void Dump() const;
    void DumpModelOnly() const;

private:

    int numSrcPops;
    //vector<int> numSites;
    vector<vector<vector<double> > > listAlleleFreqsPops;
    vector<vector<double> > listRecFractions;
    int numMixGens;
    double ratioMix;
};

//************************************************************************
//  Pedigree likelihood computing interface

class MixLikelihood
{
public:
    virtual double Compute(const PopMixingModel &model, int locus, const PedigreeMixHaplotype &haplotype1, const PedigreeMixHaplotype &haplotype2, vector<double> PP,double phopara,int ll,double probswitch) const = 0;
    //virtual double ComputeBackwards(const PopMixingModel &model, int locus, const PedigreeMixHaplotype &haplotype) const = 0;
    virtual int GetMaxNumPedigreeGen() const = 0;
};

// perfect tree based likelihood computation
class PedigreeMixLikelihood : public MixLikelihood
{
public:
    PedigreeMixLikelihood();
    //void AncestryInfer( const PopMixingModel &model, int locus, const PedigreeMixHaplotype &haplotype, vector<int> &listMPStates, vector<vector<int> > &listFounderPops ) const;
    virtual double Compute(const PopMixingModel &model, int locus, const PedigreeMixHaplotype &haplotype1, const PedigreeMixHaplotype &haplotype2, vector<double> PP,double phopara, int ll,double probswitch) const;
    //virtual double ComputeBackwards(const PopMixingModel &model, int locus, const PedigreeMixHaplotype &haplotype) const;
    void UseFFT(bool f) { fFFTMode = f; fRecursive = false; }
    void UseRecursive(bool f) { fRecursive = f; fFFTMode = false; }
    void SetLogMode(bool f) { fLogMode = f; }
    virtual int GetMaxNumPedigreeGen() const { return 3; }

private:
    int GetNumACs() const;
    int GetNumPerfectPedLeaves() const;
    int GetNumInternalNodes() const { return GetNumPerfectPedLeaves()-1; }
    int GetLengthIndexAC() const { return 2*(GetNumPerfectPedLeaves()+GetNumInternalNodes())+1; }
    void CalcProbsFromPrev( const vector<double> &vecProbsPrev, int siteCur, const PopMixingModel &model, int locus, int alleleCur1, int alleleCur2, vector<double> &vecProbsCur,vector<double> PP, double phopara, int ll) const;
    void CalcProbsFromPrevBackwards( const vector<double> &vecProbsLater, int siteCur, const PopMixingModel &model, int locus, int alleleNext, vector<double> &vecProbsCur ) const;
    void CalcProbsFromPrevFFT(const vector<double> &vecProbsPrev, int siteCur, const PopMixingModel &model, int locus, int alleleCur, vector<double> &vecProbsCur) const;
    void CalcProbsFromPrevRecur(const vector<double> &vecProbsPrev, int siteCur, const PopMixingModel &model, int locus, int alleleCur1, int alleleCur2, vector<double> &vecProbsCur, vector<double> PP,double phopara, int ll, double probswitch) const;
    //void CalcProbsFromPrevRecurBackwards(const vector<double> &vecProbsPrev, int siteCur, const PopMixingModel &model, int locus, int alleleCur, vector<double> &vecProbsCur) const;
    double CalcProbsFromPrevForIndexSingle( const vector<double> &vecProbsPrev, int siteCur, const PopMixingModel &model, int locus, int indexAC,vector<double> PP, double phopara, int ll) const;
    double CalcProbsFromPrevForIndexSingleBackwards( const vector<double> &vecProbsLater, int siteCur, const PopMixingModel &model, int locus, int alleleNext, int indexAC) const;
    double CalcTransProb( int indexPrev, int indexCur, int siteCur, const PopMixingModel & model, int locus, vector<double> PP, double phopara, int ll) const;
    int GetPopSrcFromACIndex(int indexAC, int pedLeafIndex,int LorR) const;
    int GetRecChoiceofACIndex(int indexAC, int indexNode) const;
    int GetPopSrcACIndex(int indexAC, int indexNode) const;
    double CalcProbObserveCurSite( int indexAC, int siteCur, const PopMixingModel &model, int locus, int alleleCur,int LorR) const;
    int GetPedLeafReachedIn(int indexAC, int LorR) const;
    double GetOverallProb(const vector<double> &listProbs) const;
    double GetACPrior() const;
    void FormTransProbFFT( vector<double> &listTransFFT, int siteCur, const PopMixingModel &model, int locus  ) const;
    void ConsTransProbVec(vector<double> &listTransFFT, int siteCur, const PopMixingModel &model, int locus ) const;
    void CalcProbsFromPrevRecurRoutine(bool fForward, int bitpos, const vector<double> &vecProbsPrev, vector<double> &vecProbCur,double probNoRecSingleGen, vector<double> PP,double recfrag,double probswitch) const;
    void SetupRecParams(bool fForward, int bitpos,vector<double> PP, double probNoRecSingleGen, vector<vector<double> > &vecParams,double recfrag,double probswitch) const;
    void GetFounderPopsForAC(int indexAC, vector<int> &listFounderPops) const;

    // number of perfect pedigree generation
    int numPerfectPedigreeGens;
    bool fFFTMode;
    bool fRecursive;
    bool fLogMode;
};



#endif /* defined(____PedigreeMixLikelihood__) */
