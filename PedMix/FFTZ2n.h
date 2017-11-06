//
//  FFTZ2n.h
//
//
//  Created by Yufeng Wu on 12/6/14.
//  Utilities for performing FFT (forward/backward) in Z2n domain
//  Refer to the paper by KRUGLYAK and LANDER
//  Caution: the length of the given FFT values must be 2^n

#ifndef ____FFTZ2n__
#define ____FFTZ2n__

#include "Utils3.h"

//************************************************************************
//  FFT utilities on Z2^n domain

class FFTZ2nUtils
{
public:
    FFTZ2nUtils(int numBits);
    void ForwardFFT(const vector<double> &listValsOnZ2n, vector<double> &listValsOnZ2nFFT) const;
    void PointwiseMulti( vector<double> &listValsOnZ2n, const vector<double> &listValsOnZ2nMulti ) const;
    void BackwardFFT( const vector<double> &listValsOnZ2nFFT, vector<double> &listValsOnZ2n ) const;

private:
    int GetNumElems() const { return 0x1 << numBitsZ2n; }

    int numBitsZ2n;
};



#endif /* defined(____FFTZ2n__) */
