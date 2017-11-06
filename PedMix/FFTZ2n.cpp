//
//  FFTZ2n.cpp
//
//
//  Created by Yufeng Wu on 12/6/14.
//
//

#include "FFTZ2n.h"
#include "Utils4.h"

//************************************************************************
//  FFT utilities on Z2^n domain


FFTZ2nUtils :: FFTZ2nUtils(int numBits) : numBitsZ2n(numBits)
{
    //
}

void FFTZ2nUtils :: ForwardFFT(const vector<double> &listValsOnZ2n, vector<double> &listValsOnZ2nFFT) const
{
    if(  (int)listValsOnZ2n.size() != GetNumElems() )
    {
        //
        cout << "listValsOnZ2n: " << listValsOnZ2n.size() << ", GetNumElemens() = " << GetNumElems() << endl;
    }
    YW_ASSERT_INFO( (int)listValsOnZ2n.size() == GetNumElems(), "ForwardFFT: Size is wrong" );
    // init
    listValsOnZ2nFFT.clear();
    for(int i=0; i<GetNumElems(); ++i)
    {
        listValsOnZ2nFFT.push_back( listValsOnZ2n[i] );
    }
    // iterate n times
    for(int k=0; k<numBitsZ2n; ++k)
    {
        vector<double> listValsOnZ2nFFTNext;
        // compute each term
        for( int i=0; i<GetNumElems(); ++i )
        {
            //
            int fac1 = 1;
            if( IsBitSetInt( i, k ) == true )
            {
                fac1 = -1;
            }
            int valIndex2 = ToggleBitInt(i,k);
//cout <<"k=" << k << ", i=" << i << ", valIndex2 = " << valIndex2 << endl;
            double vnext = fac1*listValsOnZ2nFFT[i] + listValsOnZ2nFFT[valIndex2];
            listValsOnZ2nFFTNext.push_back(vnext);
        }

        listValsOnZ2nFFT = listValsOnZ2nFFTNext;
    }
}

void FFTZ2nUtils :: PointwiseMulti( vector<double> &listValsOnZ2n, const vector<double> &listValsOnZ2nMulti ) const
{
    // just multiple point by point and save into the first vector
    YW_ASSERT_INFO(listValsOnZ2n.size() == listValsOnZ2nMulti.size(), "PointwiseMulti: size mismatch");
    for(int i=0; i<(int)listValsOnZ2n.size(); ++i)
    {
        listValsOnZ2n[i] *= listValsOnZ2nMulti[i];
    }
}

void FFTZ2nUtils :: BackwardFFT( const vector<double> &listValsOnZ2nFFT, vector<double> &listValsOnZ2n ) const
{
    // perform backwards FFT. This is simple: just do as in the forward
    // the only difference is multiply a factor
    ForwardFFT(listValsOnZ2nFFT, listValsOnZ2n);
    int numItems = GetNumElems();
    for(int i=0; i<(int)listValsOnZ2n.size(); ++i)
    {
        listValsOnZ2n[i] = listValsOnZ2n[i]/numItems;
    }
}

