//
//  Utils4.cpp
//
//
//  Created by Yufeng Wu on 9/30/13.
//
//

#include "Utils4.h"
#include <cmath>

// template functios goes to header file
#if 0

// utilities
template<class TYPE1, class TYPE2>
void CreateMapFromTwoVec( const vector<TYPE1> &vecKey, const vector<TYPE2> &vecval, map<TYPE1,TYPE2> &mapCreated  )
{
    //
    YW_ASSERT_INFO( vecKey.size() == vecval.size(), "veckey has different size as vecval" );
    mapCreated.clear();
    for( int i=0; i<(int)vecKey.size(); ++i )
    {
        //
        mapCreated.insert(map<TYPE1,TYPE2> :: value_type(vecKey[i], vecval[i] ) );
    }
}

template<class TYPE1,class TYPE2>
void KeepCommonInMaps( map<TYPE1,TYPE2> &mapSubtracted, const map<TYPE1,TYPE2> &mapToSub )
{
    // only keep those that is also in the second map
    map<TYPE1, TYPE2> mapNew;
    for( typename map<TYPE1,TYPE2> :: iterator it = mapSubtracted.begin(); it != mapSubtracted.end(); ++it )
    {
        //
        if( mapToSub.find(it->first) != mapToSub.end() )
        {
            // appear in second map so keep
            mapNew.insert( map<TYPE1,TYPE2> :: value_type(it->first, it->second) );
        }
    }
    mapSubtracted = mapNew;
}

template<class TYPE1,class TYPE2>
void KeepCommonInMapsSet( map<TYPE1,TYPE2> &mapSubtracted, const set<TYPE1> &setKept )
{
    // only keep those that is also in the second map
    map<TYPE1, TYPE2> mapNew;
    for( typename map<TYPE1,TYPE2> :: iterator it = mapSubtracted.begin(); it != mapSubtracted.end(); ++it )
    {
        //
        if( setKept.find(it->first) != setKept.end() )
        {
            // appear in second map so keep
            mapNew.insert( map<TYPE1,TYPE2> :: value_type(it->first, it->second) );
        }
    }
    mapSubtracted = mapNew;

}

template<class TYPE1, class TYPE2>
void CreateTwoVecFromMap(const map<TYPE1,TYPE2> &mapIn, vector<TYPE1> &vecKey, vector<TYPE2> &vecval )
{
    vecKey.clear();
    vecval.clear();
    for( typename map<TYPE1,TYPE2> :: iterator it = mapIn.begin(); it != mapIn.end(); ++it )
    {
        vecKey.push_back( it->first );
        vecval.push_back( it->second );
    }
}

#endif

void printvec(vector<double> VEC)
{
	cout<<"[";
	for (int i=0;i<VEC.size();i++)
	{
		cout<<VEC[i];
		if (i<VEC.size()-1)
		{
			cout<<",";
		}
	}
	cout<<"]"<<endl;
}

int GetZeroOneDiff(int x, int y)
{
    if( x == y )
    {
        return 0;
    }
    else
    {
        return 1;
    }
}

void GetMatchingPosIntVec(const int val, const vector<int> &listVals, vector<int> &listPos)
{
    listPos.clear();
    for(int i=0; i<(int)listVals.size(); ++i)
    {
        if( val == listVals[i])
        {
            listPos.push_back(i);
        }
    }
}

void FormUnitVector(int numItems, int posUnit, vector<int> &vecUnit)
{
    //
    YW_ASSERT_INFO(posUnit<numItems, "Wrong");
    vecUnit.clear();
    for(int i=0; i<numItems; ++i)
    {
        //
        vecUnit.push_back(0);
    }
    vecUnit[posUnit] = 1;
}

void FormZeroVector(int numItems, vector<int> &vecZero)
{
    //
    vecZero.clear();
    for(int i=0; i<numItems; ++i)
    {
        //
        vecZero.push_back(0);
    }
}

bool AreTwoSetsCompatible( const set<int> &set1, const set<int> &set2)
{
    // are two sets either disjoint or one contins another
    set<int> sint;
    JoinSets(set1, set2, sint);
    if( sint.size() == 0 || sint.size() == set1.size() || sint.size() == set2.size() )
    {
        return true;
    }
    return false;
}

bool IsSetCompatibleWithSets( const set<int> &set1, const set< set<int> > &setSets)
{
    bool res = true;
    for( set<set<int> > :: const_iterator it = setSets.begin(); it != setSets.end(); ++it )
    {
        if( AreTwoSetsCompatible(set1, *it) == false )
        {
            res = false;
            break;
        }
    }
    return res;
}

bool IsSignificantFraction(int totNum, int numTypes, int numOneType, double minFrac )
{
    // test whether the num of one type occupies a siinficant portion of the totNum (composed of numTypes types)
    if( minFrac >=0.0)
    {
        return numOneType >= totNum*minFrac;
    }
    // if not specific fraction is givn, then use the following rule based on the number of types
    // basicallly require appearing at least two times
    return numOneType >= 2;
}

void IncAllNumInSet( set<int> &sint )
{
    //
    set<int> res;
    for(set<int> :: iterator it = sint.begin(); it != sint.end(); ++it)
    {
        res.insert( *it+1);
    }
    sint = res;
}

void IncAllNumInSets( set< set<int> > &setInts)
{
    //
    set<set<int> > res;
    for( set<set<int> > :: iterator it = setInts.begin(); it != setInts.end(); ++it )
    {
        set<int> sint = *it;
        IncAllNumInSet(sint);
        res.insert( sint);
    }
    setInts = res;
}

void GetNonZeroPosofVec( const vector<int> &vec, set<int> &setpos)
{
    //
    setpos.clear();
    for( int i=0; i<(int)vec.size(); ++i )
    {
        if( vec[i] != 0 )
        {
            setpos.insert(i);
        }
    }
}

int GetSegIndex(int val, const vector<int> &listSegSizes)
{
    //
    int res = -1;
    int szSoFar = 0;
    while( val >= szSoFar  && res < (int)listSegSizes.size() )
    {
        ++res;
        szSoFar += listSegSizes[res];
    }
    return res;
}

// Prob related utilties
double CalcPoisonProb(double rate, int numEvts)
{
    //
    double res = exp(-1.0*rate);
    for(int i=1; i<=numEvts; ++i)
    {
        res *= rate/i;
    }
    return res;
}

void GetDiffPosOfTwoVec(const vector<int> &vec1, const vector<int> &vec2, set<int> &setpos)
{
    //
    YW_ASSERT_INFO(vec1.size() == vec2.size(), "Size: mismatch");
    setpos.clear();
    for(int i=0; i<(int)vec1.size(); ++i)
    {
        if( vec1[i] != vec2[i])
        {
            setpos.insert( i );
        }
    }
}

void ComplementBoolVec(vector<bool> &listVals)
{
    // T->F and vice versa
    for(int i=0; i<(int)listVals.size(); ++i)
    {
        if( listVals[i] == true)
        {
            listVals[i] = false;
        }
        else
        {
            listVals[i] = true;
        }
    }
}

void GetAllGridPoints( int gridLB, int gridUB, int dimGrid, set< vector<int> > &setGridPts )
{
    //  get all grid points whose num is within the range [lb,ub]
    YW_ASSERT_INFO( gridLB <= gridUB, "Bounds wrong");
    YW_ASSERT_INFO( dimGrid >=1, "Dimension must be positive");
    // apply recurrence
    setGridPts.clear();
    if( dimGrid == 1 )
    {
        for( int v=gridLB; v<= gridUB; ++v )
        {
            //
            vector<int> vec;
            vec.push_back( v );
            setGridPts.insert(vec);
        }
    }
    else
    {
        //
        set< vector<int> > setGridPtsSmall;
        GetAllGridPoints( gridLB, gridUB, dimGrid-1, setGridPtsSmall );
        for( set< vector<int> > :: iterator it = setGridPtsSmall.begin(); it != setGridPtsSmall.end(); ++it)
        {
            //
            for( int v=gridLB; v<= gridUB; ++v )
            {
                //
                vector<int> vec = *it;
                vec.push_back( v );
                setGridPts.insert(vec);
            }
        }
    }
}

void MapIntListToAnother( const vector<int> &vec1, const vector<int> &vec2, map<int,int> &mapVec1IndexToVec2)
{
    // given two vectors, e.g. vec1 = [2,1,3] and vec2 = [3,2,1]. Create a map from vec1's index to vec2
    // map = [0,1], [1,2], [2,0]
    // we assume there is no dupllicate for now
//cout << "MapIntListToAnother: vec1: ";
//DumpIntVec(vec1);
//cout << "vec2: ";
//DumpIntVec(vec2);
    mapVec1IndexToVec2.clear();
    YW_ASSERT_INFO( vec1.size() == vec2.size(), "size: mismatch");
    map<int,int> mapValToIndex1, mapValToIndex2;
    for( int i = 0; i<(int)vec1.size(); ++i )
    {
        //
        YW_ASSERT_INFO( mapValToIndex1.find( vec1[i] ) == mapValToIndex1.end(), "Duplicate found" );
        mapValToIndex1.insert( map<int,int> :: value_type( vec1[i], i ) );
//cout << "mapValToIndex1: " << vec1[i] << ", " << i << endl;
    }
    for( int i = 0; i<(int)vec2.size(); ++i )
    {
        //
        YW_ASSERT_INFO( mapValToIndex2.find( vec2[i] ) == mapValToIndex2.end(), "Duplicate found" );
        mapValToIndex2.insert( map<int,int> :: value_type( vec2[i], i ) );
//cout << "mapValToIndex12 " << vec2[i] << ", " << i << endl;
    }
    for( map<int,int> :: iterator it = mapValToIndex1.begin(); it != mapValToIndex1.end(); ++it )
    {
        YW_ASSERT_INFO( mapValToIndex2.find(it->first) != mapVec1IndexToVec2.end(), "Two lists: not idential" );
        mapVec1IndexToVec2.insert( map<int,int> :: value_type( it->second, mapVec1IndexToVec2[it->first] ) );
    }
}

void FindEvenDistriPoints(double valMin, double valMax, double valResolution, int maxNumPoints, vector<double> &listChosenVals)
{
    // pick uniformly some number (<= maxNumPoints) of points within [valMin, valMax}, with distance no more than resolution
    // first figure out spacing
    double valSpacing = (valMax-valMin)/maxNumPoints;
    if( valSpacing < valResolution )
    {
        valSpacing = valResolution;
    }
    for(int i=0; i<(int) (valMax-valMin)/valSpacing ; ++i)
    {
        //
        double val = (i+0.5)*valSpacing;
        listChosenVals.push_back(val);
    }
}

// bits operation
bool IsBitSetInt(int val, int posBit)
{
    //
    // for an index of AC, which src populaiton is a leave
    int mask = (0x1 << posBit);
    // assume only two populaitons for now
    bool res = false;
    if( ( val & mask ) != 0 )
    {
        res = true;
    }
    return res;
}

int ToggleBitInt(int val, int posBit)
{
    //
    return val ^ (1 << posBit);
}


