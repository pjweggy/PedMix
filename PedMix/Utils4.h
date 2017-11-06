//
//  Utils4.h
//
//
//  Created by Yufeng Wu on 9/30/13.
//  More utilities that may be useful
//

#ifndef ____Utils4__
#define ____Utils4__

#include "Utils3.h"
#include <cmath>
#include <algorithm>
#include <sstream>
#include <string>
#include <queue>

#define YW_VERY_SMALL_FRACTION   0.000000000001

#if 0
// utilities
template<class TYPE1, class TYPE2>
void CreateMapFromTwoVec( const vector<TYPE1> &vecKey, const vector<TYPE2> &vecval, map<TYPE1,TYPE2> &mapCreated  );

template<class TYPE1,class TYPE2>
void KeepCommonInMaps( map<TYPE1,TYPE2> &mapSubtracted, const map<TYPE1,TYPE2> &mapToSub );

template<class TYPE1,class TYPE2>
void KeepCommonInMapsSet( map<TYPE1,TYPE2> &mapSubtracted, const set<TYPE1> &setKept );

template<class TYPE1, class TYPE2>
void CreateTwoVecFromMap(const map<TYPE1,TYPE2> &mapIn, vector<TYPE1> &vecKey, vector<TYPE2> &vecval);
#endif

// a list of templates
template<class TYPE1, class TYPE2>
void CreateMapFromTwoVec( const vector<TYPE1> &vecKey, const vector<TYPE2> &vecval, map<TYPE1,TYPE2> &mapCreated  )
{
    //
    YW_ASSERT_INFO( vecKey.size() == vecval.size(), "veckey has different size as vecval" );
    mapCreated.clear();
    for( int i=0; i<(int)vecKey.size(); ++i )
    {
        //
        mapCreated.insert(typename map<TYPE1,TYPE2> :: value_type(vecKey[i], vecval[i] ) );
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
            mapNew.insert( typename map<TYPE1,TYPE2> :: value_type(it->first, it->second) );
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
            mapNew.insert( typename map<TYPE1,TYPE2> :: value_type(it->first, it->second) );
        }
    }
    mapSubtracted = mapNew;

}

template<class TYPE1, class TYPE2>
void CreateTwoVecFromMap(const map<TYPE1,TYPE2> &mapIn, vector<TYPE1> &vecKey, vector<TYPE2> &vecval )
{
    vecKey.clear();
    vecval.clear();
    for( typename map<TYPE1,TYPE2> :: const_iterator it = mapIn.begin(); it != mapIn.end(); ++it )
    {
        vecKey.push_back( it->first );
        vecval.push_back( it->second );
    }
}

template<class TYPE>
TYPE GetSumOfVecElements(const vector<TYPE> &listVals)
{
    TYPE sum = 0;
    for(int i=0; i<(int)listVals.size(); ++i)
    {
        sum += listVals[i];
    }
    return sum;
}

template<class TYPE>
double MyCalcStdError(const vector<TYPE> &listVals, const vector<TYPE> &listValsRef)
{
    YW_ASSERT_INFO(listVals.size() == listValsRef.size(), "CalcStdError: Size mismatch");
    double sum = 0.0;
    if( listValsRef.size() == 0)
    {
        return 0.0;
    }
    for(int i=0; i<(int)listVals.size(); ++i)
    {
        double diff = listVals[i] - listValsRef[i];
        sum += diff*diff;
    }
    return sqrt(sum/listValsRef.size());
}

template<class TYPE>
void GetPositionsOverThres(const vector<TYPE> &listVals, const TYPE& val, int maxNum, set<int> &listPoses)
{
    // get positions that are either over or at, as long as the total number is not over
    listPoses.clear();
    for(int p=0; p<(int)listVals.size(); ++p)
    {
        //
        if(  val < listVals[p] && (int)listPoses.size() < maxNum  )
        {
            //
            listPoses.insert(p);
        }
    }
    // now also check for those equal
    for(int p=0; p<(int)listVals.size(); ++p)
    {
        //
        if(  val == listVals[p] && (int)listPoses.size() < maxNum  )
        {
            //
            listPoses.insert(p);
        }
    }
}

template<class TYPE>
void AddVecTo(vector<TYPE> &vecAdded, const vector<TYPE> &vecAdding)
{
    //
    YW_ASSERT_INFO(vecAdded.size() == vecAdding.size(), "Size mismatch");
    for(int i=0; i<(int)vecAdding.size(); ++i)
    {
        vecAdded[i] += vecAdding[i];
    }
}

template<class TYPE>
int FindMajorityElemVal(const vector<TYPE> &listItems, double valThres)
{
    // find out whether there is an item that is over some percentage of all the sum (say 50%)
    //TYPE sum = GetSumOfVecElements(listItems);
    //YW_ASSERT_INFO(sum > 0, "Can not only have zero");
    for(int i=0; i<(int)listItems.size(); ++i)
    {
        if( (double)(listItems[i]) > valThres )
           {
               return i;
           }
    }
    // no majority item
    return -1;
}

template<class TYPE>
int FindMajorityElem(const vector<TYPE> &listItems, double fracMaj)
{
    TYPE sum = GetSumOfVecElements(listItems);
    return FindMajorityElemVal(listItems, fracMaj*sum);
}

template<class TYPE>
void FindMajorityMultiElemVal(const vector<TYPE> &listItems, double valThres, int maxNum, set<int> &listChosenPos)
{
    // find out whether there is an item that is over some percentage of all the sum (say 50%)
    listChosenPos.clear();
    //TYPE sum = GetSumOfVecElements(listItems);
    //YW_ASSERT_INFO(sum > 0, "Can not only have zero");
    vector<TYPE> listItemsSort = listItems;
    YWSort(listItemsSort);
    TYPE sumSoFar = 0;
    int numAdd = 0;
    int indexPicked = -1;
    for(int i=(int)listItemsSort.size()-1; i>=0;--i)
    {
        ++numAdd;
        if( numAdd > maxNum )
        {
            break;
        }
        sumSoFar += listItemsSort[i];

        if( (double)(sumSoFar) > valThres )
        {
            indexPicked = i;
            break;
        }
    }
    // no majority item
    if( indexPicked < 0 )
    {
        return;
    }
    // find the set of items that is at least that much
    // get those over first
    for( int i=0; i<(int)listItems.size(); ++i)
    {
        if( listItems[i] > listItemsSort[indexPicked])
        {
            listChosenPos.insert(i);
        }
    }
    for( int i=0; i<(int)listItems.size(); ++i)
    {
        if( listItems[i] == listItemsSort[indexPicked])
        {
            if( (int)listChosenPos.size() < numAdd )
            {
                listChosenPos.insert(i);
            }
            else
            {
                break;
            }
        }
    }
}

template<class TYPE>
void FindMajorityMultiElem(const vector<TYPE> &listItems, double fracMaj, int maxNum, set<int> &listChosenPos)
{
    // find out whether there is an item that is over some percentage of all the sum (say 50%)
    listChosenPos.clear();
    TYPE sum = GetSumOfVecElements(listItems);
    YW_ASSERT_INFO(sum > 0, "Can not only have zero");
    FindMajorityMultiElemVal( listItems, fracMaj*sum, maxNum, listChosenPos );
}


template<class TYPE>
bool YWSortCmpFunc (TYPE i,TYPE j) { return (i<j); }

template<class TYPE>
void YWSort(vector<TYPE> &vecIn)
{
    //
    typedef bool (*comparer_t)(const TYPE, const TYPE);
    comparer_t cmp = &YWSortCmpFunc;
    std::sort (vecIn.begin(), vecIn.end(), cmp);
}

template<class TYPE>
void DumpPair(const pair<TYPE,TYPE> &pp)
{
    cout << "[" << pp.first << "," << pp.second << "]";
}

template<class TYPE>
void GetSubsetItem( const vector<TYPE> &listItems, const vector<int> &vecpos, set<TYPE> &subsetItems )
{
    //
    subsetItems.clear();
    for(int i=0; i<(int)vecpos.size(); ++i)
    {
        YW_ASSERT_INFO( vecpos[i] <(int)listItems.size(), "Fail");
        subsetItems.insert( listItems[ vecpos[i] ]);
    }
}

template<class TYPE>
void GetSubsetSets( const vector<set<TYPE> > &listItems, const vector<int> &vecpos, set<TYPE> &subsetItems )
{
    //
    subsetItems.clear();
    for(int i=0; i<(int)vecpos.size(); ++i)
    {
        YW_ASSERT_INFO( vecpos[i] <(int)listItems.size(), "Fail");
        UnionSetsGen( subsetItems, listItems[ vecpos[i] ] );
    }
}

template<class TYPE>
int FindMaxValPositionFromList(const vector<TYPE> &listItems)
{
    // find out whether there is an item that is over some percentage of all the sum (say 50%)
    //TYPE sum = GetSumOfVecElements(listItems);
    //YW_ASSERT_INFO(sum > 0, "Can not only have zero");
    YW_ASSERT_INFO( listItems.size() > 0, "Must have at least one");
    TYPE valMaxCur = listItems[0];
    int posMaxCur = 0;
    for(int i=1; i<(int)listItems.size(); ++i)
    {
        if( valMaxCur < listItems[i] )
        {
            posMaxCur = i;
            valMaxCur = listItems[i];
        }
    }
    // no majority item
    return posMaxCur;
}

template<class TYPE>
bool IsSetContainerGen(const set<TYPE> &container, const set<TYPE> &contained)
{
    //
    for(typename set<TYPE> ::iterator it = contained.begin(); it != contained.end(); ++it)
    {
        if(  container.find( *it ) == container.end()  )
        {
            return false;
        }
    }
    return true;
}

template<class TYPE>
bool AreSetsIntersecting(const set<TYPE> &s1In, const set<TYPE> &s2In)
{
    //
    const set<TYPE> *ptrSet1 = &s1In;
    const set<TYPE> *ptrSet2 = &s2In;
    if( s1In.size() > s2In.size() )
    {
        ptrSet1 = &s2In;
        ptrSet2 = &s1In;
    }
    for(typename set<TYPE> ::iterator it = ptrSet1->begin(); it != ptrSet1->end(); ++it)
    {
        if(  ptrSet2->find( *it ) != ptrSet2->end()  )
        {
            return true;
        }
    }
    return false;
}

template<class TYPE>
string ConvToString(const TYPE &val)
{
    ostringstream convert;   // stream used for the conversion
    convert << val;      // insert the textual representation of 'Number' in the characters in the stream
    return convert.str();
}

template<class TYPE>
string ConsNewickTreeFromClades(const set<set<TYPE> > &setClades)
{
    // clade: a collection of taxa (int or string); output newick format
    // first, the set of taxa is always the outmost clade
    set<TYPE> setTaxa;
    map<set<TYPE>, set<TYPE>  > mapCladePars;
    for(typename set<set<TYPE> > :: const_iterator it = setClades.begin(); it!=setClades.end(); ++it)
    {
        // find it out the set of taxa
        for( typename set<TYPE> :: iterator itg=it->begin(); itg != it->end(); ++itg)
        {
            setTaxa.insert(*itg);
        }
    }
    set<set<TYPE> > setCladesUsed = setClades;
    setCladesUsed.insert(setTaxa);
    // also ensure single taxon is in
    for( typename set<TYPE> :: iterator it = setTaxa.begin(); it!= setTaxa.end(); ++it )
    {
        //
        set<TYPE> ss;
        ss.insert(*it);
        setCladesUsed.insert(ss);
    }
    // order the clades by size (YW: not the best implementation but hope it will work)
    map<int, set<set<TYPE> > > mapCladesSz;
    for( typename set<set<TYPE> > :: iterator it = setCladesUsed.begin(); it != setCladesUsed.end(); ++it)
    {
        if( mapCladesSz.find(it->size()) == mapCladesSz.end() )
        {
            set<set<TYPE> > ss;
            mapCladesSz.insert( typename map<int,set<set<TYPE> > > :: value_type(it->size(), ss) );
        }
        mapCladesSz[it->size()].insert(*it);
    }
    // find par of each clade
    for( typename set<set<TYPE> > :: iterator it = setCladesUsed.begin(); it != setCladesUsed.end(); ++it )
    {
        //
        for( typename set<set<TYPE> > :: iterator it2 = setCladesUsed.begin(); it2 != setCladesUsed.end(); ++it2 )
        {
            //
            if( it2 != it && IsSetContainerGen(*it2, *it) == true  )
            {
                //
                if( mapCladePars.find( *it) == mapCladePars.end() )
                {
                    mapCladePars.insert( typename map<set<TYPE>,set<TYPE> > :: value_type(*it, *it2) );
                }
                else if(mapCladePars[*it].size() > it2->size())
                {
                    mapCladePars[*it] = *it2;
                }
            }
        }
    }
    // now assign each clade a string
    map<set<TYPE>, string> mapCladeToStr;
    queue<set<TYPE> > queueToProc;
    // init leaves
    for(typename set<TYPE> :: iterator it = setTaxa.begin(); it != setTaxa.end(); ++it)
    {
        set<TYPE> ss;
        ss.insert(*it);
        string strLbl = ConvToString(*it);
        mapCladeToStr.insert( typename map<set<TYPE>, string>  :: value_type(ss, strLbl) );
        queueToProc.push(ss);
    }
    // now proc from bottom up
    for( typename map<int, set<set<TYPE> > > :: iterator it = mapCladesSz.begin(); it != mapCladesSz.end(); ++it)
    {
        for( typename set<set<TYPE> > :: iterator itg = it->second.begin(); itg != it->second.end(); ++itg)
        {
            YW_ASSERT_INFO( mapCladeToStr.find( *itg) != mapCladeToStr.end(), "Fail to find string");
            // pass it to parent
            if( mapCladePars.find( *itg) != mapCladePars.end() )
            {
                string strBase = mapCladeToStr[*itg];
                if( itg->size() > 1 )
                {
                    // add parenthsis
                    strBase = "(" + mapCladeToStr[*itg]+")";
                }

                //
                set<TYPE> sPar = mapCladePars[*itg];
                if( mapCladeToStr.find(sPar) == mapCladeToStr.end() )
                {
                    mapCladeToStr.insert( typename map<set<TYPE>, string> :: value_type( sPar, strBase ) );
                }
                else
                {
                    mapCladeToStr[sPar] = mapCladeToStr[sPar] + "," + strBase;
                }
            }
        }
    }
    // finally
    YW_ASSERT_INFO( mapCladeToStr.find( setTaxa) != mapCladeToStr.end(), "Wrong");
    string res = "(" + mapCladeToStr[setTaxa] + ")";
    return res;
}

template<class TYPE>
void FindMaximalSets( set<set<TYPE> > &setsItems )
{
    // only keep those with no super set
    set<set<TYPE> > setsItemsRes;
    for( typename set<set<TYPE> > :: iterator it = setsItems.begin(); it != setsItems.end(); ++it)
    {
        bool fSuperSet = false;
        for( typename set<set<TYPE> > :: iterator itg = setsItems.begin(); itg != setsItems.end(); ++itg)
        {
            // is itg the super set?
            if( itg->size() > it->size() && IsSetContainerGen( *itg, *it) == true )
            {
                fSuperSet = true;
                break;
            }
        }
        if( fSuperSet == false )
        {
            setsItemsRes.insert(*it);
        }
    }
    setsItems = setsItemsRes;
}

template<class TYPE>
void InitVecWithVal( vector<TYPE> &listVec, TYPE valInit, int numItems)
{
    listVec.clear();
    for(int i=0; i<numItems; ++i)
    {
        listVec.push_back(valInit);
    }
}

template<class TYPE>
void PopulateVecBySetGen( vector<TYPE> &vec, const set<TYPE> &sset)
{
    //
    vec.clear();
    for( typename set<TYPE> :: const_iterator it = sset.begin(); it != sset.end(); ++it)
    {
        vec.push_back( *it );
    }
}

template<class TYPE>
void PopulateVecBySetPtrGen( vector<const TYPE *> &vec, const set<TYPE> &sset)
{
    //
    vec.clear();
    for( typename set<TYPE> :: const_iterator it = sset.begin(); it != sset.end(); ++it)
    {
        vec.push_back( &(*it) );
    }
}

template<class TYPE>
void PopulateSetPtrBySetGen( set<const TYPE *> &sptrs, const set<TYPE> &sset)
{
    //
    sptrs.clear();
    for( typename set<TYPE> :: const_iterator it = sset.begin(); it != sset.end(); ++it)
    {
        sptrs.insert( &(*it) );
    }
}


template<class TYPE>
void PopulateSetByVecGen( set<TYPE> &sset, const vector<TYPE> &vec)
{
    //
    vec.clear();
    for( typename vector<TYPE> :: const_iterator it = vec.begin(); it != vec.end(); ++it)
    {
        sset.insert( *it );
    }
}

template<class TYPE>
void PopulateSetBySetPtrGen( set<TYPE> &sset, const set<const TYPE *> &ssetPtr)
{
    //
    sset.clear();
    for( typename set<const TYPE *> :: const_iterator it = ssetPtr.begin(); it != ssetPtr.end(); ++it)
    {
        sset.push_back( *(*it) );
    }
}


template<class TYPE1, class TYPE2>
void MergeMapGen(map<TYPE1,TYPE2> &mapCombined, const map<TYPE1,TYPE2> &mapToAdd)
{
    for( typename map<TYPE1,TYPE2> :: const_iterator it = mapToAdd.begin(); it != mapToAdd.end(); ++it)
    {
        mapCombined.insert( typename map<TYPE1,TYPE2> :: value_type(it->first, it->second) );
    }
}

template<class TYPE>
void SplitItemsBySetOfPartition( const set<TYPE> &setItems, const set<set<TYPE> > &setPartitions, vector<set<TYPE> > &vecSplitParts )
{
    // setItems: a list of items; setpartitions: parition the space of items; vecSplitParts: split setItems into unit of those partitions
    // approach, take join repeatitively
    vecSplitParts.clear();
    set<TYPE> setItemsUse = setItems;
    while(setItemsUse.size() > 0)
    {
        bool fSub = false;
        for( typename set<set<TYPE> > :: iterator it = setPartitions.begin(); it != setPartitions.end(); ++it )
        {
            //
            set<TYPE> setItemSub;
            JoinSetsGen(*it, setItemsUse, setItemSub);
            YW_ASSERT_INFO( setItemSub.size() == 0 || setItemSub.size() == it->size(), "Not a partition");
            if( setItemSub.size() == it->size() && it->size() > 0 )
            {
                vecSplitParts.push_back(*it);
                SubtractSetsGen( setItemsUse, *it  );
                fSub = true;
            }
        }
        YW_ASSERT_INFO( fSub == true || setItemsUse.size() == 0, "FATAL ERROR: not progress made in SplitItemsBySetOfPartition");
    }
}

template<class TYPE>
bool SplitItemsBySetOfPartitionTF( const set<TYPE> &setItems, const set<set<TYPE> > &setPartitions, vector<set<TYPE> > &vecSplitParts )
{
    // setItems: a list of items; setpartitions: parition the space of items; vecSplitParts: split setItems into unit of those partitions
    // approach, take join repeatitively
    vecSplitParts.clear();
    set<TYPE> setItemsUse = setItems;
    while(setItemsUse.size() > 0)
    {
        bool fSub = false;
        for( typename set<set<TYPE> > :: iterator it = setPartitions.begin(); it != setPartitions.end(); ++it )
        {
            //
            set<TYPE> setItemSub;
            JoinSetsGen(*it, setItemsUse, setItemSub);
            if( setItemSub.size() > 0 && setItemSub.size() < it->size()  )
            {
                return false;
            }
            if( setItemSub.size() == it->size() && it->size() > 0 )
            {
                vecSplitParts.push_back(*it);
                SubtractSetsGen( setItemsUse, *it  );
                fSub = true;
            }
        }
        YW_ASSERT_INFO( fSub == true || setItemsUse.size() == 0, "FATAL ERROR: not progress made in SplitItemsBySetOfPartition");
    }
    return true;
}

template<class TYPE>
void SplitItemsofVecIntoTwoParts( const vector<TYPE> &vecItems, vector<TYPE> &vecFirstPart, vector<TYPE> &vecSecondPart, int posStartof2ndPart )
{
    // caution: position is 0 based
    vecFirstPart.clear();
    vecSecondPart.clear();
    for(int i=0; i<(int)vecItems.size() && i<posStartof2ndPart; ++i )
    {
        vecFirstPart.push_back(vecItems[i]);
    }
    for(int i=posStartof2ndPart; i<(int)vecItems.size(); ++i )
    {
        vecSecondPart.push_back(vecItems[i]);
    }
}

template<class TYPE>
void MergeTwoVectorsInto( vector<TYPE> &vecItems, const vector<TYPE> &vecFirstPart, const vector<TYPE> &vecSecondPart)
{
    //
    vecItems.clear();
    for(int i=0; i<(int)vecFirstPart.size(); ++i )
    {
        vecItems.push_back(vecFirstPart[i]);
    }
    for(int i=0; i<(int)vecSecondPart.size(); ++i )
    {
        vecItems.push_back(vecSecondPart[i]);
    }
}

template<class TYPE>
void ScaleVectorValBy(vector<TYPE> &vecItems, const TYPE &factor)
{
    for(int i=0; i<(int)vecItems.size(); ++i)
    {
        vecItems[i] *= factor;
    }
}

template<class TYPE>
void OffsetVectorValBy(vector<TYPE> &vecItems, const TYPE &factor)
{
    for(int i=0; i<(int)vecItems.size(); ++i)
    {
        vecItems[i] += factor;
    }
}

template<class TYPE>
void PointwiseMultiVectorBy(vector<TYPE> &vecItems, const vector<TYPE> &vecItemsFactors)
{
    YW_ASSERT_INFO( vecItems.size() == vecItemsFactors.size(), "PointwiseMultiVectorBy: size wrong" );
    for(int i=0; i<(int)vecItems.size(); ++i)
    {
        vecItems[i] *= vecItemsFactors[i];
    }
}

template<class TYPE>
void PointwiseAddVectorBy(vector<TYPE> &vecItemsAdded, const vector<TYPE> &vecItemsAdding)
{
    YW_ASSERT_INFO( vecItemsAdded.size() == vecItemsAdding.size(), "PointwiseMultiVectorBy: size wrong" );
    for(int i=0; i<(int)vecItemsAdded.size(); ++i)
    {
        vecItemsAdded[i] += vecItemsAdding[i];
    }
}

template<class TYPE>
void CopyVecToArray( const vector<TYPE> &vecItems, TYPE *parray)
{
    // CAUTION: the array must have adequate size to avoid buffer overrun
    for(int i=0; i<(int)vecItems.size(); ++i)
    {
        parray[i] = vecItems[i];
    }
}

template<class TYPE>
void CopyArrayToVec( TYPE *parray, int sz, vector<TYPE> &vecItems)
{
    // CAUTION: the array must have adequate size to avoid buffer overrun
    vecItems.clear();
    for(int i=0; i<sz; ++i)
    {
        vecItems.push_back( parray[i] );
    }
}

template<class TYPE>
void SwapItemsInVec(vector<TYPE> &vecItems, int pos1, int pos2)
{
    YW_ASSERT_INFO(pos1<(int)vecItems.size() && pos2 <(int)vecItems.size(), "Overflow");
    TYPE tmp = vecItems[pos1];
    vecItems[pos1] = vecItems[pos2];
    vecItems[pos2] = tmp;
}

#if 0
// donot work; so DO NOT USE
// copy from a pointer to a one-dimensional array to a 2d vector
template<class TYPE>
void Copy2DArrayToVec(TYPE *parray, int nr, int nc, vector<vector<TYPE> > &vec2DArray)
{
    vec2DArray.clear();
    for(int i=0; i<nr; ++i)
    {
        vector<TYPE> row;
        for(int j=0; j<nc; ++j)
        {
            row[j] = parray[i*nc+j];
        }
cout << "Row " << i << " is done\n";
        vec2DArray.push_back(row);
    }
}
#endif

//#if 0
template<class TYPE>
void ReduceContainerSetsForSetsGen(vector<set<TYPE> > &listSets)
{
    // give a list of sets, if one set A contains another set B, then remove
    // the intersection between them from A (not B)
    // if there is non-empty intersection but neither contains one another, DO NOTHING!
    // note: there may be multiple ways for doing this; fornow, this procedure just finds a legal solution
    vector<set<TYPE> > listSetsNext;     // we ensure there is no container sets here
    // process each input set, if it contains any set in the new list, reduces it and add to the ist
    // if contained, reduce the one already in teh list (which still introduce no new container in the old list)
    for( int i=0; i<(int)listSets.size(); ++i )
    {
        //
        set<int> setToAdd = listSets[i];
        // loop until no more container is found
        bool fCont = true;
        while(fCont == true)
        {
            fCont = false;
            for(int j=0; j<(int)listSetsNext.size(); ++j)
            {
                // test whether the new set contains any of
                set<TYPE> setInt;
                JoinSetsGen(setToAdd, listSetsNext[j], setInt);
                if( setInt.size() == listSetsNext[j].size() )
                {
                    // reduce the one to add
                    SubtractSetsGen( setToAdd, setInt );
                    fCont = true;       // since we updated the one to add (so maybe new containment emerage), need to continue looping
                }
                else
                {
                    if( setInt.size() == setToAdd.size() )
                    {
                        SubtractSetsGen( listSetsNext[j], setInt );
                    }
                }
            }
        }
//cout << "Adding a new set to next set:";
//DumpIntSet(setToAdd);
        // add it
        listSetsNext.push_back(setToAdd);
    }
    // this is the updated sets that contains no containers
    listSets = listSetsNext;
//cout << "Resulting sets: ";
//for(int i=0; i<(int)listSets.size(); ++i)
//{
//DumpIntSet(listSets[i]);
//}
}
//#endif

template<class TYPE>
void RemoveVecElementAt(vector<TYPE> &listItems, int pos)
{
    // remove the item at the pos
    if( pos <(int) listItems.size() )
    {
        listItems.erase(listItems.begin() + pos);
    }
}

template<class TYPE>
void AppendItemToBoundedVec(const TYPE &item, vector<TYPE> &listItem, int posvecToAdd, int maxSize)
{
    // add an item to the position in a vector
    // if max capacity is reached, then drop the last one
    YW_ASSERT_INFO( posvecToAdd <=(int)listItem.size(), "Position: wrong" );
    if( (int)listItem.size() == maxSize && posvecToAdd == (int)listItem.size() )
    {
        // no room for it
        return;
    }
    else
    {
        // create a new list
        vector<TYPE> listItemNew;
        int pos = 0;
        for(; pos<posvecToAdd; ++pos)
        {
            listItemNew.push_back(listItem[pos]);
        }
        // add this item
        listItemNew.push_back(item);
        // add the rest if needed
        for(; pos<(int)listItem.size(); ++pos)
        {
            if( (int)listItemNew.size() >= maxSize )
            {
                // overflow, stop
                break;
            }
            else
            {
                listItemNew.push_back( listItem[pos] );
            }
        }
        listItem = listItemNew;
    }
}

// create a combined list by merging items (and then take average)
template<class TYPE>
void PutItemsInBuckets(int numBuckets, const vector<TYPE> &listItemsIn, vector<TYPE> &itemsInBuckets)
{
    // if list is empty, then dont do it
    if( listItemsIn.size() > 0 )
    {
        // here buckets contains the average items in the original list
        int stepNum = listItemsIn.size()/numBuckets;
        if( stepNum*numBuckets < (int)listItemsIn.size() )
        {
            stepNum += 1;
        }
        int pos = 0;
        for(int i=0; i<numBuckets; ++i)
        {
            //
            bool fStop = false;
            TYPE tot = 0;
            for(int j=0; j<stepNum; ++j)
            {
                if( pos >= (int)listItemsIn.size() )
                {
                    fStop = true;
                    break;
                }
                tot += listItemsIn[pos];
                ++pos;
            }
            if( fStop == false )
            {
                itemsInBuckets.push_back( tot/stepNum );
            }
        }
    }
    // fill in 0 if otherwise
    while( (int)itemsInBuckets.size() < numBuckets )
    {
        itemsInBuckets.push_back(0);
    }
}

template<class TYPE>
void ReverseVec( vector<TYPE> &vec)
{
    //cout << "Before switching: vec = ";
    //DumpIntVec( vec );
    // This function would reverse the integer vector, i.e. vec[0] = vec[n-1] and so on
    for(int i=0; i<(int)vec.size()/2; ++i)
    {
        TYPE tmp = vec[ (int)vec.size()- 1 -i ];
        vec[ (int)vec.size()- 1 -i ] = vec[i];
        vec[i] = tmp;
    }
    //cout << "After switching: vec = ";
    //DumpIntVec( vec );
}

// extract 1D array from 2D array
template<class TYPE>
void ExtractColFrom2DArray( const vector<vector<TYPE> > &array2D, int col, vector<TYPE> &vecCol )
{
    vecCol.clear();
    YW_ASSERT_INFO(  array2D.size() == 0 || col < (int)array2D[0].size(), "Overflow"  );
    for(int i=0; i<(int)array2D.size(); ++i )
    {
        vecCol.push_back( array2D[i][col] );
    }
}

// calc mean and variance
template<class TYPE>
void CalcMeanVarianceFor( const vector<TYPE> &listVals, double &valMean, double &valVar)
{
    YW_ASSERT_INFO( listVals.size() > 0, "Empty input" );

    //
    double valSum = 0.0;
    for(int i=0; i<(int)listVals.size(); ++i)
    {
        valSum += (double)listVals[i];
    }
    valMean = valSum/listVals.size();
    valVar = 0.0;
    for(int i=0; i<(int)listVals.size(); ++i)
    {
        double vdiff = listVals[i] - valMean;
        valVar += vdiff*vdiff;
    }
}

template<class TYPE1, class TYPE2>
void FindMinFromPairedListGen(const vector<pair<TYPE1,TYPE2> > &vecListInput, vector<pair<TYPE1,TYPE2> > &listMinItems)
{
    // TYPE1: value (key), TYPE2: can be anything (maybe a pointer for example)
    // there may be multiple items with value (type1) are minimum; listMinItems: contain all such items
    listMinItems.clear();
    if( vecListInput.size() == 0 )
    {
        return;
    }
    TYPE1 valMin = vecListInput[0].first;
    listMinItems.push_back( vecListInput[0] );
    for( int i=1; i<(int)vecListInput.size(); ++i )
    {
        //
        if( vecListInput[i].first < valMin )
        {
            valMin = vecListInput[i].first;
            listMinItems.clear();
            listMinItems.push_back( vecListInput[i] );
        }
        else if( vecListInput[i].first == valMin )
        {
            listMinItems.push_back( vecListInput[i] );
        }
    }
}

template<class TYPE>
void FindRangeInSortedVector( const vector<TYPE> &listSortVals, const TYPE &valLB, const TYPE &valUB, int &posLB, int &posUB)
{
    // given a sorted list, and a range [lb,ub]; want to find the range in the list that contain the list
    // if there is no such range, set as -1
    posLB = 0;
    posUB = (int)listSortVals.size()-1;
    while( listSortVals[posLB] < valLB  )
    {
        ++posLB;
    }
    while( listSortVals[posUB] > valUB  )
    {
        --posUB;
    }
    if( posLB > posUB)
    {
        posLB = -1;
        posUB = -1;
    }
}

template<class TYPE>
void DumpVecWithSpace(const vector<TYPE> &listItems)
{
    // remove the item at the pos
    for(int i=0; i <(int) listItems.size(); ++i )
    {
        cout << listItems[i];
        if( i<(int)listItems.size()-1 )
        {
            cout << "  ";
        }
    }
}

// other utilities
void printvec(vector<double> VEC);
int GetZeroOneDiff(int x, int y);
void GetMatchingPosIntVec(const int val, const vector<int> &listVals, vector<int> &listPos);
void FormUnitVector(int numItems, int posUnit, vector<int> &vecUnit);
void FormZeroVector(int numItems, vector<int> &vecUnit);
bool AreTwoSetsCompatible( const set<int> &set1, const set<int> &set2);
bool IsSetCompatibleWithSets( const set<int> &set1, const set< set<int> > &setSets);
bool IsSignificantFraction(int totNum, int numTypes, int numOneType, double minFrac = -1.0);
void IncAllNumInSet( set<int> &sint );
void IncAllNumInSets( set< set<int> > &setInts);
void GetNonZeroPosofVec( const vector<int> &vec, set<int> &setpos);
void GetDiffPosOfTwoVec(const vector<int> &vec1, const vector<int> &vec2, set<int> &setpos);
int GetSegIndex(int val, const vector<int> &listSegSizes);
void ComplementBoolVec(vector<bool> &listVals);
void GetAllGridPoints( int gridLB, int gridUB, int dimGrid, set< vector<int> > &setGridPts );
//void ReduceContainerSetsForSets(vector<set<int> > &listSets);
void MapIntListToAnother( const vector<int> &vec1, const vector<int> &vec2, map<int,int> &mapVec1IndexToVec2);
void FindEvenDistriPoints(double valMin, double valMax, double valResolution, int maxNumPoints, vector<double> &listChosenVals);

// bits operation
bool IsBitSetInt(int val, int posBit);
int ToggleBitInt(int val, int posBit);

// Prob related utilties
double CalcPoisonProb(double rate, int numEvts);


#endif /* defined(____Utils4__) */
