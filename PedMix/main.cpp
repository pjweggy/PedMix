#include <iostream>
#include <fstream>
#include <cstdio>
#include <vector>

#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include <cmath>
using namespace std;


#include <cstring>
#include "Utils3.h"
#include "PedigreeMixLikelihood.h"
#include "PedMixTest.h"
#include "PedigreeMixIO.h"
#include "MixModelExplorer.h"



///////////////////////////////////////////////////////////////////////////////////////////////
// Real options now

const char *versionMIG = "**************************************************************************\n            PedMix ver. 1.1.0 (beta), \n           released on May, 2017 \n**************************************************************************\n";
static int fileArgIndex = -1;
//static int fileArgIndexModel = -1;
static const char *fileInputFile = NULL;            // file for model
//static const char *fileHapFile = NULL;     // file for haplotype
static const char *ParFile =NULL;     //file for parameter
static double mixRatioInput = 0.5;
static int mixGenerations = 10;
static bool fCalcSingleParam = false;

static void Usage()
{
    cout << "Usage: ./PedMix <OPTIONS> <file-name>\n";
	cout << "Example: ./PedMix -g 1 -m 0.5 5 -p parfile testPedMix.inp\n";
    cout << "Options: \n";
	cout << "   -g : number of generations to trace back, for example parents is '-g 1'.\n";
        cout << "   -p parfile : input parameter file to PedMix, including phasing error rate, recombination rate, length and BFGS step size.\n";
    exit(1);
}


static bool CheckArguments(int argc, char **argv)
{
    fileArgIndex = -1;
    if( argc <= 1  )
    {
        return false;
    }
    if( argc == 2)
    {
        fileArgIndex = 1;
        fileInputFile = argv[fileArgIndex];
        //fileArgIndexModel = 1;
        //fileModelFile = argv[fileArgIndexModel];
    }
    else
    {

        for(int i = 1; i< argc; ++i)
        {
            if( argv[i][0] == '-' && argv[i][1] == 'v' )
            {
                //
                //cout << "";
                //SetVerbose(true);
            }
            else if( argv[i][0] == '-' && argv[i][1] == 'g' )
            {
                ++i;
                int numPerfGen = 2;
                sscanf(argv[i], "%d", &numPerfGen);
                SetNumofPerfectPedGens( numPerfGen );
                cout << "Setting the number of perfect generations to be " << numPerfGen << endl;
            }
/*
            else if( argv[i][0] == '-' && argv[i][1] == 'm' )
            {
                // specify two parameters: generaiton time and mix ratio
                //cout << "Turn on verbose mode\n";
                //SetVerbose(true);
                YW_ASSERT_INFO( i+2 < argc, "Not enough parameters" );
                ++i;
                float val1;
                sscanf(argv[i], "%f", &val1);
                mixGenerations = val1;
                ++i;
                float val2;
                sscanf(argv[i], "%f", &val2);
                mixRatioInput = val2;
                fCalcSingleParam = true;
                cout << "Setting the model generaiton time to " << mixGenerations << ", mixing ratio to " << mixRatioInput << endl;
            }
*/
            else if ( argv[i][0] == '-' && argv[i][1] == 'p')
            {
            	// specify parameter.file
            	++i;
            	fileArgIndex = i;
            	ParFile = argv[fileArgIndex];
            	//cout<<"NOW get Parfile"<<endl;
            }

            else if( argv[i][0] != '-')
            {
                // this is the input file of haplotypes
                fileArgIndex = i;
                fileInputFile = argv[fileArgIndex];
            }
            else
            {
                YW_ASSERT_INFO(false, "Not implemented yet");
            }
            fCalcSingleParam = true;
        }
    }
    return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// testiing
void TempTest( const vector<vector<PedigreeMixHaplotype> > &listHaps, PopMixingModel &modelPedMix, vector<double> Par )
{
    //
    //TestCellTreeCostSimple();
    YW_ASSERT_INFO( fileInputFile != NULL , "Haplotype or model files: one of is not set" );
    //TestCellTreeInfGenoFile( fileGenoFile );
//TestLikelihoodDirect();
//exit(1);
    //TestLikelihoodForHaps( listHaps, modelPedMix );

    if( fCalcSingleParam == false )
    {
        // try to search for the model
        PedigreeMixModelExplorer mixModelSearcher( modelPedMix, listHaps );
        mixModelSearcher.Explore();
        // dump out the model
        cout << "*********************************************************************************\n";
        cout << "The optimal model is: \n";
        //cout << "Number of mixing generaitons: " << modelPedMix.GetNumMixGens() << endl;
        //cout << "Fraction of mixing: " << modelPedMix.GetMixRatio(1) << endl;
        cout << "Log-likelihood is: " << mixModelSearcher.GetLogLikeli() << endl;
        //modelPedMix.Dump();
    }
    else
    {
        // just compute for the fixed model
    	//int gen=GetNumofPerfectPedGens();
    	//cout<<gen<<endl;
        TestLikelihoodForHaps( listHaps, modelPedMix, mixGenerations, mixRatioInput, Par );

    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// for ms file only


int main(int argc, char **argv)
{
	// info
	cout << versionMIG << endl;
    if( CheckArguments( argc, argv) == false)
    {
        Usage();
    }

    // read in haplotypes
    vector<vector<PedigreeMixHaplotype> > listHaplotypes;       // may from multiple loci
    const int numPops = 2;
    const int numMixingGen = 10;
    const double ratioMix = 0.5;
    PopMixingModel modelPedMix( numPops, numMixingGen, ratioMix );
    ReadPedMixInputFromFile(fileInputFile, listHaplotypes, modelPedMix);//Dump info: freqA freqB and haplotypes---->modelPedMix
    vector<double> Par=ReadParFromFile(ParFile);
/*
    cout<<"Haplotype size "<<listHaplotypes.size()<<endl;
    for (int k=0;k<(int)listHaplotypes.size();k++)
    {
    	cout<<"this is locus "<<k<<":"<<endl;
    	for (int l=0;l<(int)listHaplotypes[k].size();l++)
    	{
    		cout<<"site "<<l<<":"<<endl;
    		for (int m=0;m<4;m++)
    		{
    			cout<<listHaplotypes[k][l].GetAlleleAt(m)<<",";
    		}
    		cout<<endl;
    	}
    }
*/
    // start timing
    long tstart1 = GetCurrentTimeTick();

    //

    TempTest( listHaplotypes, modelPedMix, Par );


    cout << "Elapsed time = " << GetElapseTime( tstart1 ) << " seconds." << endl;

    // for now, do nothing
    return 0;
}
