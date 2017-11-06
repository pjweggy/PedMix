//
//  PedigreeMixIO.cpp
//
//
//  Created by Yufeng Wu on 12/7/14.
//
//  (1) File input:
//  The file should contain the following: each line refers to a single site
//  Each line should have the following:
//  <base refer to as 1> <freq of allele1 in populaiton A> <freq of allele1 in populaiton B> <Recombinaiton frac to next site> <allele of hap1> <allele of hap2>....
//  base refer to as 1: if alleles are in 0/1 form, must be 1; otherwise refers
//  to the one designated for 1-allele (and this is what the allele frequencey is for)
//  recombination fraction: genetic distance betwen this site and the one to the next; the last site can just have arbitary value (say 0.0); should be between 0 and 1
//  format of haplotypes: each column is a haplotype and there can be spaces
//  between two alleles in the same row (but no empty lines in columns)
//  Caution: assume each row has the same number of alleles
//  Allele: can be eithr 0/1 or A/T/C/G (or a/t/c/g); all other symbols (except space) are invalid and will cause the program to stop
//  For now don't allow missing alleles
//  For convenience, the internal coding uses 0/1, and we use whatever the first
//  base (if A/C/G/T) is used
//  More than one loci is allowed. Use a single line w/ "//" to separate loci
//

#include "PedigreeMixIO.h"
#include <fstream>
#include <sstream>
using namespace std;

//************************************************************************
// Input processing

const double AF_MIN_DEF = 0.0001;

static double FixExtremeVal0to1(double val)
{
    // make sure the value is not too close to 0 or 1
    if( val < AF_MIN_DEF )
    {
        return AF_MIN_DEF;
    }
    else if( val > 1.0-AF_MIN_DEF)
    {
        return 1.0-AF_MIN_DEF;
    }
    else
    {
        return val;
    }
}


static bool IsEmptyAllele(char allele)
{
    return allele == ' ';
}

static int GetBinAllele( char allele, char alleleFirstPos )
{
    // rule: if the allele is already 0/1, just take it
    // otherwise, it is 0 if it matches alleleFirstPos, and 1 otherwise
    // CAUTION: we don't check for wrong alleles (say a column w/ more than 2 allele types)
    if( allele == '0')
    {
        return 0;
    }
    else if(allele == '1')
    {
        return 1;
    }
    else if( allele == alleleFirstPos )
    {
        return 0;
    }
    return 1;   // regardless
}


vector<double> ReadParFromFile(const char *ParFile)
{
	ifstream inFile;
	inFile.open(ParFile);
	vector<double> Par;
	while(inFile.eof() == false )
	{
		char buf[10240];
		inFile.getline( buf, sizeof(buf) );

		if ( strncmp(buf,"phas",4)==0 )
		{
			char * pch[10];
			pch[0] = strtok(buf," ");
			int k=0;
			while ( pch[k] !=NULL)
			{
				pch[++k]=strtok(NULL," ");
			}
			Par.push_back(atof(pch[2]));
			//cout<<buf<<endl;
		}
		else if( strncmp(buf,"recom",5)==0 )
		{
			char * pch[10];
			pch[0] = strtok(buf," ");
			int k=0;
			while ( pch[k] !=NULL)
			{
				pch[++k]=strtok(NULL," ");
			}
			Par.push_back(atof(pch[2]));
		}
		else if( strncmp(buf,"len",3)==0 )
		{
			char * pch[10];
			pch[0] = strtok(buf," ");
			int k=0;
			while ( pch[k] !=NULL)
			{
				pch[++k]=strtok(NULL," ");
			}
			Par.push_back(atof(pch[2]));
		}
		else if( strncmp(buf,"step",4)==0 )
		{
			char * pch[10];
			pch[0] = strtok(buf," ");
			int k=0;
			while ( pch[k] !=NULL)
			{
				pch[++k]=strtok(NULL," ");
			}
			Par.push_back(atof(pch[2]));
		}

		/*
		std::stringstream bufStream(buf);

		double PPelem;
		for (int i=0;i<8;i++)
		{
			bufStream >> PPelem;
			Par.push_back(PPelem);
		}
		*/

	}
	/*
	for (int i=0;i<4;++i)
	{
		cout<<Par[i]<<endl;
	}
	*/

	inFile.close();
	return Par;
}



void ReadPedMixInputFromFile(const char *fileName, vector<vector<PedigreeMixHaplotype> > &listHaps, PopMixingModel &modelPedMix)
{
    //
    ifstream inFile;
    inFile.open( fileName );

    listHaps.clear();
    vector<PedigreeMixHaplotype> listHapsCur;

    bool fFirstLine = true;
    int lociCurr = -1;
    //vector<char> listFirstAlleles;

    // read line by line
    while( inFile.eof() == false )
    {
        char buf[10240];
        inFile.getline( buf, sizeof(buf) );
        if( strlen(buf) == 0 )
        {
            break;
        }
//cout << "Read in one line: " << buf << endl;
        // if finding a //, start a new locus
        if( strncmp( buf, "//", 2 ) == 0 )
        {
            // if there is something read, add the last locus
            if( listHapsCur.size() > 0)
            {
                listHaps.push_back( listHapsCur );
            }

            // clear out to prepare for next locus
            listHapsCur.clear();

            ++lociCurr;
//cout << "Processing loci " << lociCurr << endl;
            continue;
        }

        // if the list of haps is non-empty, then not the first line
        if( listHapsCur.size() > 0 )
        {
            fFirstLine = false;
        }
        else
        {
            fFirstLine = true;
        }

        std::stringstream bufStream(buf);
        int numAlleleRead = 0;

        // first read in header of each site
        char allele1Name;
        bufStream >> allele1Name;
        //listFirstAlleles.push_back(allele1Name);

        double freqPopA, freqPopB, recFrac;
        bufStream >> freqPopA >> freqPopB >> recFrac;
        freqPopA = FixExtremeVal0to1(freqPopA);
        freqPopB = FixExtremeVal0to1(freqPopB);
//cout << "freqPopA: " << freqPopA << ", freqPopB: " << freqPopB << ", recFrac: " << recFrac << endl;
        modelPedMix.AddAlleleFreqs( lociCurr, freqPopA, freqPopB );
        modelPedMix.AddRecFrac( lociCurr, recFrac);

        // now read in the name
        char allele;
        while(  bufStream.get(allele) )
        {
            //char allele;
            //bufStream >> allele;
//cout << "read in char: " << allele << endl;
            if( IsEmptyAllele(allele) == true)
            {
                continue;
            }
//cout << "allele: " << allele << endl;

            //cout << "Found one kmer: " << kmer << ", with freq: " << kmerFreq << endl;
            // if first line, add a new hap
            if( fFirstLine == true )
            {
                // add a new hap
                PedigreeMixHaplotype hap;
                listHapsCur.push_back(hap);
                //listFirstAlleles.push_back(allele);
            }

            // append but first ensure we have enough hapltypes
            int allele0or1 = GetBinAllele( allele, allele1Name );
//cout << "Processing " << allele0or1 <<", listHapsCur size: " << listHapsCur.size() << ", numAlleleRead = " << numAlleleRead  << endl;
            YW_ASSERT_INFO( numAlleleRead <(int)listHapsCur.size(), "FATAL ERROR: some haplotypes have missing alleles1"  );
            listHapsCur[numAlleleRead].AddAllele(allele0or1);

            ++numAlleleRead;
        }
        // ensure the number of alleles read match the num of haps
        YW_ASSERT_INFO( numAlleleRead == (int)listHapsCur.size(), "FATAL ERROR: some haplotypes have missing alleles." );

    }

    // add the last locus
    if( listHapsCur.size() > 0)
    {
        listHaps.push_back( listHapsCur );
    }

    inFile.close();

#if 0
cout << "*** Model read: ";
modelPedMix.Dump();
cout << "*** Haplotypes read: \n";
for(int i=0; i<(int)listHaps.size(); ++i)
{
cout << "----Loci: " << i << endl;
for(int j=0; j<(int)listHaps[i].size(); ++j)
{
listHaps[i][j].Dump();
}
}
#endif
}

