/*
 *  read_fasta.cpp
 *  
 *
 *  Created by Ania on 6/25/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cmath>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h> 

#include "read_fasta.h"


bool _calculated;

using namespace std;
typedef long long int LONG_INTEGER;
typedef int INTEGER;

#define NUM_ADN_BASES 4
#define NUM_ADN_OPTIONS 4
enum ADNBases {
	ADENINA, CITOSINA, TIMINA, GUANINA, NONE
};

#define MAX_CHAIN_LENGTH 100000 // max OTUs length
int MAX_NUM_SPECIES;
int SHIFT_CORRIMENT;
static int g_mask;

// static INTEGER g_valsSaved; // The chains that are not repeated...
int g_chainLength = -1;
int g_numSpecies = -1;
static char* g_chainCharsSpeciesADN = 0; // Matrix dim: g_chainLength*g_numSpecies

static LONG_INTEGER *g_chainValsADN = 0; // Matrix dim: g_chainLength*NUM_ADN_OPTIONS
static INTEGER *g_numCountsData = 0; // Matrix dim: g_chainLength*NUM_ADN_OPTIONS
std::string *g_nameSpecies = 0;      // Vector dim: g_numSpecies

static std::string changeADNBases[1][6] = { { "(at)(cg)", "(ac)(gt)",
	"(ag)(ct)", "(ag)", "(ac)", "(at)" } };
std::vector<LONG_INTEGER> _orbitals;
std::vector<INTEGER> _counts;



// Free System Memory
void FreeMemory() {
	if (g_chainCharsSpeciesADN) {
		delete[] g_chainCharsSpeciesADN;
		g_chainCharsSpeciesADN = 0;
	}
	if (g_chainValsADN) {
		delete[] g_chainValsADN;
		g_chainValsADN = 0;
	}
	if (g_numCountsData) {
		delete[] g_numCountsData;
		g_numCountsData = 0;
	}
	if (g_nameSpecies) {
		delete[] g_nameSpecies;
		g_nameSpecies = 0;
	}
}

// Reserve Memory
void ReserveMemory(int numSpecies, int chainLength) { 
	g_chainCharsSpeciesADN = new char[numSpecies * chainLength];
	g_chainValsADN = new LONG_INTEGER[chainLength * NUM_ADN_OPTIONS];
	g_numCountsData = new INTEGER[chainLength * NUM_ADN_OPTIONS];
	g_nameSpecies = new std::string[numSpecies];
}

int MaxNumSpecies() {
	int numBits = sizeof(LONG_INTEGER) * 8;
	int val = -1;
	
	if (NUM_ADN_BASES + 1 <= 0) {
		std::cout << "Num of ADN Bases cannot be negative or zero."
		<< std::endl;
	} else if (NUM_ADN_BASES + 1 <= 2) {
		val = numBits;
	} else if (NUM_ADN_BASES + 1 <= 4) {
		val = numBits / 2;
	} else if (NUM_ADN_BASES + 1 <= 8) {
		val = numBits / 3;
	} else if (NUM_ADN_BASES + 1 <= 16) {
		val = numBits / 4;
	} else {
		std::cout
		<< "Num of ADN Bases superior to 16. Contact the program administrator."
		<< std::endl;
	}
	
	return val;
}



// Get SHIFT corriment info for bit manipulation
int GetShiftCorriment() {
  //	int numBits = sizeof(LONG_INTEGER) * 8;
	int val = -1;
	
	if (NUM_ADN_BASES + 1 <= 0) {
		std::cout << "Num of DNA Bases cannot be negative or zero."
		<< std::endl;
	} else if (NUM_ADN_BASES + 1 <= 2) {
		val = 1;
	} else if (NUM_ADN_BASES + 1 <= 4) {
		val = 2;
	} else if (NUM_ADN_BASES + 1 <= 8) {
		val = 3;
	} else if (NUM_ADN_BASES + 1 <= 16) {
		val = 4;
	} else {
		std::cout
		<< "Num of DNA Bases superior to 16. Contact the program administrator."
		<< std::endl;
	}
	
	return val;
}

void SetGlobalMask() {
	if (NUM_ADN_BASES + 1 <= 0) {
		std::cout << "Num of DNA Bases cannot be negative or zero." << std::endl;
	} else if (NUM_ADN_BASES + 1 <= 2) {
		g_mask = 1;
	} else if (NUM_ADN_BASES + 1 <= 4) {
		g_mask = 3;
	} else if (NUM_ADN_BASES + 1 <= 8) {
		g_mask = 7;
	} else if (NUM_ADN_BASES + 1 <= 16) {
		g_mask = 15;
	} else {
	}
}


inline char ADNBaseValToChar(int adnVal) {
	if (adnVal == 0)
		return 'a';
	else if (adnVal == 1)
		return 'c';
	else if (adnVal == 2)
		return 'g';
	else if (adnVal == 3)
		return 't';
		return 0;
}

inline int ADNBaseCharToVal(int adnChar) {
	if (adnChar == 'a')
		return 0;
	else if (adnChar == 'c')
		return 1;
	else if (adnChar == 'g')
		return 2;
	else if (adnChar == 't')
		return 3;
		return -1;
}

LONG_INTEGER change_adn_bases(LONG_INTEGER value, int col, int fil) {
	std::string change = changeADNBases[fil][col];
	LONG_INTEGER val1 = ADNBaseCharToVal(change[1]);
	LONG_INTEGER val2 = ADNBaseCharToVal(change[2]);
	LONG_INTEGER val3 = ADNBaseCharToVal(change[5]);
	LONG_INTEGER val4 = ADNBaseCharToVal(change[6]);
	
	LONG_INTEGER new_value = 0;
	for (int i = 0; i < g_numSpecies; i++) {
		LONG_INTEGER val = value >> i * SHIFT_CORRIMENT & g_mask;
		if (val == val1) {
			new_value |= (val2 << i * SHIFT_CORRIMENT);
		} else if (val == val2) {
			new_value |= (val1 << i * SHIFT_CORRIMENT);
		} else if (val == val3) {
			new_value |= (val4 << i * SHIFT_CORRIMENT);
		} else if (val == val4) {
			new_value |= (val3 << i * SHIFT_CORRIMENT);
		}
	}
	
	return new_value;
}

std::string transform_adn_chain_val_to_string(LONG_INTEGER val) {
	std::string adnChain;

	for (int i = g_numSpecies - 1; i >= 0; i--) {
		
		LONG_INTEGER value = val >> i * SHIFT_CORRIMENT & g_mask;

		adnChain.push_back(ADNBaseValToChar(value));
	}
	
	return adnChain;
}

bool IsCalculated() {
	return _calculated;
}

// Save Information--

void SaveOrbital(LONG_INTEGER &newOrbital) {
        unsigned int i;
	bool orbFound = false;
	int orbPosCol = -1;
	for (i = 0; i < _orbitals.size() && !orbFound; i++) {
		
		if (_orbitals[i] == newOrbital) {
			orbPosCol = i;
			orbFound = true;
		}
	}
	if (orbFound) {
		_counts[orbPosCol] += 1;
	} else {
		_counts.push_back(1);
		_orbitals.push_back(newOrbital);		
	}
} 


void CountPatterns() {
	
	int _check = 0;
	
	for (int i = 0; i < g_chainLength; i++) {
		LONG_INTEGER valOrbital = 0;
		for (int j = 0; j < g_numSpecies; j++) {
			int val = ADNBaseCharToVal(g_chainCharsSpeciesADN[j * g_chainLength + i]); // returns -1 if the sumbol is not correct (one of the nt)
			
			if(val==-1){
				_check=1; // skip this site in the alignment, it contains signs that are not allowed
				
				j=g_numSpecies; // no need to read further across the species, go to the next site
			}
			else{ 
				valOrbital = (valOrbital << SHIFT_CORRIMENT) | val;
			}
			
		}
		// Save the val (only if the signs were ok)
		if(_check==0){   
			
			SaveOrbital(valOrbital);
		}
		else
		{
			_check=0; 
		};
		
	}
	_calculated = true;
	if(_orbitals.size()==0)
	{
 		std::cout << "Every site in the multiple DNA sequence alignment contains undefined signs. Only bases A, C, G, T are allowed.\n" << std::endl;
		exit(0);	}
	
	
}


// Reads an alignment from a fasta file and returns a structure:
// that contains a vector of patterns and relative frequencies.
void read_alignment(std::string fileName) {
    int pos;
	if ((pos = fileName.find(".fa")) == -1)  // checking the extension .fa
	{
		std::cout << "Invalid format: fasta file should end with \".fa\"" << std::endl;
		exit(-1);
		//fileName.replace(pos, 4, ".fa");
	}
	
	std::ifstream myfile;
	
	myfile.open(fileName.c_str(), ios::in);
	
	INTEGER numSpecies=0;
	//	INTEGER chainLength;
	bool read = true;
	INTEGER min_seq_length = 0; // length of the shortest sequence in the alignment
	INTEGER seqLengthRead = 0; // we have to get the length of the shortest OTU
	
	std::string data;
	
	if (!myfile) {
		
		std::cout <<  "Unable to open the fasta file "<< std::endl;
		exit(1);   // call system to stop
	}
	else
	{
		
		
		FreeMemory();// !!!!!!!!!!!!!
		
		// first pass to count the number of species and the length of the chainLength
		// Load Info
		MAX_NUM_SPECIES = MaxNumSpecies();
		SHIFT_CORRIMENT = GetShiftCorriment();
		SetGlobalMask();
		
		while (read) {
			myfile >> data;
			if (data.size() > 0) {
				if (data[0] != '>' && numSpecies == 0) 
				{ 
					std::cout << "Invalid fasta format: first line should start with \">\" " << std::endl;
					read = false;
					exit(-1);
				}
				if (data[0] == '>') {
					numSpecies++;  // adding species
					
					if (seqLengthRead < min_seq_length )  // if the newly added sequence is shorter than the one before
					{
				        min_seq_length=seqLengthRead; 
					} 
					seqLengthRead=0;  // reset the counter
					
				} else if (myfile.eof()) {
					read = false;
				} else {
					for (unsigned int i = 0; i < data.size(); i++) {
						seqLengthRead++;
					
					}
				}
				
				if ( numSpecies == 1)  // if we there was only a single sequence
				{
					min_seq_length=seqLengthRead;
				}	
			}	 
		}
		
	}
	g_chainLength = min_seq_length;
	std::cout << "Mulitple sequence alignment of length " << g_chainLength << "bp on "<< numSpecies << " taxa." << std::endl;
	////////////////////////////////////////////////// end first pass
	myfile.clear();              // forget we hit the end of file
	
	////////////  second round - goal: save the data, now knowing the number of OTUs and the length of the min.one
	
	
	myfile.seekg(0, ios::beg);   // read the file from the start	
		int trunc_num = 0;
	if (!myfile) {
		
		cerr << "Unable to open file datafile " << std::endl;
		exit(1);   // call system to stop
	}
	else
	{
		// Reserve Memory
		g_numSpecies = numSpecies;
		
		ReserveMemory(numSpecies, g_chainLength);
		// Load Info
	
		read = true;
		INTEGER specie = -1;
		INTEGER valSavedForSpecie = 0;
			
		while (read) {
			
			myfile >> data;  // reads line by line 
			
			if (data.size() > 0) {
				if (data[0] == '>') {
					specie++;
					trunc_num = 0;
					valSavedForSpecie = 0;
					g_nameSpecies[specie] = data;
					std::cout << std::endl << "Reading species " << specie +1 << " " << g_nameSpecies[specie] << std::endl;
				} else if (myfile.eof()) {
					read = false;
					
				} else {
					for (unsigned int i = 0; i < data.size(); i++) {
						if (valSavedForSpecie < g_chainLength) {
							
							g_chainCharsSpeciesADN[specie * g_chainLength
												   + valSavedForSpecie] = tolower(data[i]);

							valSavedForSpecie++;
							
						} else {
							trunc_num++;
							if (trunc_num == 1)
							{
								std::cout << " sequence truncated for " << data[i]  << std::endl;

							}
							else
							{
								std::cout << ", " << data[i] ;
							}
						}
					}
				}
			}
			
		}
		memset(g_numCountsData, 0, NUM_ADN_OPTIONS * g_chainLength
			   * sizeof(INTEGER));
		memset(g_chainValsADN, 0, NUM_ADN_OPTIONS * g_chainLength
			   * sizeof(LONG_INTEGER));
		myfile.close();
	}
	
	CountPatterns();
	
}


// Prints the data: patterns and frequencies.
void print_data() {
        unsigned int i;
	if (!g_chainCharsSpeciesADN) {
		std::cout << "Error: No file Loaded." << std::endl;
		return;
	}
	if (!_calculated) {
		std::cout << "No computation done." << std::endl;
	}
	
	int option = 0;
	
	if (option == 0) {
		double expected = pow( 4, double(g_numSpecies));		
		std::cout << "Observed site patterns: " << _orbitals.size() << " out of " << expected << " possible ones" << std::endl;
				for (i = 0; i < _orbitals.size(); i++) {
			std::cout
			<< "------------------------------------------------------------------"
			<< std::endl;

			std::cout << "P(" << transform_adn_chain_val_to_string(_orbitals[i]) << ") = " << _counts[i] << std::endl;
	
	}

	}	
}

