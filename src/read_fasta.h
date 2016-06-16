/*
 *  read_fasta.h
 *  
 *
 *  Created by Ania on 6/25/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __READ_FASTA_H__
#define __READ_FASTA_H__

// Data types used outside read_fasta.
typedef long long int LONG_INTEGER;
typedef int INTEGER;

// Global variables with the data. Used outside read_fasta.
extern std::vector<LONG_INTEGER> _orbitals;
extern std::vector<INTEGER> _counts;
extern int g_chainLength;
extern int g_numSpecies;
extern std::string *g_nameSpecies;

void print_data();
void read_alignment(std::string fname);
std::string transform_adn_chain_val_to_string(LONG_INTEGER val);

#endif
